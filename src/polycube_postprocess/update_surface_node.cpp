#include <fstream>
#include <boost/property_tree/ptree.hpp>

#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/io.h>
#include <jtflib/util/util.h>

#include "../foldfree/lscm_ipopt_foldfree.h"
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"
#include "../hex_param/io.h"
#include "../common/util.h"
#include "../common/vtk.h"
#include "util.h"
#include "io.h"

using namespace std;
using namespace zjucad::matrix;
using boost::property_tree::ptree;

int uv_to_3dcoord(const char * orig_uv,
                  const char * output_uv,
                  boost::unordered_map<size_t, matrix<double> > & update_node)
{
  ifstream ifs_input_uv(orig_uv);
  ifstream ifs_output_uv(output_uv);

  update_node.clear();

  matrix<double> uv_basis = zeros<double>(3,2);
  matrix<double> points(3,3);

  string temp;size_t idx;
  for(size_t pi = 0; pi < 3; ++pi){
    ifs_input_uv >> temp >> idx >> points(0,pi) >> points(1,pi) >> points(2,pi);
    if(temp != "#"){
      cerr << "# [error] wrong uv format." << endl;
      return __LINE__;
    }
  }

  uv_basis(colon(),0) = points(colon(),1) - points(colon(),0);
  uv_basis(colon(),1) = points(colon(),2) - points(colon(),0);

  double len_0 = norm(uv_basis(colon(),0));
  double len_1 = norm(uv_basis(colon(),1));
  if(len_0 < 1e-6) len_0 = 1.0;
  if(len_1 < 1e-6) len_1 = 1.0;

  uv_basis(colon(),0) /= len_0;
  uv_basis(colon(),1) /= len_1;

  matrix<double> uvw(3,1);
  matrix<double> new_point = zeros<double>(3,1);
  while(!ifs_output_uv.eof()){
    ifs_output_uv >> idx >> uvw[0] >> uvw[1] >> uvw[2];
    if(ifs_output_uv.eof()) break;
    new_point = points(colon(),0) + uvw[0] * uv_basis(colon(),0) +
                uvw[1] * uv_basis(colon(),1);
    update_node[idx] = new_point;
  }
  return 0;
}


int find_uv_basis_point(
    const deque<pair<size_t,size_t> > & boundary_edges,
    const matrix<double> & node,
    vector<size_t> & uv_basis_point)
{
  uv_basis_point.clear();
  vector<double> cos_angle;
  for(size_t ei = 0; ei +1 < boundary_edges.size(); ++ei){
    const pair<size_t,size_t> & current_edge = boundary_edges[ei];
    const pair<size_t,size_t> & next_edge = boundary_edges[ei+1];
    if(current_edge.second != next_edge.first){
      cerr << "# [error] this chain is not consistent." << endl;
      return __LINE__;
    }
    matrix<double> e1 = node(colon(),current_edge.second)
                        - node(colon(), current_edge.first);
    matrix<double> e2 = node(colon(),next_edge.second) -
                        node(colon(), next_edge.first);

    double len_1 = norm(e1); if(len_1 < 1e-5) len_1 = 1;
    double len_2 = norm(e2); if(len_2 < 1e-5) len_2 = 1;
    e1 /= len_1; e2 /= len_2;
    cos_angle.push_back(fabs(dot(e1,e2)));
    if(fabs(dot(e1, e2)) < 1e-1){
      uv_basis_point.push_back(current_edge.second);
      uv_basis_point.push_back(current_edge.first);
      uv_basis_point.push_back(next_edge.second);
      break;
    }
  }
  if(uv_basis_point.empty()){
    cerr << "# [error] strange can not find orthogonal edges." << endl;
    {// vis
      ofstream ofs("no_orthogonal_edge_boundary.vtk");
      vector<size_t> lines;
      for(size_t ei = 0; ei < boundary_edges.size(); ++ei){
        const pair<size_t,size_t> & one_edge = boundary_edges[ei];
        lines.push_back(one_edge.first);
        lines.push_back(one_edge.second);
      }
      line2vtk(ofs, &node[0], node.size(2), &lines[0], lines.size()/2);
      for(size_t i = 0; i < cos_angle.size(); ++i)
        cerr << cos_angle[i] << endl;
    }
    return __LINE__;
  }
  return 0;
}



int dump_flipped_patches(
    const vector<matrix<size_t> > &surface_patches,
    const matrix<double> &orig_node,
    const matrix<double> &polycube_node,
    const vector<size_t> &flipped_patches,
    const vector<vector<size_t> > &boundary_nodes,
    const vector<vector<size_t> > &uv_basis_point)
{
  if(orig_node.size() != polycube_node.size()){
    cerr << "# [error] input orig node size does not equals polycube." << endl;
    return __LINE__;
  }
  if(flipped_patches.size() != boundary_nodes.size() ||
     boundary_nodes.size() != uv_basis_point.size()){
    cerr << "# [error] input configuration is not compatiable." << endl;
    return __LINE__;
  }


  for(size_t i = 0; i < flipped_patches.size(); ++i){
    ostringstream os;
    os << i;
    string name = "patch_" + os.str();
    jtf::mesh::save_obj((name+".obj").c_str(), surface_patches[flipped_patches[i]],orig_node);
    jtf::mesh::save_obj((name+".init.obj").c_str(), surface_patches[flipped_patches[i]], polycube_node);
    save_to_uv((name+".uv").c_str(),uv_basis_point[i],boundary_nodes[i], polycube_node);
  }
  return 0;
}

int check_flipped_patches(
    const vector<matrix<size_t> > &surface_patches,
    const matrix<double> &node,
    vector<size_t> &flipped_patches,
    vector<vector<size_t> > & boundary_nodes,
    vector<vector<size_t> > & uv_basis_point)
{
  flipped_patches.clear();
  boundary_nodes.clear();
  boost::unordered_set<size_t> one_boundary;
  matrix<double> face_normal;
  for(size_t pi = 0; pi < surface_patches.size(); ++pi){
    const matrix<size_t> & one_patch = surface_patches[pi];
    unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
          jtf::mesh::edge2cell_adjacent::create(one_patch));
    if(!ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }
    jtf::mesh::cal_face_normal(one_patch, node, face_normal);
    bool is_flipped_patch = false;

    vector<pair<size_t,size_t> > boundary_edges;

    for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
      const pair<size_t,size_t> & one_edge = ea->edges_[ei];
      const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
      if(ea->is_boundary_edge(tri_pair)){
        boundary_edges.push_back(one_edge);
        one_boundary.insert(one_edge.first);
        one_boundary.insert(one_edge.second);
      }else{
        if(dot(face_normal(colon(), tri_pair.first),
               face_normal(colon(), tri_pair.second)) < 0)
          is_flipped_patch = true;
      }
    }
    if(boundary_edges.empty()){
      cerr << "# [error] strange, can not find boundary edges." << endl;
      //return __LINE__;
      continue;
    }
    if(is_flipped_patch) {
      flipped_patches.push_back(pi);
      vector<deque<pair<size_t,size_t> > > chains;
      jtf::util::extract_chain_from_edges(boundary_edges, chains);
      // since one patch may be genus != 0
//      if(chains.size() != 1){
//        cerr << "# [error] strange, chains of one patch should be only one." << endl;
//        {// vis
//          ofstream ofs("strange_patch.vtk");
//          tri2vtk(ofs, &node[0], node.size(2), &surface_patches[pi][0], surface_patches[pi].size(2));
//        }
//        return __LINE__;
//      }
      vector<size_t> uv_basis;
      int rtn = find_uv_basis_point(chains[0], node, uv_basis);
      if(rtn)    return __LINE__;
      vector<size_t> one_boundary_vec(one_boundary.begin(), one_boundary.end());
      boundary_nodes.push_back(one_boundary_vec);
      uv_basis_point.push_back(uv_basis);
    }
    one_boundary.clear();
  }
  return 0;
}

int recovery_from_optimized_patch(
    const vector<matrix<size_t> > &surface_patches,
    const vector<size_t> &flipped_patches,
    matrix<double> &polycube_node)
{
  boost::unordered_map<size_t, matrix<double> > update_node;
  for(size_t i = 0; i < flipped_patches.size(); ++i){
    //if(i > 29) continue;
    cerr << i << endl;
    ostringstream os;
    os << i;
    string patch_name = "patch_" + os.str();
    bool rtn = jy::ipopt_foldfree_parameterization(
                (patch_name + ".obj").c_str(), (patch_name + ".uv").c_str(), "temp.uv");
    if(rtn){ // success
      uv_to_3dcoord((patch_name + ".uv").c_str(), "temp.uv", update_node);
      for(boost::unordered_map<size_t,matrix<double> >::const_iterator cit =
          update_node.begin(); cit != update_node.end(); ++cit){
        polycube_node(colon(), cit->first) = cit->second;
       }
    }
  }
  return 0;
}


int update_surface_node(ptree & pt)
{
  jtf::mesh::meshes orig_tm, polycube_tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("orig_tet.value").c_str(),
                               &orig_tm.node_, & orig_tm.mesh_))
    return __LINE__;

  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("polycube_tet.value").c_str(),
                               &polycube_tm.node_, & polycube_tm.mesh_))
    return __LINE__;

  if(orig_tm.mesh_.size() != polycube_tm.mesh_.size()){
    cerr << "# [error] orig tet size is not compatiable with polycube tet." << endl;
    return __LINE__;
  }

  boost::unordered_map<size_t,size_t> surface_type;
  if(load_surface_restricted_type_static(
       pt.get<string>("surface_type.value").c_str(),
       surface_type))
    return __LINE__;

  vector<matrix<size_t> > surface_patches;
  boost::unordered_set<pair<size_t,size_t> > patch_boundary;

  convert_surface_type_to_surface_patches(
        orig_tm.mesh_, surface_type, surface_patches, patch_boundary);

  const size_t iter_b = pt.get<size_t>("iter_b.value", 1);
  smooth_boundary(patch_boundary, orig_tm.node_, polycube_tm.node_, iter_b);
  //smooth_patch(surface_patches, patch_boundary, polycube_tm.node_, 10);
  vector<size_t> flipped_patches;
  vector<vector<size_t> >  boundary_nodes, uv_basis_point;

  pt.put("flip_or_all.value", "dump out flipped patch or all, [0/1]");
  const size_t only_flipped_patch = pt.get<size_t>("flip_or_all.value",0);

  if(only_flipped_patch == 0){
    if(check_flipped_patches(surface_patches, polycube_tm.node_, flipped_patches,
                          boundary_nodes, uv_basis_point))
      return __LINE__;
  }else if(only_flipped_patch == 1){
    flipped_patches.resize(surface_patches.size());
    for(size_t pi = 0; pi < flipped_patches.size(); ++pi){
      flipped_patches[pi] = pi;
    }
  }else {
    cerr << "# [error] wrong input for flip_or_all." << endl;
    return __LINE__;
  }

  dump_flipped_patches(surface_patches, orig_tm.node_, polycube_tm.node_,
                       flipped_patches, boundary_nodes, uv_basis_point);

  matrix<double> node = polycube_tm.node_;
  const size_t limit_patch_number = pt.get<size_t>("limit_patch.value",-1);

  vector<size_t > limit_patches;
  if(limit_patch_number == -1 || limit_patch_number >= surface_patches.size())

    limit_patches = flipped_patches;
  else {
    limit_patches.resize(limit_patch_number);
    for(size_t pi = 0; pi < limit_patch_number; ++pi){
      limit_patches[pi] = flipped_patches[pi];
    }
  }

  recovery_from_optimized_patch(surface_patches, limit_patches, node);

  cerr << "# [info] node diff " << norm(node - polycube_tm.node_) << endl;
  jtf::mesh::tet_mesh_write_to_zjumat("optimized_polycube_tet.tet", &node,
                           &polycube_tm.mesh_);
  return 0;
}
