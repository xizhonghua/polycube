#include <boost/property_tree/ptree.hpp>

#include <string>
#include <fstream>
#include <numeric>
#include <hjlib/function/func_aux.h>
#include <hjlib/math/polar.h>
#include <stack>

#include <jtflib/mesh/io.h>
#include <jtflib/optimizer/optimizer.h>

#include <zjucad/optimizer/optimizer.h>
#include <zjucad/ptree/ptree.h>

#include "../tetmesh/util.h"
#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"
#include "../common/util.h"
#include "../common/IO.h"
#include <jtflib/mesh/mesh.h>

#include "util.h"
#include "io.h"

#include "polycube_surface_func.h"

#include "../mesh_func/tri-area.h"
#include "tet_func.h"
#include "../mesh_func/tri-normal.h"

#include "util.h"

#include "../tet_mesh_sxx/tet_mesh_sxx.h"
#include "adaptive.h"
#include "../tetmesh/subdivide_tet.h"

using namespace std;
using boost::property_tree::ptree;
using namespace hj::function;
using namespace zjucad::matrix;

double arap_w_2;
double fix_w_2;
double fl_w_2;
double patch_w_2;
double group_w_2;
double normal_fix_w;
double vol_w_2;

size_t split_percent;

hj::function::function_t<double, int32_t> *
build_polycube_patch_function(
    const matrix<double> &node,
    const matrix<size_t> &tet,
    const vector<matrix<size_t> > &surface_patches,
    const vector<matrix<double> > & patch_normal,
    const boost::unordered_set<pair<size_t,size_t> > &patch_boundary,
    const vector<boost::unordered_set<size_t> > & fnode_group,
    const vector<vector<size_t> > * select_flipped_face_ptr = 0)
{
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
    cerr  << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return 0;
  }

  matrix<size_t> faces;
  get_outside_face(*fa, faces, true);

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(faces));
  if(!ea.get()){
    cerr  << "# [error] can not build edge2triangle adjacent." << endl;
    return 0;
  }

  //jtf::mesh::save_to_obj("outside_face.obj", node, faces);

  //////////////////////////////////////////////////////////////////
  ////////////////  add function term  /////////////////////////////
  boost::shared_ptr<function_t<double, int32_t> > tetmesh_arap(
        build_tetmesh_arap_func(node, tet, arap_w_2));

  boost::shared_ptr<function_t<double, int32_t> > fix_pos(
        build_fix_zero_node_func(node, node(colon(),0), fix_w_2));

  boost::shared_ptr<function_t<double, int32_t> > feature_line_func(
        build_feature_line_func(node, faces, *ea, patch_boundary, fl_w_2));

  boost::shared_ptr<function_t<double, int32_t> > surface_patch_algin_func(
        build_surface_patch_align_func(node, surface_patches, patch_w_2));

  boost::shared_ptr<function_t<double, int32_t> > group_func(
        build_group_equation_func(node, fnode_group, group_w_2));

  boost::shared_ptr<function_t<double, int32_t> > patch_normal_fix_func(
        build_patch_normal_fix_func(node, patch_normal, surface_patches, normal_fix_w));


  boost::shared_ptr<vector<boost::shared_ptr<function_t<double,int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double,int32_t> > > );

  if(tetmesh_arap.get() && non_zero(arap_w_2)){
    funcs->push_back(tetmesh_arap);
    cerr << "# [info] add jtf::mesh::mesh arap function, weight " << arap_w_2 << endl;
  }

  if(fix_pos.get() && non_zero(fix_w_2)){
    funcs->push_back(fix_pos);
    cerr << "# [info] add fix first node zero function, weight "
         << fix_w_2 << endl;
  }

  if(feature_line_func.get() && non_zero(fl_w_2)){
    funcs->push_back(feature_line_func);
    cerr << "# [info] add feature line align function, weight " << fl_w_2 << endl;
  }

  if(surface_patch_algin_func.get() && non_zero(patch_w_2)){
    funcs->push_back(surface_patch_algin_func);
    cerr << "# [info] add surface_patch_algin_function, weight "
         << patch_w_2 << endl;
  }

  if(group_func.get() && non_zero(group_w_2)){
    funcs->push_back(group_func);
    cerr << "# [info] add group equation function, weight " << group_w_2 << endl;
  }

  if(patch_normal_fix_func.get() && non_zero(normal_fix_w)){
    funcs->push_back(patch_normal_fix_func);
    cerr << "# [info] add normal fix function, weight " << normal_fix_w << endl;
  }

  return new_catenated_function<double, int32_t>(funcs);
}


int check_patch_and_boundary(
    const matrix<double> & node,
    const jtf::mesh::face2tet_adjacent &fa,
    const matrix<size_t> & outside_face,
    boost::unordered_set<pair<size_t,size_t> > & boundary,
    const size_t idx = 0)
{
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  matrix<double> face_normal(3, outside_face.size(2));

  jtf::mesh::cal_face_normal(outside_face, node, face_normal);

  matrix<double> edge_error = zeros<double>(ea->edges_.size(),1);

  vector<pair<size_t,size_t> > err_edge;
  vector<double> diff_vec;
  for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
    const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
    pair<size_t,size_t> one_edge = ea->edges_[ei];
    const double diff = dot(face_normal(colon(), tri_pair.first),
                            face_normal(colon(), tri_pair.second));
    if(one_edge.first > one_edge.second)
      swap(one_edge.first, one_edge.second);
    boost::unordered_set<pair<size_t,size_t> >::const_iterator cit =
        boundary.find(one_edge);
    if(cit == boundary.end()){
      if(fabs(diff-1) > 1e-2){
        cerr << "# [error] dot of edge " << one_edge.first << " "
             << one_edge.second << " should be 1, but " << diff << endl;
        err_edge.push_back(tri_pair);
        diff_vec.push_back(fabs(diff-1.0));
        diff_vec.push_back(fabs(diff-1.0));
      }
    }else{
      if(fabs(diff) > 1e-2){
        cerr << "# [error] dot of edge " << one_edge.first << " "
             << one_edge.second << " should be 0, but " << diff << endl;
        err_edge.push_back(tri_pair);
        diff_vec.push_back(fabs(diff));
        diff_vec.push_back(fabs(diff));
      }
    }
  }

  {
    ostringstream os;
    os << idx;
    string name = "large_error_face_pair_" + os.str() + ".vtk" ;
    ofstream ofs(name.c_str());
    vector<size_t> faces;
    for(size_t fi = 0; fi < err_edge.size(); ++fi){
      const pair<size_t,size_t> & tri_pair = err_edge[fi];
      faces.insert(faces.end(), outside_face(colon(), tri_pair.first).begin(),
                   outside_face(colon(), tri_pair.first).end());
      faces.insert(faces.end(), outside_face(colon(), tri_pair.second).begin(),
                   outside_face(colon(), tri_pair.second).end());
    }
    tri2vtk(ofs, &node[0], node.size(2), &faces[0], faces.size()/3);
    cell_data(ofs,&diff_vec[0], diff_vec.size(), "err");
  }

  return 0;
}

int update_fnode_group(const matrix<size_t> &new_node_parent_id,
                       const size_t &orig_node_num,
                       vector<boost::unordered_set<size_t> > &fnode_group)
{
  if(new_node_parent_id.size() == 0) return 0;

  for(size_t ei = 0; ei < new_node_parent_id.size(1); ++ei){
    const pair<size_t,size_t> one_edge(new_node_parent_id(ei,0),
                                       new_node_parent_id(ei,1));
    for(size_t di = 0; di < 3; ++di){
      for(size_t gi = 0; gi < fnode_group.size(); ++gi){
        boost::unordered_set<size_t> & one_group = fnode_group[gi];
        if(one_group.find(3 * one_edge.first + di) != one_group.end() &&
           one_group.find(3 * one_edge.second + di) != one_group.end()){
          one_group.insert(3 * (orig_node_num + ei) + di);
          break;
        }
      }
    }
  }
  return 0;
}

int update_patch_boundary(const matrix<size_t> &new_node_parent_id,
                          const size_t & orig_node_number,
                          boost::unordered_set<pair<size_t,size_t> > &patch_boundary)
{
  if(new_node_parent_id.size() == 0) return 0;

  boost::unordered_map<pair<size_t,size_t>,size_t> edge2idx;
  for(size_t ei = 0; ei < new_node_parent_id.size(1); ++ei){
    pair<size_t,size_t> one_edge(new_node_parent_id(ei,0), new_node_parent_id(ei,1));
    if(one_edge.first > one_edge.second)
      swap(one_edge.first, one_edge.second);
    edge2idx[one_edge] = ei;
  }

  boost::unordered_set<pair<size_t,size_t> > new_boundary;
  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit
      = patch_boundary.begin(); cit != patch_boundary.end(); ++cit){
    pair<size_t,size_t> one_edge = *cit;
    if(one_edge.first > one_edge.second)
      swap(one_edge.first, one_edge.second);
    boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator ccit
        = edge2idx.find(one_edge);
    if(ccit == edge2idx.end())
      new_boundary.insert(one_edge);
    else{
      new_boundary.insert(make_pair(one_edge.first, ccit->second + orig_node_number));
      new_boundary.insert(make_pair(one_edge.second, ccit->second + orig_node_number));
    }
  }
  patch_boundary = new_boundary;
  return 0;
}

//! This function takes such an assumption that all new_points are insert behind original nodes
int update_surface_patch(const matrix<size_t> &tet,
                         const matrix<size_t> &new_node_parent_id,
                         const size_t orig_node_number,
                         const jtf::mesh::face2tet_adjacent & fa_orig,
                         vector<matrix<size_t> >  &surface_patches)
{
  if(new_node_parent_id.size() == 0) return 0;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  matrix<size_t> outside_face;
  get_outside_face(*fa, outside_face);

  boost::unordered_map<vector<size_t>, vector<size_t> > orig_face2_new_face;

  for(size_t fi = 0; fi < outside_face.size(2); ++fi){
    set<size_t> one_face;
    for(size_t pi = 0; pi < outside_face.size(1); ++pi){
      const size_t &point_idx = outside_face(pi, fi);
      if(point_idx >= orig_node_number)
        one_face.insert(new_node_parent_id(point_idx - orig_node_number,colon()).begin(),
                        new_node_parent_id(point_idx - orig_node_number,colon()).end());
      else
        one_face.insert(point_idx);
    }
    if(one_face.size() != 3){
      cerr << "# [error] strange, this face contains " << one_face.size() << endl;
    }
    vector<size_t> one_face_vec(one_face.begin(), one_face.end());

    // recheck
    const size_t face_idx = fa_orig.get_face_idx(&one_face_vec[0]);
    if(face_idx == -1){
      cerr << "# [error] " << one_face_vec[0] << " " << one_face_vec[1]
           << " " << one_face_vec[2] << " does not belong to original tet." << endl;
      return __LINE__;
    }
    orig_face2_new_face[one_face_vec].push_back(fi);
  }

  // after build orig_face2new_face mapping, we should update the output surface_patches

  for(size_t pi = 0; pi < surface_patches.size(); ++pi){
    matrix<size_t> & one_patch = surface_patches[pi];
    vector<size_t> new_patch;
    for(size_t fi = 0; fi < one_patch.size(2); ++fi){
      vector<size_t> one_face(one_patch(colon(),fi).begin(), one_patch(colon(),fi).end());
      sort(one_face.begin(), one_face.end());
      boost::unordered_map<vector<size_t>, vector<size_t> >::const_iterator cit
          = orig_face2_new_face.find(one_face);
      if(cit == orig_face2_new_face.end()){
        cerr << "# [error] can not find face mapping. " << endl;
        return __LINE__;
      }
      const vector<size_t> & split_face_idx = cit->second;
      for(size_t sfi = 0; sfi < split_face_idx.size(); ++sfi){
        new_patch.insert(new_patch.end(), outside_face(colon(), split_face_idx[sfi]).begin(),
                         outside_face(colon(), split_face_idx[sfi]).end());
      }
    }
    one_patch.resize(3, new_patch.size()/3);
    copy(new_patch.begin(), new_patch.end(), one_patch.begin());
  }

  return 0;
}

int polycube_with_type(ptree &pt)
{
  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm.node_, &tm.mesh_))
    return __LINE__;

  cerr << "# read in tet " << tm.node_.size(2) << " " << tm.mesh_.size(2) << endl;

  pt.put("arap_w.desc", "arap_weighting, default is 1");
  pt.put("fix_w.desc", "fix zero node weighting, default is 1");
  pt.put("fl_w.desc", "align patch boundary, default is 0");
  pt.put("patch_w.desc", "patch align weight, defaut is 0");
  pt.put("iter_w.desc", "iteration number, default is 8");
  pt.put("init.desc", "init_tet");
  pt.put("group_w.desc", "equation group node weight, default is 1");
  pt.put("surface_type.desc", "surface type file.");
  pt.put("normal_fix_w.desc", "normal_fix_w, default is 0");

  boost::unordered_map<size_t,size_t> surface_type;
  if(load_surface_restricted_type_static(
       pt.get<string>("surface_type.value").c_str(), surface_type))
    return __LINE__;

  vector<boost::unordered_set<size_t> >  fnode_group;
  if(zjucad::has("group_file.value",pt)){
    if(load_fnode_group_static(
         pt.get<string>("group_file.value").c_str(), tm.mesh_, fnode_group))
      return __LINE__;
  }

  vector<matrix<size_t> > surface_patches;
  boost::unordered_set<pair<size_t,size_t> > patch_boundary;
  convert_surface_type_to_surface_patches(
        tm.mesh_, surface_type, surface_patches, patch_boundary);


  //#define VISUAL 0
#ifdef VISUAL
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh));
  matrixst outside_face;
  get_outside_face(*fa, outside_face);
#endif

  vector<matrix<double> > patch_normal;
  matrixd node = tm.node_;

  if(zjucad::has("init.value",pt)){
    jtf::mesh::meshes tm_init;
    if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("init.value").c_str(),
                                 &tm_init.node_, &tm_init.mesh_))
      return __LINE__;
    if(tm_init.mesh_.size() != tm.mesh_.size()){
      cerr << "# [error] init tet is not compatiable with orig tet." << endl;
      return __LINE__;
    }

    matrix<double> eye_mat = eye<double>(3);
    patch_normal.resize(surface_patches.size());
    matrix<double> normal = zeros<double>(3,1), avg_normal = zeros<double>(3,1);
    vector<size_t> surface_with_normal_dir;
    vector<size_t> surface_type_with_normal_dir;

    for(size_t pi = 0; pi < surface_patches.size(); ++pi){

      avg_normal = zeros<double>(3,1);
      vector<pair<double, size_t> > normal_diff(6);
      surface_with_normal_dir.insert(
            surface_with_normal_dir.end(), surface_patches[pi].begin(),
            surface_patches[pi].end());

      for(size_t fi = 0; fi < surface_patches[pi].size(2); ++fi){
        const matrix<size_t> & one_patch = surface_patches[pi];
        zjucad::matrix::matrix<double> node_ = tm_init.node_(colon(), one_patch(colon(),fi));
        calc_tri_normal_(&normal[0], &node_[0]);
        if(!isfinite(normal[0]) || !isfinite(normal[1]) || !isfinite(normal[2]))
          continue;
        avg_normal += normal;
      }

      double len = norm(avg_normal);
      if(len < 1e-6) len = 1;
      avg_normal /= len;

      for(size_t di = 0; di < 6; ++di){
        normal_diff[di] = make_pair(dot(avg_normal,  eye_mat(colon(),di/2) * (di%2==0?1.0:-1.0)),di);
      }
      sort(normal_diff.begin(), normal_diff.end());
      const size_t nearest_axis = normal_diff.back().second;
      patch_normal[pi] = eye_mat(colon(), nearest_axis/2) * (nearest_axis%2==0?1.0:-1.0);

      for(size_t fi = 0; fi < surface_patches[pi].size(2); ++fi)
        surface_type_with_normal_dir.push_back(nearest_axis);
    }

    {
      ofstream ofs("surface_with_normal_dir.vtk");
      tri2vtk(ofs, &tm_init.node_[0], tm_init.node_.size(2),
              &surface_with_normal_dir[0], surface_with_normal_dir.size()/3);
      cell_data(ofs, &surface_type_with_normal_dir[0], surface_type_with_normal_dir.size(), "normal_dir");
    }

    {// visual
      ostringstream os;
      os << "surface_patch_type.vtk";
      ofstream ofs(os.str().c_str());
      vector<size_t> surface_patch_vec;
      vector<size_t> surface_patch_type_vec;
      for(size_t pi = 0; pi < surface_patches.size(); ++pi){
        const matrix<size_t> & one_patch = surface_patches[pi];
        surface_patch_vec.insert(surface_patch_vec.end(), one_patch.begin(), one_patch.end());
        for(size_t fi = 0; fi < one_patch.size(2); ++fi){
          surface_patch_type_vec.push_back(pi);
        }
      }
      tri2vtk(ofs, &tm_init.node_[0], tm_init.node_.size(2), &surface_patch_vec[0], surface_patch_vec.size()/3);
      cell_data(ofs, &surface_patch_type_vec[0], surface_patch_type_vec.size(), "type");
    }

    pt.put("use_init.desc","whether use polycube node as init value [y/n]");
    string use_init_node = pt.get<string>("use_init.value");
    if(use_init_node == "y" || use_init_node == "Y"){
      cerr << "# [info] use init node." << endl;
      node = tm_init.node_;
    }

    if(!fnode_group.empty()){
      for(size_t gi = 0; gi < fnode_group.size(); ++gi){
        const boost::unordered_set<size_t> & one_group = fnode_group[gi];
        double avg = 0;
        for(boost::unordered_set<size_t>::const_iterator cit = one_group.begin();
            cit != one_group.end(); ++cit){
          avg += node[*cit];
        }
        avg /= one_group.size();
        for(boost::unordered_set<size_t>::const_iterator cit = one_group.begin();
            cit != one_group.end(); ++cit){
          node[*cit] = avg;
        }
      }
    }

  }

#ifdef VISUAL
  { // visual
    {
      ofstream ofs("feature_line.vtk");
      vector<size_t> edges;

      for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit =
          patch_boundary.begin(); cit != patch_boundary.end(); ++cit){
        edges.push_back(cit->first);
        edges.push_back(cit->second);

      }
      line2vtk(ofs, &node[0], node.size(2), &edges[0], edges.size()/2);
    }

    {
      ofstream ofs("patch.vtk");
      vector<size_t> patch_faces_vis;
      vector<size_t> patch_type;
      for(size_t pi = 0; pi < surface_patches.size(); ++pi){
        patch_faces_vis.insert(patch_faces_vis.end(), surface_patches[pi].begin(),
                               surface_patches[pi].end());
        patch_type.insert(patch_type.end(), surface_patches[pi].size(2), pi);
      }
      tri2vtk(ofs, &node[0], node.size(2), &patch_faces_vis[0], patch_faces_vis.size()/3);
      cell_data(ofs, &patch_type[0], patch_type.size(), "patch_idx");
    }
  }
#endif

  vol_w_2 = pt.get<double>("vol_w.value",0.0);
  arap_w_2 = pt.get<double>("arap_w.value",1);
  fix_w_2 = pt.get<double>("fix_w.value",1);
  fl_w_2 = pt.get<double>("fl_w.value",0);
  patch_w_2 = pt.get<double>("patch_w.value",0);
  group_w_2 = pt.get<double>("group_w.value",1);
  normal_fix_w = pt.get<double>("normal_fix_w.value",0);
  split_percent = pt.get<size_t>("split_percent_w.value",1);

  const size_t iter_w = pt.get<size_t>("iter_w.value",8);
  //matrix<double> R;

  //  if(pt.get_child_optional("subdivide.value")) {
  //    cerr << "# beg subdivide: " << tm.node.size(2) << " " << tm.mesh.size(2) << endl;
  //    matrix<double> *p_node[2] = {&tm.node, &node};
  //    matrix<size_t> c_tet, node_parent;
  //    matrix<double> c_node;
  //    time_t beg = clock();
  //    boundary_uniform_subdivide(tm.node, tm.mesh, c_node, c_tet, node_parent);
  //    cout << "# subd time: " << (clock()-beg)/double(CLOCKS_PER_SEC) << endl;
  //    tm.mesh = c_tet;
  //    tm.node = c_node;

  //    subdivide_top2geo(node, node_parent, c_node);
  //    node = c_node;
  //    ofstream ofs("subdivide.vtk");
  //    tet2vtk(ofs, &tm.node[0], tm.node.size(2), &tm.mesh[0], tm.mesh.size(2));
  //    jtf::mesh::tet_mesh_write_to_zjumat("subdivide.tet", &tm.node, &tm.mesh);
  //    cerr << "# end subdivide: " << tm.node.size(2) << " " << tm.mesh.size(2) << endl;
//    }chr


  matrix<size_t> new_node_parent_id;
  pt.put("split.desc", "split jtf::mesh::mesh while deformation, [0/1] no_split/split");
  const size_t split_or_not = pt.get<size_t>("split.value",0);



  for(size_t i = 0; i < iter_w; ++i){
    normal_fix_w *= 2;
    group_w_2 *= 2;

    if(!fnode_group.empty()){
      for(size_t gi = 0; gi < fnode_group.size(); ++gi){
        const boost::unordered_set<size_t> & one_group = fnode_group[gi];
        double avg = 0;
        for(boost::unordered_set<size_t>::const_iterator cit = one_group.begin();
            cit != one_group.end(); ++cit){
          avg += node[*cit];
        }
        avg /= one_group.size();
        for(boost::unordered_set<size_t>::const_iterator cit = one_group.begin();
            cit != one_group.end(); ++cit){
          node[*cit] = avg;
        }
      }
    }

#if 1 // need subdivided
    {
      if(0){
        const size_t orig_node_num = tm.node_.size(2);
        unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
        if(!fa.get()){
          cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
          return __LINE__;
        }

        //local_remesh_tet(tm.mesh, tm.node, node, fnode_group);

        //        cerr << "# [info] split_percent, " << split_percent << endl;
        //        split_tet_at_large_arap_distortion(tm.mesh, tm.node, node,new_node_parent_id);
        //        update_surface_patch(tm.mesh, new_node_parent_id, orig_node_num, *fa, surface_patches);
        //        update_patch_boundary(new_node_parent_id, orig_node_num, patch_boundary);
        //        update_fnode_group(new_node_parent_id, orig_node_num, fnode_group);

        {// visual
          ostringstream os;
          os << "surface_patch_" << i << ".vtk";
          ofstream ofs(os.str().c_str());
          vector<size_t> surface_patch_vec;
          vector<size_t> surface_patch_type_vec;
          for(size_t pi = 0; pi < surface_patches.size(); ++pi){
            const matrix<size_t> & one_patch = surface_patches[pi];
            surface_patch_vec.insert(surface_patch_vec.end(), one_patch.begin(), one_patch.end());
            for(size_t fi = 0; fi < one_patch.size(2); ++fi){
              surface_patch_type_vec.push_back(pi);
            }
          }
          tri2vtk(ofs, &node[0], node.size(2), &surface_patch_vec[0], surface_patch_vec.size()/3);
          cell_data(ofs, &surface_patch_type_vec[0], surface_patch_type_vec.size(), "type");
        }
      }
    }
#endif

    shared_ptr<function_t<double, int32_t> > func(
          build_polycube_patch_function(tm.node_, tm.mesh_,surface_patches, patch_normal,
                                        patch_boundary, fnode_group));

    shared_ptr<jtf::function::functionN1_t<double,int32_t> >  lsgn(jtf::function::least_square_warpper(func));
    jtf::optimize(*lsgn, node, pt, nullptr, nullptr, nullptr);

    { // visualize
      ostringstream os_i;
      os_i << i;
      string output = "output-refine-0";
      output +=  os_i.str() + ".vtk";
      ofstream ofs(output.c_str());
      tet2vtk(ofs, &node[0], node.size(2), &tm.mesh_[0], tm.mesh_.size(2));
    }

#ifdef VISUAL
    {
      ostringstream os;
      os << i;
      {
        string patch_str = "patch_" + os.str() + ".vtk";
        ofstream ofs(patch_str.c_str());
        vector<size_t> patch_faces_vis;
        vector<size_t> patch_type;
        //for(size_t pi = 0; pi < surface_patches.size(); ++pi){
        for(boost::unordered_map<size_t,size_t>::const_iterator cit = surface_type.begin();
            cit != surface_type.end(); ++cit){
          const vector<size_t> & one_face = fa->faces_[cit->first];
          patch_faces_vis.insert(patch_faces_vis.end(), one_face.begin(),
                                 one_face.end());
          patch_type.push_back(cit->second);
        }
        tri2vtk(ofs, &node[0], node.size(2), &patch_faces_vis[0], patch_faces_vis.size()/3);
        cell_data(ofs, &patch_type[0], patch_type.size(), "patch_idx");
      }
      {
        string feature = "feature_" + os.str() + ".vtk";
        ofstream ofs(feature.c_str());
        vector<size_t> edges;

        for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit =
            patch_boundary.begin(); cit != patch_boundary.end(); ++cit){
          edges.push_back(cit->first);
          edges.push_back(cit->second);

        }
        line2vtk(ofs, &node[0], node.size(2), &edges[0], edges.size()/2);
      }
    }
#endif
  }


  jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("output.value").c_str(), &node, &tm.mesh_);

  cerr << "success." << endl;
  return 0;
}
