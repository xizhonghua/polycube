#include <zjucad/ptree/ptree.h>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/io.h>

#include "../hex_param/cut_tet.h"
#include "../equation_graph/cycle_detection.h"
#include "../common/vtk.h"
#include "../common/IO.h"
#include "../hex_param/io.h"
#include "../equation_graph/equation_graph.h"
#include "../equation_graph/util.h"
#include "../hex_param/topology_analysis.h"

using namespace std;
using boost::property_tree::ptree;
using namespace zjucad::matrix;



int convert_configuration_to_group(
    const matrix<size_t> & tetmesh,
    const matrix<double> & node,
    matrix<size_t> &cut_tet,
    const boost::unordered_map<pair<size_t,size_t>, size_t> & inner_type,
    const boost::unordered_map<size_t,size_t> & surface_type,
    const bool is_restricted_type,
    const jtf::mesh::face2tet_adjacent &fa,
    vector<vector<size_t> > & node_group_vec,
    vector<bool> & integer_group_flag,
    bool is_cut_tet,
    matrix<matrix<double> > * frame_ptr = 0);


int convert_configuration_to_group_wrapper(
    const matrix<size_t> & tetmesh,
    const matrix<double> & node,
    const char * inner_type_file,
    const char * surface_type_file,
    vector<vector<size_t> > & node_group_vec,
    vector<bool> & integer_group_flas,
    matrix<size_t> &cut_tet,
    bool is_cut_tet,
    const jtf::mesh::face2tet_adjacent *fa_ptr = 0,
    matrix<matrix<double> > * frame_ptr = 0)
{
  boost::unordered_map<pair<size_t,size_t>,size_t> inner_type;
  boost::unordered_map<size_t,size_t> surface_type;
  load_inner_face_jump_type(inner_type_file, inner_type);

  unique_ptr<jtf::mesh::face2tet_adjacent> fa_p;
  if(fa_ptr == 0){
      fa_p.reset(jtf::mesh::face2tet_adjacent::create(tetmesh));
      if(!fa_p.get()){
          cerr << "# [error] can not build face2tet_adjacent." << endl;
          return __LINE__;
        }
    }

  int restricted_type_or_not = load_surface_type(surface_type_file, surface_type,
                                                 (fa_ptr?fa_ptr:fa_p.get()));

  {
    if(restricted_type_or_not != 0)
      dump_surface_restricted_type_to_vtk("surface_type_before_global_aligned.vtk", "normal_type", node, *fa_p, surface_type);
  }
  convert_configuration_to_group(tetmesh, node, cut_tet, inner_type, surface_type,
                                 (restricted_type_or_not==0?true:false),
                                 (fa_ptr?*fa_ptr:*fa_p), node_group_vec,
                                 integer_group_flas, is_cut_tet, frame_ptr);

  {// save node_groups
    ofstream ofs("node_group");
    assert(node_group_vec.size() == integer_group_flas.size());
    for(size_t i = 0 ; i < node_group_vec.size(); ++i){
        ofs << "g " << i << " " << node_group_vec[i].size() << " " << integer_group_flas[i] << endl;
        for(size_t vi = 0; vi < node_group_vec[i].size(); ++vi){
            ofs << node_group_vec[i][vi] << " ";
          }
        ofs << endl;
      }
  }

  return 0;
}

static int load_data(ptree &pt, jtf::mesh::meshes &tm,
                     matrix<size_t> & cut_tet,
                     vector<vector<size_t> > & node_groups,
                     vector<bool> & integer_group_flag,
                     const jtf::mesh::face2tet_adjacent * fa = 0)
{
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm.node_, &tm.mesh_))
    return __LINE__;

  if(zjucad::has("group.value", pt)){
      if(load_group_file(pt.get<string>("group.value").c_str(), node_groups, &integer_group_flag))
        return __LINE__;
      cut_tet = tm.mesh_;
    }else{
      unique_ptr<matrix<matrixd> > frame_ptr = nullptr;
      pt.put("zyz.desc", "tet zyz frame file");

      if(zjucad::has("zyz.value", pt)){
          matrixd zyz;
          if(read_zyz(pt.get<string>("zyz.value").c_str(), zyz))
            return __LINE__;
          frame_ptr.reset(new matrix<matrixd>);
          frame_ptr->resize(zyz.size(2),1);
          zyz2frame(zyz, *frame_ptr);
        }

      bool is_cut_tet = false;

      pt.put("input/cut_tet.desc", "input cut_tet");

      if(zjucad::has("input/cut_tet.value",pt)){
          matrix<double> cut_node;
          if(jtf::mesh::tet_mesh_read_from_zjumat(
               pt.get<string>("input/cut_tet.value").c_str(), &cut_node, &cut_tet))
            return __LINE__;
          is_cut_tet = true;
        }
      convert_configuration_to_group_wrapper(
            tm.mesh_,tm.node_, pt.get<string>("inner_type.value").c_str(),
            pt.get<string>("surface_type.value","").c_str(), node_groups,
            integer_group_flag,cut_tet, is_cut_tet, fa, frame_ptr.get());

      if(frame_ptr != nullptr){
          matrixd zyz(3, frame_ptr->size());
          frame2zyz(*frame_ptr, zyz);
          write_zyz("aligned.zyz", zyz);
        }
    }
  return 0;
}

static void save_surface_type(
    const char * filename,
    const boost::unordered_map<size_t,size_t> & surface_type)
{
  ofstream ofs(filename);
  if(ofs.fail())
    throw std::logic_error("# [error] can not open surface restricted type.");

  for(boost::unordered_map<size_t,size_t>::const_iterator
      cit = surface_type.begin();
      cit != surface_type.end(); ++cit){
      ofs << cit->first << " " << cit->second << endl;
    }
  return;
}

int get_signed_surface_type(ptree &pt)
{
  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),
                                          &tm.node_, &tm.mesh_))
    return __LINE__;

  boost::unordered_map<size_t,size_t> surface_signed_type;
  load_surface_type(pt.get<string>("surface_type.value").c_str(),
                    surface_signed_type);

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));

  if(!fa.get()){
      throw std::logic_error("# [error] can not buildjtf::mesh::face2tet_adjacent.");
    }

  matrix<double> face_normal(3,1);
  const matrix<double> eye_mat = zjucad::matrix::eye<double>(3);

  for(auto& one_face: surface_signed_type){
      assert(one_face.first < fa->faces_.size());
      const vector<size_t> & one_face_vec = fa->faces_[one_face.first];
      itr_matrix<const size_t*> face_mat(3,1, &one_face_vec[0]);
      jtf::mesh::cal_face_normal(face_mat, tm.node_, face_normal, true);
      if(dot(face_normal, eye_mat(colon(), one_face.second)) > 0)
        one_face.second = one_face.second * 2;
      else
        one_face.second = one_face.second * 2 + 1;
    }

  save_surface_type(pt.get<string>("surface_signed_type.value").c_str(),
                    surface_signed_type);

  matrix<size_t> surface_faces(3, surface_signed_type.size());
  vector<size_t> surface_type;
  surface_type.reserve(surface_signed_type.size());

  size_t fi = 0;
  for(boost::unordered_map<size_t,size_t>::const_iterator
      cit = surface_signed_type.begin();
      cit != surface_signed_type.end(); ++cit){
      surface_type.push_back(cit->second);
      std::copy(fa->faces_[cit->first].begin(), fa->faces_[cit->first].end(),
          surface_faces(colon(),fi++).begin());
    }
  ofstream ofs("sigend_surface_type.vtk");
  tri2vtk(ofs, &tm.node_[0], tm.node_.size(2), &surface_faces[0], surface_faces.size(2));
  cell_data(ofs, &surface_type[0], surface_type.size(), "signed_type");
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

/// @brief it converts a tetmesh with inner/outside transition to group information
/// @param tetmesh
/// @param inner_type tet pair and its type
/// @param surface_type input surface type
/// @return
int convert_configuration_to_group(
    const matrix<size_t> & tetmesh,
    const matrix<double> & node,
    matrix<size_t> &cut_mesh,
    const boost::unordered_map<pair<size_t,size_t>, size_t> & inner_type,
    const boost::unordered_map<size_t,size_t> & surface_type,
    const bool is_restricted_type,
    const jtf::mesh::face2tet_adjacent &fa,
    vector<vector<size_t> > & node_group_vec,
    vector<bool> &integer_group_flag ,
    bool is_cut_tet,
    matrix<matrix<double> > * frame_ptr)
{
  // this function takes several steps:
  // 1. cut mesh
  // 2. extract node_group

  // note that inner_type_after_cut/surface_type_after_cut are defined on original mesh
  boost::unordered_map<pair<size_t,size_t>, size_t> inner_type_after_cut = inner_type;
  boost::unordered_map<size_t,size_t> surface_type_after_cut = surface_type;

  if(inner_type.empty()){// it's a polycube mesh.
      cut_mesh = tetmesh;
    }else{ // with inner singularity
      if(!is_cut_tet)
        cut_tetmesh(tetmesh, node, fa, inner_type_after_cut,
                    surface_type_after_cut, false,cut_mesh, false, frame_ptr);
    }

  // build cut_tet to tet mapping.
  matrix<size_t> cut_tet2tet(max(cut_mesh)+1);
  cut_tet2tet(cut_mesh) = tetmesh(colon());

  unique_ptr<transition_elimination> tre(
        transition_elimination::create(
          tetmesh,cut_mesh, cut_tet2tet, node,inner_type_after_cut,
          surface_type_after_cut,is_restricted_type ));
  if(!tre.get()){
      cerr << "# [error] can not get transition_elimination." << endl;
      return __LINE__;
    }

  vector<group<size_t> > node_group = tre->out();
  vector<size_t> integer_group_idx = tre->get_integer_group_idx();

  cerr << "# [info] integer group number " << integer_group_idx.size() << endl;
  vector<size_t> one_group_vec;
  //for(const auto & one_group : node_group){
  for(size_t gi = 0; gi < node_group.size(); ++gi){
      const auto & one_group = node_group[gi];
      if(one_group.empty()) continue;
      one_group_vec.resize(one_group.size());
      copy(one_group.begin(), one_group.end(), one_group_vec.begin());
      node_group_vec.push_back(one_group_vec);

      if(find(integer_group_idx.begin(), integer_group_idx.end(), gi)
         != integer_group_idx.end())
        integer_group_flag.push_back(true);
      else
        integer_group_flag.push_back(false);
    }

  dump_inner_face_jump_type("inner_type_cut",inner_type_after_cut);
  dump_inner_face_jump_type2vtk("inner_type.vtk", node, tetmesh, inner_type_after_cut);
  if(is_restricted_type){
      dump_surface_restricted_type_to_vtk("surface_type.vtk", "restricted_type", node, fa, surface_type_after_cut);
      dump_surface_restricted_type("surface_type_cut",surface_type_after_cut);
    }
  else{
      dump_surface_normal_align_type2vtk("surface_type.vtk", node, fa, surface_type_after_cut);
      matrix<size_t> outside_face, outside_face_idx;
      jtf::mesh::get_outside_face(fa, outside_face);
      jtf::mesh::get_outside_face_idx(fa, outside_face_idx);
      dump_surface_normal_align_type_map("surface_type_cut", outside_face,
                                         outside_face_idx, surface_type_after_cut);
    }


  return 0;
}

static void get_edges_from_tet(
    const matrix<size_t> & tet, vector<pair<size_t,size_t> > & edges)
{
  set<pair<size_t,size_t> > edges_set;
  for(size_t ti = 0; ti < tet.size(2); ++ti){
      for(size_t pi = 0; pi < tet.size(1); ++pi){
          for(size_t pj = pi + 1; pj < tet.size(1); ++pj){
              pair<size_t,size_t> one_edge(tet(pi,ti), tet(pj,ti));
              if(one_edge.first > one_edge.second)
                swap(one_edge.first, one_edge.second);
              edges_set.insert(one_edge);
            }
        }
    }
  edges.resize(edges_set.size());
  copy(edges_set.begin(), edges_set.end(), edges.begin());
}

int equation_graph_check(ptree &pt)
{
  jtf::mesh::meshes tm;
  vector<vector<size_t> > node_groups;
  vector<bool> integer_group_flag;

  zjucad::matrix::matrix<size_t> cut_tet;
  // all equality constraints are formulated in node group.
  if(load_data(pt, tm, cut_tet, node_groups,integer_group_flag))
    return __LINE__;

  matrix<size_t> cut_tet2tet(max(cut_tet)+1);
  cut_tet2tet(cut_tet) = tm.mesh_(colon());
  matrix<double> cut_node(3, max(cut_tet)+1);
  if(zjucad::has("input/cut_tet.value",pt)){
      if(jtf::mesh::tet_mesh_read_from_zjumat(
           pt.get<string>("input/cut_tet.value").c_str(), &cut_node, &cut_tet))
        return __LINE__;

    }else{
      for(size_t pi = 0; pi < cut_node.size(2); ++pi){
          cut_node(colon(), pi) = tm.node_(colon(), cut_tet2tet[pi]);
        }
    }

  {// debug
    if(zjucad::has("output/cut_tet.value",pt)){
        if(jtf::mesh::tet_mesh_write_to_zjumat(
             pt.get<string>("output/cut_tet.value").c_str(),
             &cut_node, &cut_tet))
          return __LINE__;
      }
  }

  equation_graph eg(cut_node.size());

  eg.bind_node(&cut_node);

  vector<pair<size_t,size_t> > edges;
  get_edges_from_tet(cut_tet, edges);
  eg.assemble_equation_graph_generally(node_groups, edges);

  bool rtn = eg.check_valid(true);
  if(!rtn){
      cerr << "# [error] wrong graph." << endl;
      return 0;
    }

  cerr << "# [info] well graph." << endl;
  return 0;
}
