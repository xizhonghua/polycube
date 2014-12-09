#include <zjucad/ptree/ptree.h>
#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include <string>
#include "../equation_graph/equation_graph.h"
#include "../common/vtk.h"
#include "../hex_param/io.h"

using namespace std;
using namespace jtf::mesh;
using namespace zjucad::matrix;
using boost::property_tree::ptree;


void construct_original_face_in_cut_mesh(
    const matrix<size_t> & uncut_tet,
    const matrix<size_t> & cut_tet,
    const jtf::mesh::face2tet_adjacent & uncut_fa,
    const jtf::mesh::face2tet_adjacent & cut_fa,
    const matrix<size_t> &cut_tet2tet,
    matrix<size_t> &orig_face_in_cut)
{
  vector<size_t> orig_face_idx_in_cut; // it store original surface in cut mesh
  for(size_t fi = 0 ; fi < cut_fa.faces_.size(); ++fi){
      const vector<size_t> & one_face = cut_fa.faces_[fi];
      const size_t face_idx  = uncut_fa.get_face_idx(
            cut_tet2tet[one_face[0]],  cut_tet2tet[one_face[1]],
          cut_tet2tet[one_face[2]]);
      if(face_idx == -1) {
          cerr << "# [error] can not find face "
               << one_face[0] << " " << one_face[1] << " " << one_face[2] << endl;
          return ;
        }
      if(uncut_fa.is_outside_face(uncut_fa.face2tet_[face_idx]))
        orig_face_idx_in_cut.push_back(fi);
    }

  orig_face_in_cut.resize(3, orig_face_idx_in_cut.size());
  for(size_t fi = 0; fi < orig_face_idx_in_cut.size(); ++fi){
      std::copy(cut_fa.faces_[orig_face_idx_in_cut[fi]].begin(),
          cut_fa.faces_[orig_face_idx_in_cut[fi]].end(),
          &orig_face_in_cut(0,fi));
    }
}

void construct_original_face_in_cut_mesh(
    const boost::property_tree::ptree &pt,
    const zjucad::matrix::matrix<size_t> & cut_tet,
    const jtf::mesh::face2tet_adjacent & cut_fa,
    zjucad::matrix::matrix<size_t> &orig_face_in_cut)
{
  matrix<size_t> uncut_tet ;
  matrix<double> uncut_node;
  if(jtf::mesh::tet_mesh_read_from_zjumat(
       pt.get<string>("input/uncut_tet.value").c_str(), &uncut_node,
       &uncut_tet))
    throw std::invalid_argument("can not open uncut_tet.");
  unique_ptr<jtf::mesh::face2tet_adjacent> uncut_fa(
        jtf::mesh::face2tet_adjacent::create(uncut_tet));
  if(!uncut_fa.get()){
      throw std::invalid_argument("can not build face2tet_adjacent.");
    }

  matrix<size_t> cut_tet2tet(max(cut_tet) + 1);
  cut_tet2tet(cut_tet) = uncut_tet(colon());

  construct_original_face_in_cut_mesh(
        uncut_tet, cut_tet, *uncut_fa, cut_fa, cut_tet2tet, orig_face_in_cut);
}


static void pt_description(ptree &pt)
{
  pt.put("input/tet.desc", "input cut tet");
  pt.put("input/node_group.desc", "node_group contains only groups come frome inner transition, need to be updated");
}


static void merge_node(const size_t &from, const size_t &to,
                       vector<size_t> & node2group,
                       vector<group<size_t> > & node_group)
{
  if(node2group[from] == node2group[to]) return;
  group<size_t> & from_group = node_group[node2group[from]];
  for(const auto & idx : from_group) node2group[idx] = node2group[to];

  node_group[node2group[to]].merge(from_group);
}

int label_surface_type_after_param(ptree &pt)
{
  cerr << "# [WARNING!!!] This function should be called after parameterization." << endl;
  pt_description(pt);

  meshes cut_tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(
       pt.get<string>("input/tet.value").c_str(),
       &cut_tm.node_, &cut_tm.mesh_))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> cut_fa(jtf::mesh::face2tet_adjacent::create(cut_tm.mesh_));
  if(!cut_fa.get()){
      cerr << "# [error] can not build face2tet_adjacent." << endl;
      return __LINE__;
    }

  matrix<size_t> uncut_tet ;
  matrix<double> uncut_node;
  if(jtf::mesh::tet_mesh_read_from_zjumat(
       pt.get<string>("input/uncut_tet.value").c_str(), &uncut_node,
       &uncut_tet))
    throw std::invalid_argument("can not open uncut_tet.");
  unique_ptr<jtf::mesh::face2tet_adjacent> uncut_fa(
        jtf::mesh::face2tet_adjacent::create(uncut_tet));
  if(!uncut_fa.get()){
      throw std::invalid_argument("can not build face2tet_adjacent.");
    }

  matrix<size_t> cut_tet2tet(max(cut_tm.mesh_) + 1);
  cut_tet2tet(cut_tm.mesh_) = uncut_tet(colon());

  matrix<size_t> orig_face_in_cut;
  construct_original_face_in_cut_mesh(
        uncut_tet, cut_tm.mesh_, *uncut_fa, *cut_fa, cut_tet2tet, orig_face_in_cut);

  matrix<double> normal(3,1);
  const matrix<double> axes = zjucad::matrix::eye<double>(3);
  matrix<size_t> orig_face_in_cut_type(orig_face_in_cut.size(2),1);

  vector<pair<double,int> > normal_diff(6);
  for(size_t fi = 0; fi < orig_face_in_cut.size(2); ++fi){
      jtf::mesh::cal_face_normal(orig_face_in_cut(colon(),fi), cut_tm.node_, normal);
      for(size_t di = 0; di < 3; ++di){
          normal_diff[di*2+0] = make_pair(dot(normal, axes(colon(),di)), di*2+0);
          normal_diff[di*2+1] = make_pair(dot(normal, -1*axes(colon(),di)), di*2+1);
        }
      sort(normal_diff.begin(), normal_diff.end());
      orig_face_in_cut_type[fi] = normal_diff.back().second/2;
    }

  matrix<size_t> orig_face_in_orig = orig_face_in_cut;
  for(size_t i = 0 ; i < orig_face_in_orig.size(); ++i)
    orig_face_in_orig[i] = cut_tet2tet[orig_face_in_orig[i]];

  {// visualize
    ofstream ofs("surface_label_in_cut.vtk");
    tri2vtk(ofs, &cut_tm.node_[0], cut_tm.node_.size(2), &orig_face_in_cut[0], orig_face_in_cut.size(2));
    cell_data(ofs, &orig_face_in_cut_type[0], orig_face_in_cut_type.size(), "label");

    ofstream ofs_uncut("surface_label_in_orig.vtk");
    tri2vtk(ofs_uncut, &uncut_node[0], uncut_node.size(2), &orig_face_in_orig[0], orig_face_in_orig.size(2));
    cell_data(ofs_uncut, &orig_face_in_cut_type[0], orig_face_in_cut_type.size(), "label");
  }

  ofstream ofs_surface_type("surface_restricted_type_uncut");
  for(size_t fi = 0; fi < orig_face_in_orig.size(2); ++fi){
      const size_t face_idx = uncut_fa->get_face_idx(&orig_face_in_orig(0,fi));
      assert(face_idx != -1);
      ofs_surface_type << face_idx << " " << orig_face_in_cut_type[fi] << endl;
    }

  ofstream ofs_surface_type_cut("surface_restricted_type_cut");
  for(size_t fi = 0; fi < orig_face_in_cut.size(2); ++fi){
      const size_t face_idx = cut_fa->get_face_idx(&orig_face_in_cut(0,fi));
      assert(face_idx != -1);
      ofs_surface_type_cut << face_idx << " " << orig_face_in_cut_type[fi] << endl;
    }

  if(zjucad::has("input/node_group.value",pt)){
      vector<size_t> node2group(cut_tm.node_.size());
      vector<group<size_t> > node_groups(cut_tm.node_.size());

      // init node groups
      for(size_t vi = 0; vi < node2group.size(); ++vi){
          node2group[vi] = vi;
          node_groups[vi] << vi;
        }

      vector<vector<size_t> > groups;
      if(load_group_file(pt.get<string>("input/node_group.value").c_str(),
                         groups))
        return __LINE__;

      for(size_t gi = 0; gi < groups.size(); ++gi){
          const vector<size_t> & one_group = groups[gi];

          for(size_t vi = 1; vi < one_group.size(); ++vi){
              merge_node(one_group[vi],one_group.front(), node2group, node_groups);
            }
        }

//      for(size_t fi = 0; fi < orig_face_in_cut.size(2); ++fi){
//          const size_t axis_type = orig_face_in_cut_type[fi];
//          assert(axis_type == 0 || axis_type == 1 || axis_type == 2);
//          const size_t dim = orig_face_in_cut.size(1);
//          for(size_t pi = 1; pi < dim; ++pi){
//              merge_node(dim*orig_face_in_cut(pi,fi)+axis_type,
//                         dim*orig_face_in_cut(0,fi)+axis_type,
//                         node2group, node_groups);
//            }
//        }


      {// save node_groups
        ofstream ofs("update_node_group_no_integer_flag");
        size_t gi = 0;
        for(size_t i = 0 ; i < node_groups.size(); ++i){
            if(node_groups[i].empty()) continue;

            ofs << "g " << gi++ << " " << node_groups[i].size() << "  0" << endl;
            for(const auto & idx : node_groups[i])
              ofs << idx << " " ;
            ofs << endl;
          }
      }

    }

  return 0;
}
