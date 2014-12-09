#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <iostream>

#include "../hex_param/cut_tet.h"
#include "../common/vtk.h"
#include "../equation_graph/equation_graph.h"
#include "../hex_param/io.h"
using namespace std;
using namespace zjucad::matrix;

void need_merged(const zjucad::matrix::matrix<size_t> & two_faces,
                 vector<size_t> & idx_need_to_be_merged)
{
  idx_need_to_be_merged.clear();
  assert(two_faces.size(1) == 3 && two_faces.size(2) == 2);
  for(size_t pi = 0; pi < two_faces.size(1); ++pi){
      if(two_faces(pi, 0) != two_faces(pi, 1)){
          idx_need_to_be_merged.push_back(pi);
        }
    }
  if(idx_need_to_be_merged.size() == 3){ // all points need to be merged? CAN NOT apply this merging,
      // or topology will be changed.
      idx_need_to_be_merged.clear();
    }
}

void merge_faces(const size_t idx0,const size_t idx1,
                 vector<matrix<size_t> > & faces)
{
  for(size_t fi = 0; fi < faces.size(); ++fi){
      for(size_t pj = 0; pj < faces[fi].size(); ++pj){
          if(faces[fi][pj] == idx0)
            faces[fi][pj] = idx1;
        }
    }
}
int eliminate_gaps(
    const zjucad::matrix::matrix<size_t> &cut_tet,
    const zjucad::matrix::matrix<size_t> &uncut_tet,
    const zjucad::matrix::matrix<size_t> &cut_tet2tet,
    zjucad::matrix::matrix<double> & cut_node,
    const zjucad::matrix::matrix<double> & uncut_node,
    const boost::unordered_map<std::pair<size_t, size_t>, size_t> &inner_type,
    zjucad::matrix::matrix<size_t> &cut_tet_new)
{// this function will add all gap node, and eliminate them, after that,
  // all gaps with zero value and identity rotation can be remvoed

  boost::unordered_map<size_t,size_t> surface_type;
  bool is_restricted_type = false; // empty surface type do not need this.
  unique_ptr<transition_elimination> te(
        transition_elimination::create(
          uncut_tet, cut_tet,cut_tet2tet,uncut_node, inner_type, surface_type,
          is_restricted_type));
  if(!te.get()){
      cerr << "# [error] can not build transition elimination." << endl;
      return __LINE__;
    }

  vector<matrix<size_t> > trivial_cut_face_pairs;
  te->get_trivial_cut_face_pair(inner_type, cut_tet, cut_tet2tet, trivial_cut_face_pairs);

  {
    vector<size_t> trivial_faces;
    vector<size_t> type;
    for(size_t fi = 0; fi < trivial_cut_face_pairs.size(); ++fi){
        trivial_faces.insert(trivial_faces.end(), trivial_cut_face_pairs[fi](colon(),0).begin(),
            trivial_cut_face_pairs[fi](colon(),0).end());
        trivial_faces.insert(trivial_faces.end(), trivial_cut_face_pairs[fi](colon(),1).begin(),
            trivial_cut_face_pairs[fi](colon(),1).end());
        type.push_back(fi);
        type.push_back(fi);
      }
    ofstream ofs("trivial_cut_face.vtk");
    tri2vtk(ofs, &cut_node[0], cut_node.size(2), &trivial_faces[0], trivial_faces.size()/3);
    cell_data(ofs, &type[0], type.size(), "type");
  }

  cut_tet_new = cut_tet;

  // checking whether non-manifold or genus is changed after two point merging
  map<size_t,size_t> node_mapping;
  {
    vector<matrix<size_t> > trivial_cut_face_pairs_merged = trivial_cut_face_pairs;
    vector<bool> face_merged(trivial_cut_face_pairs.size(), false);
    vector<size_t> idx_need_to_be_merged;
    while(1){
        bool merged = false;
        for(size_t i = 0; i < face_merged.size(); ++i){
            if(face_merged[i]) continue;
            need_merged(trivial_cut_face_pairs_merged[i],idx_need_to_be_merged);
            if(!idx_need_to_be_merged.empty()){ // merge two points at idx position
                for(const auto & idx : idx_need_to_be_merged){
                    merge_faces(trivial_cut_face_pairs_merged[i](idx,0),
                        trivial_cut_face_pairs_merged[i](idx,1),
                        trivial_cut_face_pairs_merged);
                  }
                merged = true;
                face_merged[i] = true;
              }
          }
        if(merged == false) break;
      }
    cerr << "# [info] merged " << std::count(face_merged.begin(), face_merged.end(), true)
         << " faces" << endl;
    for(size_t fi = 0; fi < trivial_cut_face_pairs_merged.size(); ++fi){
        for(size_t pi = 0 ; pi < trivial_cut_face_pairs_merged[fi].size(); ++pi){
            if(trivial_cut_face_pairs_merged[fi][pi] !=
               trivial_cut_face_pairs[fi][pi]){
                node_mapping[trivial_cut_face_pairs[fi][pi]] = trivial_cut_face_pairs_merged[fi][pi];
              }
          }
      }
  }

  for(const auto & one_mapping : node_mapping){
      auto jt = find(cut_tet_new.begin(), cut_tet_new.end(), one_mapping.first);
      while(jt != cut_tet_new.end()){
          *jt = one_mapping.second;
          jt = find(cut_tet_new.begin(), cut_tet_new.end(), one_mapping.first);
        }
    }

  //  {
  //    map<size_t, set<size_t> > node_mapping;
  //    for(size_t fi = 0; fi < trivial_cut_face_pairs.size(); ++fi){
  //        for(size_t pi = 0; pi < trivial_cut_face_pairs[fi].size(1); ++pi){
  //            node_mapping[cut_tet2tet[trivial_cut_face_pairs[fi](pi,0)]].insert(
  //                  trivial_cut_face_pairs[fi](pi,0));
  //            node_mapping[cut_tet2tet[trivial_cut_face_pairs[fi](pi,1)]].insert(
  //                  trivial_cut_face_pairs[fi](pi,1));
  //          }
  //      }
  //    for(const auto & one_node_mapping : node_mapping){
  //        const set<size_t> & mapping = one_node_mapping.second;
  //        if(mapping.size() == 1) continue;
  //        const size_t target = *mapping.begin();
  //        auto it = mapping.begin();
  //        ++it;
  //        for(; it != mapping.end(); ++it){
  //            auto jt = find(cut_tet_new.begin(), cut_tet_new.end(), *it);
  //            while(jt != cut_tet_new.end()){
  //                *jt = target;
  //                jt = find(cut_tet_new.begin(), cut_tet_new.end(), *it);
  //              }
  //          }
  //      }
  //  }


  remove_extra_node(cut_tet_new, cut_node);

  return 0;
}


int simply_cut_tet(int argc, char * argv[])
{
  if(argc != 5){
      cerr << "# [usage] simply_cut_tet uncut_tet cut_tet inner_type output_tet" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes cut_mesh, uncut_mesh;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &uncut_mesh.node_, &uncut_mesh.mesh_))
    return __LINE__;

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &cut_mesh.node_, &cut_mesh.mesh_))
    return __LINE__;

  boost::unordered_map<pair<size_t,size_t>,size_t> inner_type;
  if(load_inner_face_jump_type(argv[3], inner_type))
    return __LINE__;

  matrix<size_t> cut_tet_new;
  matrix<size_t> cut_tet2tet(max(cut_mesh.mesh_)+1);
  cut_tet2tet(cut_mesh.mesh_) = uncut_mesh.mesh_(colon());
  eliminate_gaps(cut_mesh.mesh_, uncut_mesh.mesh_, cut_tet2tet, cut_mesh.node_,
                 uncut_mesh.node_, inner_type, cut_tet_new);

  ofstream ofs("output_simpl.vtk");
  tet2vtk(ofs, &cut_mesh.node_[0], cut_mesh.node_.size(2), &cut_tet_new[0], cut_tet_new.size(2));

  if(jtf::mesh::tet_mesh_write_to_zjumat(argv[4], &cut_mesh.node_, &cut_tet_new))
    return __LINE__;
  return 0;
}
