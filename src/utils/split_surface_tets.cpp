#include <boost/unordered_map.hpp>
#include <numeric>
#include <jtflib/mesh/io.h>

#include "../tetmesh/hex_io.h"
#include "../tetmesh/tetmesh.h"
#include "../tet_mesh_sxx/tet_mesh_sxx.h"
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"


using namespace std;
using namespace zjucad::matrix;

//! @brief assume face1,face2 are adjacent tri, and this prog find the non-connected
//!  edge
//! @param face1 input face1
//! @param face2 input face2
//! @param edge output edge
//! @return 0 if run correct, otherwise non-zeros
int find_non_connect_edge(const vector<size_t> & face1,
                          const vector<size_t> & face2,
                          pair<size_t,size_t> & edge)
{
  set<size_t> collect_points;
  collect_points.insert(face1.begin(), face1.end());
  collect_points.insert(face2.begin(), face2.end());
  if(collect_points.size() != 4) {
      cerr << "# [error] input face1 face2 are not adjacent triangles." << endl;
      return __LINE__;
    }
  const size_t sum = std::accumulate(collect_points.begin(), collect_points.end(),0);
  const size_t sum_face1 = std::accumulate(face1.begin(), face1.end(),0);
  const size_t sum_face2 = std::accumulate(face2.begin(), face2.end(),0);
  edge.first = sum - sum_face1;
  edge.second = sum - sum_face2;
  return 0;
}

//! @brief check whether a tet have multi surface triangle,
//! if so split opposite edge to one pair of such adjacent faces
//! @param node input/output node
//! @param tet  input/output tet
int split_multi_surface_tets(matrix<double> & node,
                             matrix<size_t> & tet,
                             const size_t check_or_split)
{
  assert(check_or_split == 0 || check_or_split == 1);
  matrix<size_t> outside_face_idx;

  std::shared_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }

  set<pair<size_t,size_t> > edges_need_split;
  get_outside_face_idx(*fa, outside_face_idx);

  map<size_t, vector<size_t> > tet_with_face_idx;
  for(size_t fi = 0; fi < outside_face_idx.size(); ++fi){
      const size_t &face_idx = outside_face_idx[fi];
      const pair<size_t,size_t> & tet_pair = fa->face2tet_[face_idx];
      assert(tet_pair.first == -1 || tet_pair.second == -1);
      const size_t tet_idx = (tet_pair.first == -1?tet_pair.second:tet_pair.first);
      tet_with_face_idx[tet_idx].push_back(face_idx);
    }
  bool find_split_edge = false;
  pair<size_t,size_t> one_split_edge;
  for(map<size_t, vector<size_t> >::const_iterator cit = tet_with_face_idx.begin();
      cit != tet_with_face_idx.end(); ++cit){
      if(cit->second.size() < 2) continue;
      assert(cit->second.size() == 2 || cit->second.size() == 3);
      find_split_edge = true;
      const vector<size_t> & faces = cit->second;
      for(size_t i = 0; i < faces.size()-1; ++i){
          find_non_connect_edge(fa->faces_[faces[i]], fa->faces_[faces[i+1]],
              one_split_edge);
          edges_need_split.insert(one_split_edge);
        }
    }

  if(check_or_split == 0 ){ // only detection
      if(find_split_edge)
        cerr << "# [error] find tet with multi-face on surface." << endl;
      else
        cerr << "# [info] all tets have only one face on surface." << endl;
      return 0;
    }
  sxx::tet_mesh stm;
  stm.create_tetmesh(node, tet);

  for(set<pair<size_t,size_t> >::const_iterator cit = edges_need_split.begin();
      cit != edges_need_split.end(); ++cit)
    stm.split_edge(*cit);

  stm.write_tetmesh_to_matrix(node, tet);

  return 0;
}

//! @brief split inner edge whose both ends touch surface
//! @param node input/output node
//! @param tet input/output tet
//! @param flag to control only check or real splitting
//! @return 0: run ok, non-zeros: meet problem
int split_inner_edge_touch_surface(matrix<double> & node,
                                   matrix<size_t> & tet,
                                   const size_t check_or_split)
{
  assert(check_or_split == 0 || check_or_split == 1);
  std::shared_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }

  matrixst outside_face_idx, outside_face;
  get_outside_face_idx(*fa, outside_face_idx);
  get_outside_face(*fa, outside_face);

  boost::unordered_set<size_t> surface_node(outside_face.begin(), outside_face.end());

  jtf::mesh::one_ring_tet_at_edge ortae;
  ortae.add_tets(tet, *fa);
  ortae.sort_into_loop(tet, node);

  boost::unordered_set<pair<size_t,size_t> > edges_need_to_split;
  bool meet_case = false;
  for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator cit =
      ortae.e2t_.begin(); cit != ortae.e2t_.end(); ++cit){
      const pair<size_t,size_t> & one_edge = cit->first;
      const vector<size_t> & loop = cit->second;
      if(ortae.is_inner_edge(loop)){
          if(surface_node.find(one_edge.first) != surface_node.end() &&
             surface_node.find(one_edge.second) != surface_node.end()){
              if(check_or_split == 0) {
                  cerr << "# [error] meet inner edge has surface points." << endl;
                  return 0;
                }
              meet_case = true;
              if(one_edge.first > one_edge.second)
                edges_need_to_split.insert(make_pair(one_edge.second, one_edge.first));
              else
                edges_need_to_split.insert(one_edge);
            }
        }
    }

  if(check_or_split == 0){ // only detection
      cerr << "# [info] all inner edge is ok." << endl;
      return 0;
    }

  sxx::tet_mesh stm;
  stm.create_tetmesh(node, tet);

  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit =
      edges_need_to_split.begin(); cit != edges_need_to_split.end(); ++cit){
      stm.split_edge(*cit);
    }

  stm.write_tetmesh_to_matrix(node, tet);
  return 0;
}

int split_surface_tets(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] split_surface_tets tet output_tet 0/1[only_check/split]." << endl;
      return __LINE__;
    }

  int check_or_split = atoi(argv[3]); // 0: only check, 1: split

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  split_multi_surface_tets(tm.node_, tm.mesh_,check_or_split);
  split_inner_edge_touch_surface(tm.node_, tm.mesh_,check_or_split);


  if(jtf::mesh::tet_mesh_write_to_zjumat(argv[2], &tm.node_, &tm.mesh_)){
      cerr << "# [error] fail to write tet." << endl;
      return __LINE__;
    }
  return 0;
}
