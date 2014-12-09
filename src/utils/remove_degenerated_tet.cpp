#include <iostream>
#include <fstream>
#include <numeric>

#include <jtflib/mesh/util.h>
#include <jtflib/mesh/io.h>

#include "../common/util.h"
#include "../common/vtk.h"
#include "../tet_mesh_sxx/tet_mesh_sxx.h"
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"
#include "../numeric/util.h"

using namespace std;
using namespace zjucad::matrix;

int check_vol_zero_tet(
    const jtf::mesh::meshes & tm,
    const jtf::mesh::meshes & polycube_tm)
{
  if(fabs(norm(tm.mesh_ - polycube_tm.mesh_)) > 1e-6){
    cerr << "# [error] orig_tet is not_compatiable with polycube_tet "
            << fabs(norm(tm.mesh_ - polycube_tm.mesh_)) << endl;
    return __LINE__;
  }

  const double threshold = 1e-5;
  vector<size_t> degenerated_tet;
  for(size_t ti = 0; ti < polycube_tm.mesh_.size(2); ++ti){
    const double vol =
        jtf::mesh::cal_tet_vol(polycube_tm.node_(colon(), polycube_tm.mesh_(colon(),ti)));
    if(fabs(vol) < threshold)
      degenerated_tet.push_back(ti);
  }

  cerr << "# [info] no degenerated tet (" << threshold << ")." << endl;

  if(degenerated_tet.empty()){
    return 0;
  }

  {// visualize
    ofstream ofs("degenerated_tet_in_orig_tet.vtk");
    ofstream ofs_polycube("degenerated_tet_in_orig_tet.vtk");
    tet2vtk(ofs, &tm.node_[0], tm.node_.size(2), &tm.mesh_[0], tm.mesh_.size(2));
    tet2vtk(ofs_polycube, &polycube_tm.node_[0], polycube_tm.node_.size(2),
            &polycube_tm.mesh_[0], polycube_tm.mesh_.size(2));
  }
  return 0;
}

int remove_face_zero_tet(
    jtf::mesh::meshes & tm,
    jtf::mesh::meshes & polycube_tm)
{
  if(fabs(norm(tm.mesh_ - polycube_tm.mesh_)) > 1e-6){
    cerr << "# [error] orig_tet is not_compatiable with polycube_tet "
            << fabs(norm(tm.mesh_ - polycube_tm.mesh_)) << endl;
    return __LINE__;
  }

  sxx::tet_mesh stm_orig, stm_polycube;
  stm_orig.create_tetmesh(tm.node_, tm.mesh_);
  stm_polycube.create_tetmesh(polycube_tm.node_, polycube_tm.mesh_);

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return  __LINE__;
  }

  matrixd face_area = zeros<double>(fa->faces_.size(),1);

  for(size_t fi = 0; fi < fa->faces_.size(); ++fi){
    const vector<size_t> & one_face = fa->faces_[fi];
    face_area[fi] = jtf::mesh::cal_face_area(&one_face[0], one_face.size(), polycube_tm.node_);
  }
  double average_area =
      std::accumulate(face_area.begin(), face_area.end(),0.0)/face_area.size();

  const double threshold = 1e-2;
  const double angle_threshold = cos(60.0/180.0 * My_PI());
  matrixd angle = zeros<double>(3,1); // (p1-p0)*(p2-p0),(p2-p1)*(p0-p1),(p1-p2)*(p0-p2)

  // each face is stored as a,b,c. ab is the edge which should be split with a new vertex d
  // and then remove dc
  boost::unordered_set<vector<size_t> > face_need_to_address;

  for(size_t fi = 0; fi < face_area.size(); ++fi){
    const vector<size_t> & one_face = fa->faces_[fi];
    double len_p0p1 =
        norm(polycube_tm.node_(colon(), one_face[1]) - polycube_tm.node_(colon(), one_face[0]));
    double len_p0p2 =
        norm(polycube_tm.node_(colon(), one_face[2]) - polycube_tm.node_(colon(), one_face[0]));
    double len_p1p2 =
        norm(polycube_tm.node_(colon(), one_face[1]) - polycube_tm.node_(colon(), one_face[2]));

    if(len_p0p1 < 1e-6) len_p0p1 = 1;
    if(len_p0p2 < 1e-6) len_p0p2 = 1;
    if(len_p1p2 < 1e-6) len_p1p2 = 1;

    if(face_area[fi] / average_area < threshold){
      angle[0] =
          dot(polycube_tm.node_(colon(), one_face[1]) - polycube_tm.node_(colon(), one_face[0]),
              polycube_tm.node_(colon(), one_face[2]) - polycube_tm.node_(colon(), one_face[0]))
          / (len_p0p1 * len_p0p2);
      angle[1] =
          dot(polycube_tm.node_(colon(), one_face[2]) - polycube_tm.node_(colon(), one_face[1]),
              polycube_tm.node_(colon(), one_face[0]) - polycube_tm.node_(colon(), one_face[1]))
          / (len_p1p2 * len_p0p1);
      angle[2] =
          dot(polycube_tm.node_(colon(), one_face[1]) - polycube_tm.node_(colon(), one_face[2]),
              polycube_tm.node_(colon(), one_face[0]) - polycube_tm.node_(colon(), one_face[2]))
          / (len_p1p2 * len_p0p2);

      vector<size_t>  obtuse_angle;
      for(size_t mi = 0; mi < angle.size(); ++mi)
        if(angle[mi] < angle_threshold) obtuse_angle.push_back(mi);
      if(obtuse_angle.size() > 1) continue; //

      // move the obtuse angle point to the last
      vector<size_t> face_point_order = one_face;
      swap(face_point_order[obtuse_angle.front()], face_point_order.back());
      face_need_to_address.insert(face_point_order);
    }
  }

  cerr << "# [info] there are " << face_need_to_address.size()
       << " faces need to address." << endl;
  if(face_need_to_address.empty()) return 1;

  { // visualliation
    vector<size_t> faces;
    faces.reserve((face_need_to_address.size() * 3));
    for(boost::unordered_set<vector<size_t> >::const_iterator cit =
        face_need_to_address.begin(); cit != face_need_to_address.end(); ++cit){
      const vector<size_t> & one_face = *cit;
      faces.insert(faces.end(), one_face.begin(), one_face.end());
    }
    ofstream ofs_polycube("face_need_to_address_on_polycube.vtk");
    ofstream ofs_orig("face_need_to_address_on_orig.vtk");
    tri2vtk(ofs_orig, &tm.node_[0], tm.node_.size(2), &faces[0], faces.size()/3);
    tri2vtk(ofs_polycube, &polycube_tm.node_[0], polycube_tm.node_.size(2), &faces[0], faces.size()/3);
  }


  for(boost::unordered_set<vector<size_t> >::const_iterator fit =
      face_need_to_address.begin(); fit != face_need_to_address.end(); ++fit){
    const vector<size_t> & one_face = *fit;
    const size_t new_point_idx =
        stm_orig.split_edge(make_pair(one_face[0], one_face[1]));
    if(new_point_idx == -1) continue;
    const size_t new_point_idx_2 =
        stm_polycube.split_edge(make_pair(one_face[0], one_face[1]));
    if(new_point_idx_2 != new_point_idx) {
      cerr << "# [error] the new point in polycube is not the same as orig_tet."
           << endl;
      return __LINE__;
    }
    int rtn = stm_orig.collapse_edge(make_pair(new_point_idx, one_face.back()),true);
    if(rtn == 1) return 1;

    int rtn2 = stm_polycube.collapse_edge(make_pair(new_point_idx_2, one_face.back()), true);
    if(rtn2 == 1) return 1;
  }

  stm_orig.write_tetmesh_to_matrix(tm.node_, tm.mesh_);
  stm_polycube.write_tetmesh_to_matrix(polycube_tm.node_, polycube_tm.mesh_);

  assert(polycube_tm.mesh_.size() == tm.mesh_.size());
  polycube_tm.mesh_ = tm.mesh_;

  return 0;
}

int remove_edge_zero_tet(
    jtf::mesh::meshes & tm,
    jtf::mesh::meshes & polycube_tm)
{
  if(fabs(norm(tm.mesh_ - polycube_tm.mesh_)) > 1e-6){
    cerr << "# [error] orig_tet is not_compatiable with polycube_tet "
            << fabs(norm(tm.mesh_ - polycube_tm.mesh_)) << endl;
    return __LINE__;
  }

  sxx::tet_mesh stm_orig, stm_polycube;
  stm_orig.create_tetmesh(tm.node_, tm.mesh_);
  stm_polycube.create_tetmesh(polycube_tm.node_, polycube_tm.mesh_);

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  jtf::mesh::one_ring_tet_at_edge ortae;
  ortae.add_tets(tm.mesh_, *fa);

  double average_edge_len = 0;

  for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator cit =
      ortae.e2t_.begin(); cit != ortae.e2t_.end(); ++cit){
    const pair<size_t,size_t> & one_edge = cit->first;
    average_edge_len += norm(polycube_tm.node_(colon(), one_edge.first) -
                             polycube_tm.node_(colon(), one_edge.second));
  }
  average_edge_len /= ortae.e2t_.size();

  const double threshold = 1e-2;
  double edge_len = 0;

  boost::unordered_set<pair<size_t,size_t> > edges_need_to_collapse;

  for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator cit =
      ortae.e2t_.begin(); cit != ortae.e2t_.end(); ++cit){
    const pair<size_t,size_t> & one_edge = cit->first;

    edge_len = norm(polycube_tm.node_(colon(), one_edge.first) -
                    polycube_tm.node_(colon(), one_edge.second));
    if(edge_len / average_edge_len < threshold)
      edges_need_to_collapse.insert(cit->first);
  }

  cerr << "# [info] there are " << edges_need_to_collapse.size()
       << " edges need to be collapsed." << endl;

  if(edges_need_to_collapse.empty())
    return 1;

  {
    ofstream ofs("edges_need_to_collapse.vtk");
    vector<size_t> edges;
    for(boost::unordered_set<pair<size_t,size_t> >::const_iterator eit =
        edges_need_to_collapse.begin(); eit != edges_need_to_collapse.end(); ++eit){
      edges.push_back(eit->first);
      edges.push_back(eit->second);
    }
    line2vtk(ofs, &polycube_tm.node_[0], polycube_tm.node_.size(2),
             &edges[0], edges.size()/2);
  }

  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator eit =
      edges_need_to_collapse.begin(); eit != edges_need_to_collapse.end(); ++eit){
    stm_orig.collapse_edge(*eit);
    stm_polycube.collapse_edge(*eit);
  }

  stm_orig.write_tetmesh_to_matrix(tm.node_, tm.mesh_);
  stm_polycube.write_tetmesh_to_matrix(polycube_tm.node_, polycube_tm.mesh_);

  assert(polycube_tm.mesh_.size() == tm.mesh_.size());
  polycube_tm.mesh_ = tm.mesh_;

  return 0;
}

int remove_degenerated_tet(int argc , char * argv[])
{
  if(argc != 5){
    cerr << "# [usage] remove_degenerated_tet orig_tet polycube_tet "
         << "output_orig_tet output_polycube_tet."  << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm, polycube_tm;

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &polycube_tm.node_, &polycube_tm.mesh_))
    return __LINE__;

  if(tm.mesh_.size() != polycube_tm.mesh_.size()){
    cerr << "# [error] orig_tet has different topology with polycube_tet." << endl;
    return __LINE__;
  }

  for(size_t try_i = 0; try_i < 10; ++try_i){
    int rtn1 = remove_edge_zero_tet(tm, polycube_tm);
    int rtn2 = remove_face_zero_tet(tm, polycube_tm);
    if(rtn1 == 1 && rtn2 == 1) break;
  }

  check_vol_zero_tet(tm, polycube_tm);

  if(jtf::mesh::tet_mesh_write_to_zjumat(argv[3], &tm.node_, &tm.mesh_))
    return __LINE__;

  if(jtf::mesh::tet_mesh_write_to_zjumat(argv[4], &polycube_tm.node_, &polycube_tm.mesh_))
    return __LINE__;

  return 0;
}
