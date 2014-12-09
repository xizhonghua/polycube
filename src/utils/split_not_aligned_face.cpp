#include <fstream>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>

#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"
#include "../tet_mesh_sxx/tet_mesh_sxx.h"
#include "../common/vtk.h"

using namespace std;
using namespace zjucad::matrix;

void split_face( const jtf::mesh::meshes & polycube_tm,
                 sxx::tet_mesh & stm,
                 matrixst & outside_face,
                 matrixd & face_normal)
{
  boost::unordered_set<size_t> face_need_to_split;
  for(size_t fi = 0; fi < outside_face.size(2); ++fi){
    double sum_fabs = 0;

    for(size_t pi = 0; pi < face_normal.size(1); ++pi)
      sum_fabs += fabs(face_normal(pi,fi));

    if(fabs(sum_fabs-1) > 1e-2)
      face_need_to_split.insert(fi);
  }

  { // visual
    vector<size_t> split_faces;
    for(boost::unordered_set<size_t>::const_iterator cit =
        face_need_to_split.begin(); cit != face_need_to_split.end(); ++cit){
      split_faces.insert(split_faces.end(),
                         outside_face(zjucad::matrix::colon(), *cit).begin(),
                         outside_face(zjucad::matrix::colon(), *cit).end());
    }
    ofstream ofs("split_face.vtk");
    tri2vtk(ofs, &polycube_tm.node_[0], polycube_tm.node_.size(2),
            &split_faces[0], split_faces.size()/3);
  }

  for(boost::unordered_set<size_t>::const_iterator cit = face_need_to_split.begin();
      cit != face_need_to_split.end(); ++cit){
    const size_t & face_idx_in_outside_face = *cit;
    stm.split_face(&outside_face(0,face_idx_in_outside_face));
  }
}

void split_edges(const jtf::mesh::meshes & polycube_tm,
                 sxx::tet_mesh & stm,
                 matrixst & outside_face,
                 matrixd & face_normal)
{
  boost::unordered_set<pair<size_t,size_t> > edges_need_to_split;
  boost::unordered_set<size_t> face_need_to_split;
  for(size_t fi = 0; fi < outside_face.size(2); ++fi){
    double sum_fabs = 0;

    for(size_t pi = 0; pi < face_normal.size(1); ++pi)
      sum_fabs += fabs(face_normal(pi,fi));

    if(fabs(sum_fabs-1) > 1e-2)
      face_need_to_split.insert(fi);
  }
  for(boost::unordered_set<size_t>::const_iterator cit =
      face_need_to_split.begin(); cit != face_need_to_split.end(); ++cit){
    const size_t & face_idx_in_outside = *cit;
    for(size_t pi = 0; pi < outside_face.size(1); ++pi){
      pair<size_t,size_t> edge(
            outside_face(pi, face_idx_in_outside),
            outside_face((pi+1)%outside_face.size(1), face_idx_in_outside));
      if(edge.first > edge.second)
        swap(edge.first, edge.second);
      edges_need_to_split.insert(edge);
    }
  }

  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator eit =
      edges_need_to_split.begin(); eit != edges_need_to_split.end(); ++eit){
    stm.split_edge(*eit);
  }
}

void collapse_degenerated_edges(
    sxx::tet_mesh &stm,
    const matrixst &outside_face,
    const matrixd &polycube_node)
{
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return ;
  }

  const double threshold = 1e-2;
  double average_len = 0;
  for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
    average_len += norm(polycube_node(colon(), ea->edges_[ei].first)
                   - polycube_node(colon(), ea->edges_[ei].second));
  }
  average_len /= ea->edges_.size();

  boost::unordered_set<pair<size_t,size_t> > edges_need_collapse;

  for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
    const double len = norm(polycube_node(colon(), ea->edges_[ei].first)
                            - polycube_node(colon(), ea->edges_[ei].second));
    if(len / average_len < threshold){
      pair<size_t,size_t> one_edge = ea->edges_[ei];
      if(one_edge.first > one_edge.second)
        swap(one_edge.first, one_edge.second);
      edges_need_collapse.insert(one_edge);
    }
  }

  cerr << "# [info] edges_need_collapse_num: " << edges_need_collapse.size() << endl;
  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit =
      edges_need_collapse.begin(); cit != edges_need_collapse.end(); ++cit){
    cerr << "# [info] ------ " << cit->first << " "  << cit->second << endl;
    int rtn = stm.collapse_edge(*cit);
    if(rtn == -1)
      cerr << "# [info] ------ collapse fail." << endl;
  }
}

int split_not_aligned_face(int argc, char * argv[])
{
  if(argc != 4){
    cerr << "# [usage] split_not_aligned_face tet polycube_tet output_tet." << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm, polycube_tm, output_tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_)){
    return __LINE__;
  }

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &polycube_tm.node_, &polycube_tm.mesh_)){
    return __LINE__;
  }

  if(fabs(norm(tm.mesh_ - polycube_tm.mesh_)) > 1e-6){
    cerr << "# [error] original tet is not the same as polycube_tet in topology." << endl;
    return __LINE__;
  }

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  matrixst outside_face;
  matrixd face_normal;
  get_outside_face(*fa, outside_face);

  jtf::mesh::cal_face_normal(outside_face, polycube_tm.node_, face_normal);

  sxx::tet_mesh stm;
  stm.create_tetmesh(tm.node_, tm.mesh_);

  //split_face(stm, outside_face, face_normal);
  //split_edges(stm, outside_face, face_normal);
  collapse_degenerated_edges(stm, outside_face, polycube_tm.node_);
  //collapse_degenerated_faces(stm, outside_face, )
  stm.write_tetmesh_to_file(argv[3]);

//  stm.write_tetmesh_to_matrix(tm.node_, tm.mesh_);
//  matrixd node_ = tm.node;
//  jtf::tetmesh::tetmesh_quality_improver tqi(tm.mesh_, node_);
//  tqi.improve();

//  cerr << "# [info] difference " <<  norm(node_ - tm.node_) << endl;


  return 0;
}
