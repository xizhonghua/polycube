#include "../tetmesh/tetmesh.h"
#include "../hex_param/io.h"
#include "../hex_param/cut_tet.h"
#include "../hex_param/topology_analysis.h"
#include <fstream>
#include <boost/unordered_map.hpp>

using namespace std;

int check_group_file(int argc, char * argv[])
{
  return __LINE__;
//  if(argc != 6){
//    cerr << "# [usage] check_group_file tet cut_tet group_file inner_face_jump_type"
//            " g_unknown_face_pair." << endl;
//    return __LINE__;
//  }

//  jtf::mesh::meshes tm, cut_tm;
//  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_)){
//    return __LINE__;
//  }

//  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &cut_tm.node_, &cut_tm.mesh_)){
//    return __LINE__;
//  }

//  vector<size_t> fnode;
//  if(load_fnode_group(argv[3],cut_tm.mesh_, fnode)){
//    return __LINE__;
//  }

//  if(fnode.size() != cut_tm.node_.size()){
//    cerr << "# [error] strange cut_tm node size is " << cut_tm.node_.size()
//         << " fnode size is " << fnode.size() << endl;
//    return __LINE__;
//  }

//  zjucad::matrix::itr_matrix<const size_t*> fnode_mat(3, fnode.size() /3,&fnode[0]);

//  boost::unordered_map<size_t,size_t> inner_face_jump_type;
//  if(load_inner_face_jump_type(argv[4], inner_face_jump_type)){
//    return __LINE__;
//  }

//  vector<std::pair<size_t,size_t> > g_unknown_face_pair;
//  if(load_g_unknown_face(argv[5], g_unknown_face_pair)){
//    return __LINE__;
//  }

//  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
//  unique_ptr<jtf::mesh::face2tet_adjacent> fa_cut(jtf::mesh::face2tet_adjacent::create(cut_tm.mesh_));

//  if(!fa.get() || !fa_cut.get()){
//    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
//    return __LINE__;
//  }

////  boost::unordered_set<size_t> g_unknown_face
////  matrixst face_pair, outside_face_idx;
////  get_outside_face_idx(*fa_cut, outside_face_idx);

//  analysis_transition_raw(tm.mesh_,tm.node_,*fa,*fa_cut,cut_tm.mesh_,face_pair);

  return 0;
}
