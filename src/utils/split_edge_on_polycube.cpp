#include <iostream>
#include <fstream>

#include <jtflib/mesh/io.h>

#include "../tet_mesh_sxx/tet_mesh_sxx.h"
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
using namespace std;

int load_split_edge_file(const char * file,
                         vector<pair<size_t,size_t> > & edges_need_split)
{
  ifstream ifs(file);
  if(ifs.fail()){
    cerr << "# [error] can not open split edge file." << endl;
    return __LINE__;
  }

  size_t edge_num;

  ifs >> edge_num;
  edges_need_split.resize(edge_num);

  for(size_t ei = 0; ei < edge_num; ++ei){
    ifs >> edges_need_split[ei].first >> edges_need_split[ei].second;
  }

  return 0;
}

int  split_edge_on_polycube(int argc, char * argv[])
{
  if(argc != 4){
    cerr << "# [usage] split_edge_on_polycube polycube_tet "
         << "orig_split_tet split_edge_file." << endl;
    return __LINE__;
  }

  jtf::mesh::meshes polycube_tm, orig_tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &polycube_tm.node_, &polycube_tm.mesh_))
    return __LINE__;

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &orig_tm.node_, &orig_tm.mesh_))
    return __LINE__;


  vector<pair<size_t,size_t> > edges_need_to_split;
  if(load_split_edge_file(argv[3], edges_need_to_split))
    return __LINE__;

  if(edges_need_to_split.size() != 0 &&
     (orig_tm.mesh_.size()  == polycube_tm.mesh_.size())){
    cerr << "# [error] polycube tet size is the same as orig splitted tet." << endl;
    return __LINE__;
  }

  cerr << "# [info] " << edges_need_to_split.size() << " edges need to split." << endl;
  sxx::tet_mesh stm;
  stm.create_tetmesh(polycube_tm.node_, polycube_tm.mesh_);
  for(size_t ei = 0; ei < edges_need_to_split.size(); ++ei){
    stm.split_edge(edges_need_to_split[ei]);
  }

  stm.write_tetmesh_to_matrix(polycube_tm.node_, polycube_tm.mesh_);

  if(polycube_tm.mesh_.size() != orig_tm.mesh_.size()){
    cerr << "# [error] wrong split." << endl;
    return __LINE__;
  }

  orient_tet(orig_tm.node_, polycube_tm.mesh_);

  jtf::mesh::tet_mesh_write_to_zjumat("polycube_after_split.tet", &polycube_tm.node_, &polycube_tm.mesh_);
  cerr << "# [info] success." << endl;
  return 0;
}
