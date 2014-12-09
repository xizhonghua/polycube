#include <iostream>
#include <vector>
#include <deque>
#include "../common/visualize_tool.h"
using namespace std;

int dump_sc_to_vtk(int argc, char * argv[])
{
  if(argc != 4){
    std::cerr << "# [usage] sc2vtk tet_file sc_file_name vtk_file_name" << std::endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  vector<deque<pair<size_t,size_t> > > singularity_edges;
  vector<deque<size_t> > singularity_type;
  if(load_singularity_chain_new(argv[2], singularity_edges, singularity_type))
    return __LINE__;

  if(dump_singularity_chain_to_vtk_2(argv[3], tm.node_,
                                     singularity_edges, singularity_type))
    return __LINE__;
  return 0;
}
