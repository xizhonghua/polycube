#include <iostream>
#include <set>
#include <fstream>

#include <jtflib/mesh/io.h>

#include "../common/vtk.h"
#include "../tetmesh/hex_io.h"
#include "../common/def.h"
#include "../hex_param/io.h"
using namespace std;

int check_type_setting_files(int argc, char * argv[])
{
  if(argc != 4 ){
    cerr << "# [usage] check_type_setting_files tet group_file chain_file "
         << endl;
    return __LINE__;
  }

  matrixst tet;
  matrixd node;



  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &node, &tet)){
    cerr << "# [error] can not open tet file." << endl;
    return __LINE__;
  }

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
  }

  matrixst outside_face;
  get_outside_face(*fa, outside_face);

  set<size_t> surface_points_set;
  for(size_t i = 0; i < outside_face.size(); ++i)
      surface_points_set.insert(outside_face[i]);

  set<size_t> surface_points_in_group;

  vector<vector<size_t> > groups,chains;
  if(load_group_file(argv[2], groups))
    return __LINE__;

  if(load_chain_file(argv[3], chains))
    return __LINE__;

  for(size_t gi = 0; gi < groups.size(); ++gi){
    const vector<size_t> & one_group = groups[gi];
    for(size_t i = 0; i < one_group.size(); ++i){
      const size_t point_idx = one_group[i]/3;
      const size_t uvw = one_group[i]%3;
      node(uvw, point_idx) = node(one_group.back()%3, one_group.back()/3);
      surface_points_in_group.insert(point_idx);
    }
  }

  if(surface_points_in_group.size() != surface_points_set.size()){
      cerr << "# [error] wrong group, since group node number do not equal surface point number." << endl;
      cerr << "# [error] If this tet is not polycube, jsut ignore it." << endl;
  }else{
      cerr << "# [info] group points equals surface points." << endl;
  }

  ofstream ofs_tet("grouped_tet.vtk");
  tet2vtk(ofs_tet, &node[0], node.size(2), &tet[0], tet.size(2));

  vector<size_t> chain2vtk;
  vector<size_t> chainidx;
  {
    for(size_t ci = 0; ci < chains.size(); ++ci){
      const vector<size_t> & one_chain = chains[ci];
      for(size_t i = 0; i < one_chain.size()-1; ++i){
        chain2vtk.push_back(one_chain[i]/3);
        chain2vtk.push_back(one_chain[i+1]/3);
        chainidx.push_back(ci);
      }
    }
    ofstream ofs_chain("grouped_chain.vtk");
    line2vtk(ofs_chain, &node[0], node.size(2),
             &chain2vtk[0],chain2vtk.size()/2);
    cell_data(ofs_chain, &chainidx[0], chainidx.size(), "idx");
  }

  {
    ofstream ofs_g_node("grouped_nodes.vtk");
    vector<size_t> group_nodes, group_node_idx;
    for(size_t gi = 0; gi < groups.size(); ++gi){
      const vector<size_t> & one_group = groups[gi];
      for(size_t i = 0; i < one_group.size(); ++i){
        group_nodes.push_back(one_group[i]/3);
        group_node_idx.push_back(gi);
      }
    }
    point2vtk(ofs_g_node, &node[0], node.size(2), &group_nodes[0], group_nodes.size());
    cell_data(ofs_g_node, &group_node_idx[0], group_node_idx.size(), "idx");
  }
  return 0;
}
