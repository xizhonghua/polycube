#include <set>
#include "../common/IO.h"
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>


using namespace std;
using namespace jtf::mesh;
using namespace zjucad::matrix;

int compress_obj(int argc, char  * argv[])
{
  if(argc != 4){
    cerr << "# [usage] compress_obj input_obj output_obj s2v_mat." << endl;
    return __LINE__;
  }

  meshes trm;
  if(load_obj(argv[1], trm.mesh_, trm.node_))
    return __LINE__;

  set<size_t> used_vertex(trm.mesh_.begin(), trm.mesh_.end());
  matrix<int32_t> used_vertex_idx(used_vertex.size(),1);
  std::copy(used_vertex.begin(), used_vertex.end(), used_vertex_idx.begin());

  map<size_t,size_t> p2p;
  for(size_t pi = 0; pi < used_vertex_idx.size(); ++pi){
    p2p[used_vertex_idx[pi]] = pi;
  }

  // clean the node
  matrix<double> new_node(3, used_vertex_idx.size());
  for(size_t pi = 0; pi < new_node.size(2); ++pi){
    new_node(colon(),pi) = trm.node_(colon(), used_vertex_idx[pi]);
  }
  trm.node_ = new_node;

  // clean the tri
  for(size_t pi = 0; pi < trm.mesh_.size(); ++pi){
    trm.mesh_[pi] = p2p[trm.mesh_[pi]];
  }

  save_obj(argv[2], trm.mesh_, trm.node_);
  jtf::mesh::write_matrix(argv[3], used_vertex_idx);

  // dump out the s2v_mat;

  return 0;
}
