#include <iostream>
#include "../tetmesh/hex_io.h"
#include "../tetmesh/tetmesh.h"
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>

#include "../common/vtk.h"
#include "../common/IO.h"

using namespace std;
using namespace zjucad::matrix;

int obj2faketet(int argc, char *argv[] )
{
  if(argc != 5 && argc != 3){
      cerr << "# [usage1] obj2faketet obj s2v orig_tet output_tet" << endl;
      cerr << "# [usage2] obj2faketet obj output_tet" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes trm;
  if(jtf::mesh::load_obj(argv[1], trm.mesh_, trm.node_))
    return __LINE__;

  if(argc == 5){
      matrix<int32_t> s2v;
      ifstream ifs(argv[2]);
      if(ifs.fail()){
          cerr << "# [error] can not open s2v file." << endl;
          return __LINE__;
        }

      if(jtf::mesh::read_matrix(ifs, s2v))
        return __LINE__;

      jtf::mesh::meshes tm;
      if(jtf::mesh::tet_mesh_read_from_zjumat(argv[3], &tm.node_, &tm.mesh_))
        return __LINE__;

      for(size_t si = 0; si < s2v.size(); ++si){
          tm.node_(colon(), s2v[si]) = trm.node_(colon(),si);
        }

      if(jtf::mesh::tet_mesh_write_to_zjumat(argv[4], &tm.node_, &tm.mesh_))
        return __LINE__;
      cerr << "# [info] success." << endl;
    }else if (argc == 3){
      matrix<double> avg_node = trm.node_ * ones<double>(trm.node_.size(2),1);
      avg_node /= trm.node_.size(2);
      jtf::mesh::meshes tm;
      tm.node_ = zeros<double>(3, trm.node_.size(2) + 1);
      tm.node_(colon(), colon(0, trm.node_.size(2)-1)) = trm.node_;
      tm.node_(colon(), trm.node_.size(2)) = avg_node;
      tm.mesh_ = ones<size_t>(4, trm.mesh_.size(2)) * trm.node_.size(2);
      tm.mesh_(colon(0,2),colon()) = trm.mesh_;
      if(jtf::mesh::tet_mesh_write_to_zjumat(argv[2], &tm.node_, &tm.mesh_))
        return __LINE__;
      cerr << "# [info] success." << endl;
    }

  return 0;
}
