#include <fstream>
#include <sstream>

#include <jtflib/mesh/io.h>
#include <zjucad/matrix/matrix.h>

#include "../tetmesh/tetmesh.h"

using namespace std;
using namespace zjucad::matrix;

int dump_mesh2tet(const char * mesh_file, const char * tet_file);

int mesh2tet(int argc, char * argv[])
{
  if(argc < 2){
    cerr << "# [useage] mesh2tet input_mesh_file output_tet_mesh" << endl;
    return __LINE__;
  }

  dump_mesh2tet(argv[1], argv[2]);
  return 0;
}

int dump_mesh2tet(const char *mesh_file, const char *tet_file)
{
  assert(mesh_file && tet_file);

  ifstream mesh_ifs(mesh_file);
  if(mesh_ifs.fail()){
    cerr << "# [error] can not open mesh file." << endl;
    return __LINE__;
  }

  matrixd node;
  matrixst tet;

  string line;
  string first_word;
  double trash;
  while(!mesh_ifs.eof()){
    getline(mesh_ifs, line);
    stringstream words(line);
    words >> first_word;
    if(first_word == "Vertices"){ // read node
      size_t node_num = 0;
      words >> node_num;
      node.resize(3, node_num);

      for(size_t pi = 0; pi < node_num; ++pi){
        mesh_ifs >> node(0, pi) >> node(1,pi) >> node(2, pi) >> trash;
      }
      cerr << "# [info] finish read vertices : " << node_num << endl;
    }
    if(first_word == "Tetrahedra"){
      size_t tet_num = 0;
      words >> tet_num;
      tet.resize(4, tet_num);
      for(size_t ti = 0; ti < tet_num; ++ti){
        mesh_ifs >> tet(0, ti) >> tet(1,ti) >> tet(2, ti) >> tet(3,ti) >> trash;
      }
      cerr << "# [info] finish read tets : " << tet_num << endl;
      break;
    }
    first_word.clear();
  }

  orient_tet(node, tet);
  jtf::mesh::tet_mesh_write_to_zjumat(tet_file, &node, &tet);
  return 0;
}
