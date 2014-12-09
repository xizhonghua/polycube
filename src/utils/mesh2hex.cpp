#include <fstream>
#include <sstream>
#include <iostream>
#include <zjucad/matrix/matrix.h>

#include "../hexmesh/io.h"

using namespace std;
using namespace zjucad::matrix;

int dump_mesh2hex(const char * mesh_file, const char * hex_file);

int mesh2hex(int argc, char * argv[])
{
  if(argc < 2){
    cerr << "# [useage] mesh2hex input_mesh_file output_hex_mesh" << endl;
    return __LINE__;
  }

  dump_mesh2hex(argv[1], argv[2]);
  return 0;
}

int dump_mesh2hex(const char *mesh_file, const char *hex_file)
{
  assert(mesh_file && hex_file);

  ifstream mesh_ifs(mesh_file);
  if(mesh_ifs.fail()){
    cerr << "# [error] can not open mesh file." << endl;
    return __LINE__;
  }

  matrixd node;
  matrixst hex;

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
    if(first_word == "Hexahedra"){
      size_t hex_num = 0;
      words >> hex_num;
      hex.resize(8, hex_num);
      for(size_t hi = 0; hi < hex_num; ++hi){
        //        for(size_t vi = 0; vi < 8; ++vi)
        //          mesh_ifs >> hex(vi, hi);
        mesh_ifs >> hex(7,hi) >> hex(5,hi)
                 >> hex(4,hi) >> hex(6,hi)
                 >> hex(3,hi) >> hex(1,hi)
                 >> hex(0,hi) >> hex(2,hi);
        mesh_ifs >> trash;
      }
      hex -= 1;
      cerr << "# [info] finish read hex : " << hex_num << endl;
      break;
    }
    first_word.clear();
  }

  jtf::hexmesh::hex_mesh_write_to_wyz(hex_file, hex, node);
  //orient_tet(node, tet);
  //hex_mesh_write_to_zjumat(hex_file, &node, &tet);
  return 0;
}
