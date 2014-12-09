#include <fstream>
#include <sstream>

#include <jtflib/mesh/io.h>
#include <zjucad/matrix/matrix.h>

#include "../tetmesh/tetmesh.h"

using namespace std;
using namespace zjucad::matrix;

int dump_mesh2obj(const char * mesh_file, const char * tet_file);

int mesh2obj(int argc, char * argv[])
{
  if(argc < 3){
    cerr << "# [useage] mesh2obj input_mesh_file output_obj" << endl;
    return __LINE__;
  }

  dump_mesh2obj(argv[1], argv[2]);
  return 0;
}

int dump_mesh2obj(const char *mesh_file, const char *obj_file)
{
  assert(mesh_file && obj_file);

  ifstream mesh_ifs(mesh_file);
  if(mesh_ifs.fail()){
    cerr << "# [error] can not open mesh file." << endl;
    return __LINE__;
  }

  matrixd node;
  matrixst faces;

  string line;
  string first_word;
  double trash;
  while(!mesh_ifs.eof()){
    getline(mesh_ifs, line);
    stringstream words(line);
    words >> first_word;
    if(first_word == "Vertices"){ // read node
      size_t node_num = 0;
      mesh_ifs >> node_num;
      node.resize(3, node_num);

      for(size_t pi = 0; pi < node_num; ++pi){
        mesh_ifs >> node(0, pi) >> node(1,pi) >> node(2, pi) >> trash;
      }
      cerr << "# [info] finish read vertices : " << node_num << endl;
    }
    if(first_word == "Triangles"){
      size_t tri_num = 0;
      mesh_ifs >> tri_num;
      faces.resize(3, tri_num);
      for(size_t ti = 0; ti < tri_num; ++ti){
        mesh_ifs >> faces(0, ti) >> faces(1,ti) >> faces(2, ti) >> trash;
      }
      cerr << "# [info] finish read Triangles : " << tri_num << endl;
      break;
    }
    first_word.clear();
  }

  faces -=1;
  jtf::mesh::save_obj(obj_file, faces, node);
  return 0;
}
