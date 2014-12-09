#include <set>
#include <iostream>
#include <sstream>

#include <zjucad/matrix/matrix.h>
#include "../hexmesh/io.h"
#include "../hexmesh/util.h"

using namespace std;
using namespace zjucad::matrix;

int convert_hex_to_obj(const matrixst & hex,
                       const matrixd & node);

int hex2obj(int argc, char * argv[])
{
  if(argc != 3){
    std::cerr << "# [usage] hex2obj hex hex_formate." << endl;
    return 0;
  }
  stringstream os(argv[2]);
  size_t hex_formate = -1;
  os >> hex_formate;

  matrixst hex;
  matrixd node;
  if(jtf::hexmesh::hex_mesh_read_from_wyz(argv[1],hex, node, hex_formate,true))
    return __LINE__;

  convert_hex_to_obj(hex, node);

  return 0;
}


int convert_hex_to_obj(const matrixst & hex,
                       const matrixd & node)
{
  set<set<size_t> > face_set;
  vector<vector<size_t> > face_original;
  set<size_t> one_face;
  vector<size_t> one_face_original;

  matrixst face_in_cube;

  size_t face_set_size;
  for(size_t i = 0; i < hex.size(2); ++i) {
    jtf::mesh::get_faces_for_one_hex(
          hex(zjucad::matrix::colon(), i), face_in_cube);
    for(size_t j = 0 ; j < face_in_cube.size(2); ++j) {
      one_face.clear();
      one_face_original.clear();
      for(size_t k = 0; k < face_in_cube.size(1); ++k) {
        one_face.insert(face_in_cube(k, j));
        one_face_original.push_back(face_in_cube(k ,j));
      }
      face_set_size = face_set.size();
      if(one_face.size() == 4){
        face_set.insert(one_face);
        if(face_set_size != face_set.size())
          face_original.push_back(one_face_original);
      }
    }
  }

  for(size_t i = 0; i < node.size(2); ++i)
    cout << "v " << node(0, i) << " " << node(1, i) << " " << node(2, i) << endl;

  for(size_t i = 0 ; i < face_original.size(); ++i)
    cout << "f " <<  face_original[i][0]+1  << " "
         <<  face_original[i][1]+1  << " "
          <<  face_original[i][2]+1  << " "
           <<  face_original[i][3]+1  << endl;
  return 0;
}
