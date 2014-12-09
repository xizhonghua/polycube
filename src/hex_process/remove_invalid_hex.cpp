#include "hex_process.h"
#include "../hexmesh/util.h"

using namespace std;
using namespace jtf::hexmesh;
using namespace zjucad::matrix;

int remove_invalid_hex(matrixst & hex,
                       const matrixd & node)
{
  return -1; // not finished
  matrixst faces;
  vector<size_t> one_face(4,0);
  typedef boost::unordered_map<std::vector<size_t>, std::vector<size_t> >
      face2adj_hex_type;
  face2adj_hex_type face_adj_hex;

  for(size_t t = 0; t < hex.size(2); ++t){
    jtf::mesh::get_faces_for_one_hex(hex(colon(),t),faces);
    for(size_t fi = 0; fi < faces.size(2); ++fi){
      copy(faces(colon(),fi).begin(), faces(colon(),fi).end(),
           one_face.begin());
      sort(one_face.begin(), one_face.end());
      face_adj_hex[one_face].push_back(t);
    }
  }

  vector<bool> hex_valid(hex.size(2), true);

  for(face2adj_hex_type::const_iterator it = face_adj_hex.begin();
      it != face_adj_hex.end(); ++it){
    if(it->second.size() > 2){

    }
  }


  return 0;
}
