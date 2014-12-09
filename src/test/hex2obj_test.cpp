#include "hex2obj_test.h"
#include "../hexmesh/io.h"
#include "../hexmesh/hexmesh.h"
#include "../hexmesh/util.h"

#include <iostream>
#include <set>
#include <vector>
using namespace std;

void hex2obj_test::setUp()
{

}

void hex2obj_test::tearDown()
{

}

void hex2obj_test::hex2obj()
{

  set<set<size_t> > face_set;
  vector<vector<size_t> > face_original;
  set<size_t> one_face;
  vector<size_t> one_face_original;
    // matrixst one_face_original;
  matrixst hex;
  matrixd node;
  matrixst face_in_cube; 
  jtf::hexmesh::hex_mesh_read_from_wyz("hex_for_wyz.hex", hex, node,1);
  size_t face_set_size;
  for(size_t i = 0; i < hex.size(2); ++i)
    {
      jtf::hexmesh::get_faces_for_one_hex(hex(zjucad::matrix::colon(), i), face_in_cube);
      for(size_t j = 0 ; j < face_in_cube.size(2); ++j)
	{
	  one_face.clear();
	  one_face_original.clear();
	  for(size_t k = 0; k < face_in_cube.size(1); ++k)
	    {
              one_face.insert(face_in_cube(k, j));
	      one_face_original.push_back(face_in_cube(k ,j));
	    }
	  face_set_size = face_set.size();
	  if(one_face.size() == 4)
	    {
              face_set.insert(one_face);
	      if(face_set_size != face_set.size())
		face_original.push_back(one_face_original);
	    }
	}
    }
  cout << node.size(2) << "  " << hex.size(2) <<endl;
  cout << face_set.size() << endl;
  std::ofstream ofs("hex2obj.obj");
  for(size_t i = 0; i < node.size(2); ++i)
    ofs << "v " << node(0, i) << " " << node(1, i) << " " << node(2, i) << endl;
  // for(face_set_it = face_set.begin(); face_set_it != face_set.end(); ++face_set_it)
  //   {
  //     ofs << "f";
  //     for(face_it = face_set_it -> begin(); face_it != face_set_it -> end(); ++ face_it)
  // 	ofs << " " << *face_it;
  //     ofs << endl;
  //   }
  for(size_t i = 0 ; i < face_original.size(); ++i)
    ofs << "f " <<  face_original[i][0]+1  << " "
                <<  face_original[i][1]+1  << " "
                <<  face_original[i][2]+1  << " "
	        <<  face_original[i][3]+1  << endl; 
}
