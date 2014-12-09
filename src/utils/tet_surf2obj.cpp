#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

#include "../common/IO.h"
#include "../tetmesh/tetmesh.h"

using namespace zjucad::matrix;

int tet_surf2obj(int argc, char *argv[])
{
  if(argc < 2) {
    cerr << "tet_surf2obj tet [s2v]. when has s2v output packed surface model." << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm;
  matrixst tri;
  if(jtf::mesh::tet_mesh_read_from_zjumat(
       argv[1],
       &tm.node_, &tm.mesh_, &tri))
    return __LINE__;

  if(tri.size(1) != 3){ // didn't get right triangle information from tetmesh
    unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
    if(!fa.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }
    get_outside_face(*fa, tri, true);
  }

  if(argc == 2) {
    for(size_t i = 0; i < tm.node_.size(2); ++i)
      cout << "v " <<  tm.node_(0, i)
           << " " <<  tm.node_(1, i)
           << " " <<  tm.node_(2, i) << "\n";
    tri += 1;
    for(size_t i = 0; i < tri.size(2); ++i)
      cout << "f " <<  tri(0, i)
           << " " <<  tri(1, i)
           << " " <<  tri(2, i) << "\n";
    return 0;
  }

  // pack
  matrix<int32_t> tri2 = tri;
  sort(tri2.begin(), tri2.end());

  matrix<int32_t>::const_iterator end = unique(tri2.begin(), tri2.end());


  matrix<int32_t> pack2unpack = tri2(colon(0, end-tri2.begin()-1));
  if(pack2unpack[pack2unpack.size()-1] == pack2unpack.size()-1)
    cerr << "# leading nodes are surface nodes." << endl;

  matrix<int32_t> unpack2pack = zeros<double>(max(tri)+1, 1)-1;
  unpack2pack(pack2unpack) = trans(colon(0, pack2unpack.size()-1));
  tri2(colon()) = unpack2pack(tri)+1;


  ofstream ofs(argv[2], ofstream::binary);
  jtf::mesh::write_matrix(ofs, pack2unpack);

  for(size_t i = 0; i < pack2unpack.size(); ++i) {
    cout << "v " <<  tm.node_(0, pack2unpack[i])
         << " " <<  tm.node_(1, pack2unpack[i])
         << " " <<  tm.node_(2, pack2unpack[i]) << "\n";
  }

  for(size_t i = 0; i < tri2.size(2); ++i) {
    cout << "f " << tri2(0, i)
         << " " << tri2(1, i)
         << " " << tri2(2, i) << "\n";
  }

  return 0;
}

int vtk_tet_surf2obj(int argc, char *argv[])
{
  if(argc < 2) {
    cerr << "vtk_tet_surf2obj tet [s2v]. when has s2v output packed surface model." << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm;
  matrixst tri;
  if(jtf::mesh::tet_mesh_read_from_vtk(
       argv[1],
       &tm.node_, &tm.mesh_))
    return __LINE__;

  if(tri.size(1) != 3){ // didn't get right triangle information from tetmesh
    unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
    if(!fa.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }
    get_outside_face(*fa, tri, true);
  }

  if(argc == 2) {
    for(size_t i = 0; i < tm.node_.size(2); ++i)
      cout << "v " <<  tm.node_(0, i)
           << " " <<  tm.node_(1, i)
           << " " <<  tm.node_(2, i) << "\n";
    tri += 1;
    for(size_t i = 0; i < tri.size(2); ++i)
      cout << "f " <<  tri(0, i)
           << " " <<  tri(1, i)
           << " " <<  tri(2, i) << "\n";
    return 0;
  }

  // pack
  matrix<int32_t> tri2 = tri;
  sort(tri2.begin(), tri2.end());

  matrix<int32_t>::const_iterator end = unique(tri2.begin(), tri2.end());


  matrix<int32_t> pack2unpack = tri2(colon(0, end-tri2.begin()-1));
  if(pack2unpack[pack2unpack.size()-1] == pack2unpack.size()-1)
    cerr << "# leading nodes are surface nodes." << endl;

  matrix<int32_t> unpack2pack = zeros<double>(max(tri)+1, 1)-1;
  unpack2pack(pack2unpack) = trans(colon(0, pack2unpack.size()-1));
  tri2(colon()) = unpack2pack(tri)+1;


  ofstream ofs(argv[2], ofstream::binary);
  jtf::mesh::write_matrix(ofs, pack2unpack);

  for(size_t i = 0; i < pack2unpack.size(); ++i) {
    cout << "v " <<  tm.node_(0, pack2unpack[i])
         << " " <<  tm.node_(1, pack2unpack[i])
         << " " <<  tm.node_(2, pack2unpack[i]) << "\n";
  }

  for(size_t i = 0; i < tri2.size(2); ++i) {
    cout << "f " << tri2(0, i)
         << " " << tri2(1, i)
         << " " << tri2(2, i) << "\n";
  }

  return 0;
}
