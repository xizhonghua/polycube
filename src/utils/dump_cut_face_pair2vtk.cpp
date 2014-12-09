#include <fstream>
#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include "../common/IO.h"
#include "../common/zyz.h"
#include "../common/vtk.h"
using namespace std;
using namespace zjucad::matrix;

int dump_cut_face_pair2vtk(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [error] dump_cut_face_pair2vtk orig_tet cut_tet" << endl;
      return __LINE__;
    }

  matrix<double> uncut_node;
  matrix<size_t> uncut_tet;

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &uncut_node, &uncut_tet))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(uncut_tet));
  if(!fa.get()){
      cerr << "# [error] can not create face2tet_adjacent."<< endl;
      return __LINE__;
    }

  matrix<size_t> outside_face;
  matrix<size_t> outside_face_idx;
  jtf::mesh::get_outside_face(*fa, outside_face, true);
  jtf::mesh::get_outside_face_idx(*fa, outside_face_idx);


  matrix<double> cut_node;
  matrix<size_t> cut_tet;

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &cut_node, &cut_tet))
    return __LINE__;


  unique_ptr<jtf::mesh::face2tet_adjacent> fa_cut(jtf::mesh::face2tet_adjacent::create(cut_tet));
  if(!fa_cut.get()){
      cerr << "# [error] can not create face2tet_adjacent."<< endl;
      return __LINE__;
    }

  matrix<size_t> outside_face_cut, outside_face_cut_idx;
  jtf::mesh::get_outside_face(*fa_cut, outside_face_cut);
  jtf::mesh::get_outside_face_idx(*fa_cut, outside_face_cut_idx);

  matrix<size_t> cut_tet2tet(max(cut_tet) + 1);
  cut_tet2tet(cut_tet) = uncut_tet(colon());

  vector<size_t> outside_face_cut_in_orig_idx;
  vector<size_t> outside_face_cut_in_orig_type;


  map<vector<size_t>, vector<size_t> > orig_face2cut_face;
  vector<size_t> one_cut_face(3);
  for(size_t fi = 0; fi < outside_face_cut.size(2); ++fi){
      for(size_t pi = 0; pi < 3; ++pi){
          one_cut_face[pi] = cut_tet2tet[outside_face_cut(pi,fi)];
        }
      sort(one_cut_face.begin(), one_cut_face.end());
      orig_face2cut_face[one_cut_face].push_back(fi);
    }

  vector<size_t> cut_faces;
  vector<size_t> cut_face_type;
  size_t type = 0;
  for(const auto & one_face: orig_face2cut_face){
      if(one_face.second.size() == 1) continue;
      assert(one_face.second.size() == 2);
      for(size_t i = 0; i < one_face.second.size(); ++i){
          cut_faces.insert(cut_faces.end(),
                           outside_face_cut(colon(),one_face.second[i]).begin(),
                           outside_face_cut(colon(),one_face.second[i]).end());
          cut_face_type.push_back(type);
        }
      ++type;
    }

  ofstream ofs("cut_face_pair.vtk");
  tri2vtk(ofs, &cut_node[0], cut_node.size(2), &cut_faces[0], cut_faces.size()/3);
  cell_data(ofs, &cut_face_type[0], cut_face_type.size(), "cut_face_pair");
  return 0;
}
