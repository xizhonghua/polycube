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

int draw_zyz_normal_error(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [error] draw_zyz_normal_error tet zyz" << endl;
      return __LINE__;
    }

  matrix<double> node;
  matrix<size_t> tet;

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &node, &tet))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
      cerr << "# [error] can not create face2tet_adjacent."<< endl;
      return __LINE__;
    }

  matrix<size_t> outside_face;
  matrix<size_t> outside_face_idx;
  jtf::mesh::get_outside_face(*fa, outside_face, true);
  jtf::mesh::get_outside_face_idx(*fa, outside_face_idx);

  matrix<double> outside_face_normal;
  jtf::mesh::cal_face_normal(outside_face, node, outside_face_normal);

  matrix<double> normal_zyz_diff(outside_face_idx.size(),1);
  matrix<int> normal_zyz_idx(outside_face_idx.size(), 1);

  matrix<double> zyz;
  if(read_zyz(argv[2], zyz)){
      cerr << "# [error] can not read zyz file." << endl;
      return __LINE__;
    }

  matrix<matrix<double> > frame;
  frame.resize(zyz.size(2),1);
  for(size_t ti = 0; ti < zyz.size(2); ++ti){
      frame[ti].resize(3,3);
      zyz_angle_2_rotation_matrix1(&zyz(0,ti), &frame[ti][0]);
    }

  vector<pair<double,int> > diff(6);
  for(size_t fi = 0; fi < outside_face_idx.size(); ++fi){
      const size_t & face_idx = outside_face_idx[fi];
      const pair<size_t,size_t> & tet_pair = fa->face2tet_[face_idx];
      assert(tet_pair.first == -1 || tet_pair.second == -1);
      const size_t tet_idx = (tet_pair.first == -1?tet_pair.second:tet_pair.first);
      for(size_t di = 0; di < 3; ++di){
          diff[di * 2 + 0] = make_pair(
                dot(frame[tet_idx](colon(),di), outside_face_normal(colon(),fi)), di * 2 + 0);
          diff[di * 2 + 1] = make_pair(
                dot(frame[tet_idx](colon(),di)*-1.0, outside_face_normal(colon(),fi)), di * 2 + 1);
        }
      sort(diff.begin(), diff.end());
      normal_zyz_diff[fi] = diff.back().first;
      normal_zyz_idx[fi] = diff.back().second;
    }

  ofstream ofs("normal_zyz_diff.vtk");
  tri2vtk(ofs, &node[0], node.size(2), &outside_face[0], outside_face.size(2));
  cell_data(ofs, &normal_zyz_diff[0], normal_zyz_diff.size(), "error");
  vtk_data(ofs, &normal_zyz_idx[0], normal_zyz_idx.size(), "idx", "new_table");

  return 0;
}
