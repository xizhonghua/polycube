
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

int map_normal_type_to_cut_tet(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [error] map_normal_type_to_cut_tet orig_tet  zyz cut_tet" << endl;
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

  matrix<double> outside_face_normal;
  jtf::mesh::cal_face_normal(outside_face, uncut_node, outside_face_normal);

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

  map<size_t,size_t> face_idx2type;
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
      normal_zyz_idx[fi] = diff.back().second;
      face_idx2type[outside_face_idx[fi]] = diff.back().second;
    }

  //  ofstream ofs("normal_zyz_diff.vtk");
  //  tri2vtk(ofs, &node[0], node.size(2), &outside_face[0], outside_face.size(2));
  //  cell_data(ofs, &normal_zyz_diff[0], normal_zyz_diff.size(), "error");
  //  vtk_data(ofs, &normal_zyz_idx[0], normal_zyz_idx.size(), "idx", "new_table");

  {
    matrix<double> cut_node;
    matrix<size_t> cut_tet;

    if(jtf::mesh::tet_mesh_read_from_zjumat(argv[3], &cut_node, &cut_tet))
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

    for(size_t fi = 0 ; fi < outside_face_cut_idx.size(); ++fi){
        //if(fa->is_outside_face(fa_cut->face2tet_[outside_face_cut_idx[fi]])){
        const size_t idx = fa->get_face_idx(
              cut_tet2tet[outside_face_cut(0,fi)],
            cut_tet2tet[outside_face_cut(1,fi)],
            cut_tet2tet[outside_face_cut(2,fi)]
            );

        assert(idx != -1);
        if(fa->is_outside_face(fa->face2tet_[idx])){
            outside_face_cut_in_orig_idx.push_back(fi);
            outside_face_cut_in_orig_type.push_back(face_idx2type[idx]);
          }
        // }
      }

    vector<size_t> dump_surface_cut;
    for(size_t fi = 0; fi < outside_face_cut_in_orig_idx.size(); ++fi){
        dump_surface_cut.insert(dump_surface_cut.end(),
                                outside_face_cut(colon(),outside_face_cut_in_orig_idx[fi]).begin(),
                                outside_face_cut(colon(),outside_face_cut_in_orig_idx[fi]).end());
      }

    ofstream ofs_cut("cut_surface_type.vtk");
    tri2vtk(ofs_cut, &cut_node[0], cut_node.size(2), &dump_surface_cut[0], dump_surface_cut.size()/3);
    cell_data(ofs_cut, &outside_face_cut_in_orig_type[0], outside_face_cut_in_orig_type.size(), "type");
  }
  return 0;
}
