#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/math/math.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <jtflib/mesh/util.h>
#include <zjucad/matrix/io.h>
#include "../common/vtk.h"

using namespace std;
using namespace zjucad::matrix;

static double signed_area(const double uv[6]){
        return (uv[2]-uv[0])*(uv[5]-uv[1]) - (uv[4]-uv[0])*(uv[3]-uv[1]);
    }

int deform_gradient_2d(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] deform_gradient_2d obj_0 obj_1 vtk" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes obj0, obj1;
  if(jtf::mesh::load_obj(argv[1], obj0.mesh_, obj0.node_))
    return __LINE__;
  if(jtf::mesh::load_obj(argv[2], obj1.mesh_, obj1.node_))
    return __LINE__;


  matrix<double> def(obj0.mesh_.size(2),1);
  matrix<double> edge0, edge1;
  matrix<double> tri0, tri1;
  matrix<double> p2e = zeros<double>(3,2);
  p2e(colon(0,1),colon(0,1)) = eye<double>(2);
  p2e(2,colon()) = ones<double>(2,1)*-1;
  matrix<double> jac,U,S,VT;

  matrix<double> area(obj0.mesh_.size(2),1);

  for(size_t fi = 0; fi < obj0.mesh_.size(2); ++fi){
      tri0 = obj0.node_(colon(0,1),obj0.mesh_(colon(),fi));
      tri1 = obj1.node_(colon(0,1),obj1.mesh_(colon(),fi));
      edge0 = tri0*p2e;
      edge1 = tri1*p2e;
      jtf::math::invert_2dmatrix(edge0);
      jac = edge1*edge0;
      svd(jac, U,S,VT);
      def[fi] = (std::max(S(0,0),S(1,1))/std::min(S(0,0),S(1,1)));

      area[fi]=signed_area(&tri1[0]);//jtf::mesh::cal_face_area(obj0.mesh_(colon(),fi), obj1.node_);
    }

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);

  cerr << "# [info] area averaged_distortion " << dot(area, def)/total_area << endl;
  ofstream ofs(argv[3]);
  tri2vtk(ofs, &obj1.node_[0], obj1.node_.size(2), &obj1.mesh_[0], obj1.mesh_.size(2));
  cell_data(ofs, &def[0], def.size(), "def");
  vtk_data(ofs, &area[0], area.size(), "area");
  return 0;
}

int test2(int argc, char * argv[])
{
  if(argc != 5){
      cerr << "# [usage] test obj0 obj1 input_mapping output_obj1" << endl;
      return __LINE__;
    }
  jtf::mesh::meshes obj0, obj1;
  if(jtf::mesh::load_obj(argv[1], obj0.mesh_, obj0.node_))
    return __LINE__;
  if(jtf::mesh::load_obj(argv[2], obj1.mesh_, obj1.node_))
    return __LINE__;

  matrix<size_t> mapping(obj0.node_.size(2),1);
  ifstream ifs(argv[3]);
  size_t a,b,c,d;
  for(size_t i = 0; i < 32; ++i){
      ifs >> a >> b >> c >> d;
      mapping[c] = a;
      mapping[d] = b;
    }

  matrix<double> new_node = obj1.node_;
  for(size_t i = 0; i < mapping.size(); ++i)
    new_node(colon(),mapping[i]) = obj1.node_(colon(), i);


  for(size_t i = 0; i < obj1.mesh_.size(); ++i)
    obj1.mesh_[i] = mapping[obj1.mesh_[i]];
  jtf::mesh::save_obj(argv[4], obj1.mesh_, new_node);

  return 0;
}

int test3(int argc, char * argv[])
{
  if(argc != 5){
      cerr << "# [usage] test3 input_obj cons source_cons_vtk target_cons_vtk" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes obj0;
  if(jtf::mesh::load_obj(argv[1], obj0.mesh_, obj0.node_))
    return __LINE__;

  map<size_t, pair<double,double> > cons_points;
  {
    ifstream ifs(argv[2]);
    if(ifs.fail()){
        cerr << "# [error] can not open cons point file." << endl;
        return __LINE__;
      }
    size_t idx; double dx, dy;
    while(!ifs.eof()){
        ifs >> idx >> dx >> dy;
        if(ifs.eof()) break;
        cons_points[idx] = make_pair(dx,dy);
      }
  }

  matrix<size_t> select_points(cons_points.size(),1);
  matrix<double> select_nodes= zeros<double>(3, cons_points.size());

  size_t pi = 0;
  for(const auto & one_p : cons_points){
      select_points[pi] = one_p.first;
      select_nodes(0,pi) = one_p.second.first;
      select_nodes(1,pi) = one_p.second.second;
      ++pi;
    }

  ofstream ofs_source(argv[3]);
  point2vtk(ofs_source, &obj0.node_[0], obj0.node_.size(2), &select_points[0], select_points.size());

  ofstream ofs_target(argv[4]);

  matrix<size_t> fake_indx = colon(0,select_points.size()-1);
  point2vtk(ofs_target, &select_nodes[0], select_nodes.size(2), &fake_indx[0], fake_indx.size());

  return 0;
}


int tri_mesh_read_from_vtk(
    const char *path,
    zjucad::matrix::matrix<double> &node,
    zjucad::matrix::matrix<size_t> &tri)
{
  ifstream ifs(path);
  if(ifs.fail()) {
      cerr << "[info] " << "can not open file" << path << endl;
      return __LINE__;
    }

  string str;
  int point_num = 0,cell_num = 0;

  vector<size_t> tri_temp;
  while(!ifs.eof()){
      ifs >> str;
      if(str == "POINTS"){
          ifs >> point_num >> str;
          node.resize(3, point_num);
          for(size_t i = 0;i < point_num; ++i){
              for(size_t j = 0;j < 3; ++j)
                ifs >> node(j, i);
            }
          continue;
        }
      if(str == "POLYGONS"){
          ifs >> cell_num >> str;
          int point_number_of_cell = 0;
          for(size_t ci = 0; ci < cell_num; ++ci){
              ifs >> point_number_of_cell;
              if(point_number_of_cell != 3){
                  for(size_t i = 0; i < point_number_of_cell; ++i)
                    ifs >> str;
                }else{
                  int p;
                  for(size_t i = 0; i < point_number_of_cell; ++i){
                      ifs >> p;
                      tri_temp.push_back(p);
                    }
                }
            }
        }
    }
  tri.resize(3, tri_temp.size()/3);
  copy(tri_temp.begin(), tri_temp.end(), tri.begin());
  return 0;
}

int test4(int argc, char *argv[])
{
  if(argc != 3){
      cerr << "# [usage] test4 input_vtk ouptut_vtk" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes tri;
  tri_mesh_read_from_vtk(argv[1],tri.node_, tri.mesh_);

  cerr << "node " << tri.node_.size(2) << endl;
  cerr << "tri " << tri.mesh_.size(2) << endl;
  matrix<double> v = ones<double>(tri.mesh_.size(2),1);

  ofstream ofs(argv[2]);
  tri2vtk(ofs, &tri.node_[0], tri.node_.size(2), &tri.mesh_[0], tri.mesh_.size(2));
  cell_data(ofs, &v[0], v.size(), "face_scalar");

  return 0;
}

int test5(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [usage] test5 tet0 tet1" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes tet0, tet1;
  jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tet0.node_, &tet0.mesh_);
  jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &tet1.node_, &tet1.mesh_);

  jtf::mesh::tet_mesh_write_to_zjumat(argv[1], &tet0.node_, &tet1.mesh_);

  return 0;
}
