#include <fstream>

#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <jtflib/mesh/io.h>

#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"

using namespace std;
using namespace zjucad::matrix;

double rigid_err(const matrix<double> & node_a,
                 const matrix<double> & node_b)
{
  assert(node_a.size(2) == 4);
  assert(node_b.size(2) == 4);

  matrix<double> Ta = node_a;
  matrix<double> Tb = node_b;
  for(size_t i =0; i < node_a.size(2); ++i){
    Ta(colon(),i) -= Ta(colon(),0);
    Tb(colon(),i) -= Tb(colon(),0);
  }

  matrix<double> r(3,3);
  // R*Ta=Tb
  // R*(Ta*Ta^T)=Tb*(Ta^T)
  // R = Tb*(Ta^T)*(Ta*Ta^T)^-1
  matrix<double> bat = Tb * trans(Ta);
  matrix<double> aat = Ta * trans(Ta);
  if(inv(aat)){
    cerr << "# [error] inverse fail." << endl;
  }
  r = temp(bat * aat);
  matrix<double> rrt = r * trans(r);
  return norm(rrt - eye<double>(3));
}

int rigid_shape_error(int argc, char *argv[])
{
  if(argc != 4){
    cerr << "# [usage] rigid_shape_error tet_a tet_b err_vtk" << endl;
    return __LINE__;
  }
  jtf::mesh::meshes a,b;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &a.node_, &a.mesh_))
    return __LINE__;

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &b.node_, &b.mesh_))
    return __LINE__;

  if(a.mesh_.size(2) != b.mesh_.size(2)){
    cerr << "# [error] tet a is not compatiable with tet b" << endl;
    return __LINE__;
  }

  matrix<double> err(a.mesh_.size(2));
  for(size_t ti = 0; ti < a.mesh_.size(2); ++ti){
    err[ti] = rigid_err(a.node_(colon(), a.mesh_(colon(),ti)),
                        b.node_(colon(), b.mesh_(colon(),ti)));
  }

  cerr << "# [info] max err " << *max_element(err.begin(), err.end()) << endl;
  const size_t max_tet_idx = max_element(err.begin(), err.end())  - err.begin();
  cerr << "# [info] max element idx " << max_tet_idx << endl;

  matrix<double> max_a = a.node_(colon(),a.mesh_(colon(), max_tet_idx));

  matrix<double> max_b = b.node_(colon(),b.mesh_(colon(), max_tet_idx));

  cerr << "# [info] max_tet_a " << a.mesh_(colon(), max_tet_idx) << endl;
  cerr << "# [info] max_tet_a " << max_a << endl;
  cerr << "# [info] max_tet_b " << b.mesh_(colon(), max_tet_idx) << endl;
  cerr << "# [info] max_tet_b " << max_b << endl;

  for(size_t ti = 0; ti < max_a.size(1); ++ti){
    max_a(colon(),ti) -= max_a(colon(),0);
    max_b(colon(),ti) -= max_b(colon(),0);
  }

  cerr << "# [info] max_a_tra " << max_a << endl;
  cerr << "# [info] max_b_tra " << max_b << endl;

//  {

//    matrix<double> r(3,3);
//    // R*Ta=Tb
//    // R*(Ta*Ta^T)=Tb*(Ta^T)
//    // R = Tb*(Ta^T)*(Ta*Ta^T)^-1
//    matrix<double> bat = max_b * trans(max_a);
//    matrix<double> aat = max_a * trans(max_a);
//    if(inv(aat) == false){
//      cerr << "# [error] inverse fail." << endl;
//    }
//    r = temp(bat * aat);
//    matrix<double> rrt = r * trans(r);
//    cerr << "# max r " << r << endl;
//  }
  cerr << "# [info] min err " << *min_element(err.begin(), err.end()) << endl;

  ofstream ofs(argv[3]);
  if(ofs.fail()){
    cerr << "# [error] can not open file." << endl;
    return __LINE__;
  }
  tet2vtk(ofs, &a.node_[0], a.node_.size(2), &a.mesh_[0], a.mesh_.size(2));
  cell_data(ofs, &err[0], err.size(), "rigid_err");

  return 0;
}

