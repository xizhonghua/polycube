#include <fstream>
#include <numeric>
#include <jtflib/mesh/io.h>

#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
#include "../common/cell_quality.h"

using namespace std;
using namespace zjucad::matrix;


int tet_jacobian(int argc, char * argv[])
{
  if(argc != 3){
    cerr << "# [usage] tet_jacobian tet len."  << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  matrix<double> jac(tm.mesh_.size(2),1);

  for(size_t ti = 0; ti < tm.mesh_.size(2); ++ti){
    jac[ti] = tet_scaled_jacobian(tm.node_(colon(), tm.mesh_(colon(),ti)));
  }

  cerr << "# [scaled_jac] min/avg/max "
       << *min_element(jac.begin(), jac.end()) << "/"
       << std::accumulate(jac.begin(), jac.end(),0.0)/jac.size() << "/"
       << *max_element(jac.begin(), jac.end()) << endl;

  size_t inverted_tet_num = 0;
  for(size_t ti = 0; ti < jac.size(); ++ti)
    if(jac[ti] < 0)
      ++inverted_tet_num;
  cerr << "# [inverted_jac] inverted_tet number/percent " << inverted_tet_num
       << "/"  << 1.0 * inverted_tet_num  / jac.size() << endl;

  int len = atoi(argv[2]);
  if(len < 1){
    cerr << "# wrong len." << endl;
    return __LINE__;
  }

  vector<double> tet_percent(len);
  vector<size_t> tet_number(len,0);


  const double delta = (1.0-(-1.0))/len;

  for(size_t ti = 0; ti < tm.mesh_.size(2); ++ti){
    const size_t idx = static_cast<size_t>((jac[ti] - (-1))/delta);
    if(idx >= tet_number.size()) continue;
    ++tet_number[idx] ;
  }

  for(size_t i = 0; i < len; ++i){
    tet_percent[i] = 1.0*tet_number[i] / tm.mesh_.size(2);
  }

  ofstream ofs("tet_jac");
  for(size_t i = 0; i < jac.size(); ++i)
    ofs << jac[i] << endl;
  return 0;
}
