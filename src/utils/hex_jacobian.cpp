#include <numeric>
#include <fstream>
#include <iostream>
#include <vector>
#include <jtflib/mesh/mesh.h>
#include "../hexmesh/io.h"
#include "../common/cell_quality.h"

using namespace std;
using namespace zjucad::matrix;

int hex_jacobian(int argc, char * argv[])
{
  if(argc != 2){
      cerr << "# [usage] hex_jac hex." << endl;
      return __LINE__;
    }

  jtf::mesh::meshes hm;
  if(jtf::hexmesh::hex_mesh_read_from_wyz(argv[1], hm.mesh_,hm.node_,1))
    return __LINE__;

  matrix<double> jac(hm.mesh_.size(2),1);
  for(size_t hi = 0; hi < hm.mesh_.size(2); ++hi){
      jac[hi] = hex_scaled_jacobian(hm.mesh_(colon(), hm.mesh_(colon(),hi)));
    }

  cerr << "# [hex num] " << hm.mesh_.size(2) << endl;
  cerr << "# [scaled_jacobian] min/avg/max "
       << *min_element(jac.begin(), jac.end()) << "/"
       << std::accumulate(jac.begin(), jac.end(), 0.0) / jac.size() << "/"
       << *max_element(jac.begin(), jac.end()) << endl;

  size_t inverted_num = 0;
  for(size_t hi = 0; hi < jac.size(); ++hi)
    if(jac[hi] < 0) ++inverted_num;

  cerr << "# [scaled_jacobian < 0] number/percent " << inverted_num << "/"
       << 1.0*inverted_num/jac.size() << endl;

  ofstream ofs("hex_jac");
  const size_t total_step = hm.mesh_.size(2)/5000;
  const double step_delta = 1.0/total_step;
  for(size_t hi = 0; hi < hm.mesh_.size(2);++hi)
      ofs << jac[hi] << endl;

  //    ofstream ofs_sort("sort_hex_jac.csv");
  //    sort(jac.begin(), jac.end());

  //    //ofs_sort << "X,Y,Z" << endl;
  //    const size_t step = hm.mesh_.size(2) * 1.0 / 100.0;
  //    for(size_t hi = 0; hi < hm.mesh_.size(2);++hi){
  //        if(hi % step == 0)
  //        ofs_sort  << jac[hi] <<  endl;
  //      }

  //  const size_t total_step = 10;
  //  const double step_delta = 1.0/total_step;
  //  matrix<size_t> his = zeros<size_t>(total_step,1);
  //  matrix<double> per(total_step,1);
  //  for(size_t hi = 0; hi < hm.mesh_.size(2); ++hi){
  //      ++his[jac[hi]/step_delta] ;
  //    }

  //  ofstream ofs_his("hex_jac_histogram.csv");
  //  ofs_his << "jac,percent" << endl;
  //  for(size_t i = 0; i < his.size(); ++i)
  //    ofs_his << i * step_delta << "," << his[i]*1.0/hm.mesh_.size(2) << endl;

  //  ofstream ofs_number("hex_jac_histogram_num");

  //  for(size_t i = 0; i < his.size(); ++i)
  //    ofs_number << his[i] << endl;
  return 0;
}
