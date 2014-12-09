#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <fstream>
#include <vector>

#include <zjucad/matrix/matrix.h>
#include <hjlib/sparse/sparse.h>
#include <jtflib/mesh/io.h>
#include <zjucad/linear_solver/linear_solver.h>

#include "../common/vtk.h"
#include "../tetmesh/tetmesh.h"
#include "../numeric/util.h"



using namespace std;
using boost::property_tree::ptree;
//using namespace hj::function;
using namespace zjucad::matrix;
using namespace hj::sparse;

typedef matrix<double> matrixd_t;
typedef hj::sparse::csc<double, int32_t> cscd_t;

int load_constraints(const char *path, size_t node_num, cscd_t &CT, matrixd_t &b)
{
  ifstream ifs(path);
  if(ifs.fail())
    return __LINE__;
  vector<pair<size_t, double> > vC;
  pair<size_t, double> entry;
  while(1) {
    ifs >> entry.first >> entry.second;
    if(entry.first >= node_num*3)
      return __LINE__;
    if(ifs.fail())
      break;
    vC.push_back(entry);
  }

  CT.resize(node_num*3, vC.size(), vC.size());
  b.resize(vC.size(), 1);
  for(size_t ci = 0; ci < CT.size(2); ++ci) {
    CT.ptr()[ci+1] = CT.ptr()[ci]+1;
    CT.idx()[CT.ptr()[ci]] = vC[ci].first;
    b[ci] = vC[ci].second;
  }
  CT.val()(colon()) = 1;
  return 0;
}

void laplace_operator(const matrixst &tet, const matrix<int> &is_boundary, cscd_t &LT)
{
  const double internal = 1e-3;
  const size_t node_num = max(tet)+1;
  hj::sparse::csc_by_vm<double> vmLT(node_num*3, node_num*3);
  for(size_t ti = 0; ti < tet.size(2); ++ti) {
    for(size_t vi = 0; vi < 4; ++vi) {
      for(size_t ni = 0; ni < 4; ++ni) {
        if(ni == vi) continue;
        for(int d = 0; d < 3; ++d) {
          const size_t i = tet(vi, ti), j = tet(ni, ti);
          if(is_boundary[i]) {
            if(is_boundary[j])
              vmLT[i*3+d][j*3+d] = -1;
            else
              vmLT[i*3+d][j*3+d] = 0;
          }
          else {
              vmLT[i*3+d][j*3+d] = -internal;
          }
        }
      }
    }
  }
  for(size_t vi = 0; vi < vmLT.size(2); ++vi)
    vmLT[vi][vi] = 0;
  convert(vmLT, LT);
  for(size_t vi = 0; vi < LT.size(2); ++vi) {
    size_t pos = -1;
    double val = 0;
    for(size_t nzi = LT.ptr()[vi]; nzi < LT.ptr()[vi+1]; ++nzi) {
      if(LT.idx()[nzi] == vi)
        pos = nzi;
      val -= LT.val()[nzi];
    }
    LT.val()[pos] = val;
  }
}

int laplace_deform(matrixd_t &x, const jtf::mesh::meshes &tm, const cscd_t &CT, const matrixd_t &b,
                   const matrix<int> &is_boundary, boost::property_tree::ptree &pt)
{
  cscd_t LT;
  laplace_operator(tm.mesh_, is_boundary, LT);
#if 0
  for(size_t i = 0; i < 10; ++i) {
    cout << LT.idx_(colon(LT.ptr()[i], LT.ptr()[i+1]-1))
         << LT.val_(colon(LT.ptr()[i], LT.ptr()[i+1]-1)) << endl;
  }
#endif

  cscd_t LTL, CTC;
  AAT<>(LT, LTL);
  AAT<>(CT, CTC);

  double alpha = 1e4;
  for(size_t ai = 0; ai < 1; ++ai) {
    cout << "# alpha: " << alpha << endl;

    cscd_t aCTC = CTC;
    aCTC.val() *= alpha;

    cscd_t H = LTL;
    acc(aCTC, H);

    matrixd_t D;
    diag(H, D);

    matrixd_t aCTb = zeros<double>(CT.size(1), 1);
    mv(false, CT, b, aCTb);
    aCTb *= alpha;

    std::unique_ptr<linear_solver> slv_;
    slv_.reset(linear_solver::create(&H.val()[0], &H.idx()[0], &H.ptr()[0],
        H.nnz(), H.size(1), H.size(2), pt));
    if(slv_.get()) {
      slv_->solve(&aCTb[0], &x[0], 1, pt);
    }
    else {
      cerr << "# solver create error." << endl;
      return __LINE__;
    }
    alpha *= 2;
  }
  return 0;
}

int deform0(ptree &pt)
{
  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm.node_, &tm.mesh_))
    return __LINE__;

  cscd_t CT; // C^T*x = b
  matrixd_t b;
  if(load_constraints(pt.get<string>("C.value").c_str(), tm.node_.size(2), CT, b))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  matrixst boundary;
  get_outside_face(*fa, boundary);
  matrix<int> is_boundary = zeros<int>(tm.node_.size(2), 1);
  is_boundary(boundary) = 1;

  matrixd_t x = tm.node_;
  if(laplace_deform(x, tm, CT, b, is_boundary, pt))
     return __LINE__;

  {
    ofstream ofs("input.vtk");
    tet2vtk(ofs, &tm.node_[0], tm.node_.size(2), &tm.mesh_[0], tm.mesh_.size(2));
    ofs << "POINT_DATA " << tm.node_.size(2) << "\n";
    matrixd_t rgba = zeros<double>(4, tm.node_.size(2));
    rgba(colon(0, 2), colon())(CT.idx()) = b;
    rgba(3, colon()) = 1;
    vtk_data_rgba(ofs, rgba.begin(), tm.node_.size(2), "Cb");
    for(size_t i = 0; i < CT.idx().size(); ++i) {
      size_t pt_idx = CT.idx()[i]/3;
      is_boundary[pt_idx] = 0;
    }
    cout << b.size() << endl;
  }
  cout << sum(is_boundary) << endl;


  {
    ofstream ofs("output.vtk");
    tet2vtk(ofs, &x[0], tm.node_.size(2), &tm.mesh_[0], tm.mesh_.size(2));
  }
  return 0;
}
