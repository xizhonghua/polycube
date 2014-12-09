#include "../tetmesh/tetmesh.h"
#include <hjlib/math_func/math_func.h>
#include <hjlib/math_func/func_aux.h>
#include <hjlib/math_func/operation.h>
#include <zjucad/ptree/ptree.h>
#include <hjlib/math/polar.h>
#include <jtflib/optimizer/optimizer.h>
#include "../common/util.h"
#include "rot.h"
using namespace std;
using namespace zjucad::matrix;

template <typename val_type, typename int_type>
class rot_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  rot_func(const double w):w_(w){}
  virtual ~rot_func(){}
  virtual size_t nx() const{
    return 9;
  }
  virtual size_t nf() const{
    return 1;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==0){
        itr_matrix<const val_type*> x0(3,3,&x[0]);
        val_type d;
        rot_matrix_(&d, &x0[0]);
        int_type c[1] = {0};
        cv[c] += d*w_;
      }
    if(k==1){
        itr_matrix<const val_type*> x0(3,3,&x[0]);
        matrix<val_type> jac(9,1);
        rot_matrix_jac_(&jac[0], &x0[0]);
        for(int_type i = 0; i < 9; ++i){
            int_type c[2] = {0,i};
            cv[c] += jac[i]*w_;
          }
      }
    if(k==2){
        itr_matrix<const val_type*> x0(3,3,&x[0]);
        matrix<val_type> hes(9,9);
        rot_matrix_hes_(&hes[0], &x0[0]);
        for(int_type i = 0; i < 9; ++i){
            for(int_type j = 0; j < 9; ++j){
                int_type c[3] = {0, i, j} ;
                cv[c] += w_*hes(i,j);
              }
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==1){
        for(int_type i = 0; i < 9; ++i){
            int_type c[2] = {0,i};
            l2g.add(cs,c);
          }
      }
    if(k==2){
        for(int_type i = 0; i < 9; ++i){
            for(int_type j = 0; j < 9; ++j){
                int_type c[3] = {0, i, j} ;
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k==0) return -1;
    if(k==1) return 9;
    if(k==2) return 81;
  }
private:
  const val_type w_;
};


template <typename val_type, typename int_type>
class rot_fit_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  rot_fit_func(const zjucad::matrix::matrix<val_type> & p0,
               const zjucad::matrix::matrix<val_type> & p1,
               const double w):w_(w),p0_(p0), p1_(p1){}
  virtual ~rot_fit_func(){}
  virtual size_t nx() const{
    return 9;
  }
  virtual size_t nf() const{
    return 3;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==0){
        itr_matrix<const val_type*> x0(3,3,&x[0]);
        matrix<val_type> diff = x0*p0_-p1_;
        for(int_type i = 0 ; i< 3; ++i){
            int_type c[1] = {i};
            cv[c] += diff[i]*w_;
          }
      }
    if(k==1){
        for(int_type i = 0 ; i< 3; ++i){
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i,j*3+i};
                cv[c] += w_*p0_[j];
              }
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==1){
        for(int_type i = 0 ; i< 3; ++i){
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i,j*3+i};
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k==0) return -1;
    if(k==1) return 9;
  }
private:
  const val_type w_;
  const zjucad::matrix::matrix<val_type> p0_;
  const zjucad::matrix::matrix<val_type> p1_;
};

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_rot_fit_func(const zjucad::matrix::matrix<double> & node0,
                   const zjucad::matrix::matrix<double> & node1,
                   const double weight)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  assert(node0.size(2) == node1.size(2));
  for(size_t ti = 0; ti < node0.size(2); ++ti){
      all_func->push_back(
            math_func_ptr(
              new rot_fit_func<double,int32_t>(node0(colon(),ti), node1(colon(),ti),
                                               sqrt(weight/node0.size(2)))));
    }

  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

void map_tets(jtf::tet_mesh &tm0, const jtf::tet_mesh &tm1,
              boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);

  double rot_w = 1;
  {
    math_func_ptr rot_f(new rot_func<double,int32_t>(rot_w));
    func->push_back( math_func_ptr(new hj::math_func::sum<double,int32_t>(rot_f)));
  }
  {
    func->push_back(math_func_ptr(new hj::math_func::sumsqr<double,int32_t>(
                                    build_rot_fit_func(tm0.tetmesh_.node_, tm1.tetmesh_.node_, 1.0))));
  }
  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  matrix<double> rot = rand<double>(3);
  jtf::optimize(*obj, rot, pt, nullptr, nullptr, nullptr);

  cerr << rot << endl;
  hj::polar3d p;
  p(rot,2);

  tm0.tetmesh_.node_ = temp(rot * tm0.tetmesh_.node_);
}

void map_tris(const zjucad::matrix::matrix<size_t> & tri,
              const zjucad::matrix::matrix<double> &tri0_node,
              const zjucad::matrix::matrix<double> &tri1_node,
              boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);

  zjucad::matrix::matrix<double> tri1_node_new = tri1_node;
  zjucad::matrix::matrix<double> tri0_node_new = tri0_node;
  matrix<double> avg_node1(3,1), avg_node0(3,1);
  cal_average_node(tri1_node_new, avg_node1);
  cal_average_node(tri0_node, avg_node0);

  tri0_node_new -= avg_node0*ones<double>(1,tri0_node.size(2));
  tri1_node_new -= avg_node1*ones<double>(1,tri1_node.size(2));

  double rot_w = 1;
  {
    math_func_ptr rot_f(new rot_func<double,int32_t>(rot_w));
    func->push_back( math_func_ptr(new hj::math_func::sum<double,int32_t>(rot_f)));
  }
  {
    func->push_back(math_func_ptr(new hj::math_func::sumsqr<double,int32_t>(
                                    build_rot_fit_func(tri0_node_new, tri1_node_new, 1.0))));
  }
  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  matrix<double> rot = rand<double>(3);
  jtf::optimize(*obj, rot, pt, nullptr, nullptr, nullptr);

  hj::polar3d p;
  p(rot,2);

  cerr << "# R*M0 + g= M1" << endl;
  cerr << "# R = " << rot << endl;
  matrix<double> g = avg_node1 - rot * avg_node0;
  cerr << "# g = " << g << endl;

  matrix<double> temp_node = tri0_node;
  temp_node = temp(rot * temp_node);
  temp_node += g * ones<double>(1, temp_node.size(2));

  jtf::mesh::save_obj("test.obj", tri, temp_node);
}
