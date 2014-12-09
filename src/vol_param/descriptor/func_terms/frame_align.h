#ifndef FRAME_ALIGN_H
#define FRAME_ALIGN_H

#include <hjlib/function/function.h>
#include <hjlib/sparse/sparse.h>
#include <jtflib/algorithm/equation.h>
#include "../def.h"
#include "../../common/util.h"


class frame_align_func: public hj::function::function_t<double,int32_t>
{
public:
  frame_align_func(const zjucad::matrix::matrix<double> & node,
                   const zjucad::matrix::matrix<size_t> & tet,
                   const zjucad::matrix::matrix<double> & frame,
                   const size_t idx, const double w,
                   const node_mapping * node_mapping_ = 0)
    :grad_op_(4, 3), tet_(tet), node_num_(node.size(2)), idx_(idx) , weight_(w),
      frame_(zjucad::matrix::trans(frame)), NM(node_mapping_){
    using namespace zjucad::matrix;
    using namespace std;

    zjucad::matrix::matrix<double> tet_node = node(colon(), tet_(colon(), idx_));
    if(calc_tet_def_grad_op(&tet_node[0], &grad_op_[0])) {
        std::cerr << "# degenerated tet" << __LINE__ << tet_node << std::endl;
      }
    assemble();
  }
  virtual ~frame_align_func(){}
public:
  virtual size_t dim_of_x(void) const{
    if(!NM)
      return node_num_ * 3;
    else
      return NM->ZT.size(1);
  }
  virtual size_t dim_of_f(void) const{
    return 9;
  }
  virtual int val(const double * x, double *f,
                  hj::function::func_ctx  *ctx=0)const{
    using namespace zjucad::matrix;
    itr_matrix<double *> f0(3, 3, f);
    if(!NM){
        const itr_matrix<const double *> T(3, node_num_, x);
        f0 = T(colon(), tet_(colon(), idx_))*grad_op_;
      }else{
        matrix<double> tet_node(3,4);
        get_cell_node(x, tet_(colon(),idx_), *NM, tet_node);
        f0 = tet_node * grad_op_;
      }
    assert(frame_.size(1) == 3 && frame_.size(2) == 3);
    f0 -= frame_;
    f0 *= weight_;

    for(size_t i = 0; i < f0.size(); ++i) jtf::math::erase_nan_inf(f0[i]);

    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t * idx = 0, hj::function::func_ctx * ctx =0) const {
    if(!NM){
        for(int c = 0; c < 3; ++c) {
            for(int r = 0; r < 3; ++r) {
                const int fi = c*3+r;
                ptr[fi+1] = ptr[fi] + tet_.size(1);
                for(int k = 0; k < tet_.size(1); ++k) {
                    idx[ptr[fi]+k] = tet_(k, idx_)*3+r;
                    val[ptr[fi]+k] = grad_op_(k, c) * weight_;
                    jtf::math::erase_nan_inf(val[ptr[fi]+k]);
                  }
              }
          }
      }else{
        assert(eqns_of_each_f.size() == 9);
        for(size_t fi = 0; fi < 9; ++fi){
            ptr[fi+1] = ptr[fi] + eqns_of_each_f[fi].e_vec_.size();
            size_t k = ptr[fi];
            for(const auto & exp : eqns_of_each_f[fi].e_vec_){
                assert(k != ptr[fi+1]);
                idx[k] = exp.index;

                val[k] = exp.coefficient * weight_;
                jtf::math::erase_nan_inf(val[k]);
                ++k;
              }
          }
      }

    return 0;
  }
  virtual size_t jac_nnz(void) const{
    if(!NM)
      return 4*dim_of_f();
    else
      return nnz_;
  }
  void assemble(){
    using  namespace zjucad::matrix;
    if(NM){
        get_node_mapping(tet_(colon(),idx_), NM->ZT, node_mapping_of_each_variabe_);
        nnz_ = 0;
        eqns_of_each_f.resize(9); // dim_of_f
        //f(i,j) = \sum_k^4 tn_i,3-k * grad_3-k, j
        for(size_t item = 0; item < 9; ++item){
            const size_t i = item % 3;
            const size_t j = item / 3;
            jtf::algorithm::equation<double> one_eqn;
            for(size_t k = 0; k < 4; ++k){
                add_to_eqn(one_eqn, node_mapping_of_each_variabe_[i + 3*(3-k)],grad_op_(3-k,j));
              }
            nnz_ += one_eqn.e_vec_.size();
            eqns_of_each_f[item] = one_eqn;
          }
      }
  }
private:
  const size_t node_num_;
  const size_t idx_;
  const double weight_;
  zjucad::matrix::matrix<double> grad_op_;
  const zjucad::matrix::matrix<size_t> &tet_;
  const zjucad::matrix::matrix<double> frame_;
  const node_mapping * NM;
  std::vector<std::vector<std::pair<size_t,double> > > node_mapping_of_each_variabe_;
  std::vector<jtf::algorithm::equation<double> > eqns_of_each_f;
  size_t nnz_ ;
}; 


//template <typename BASIC_INFO>
//hj_func * build_frame_align_func(
//    const zjucad::matrix::matrix<size_t> & mesh,
//    const zjucad::matrix::matrix<double> & node,
//    const zjucad::matrix::matrix<zjucad::matrix::matrix<double> >& frame,
//    const BASIC_INFO & bi,
//    zjucad::matrix::matrix<double> * vol_w = 0,
//    SIMPLEX sim = TET)
//{
//  if(vol_w){
//      if(vol_w->size() != mesh.size(2)){
//          std::cerr << "# [error] vol_weight size do not equal mesh.size(2) "
//                    << std::endl;
//          return nullptr;
//        }
//    }
//  using namespace std;
//  using namespace zjucad::matrix;

//  matrix<double> weight = ones<double>(mesh.size(2),1);

//  matrix<double> one_simplex_node = zeros<double>(3,mesh.size(1));
//  for(size_t ti = 0; ti < mesh.size(2); ++ti){
//      one_simplex_node = node(colon(), mesh(colon(), ti));
//      weight[ti] = (sim==TET?jtf::mesh::cal_tet_vol(one_simplex_node)
//                           :jtf::mesh::cal_face_area(one_simplex_node)) *
//          (vol_w?(*vol_w)[ti]:1.0);
//    }

//  const double total_w = std::accumulate(weight.begin(), weight.end(), 0.0);

//  std::shared_ptr<vector<std::shared_ptr<hj::function::function_t<double,int32_t> > > >
//      funcs(new vector<std::shared_ptr<hj::function::function_t<double,int32_t> > >(mesh.size(2)));

//  for(size_t ti = 0; ti < mesh.size(2); ++ti) {
//      if(sim == TET)
//        (*funcs)[ti].reset(
//            new frame_align_func<TET>( node, mesh, frame[ti], ti, sqrt(weight[ti] / total_w)));
//      else {
//          std::cerr << "# [error] unsupported simplex type." << std::endl;
//          return 0;
//        }
//    }
//  return hj::function::new_catenated_function<double, int32_t>(funcs);
//}

template <typename BASIC_INFO>
void build_frame_align_func(
    const zjucad::matrix::matrix<size_t> & mesh,
    const zjucad::matrix::matrix<double> & node,
    const zjucad::matrix::matrix<zjucad::matrix::matrix<double> >& frame,
    const BASIC_INFO & bi,
    std::vector<hj_func_cons_ptr> & func,
    zjucad::matrix::matrix<double> * vol_w = 0,
    SIMPLEX sim = TET,
    const node_mapping* node_mapping_ = 0)
{
  if(vol_w){
      if(vol_w->size() != mesh.size(2)){
          std::cerr << "# [error] vol_weight size do not equal mesh.size(2) "
                    << std::endl;
          return;
        }
    }

  if(func.size() != mesh.size(2))
    func.resize(mesh.size(2));

  using namespace std;
  using namespace zjucad::matrix;

  matrix<double> weight = ones<double>(mesh.size(2),1);

  matrix<double> one_simplex_node = zeros<double>(3,mesh.size(1));
  for(size_t ti = 0; ti < mesh.size(2); ++ti){
      one_simplex_node = node(colon(), mesh(colon(), ti));
      weight[ti] = (sim==TET?jtf::mesh::cal_tet_vol(one_simplex_node)
                           :jtf::mesh::cal_face_area(one_simplex_node)) *
          (vol_w?(*vol_w)[ti]:1.0);
    }

  const double total_w = std::accumulate(weight.begin(), weight.end(), 0.0);

  for(size_t ti = 0; ti < mesh.size(2); ++ti) {
      if(sim == TET)
        func[ti].reset(
              new frame_align_func( node, mesh, frame[ti], ti, sqrt(weight[ti] / total_w), node_mapping_));
      else {
          std::cerr << "# [error] unsupported simplex type." << std::endl;
          return;
        }
    }
}
#endif
