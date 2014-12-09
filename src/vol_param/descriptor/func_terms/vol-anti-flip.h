#ifndef VOLANTIFLIP_H
#define VOLANTIFLIP_H

#include <hjlib/function/function.h>
#include "../../common/util.h"
#include "../def.h"
#include "../../../mesh_func/tet-vol.h"
#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/util.h>

class vol_anti_flip_hj :public hj_func
{
public:
  vol_anti_flip_hj(const zjucad::matrix::matrix<size_t> &all_tets,
                   const size_t idx,
                   const double weight,
                   const node_mapping * nm = 0)
    :all_tets_(all_tets), idx_(idx), w_(weight), nm_(nm){
    if(nm)
      assemble();
  }
  virtual size_t dim_of_x() const{
    if(!nm_){
        return (zjucad::matrix::max(all_tets_)+1) * 3;
      }else{
        return nm_->ZT.size(1);
      }
  }
  virtual size_t dim_of_f() const{
    return 1;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx) const{
    using namespace zjucad::matrix;
    zjucad::matrix::matrix<double> tet_node;
    if(!nm_){
        itr_matrix<const double*> x0(3, dim_of_x()/3, x);
        tet_node = x0(colon(), all_tets_(colon(), idx_));
      }else{
        get_cell_node(x, all_tets_(colon(),idx_), *nm_, tet_node);
      }

    double v = 0;
    calc_tet_vol_(&v, &tet_node[0]);
    *f = v * w_;
    return 0;
  }
  virtual int jac(const double *x, double * val, int32_t * ptr, int32_t * idx,
                  hj::function::func_ctx *ctx)const{
    using namespace zjucad::matrix;
    zjucad::matrix::matrix<double> tet_node;
    double  g[12];
    if(!nm_){
        itr_matrix<const double*> x0(3, dim_of_x()/3, x);
        tet_node = x0(colon(), all_tets_(colon(),idx_));
      }else{
        get_cell_node(x, all_tets_(colon(),idx_), *nm_, tet_node);
      }

    calc_tet_vol_jac_(&g[0], &tet_node[0]);

    if(!nm_){
        ptr[1] = ptr[0] + 12;
        for(size_t i = 0; i < 12; ++i){
            const size_t point_idx = i / 3;
            const size_t di = i%3;
            idx[ptr[0]+i] = 3*all_tets_(point_idx,idx_) + di;
            val[ptr[0]+i] = g[i] * w_;
          }

      }else{
        assert(real_v_to_orig_v_coeff_.size());
        ptr[1] = ptr[0] + real_v_to_orig_v_coeff_.size();
        size_t idx_ = ptr[0];
        for(const auto & one_real_v : real_v_to_orig_v_coeff_){
            idx[idx_] = one_real_v.first;
            val[idx_] = 0;
            const std::vector<std::pair<size_t,double> >  & to_orig_v = one_real_v.second;
            for(const auto & one_exp : to_orig_v){
                assert(one_exp.first < 12);
                val[idx_] += g[one_exp.first] * one_exp.second * w_;
              }
            ++idx_;
          }
      }
    return 0;
  }
  virtual size_t jac_nnz() const {
    if(!nm_){
        return 3 * 4;
      }else{
        return real_v_to_orig_v_coeff_.size();
      }
  }
  void assemble(){
    using namespace zjucad::matrix;
    if(nm_){
        get_node_mapping(all_tets_(colon(),idx_), nm_->ZT, node_mapping_of_each_variable_);
        for(size_t ni = 0; ni < node_mapping_of_each_variable_.size(); ++ni){
            const std::vector<std::pair<size_t,double> > & one_eqn =
                node_mapping_of_each_variable_[ni];
            for(const auto & one_exp : one_eqn){
                real_v_to_orig_v_coeff_[one_exp.first].push_back(
                      std::make_pair(ni, one_exp.second));
              }
          }
        for(auto & one_real_v : real_v_to_orig_v_coeff_){
            std::vector<std::pair<size_t,double> > & to_real_v = one_real_v.second;
            std::sort(to_real_v.begin(), to_real_v.end());
          }
      }
  }
private:
  const zjucad::matrix::matrix<size_t> &all_tets_;
  const size_t idx_;
  const double w_;
  const node_mapping * nm_;
  std::vector<std::vector<std::pair<size_t,double> > > node_mapping_of_each_variable_;
  std::map<size_t, std::vector<std::pair<size_t,double> > > real_v_to_orig_v_coeff_;
};



inline hj::function::function_t<double, int32_t> *
build_vol_anti_flip_hj_func(const zjucad::matrix::matrix<size_t> & mesh,
                            const zjucad::matrix::matrix<double> & node,
                            zjucad::matrix::matrix<double> &weight,
                            const node_mapping * node_mapping_ = 0)
{
  using namespace std;
  using namespace zjucad::matrix;

  std::shared_ptr<vector<hj_func_ptr> >
      funcs(new vector<hj_func_ptr>);

  weight.resize(mesh.size(2));
  for(size_t ti = 0; ti < mesh.size(2); ++ti)
    weight[ti] = jtf::mesh::cal_tet_vol(node(colon(), mesh(colon(),ti)));
  const double avg_vol = std::accumulate(weight.begin(), weight.end(),0.0) / mesh.size(2);

  size_t count = 0;
  for(size_t ti = 0; ti < mesh.size(2); ++ti){
      if(weight[ti] < 0){
          funcs->push_back(hj_func_ptr(new vol_anti_flip_hj(mesh, ti,1.0,node_mapping_)));
          ++count;
        }
    }
  cerr << "# [info] add " << count << "/" << mesh.size(2) << " = " << 1.0 * count/mesh.size(2)
       << " vol-anti-flip constraints." << endl;
  if(funcs->empty()) return 0;
  return hj::function::new_catenated_function<double, int32_t>(funcs);
}

#endif // VOLANTIFLIP_H
