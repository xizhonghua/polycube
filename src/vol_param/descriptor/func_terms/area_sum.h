#ifndef AREA_SUM_H
#define AREA_SUM_H

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/sparse/sparse.h>
#include "../def.h"
#include "../../../mesh_func/tri-area.h"
#include "../../common/util.h"
class area_sum : public jtf_func
{
public:
  area_sum(size_t node_num, const zjucad::matrix::matrix<size_t> &faces,
           double ori_area,
           const node_mapping * NM = 0)
    :node_num_(node_num), faces_(faces), ori_area_(ori_area),node_mapping_(NM) {
  }
  virtual size_t dim(void) const {
    if(!node_mapping_){
        return node_num_ * 3;
      }else{
        return node_mapping_->ZT.size(1);
      }
  }
  //NOTE: add to v
  virtual int val(const double *x, double &v) {
    using namespace zjucad::matrix;

    zjucad::matrix::matrix<double> tri(3, 3);
    v = -ori_area_;
    for(size_t fi = 0; fi < faces_.size(2); ++fi) {
        if(!node_mapping_){
            const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
            tri = T(colon(), faces_(colon(), fi));
          }else{
            tri.resize(3, faces_.size(1));
            get_cell_node(x, faces_(colon(),fi), *node_mapping_, tri);
          }
        double a = 0;
        calc_tri_area_(&a, &tri[0]);
        v += a;
      }
    return 0;
  }

  virtual int gra(const double *x, double *g) {
    using namespace zjucad::matrix;

    matrix<double> tri(3, 3), g9(9);
    for(size_t fi = 0; fi < faces_.size(2); ++fi) {

        if(!node_mapping_){
            const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
            tri = T(colon(), faces_(colon(), fi));
          }else{
            tri.resize(3, faces_.size(1));
            get_cell_node(x, faces_(colon(),fi), *node_mapping_, tri);
          }

        g9(colon()) = 0;
        calc_tri_area_jac_(&g9[0], &tri[0]);
        for(size_t i = 0; i < 9; ++i) {
            if(!node_mapping_)
              g[i%3+faces_(i/3, fi)*3] += g9[i];
            else{
                const size_t orig_v = i%3+faces_(i/3, fi)*3;
                for(size_t idx = node_mapping_->ZT.ptr()[orig_v];
                    idx != node_mapping_->ZT.ptr()[orig_v+1]; ++idx){
                    g[node_mapping_->ZT.idx()[idx]] += g9[i] * node_mapping_->ZT.val()[idx];
                  }
              }
          }
      }
    return 0;
  }

  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx) {
    using namespace zjucad::matrix;

    if(g == 0 && idx == 0) {
        nnz = dim();
        return 0;
      }

    zjucad::matrix::itr_matrix<double*> g_m(dim(),1, g);
    g_m *= 0;
    gra(x, g);
    for(size_t i = 0; i < dim(); ++i) idx[i] = i;
    return 0;
  }

  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h,
                  int32_t *ptr, int32_t *idx, double alpha = 1) {
    nnz = 0;
    return 1;
  }
  virtual int hes_block(const double *x, double *h, double alpha = 1) {return -1;}
protected:
  const zjucad::matrix::matrix<size_t> &faces_;
  const size_t node_num_;
  double ori_area_;
  const node_mapping * node_mapping_;
};

#endif // AREA_SUM_H
