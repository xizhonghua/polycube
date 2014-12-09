#ifndef NORMAL_DIFF_H
#define NORMAL_DIFF_H

#include <vector>
#include <algorithm>
#include <numeric>

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <hjlib/sparse/sparse.h>
#include <jtflib/function/func_aux.h>
#include "../../../mesh_func/tri-area.h"

#include "../../common/util.h"
#include "../def.h"
#include "../../../mesh_func/tri-normal.h"


template <typename T1>
double edge_area(
    const zjucad::matrix::matrix_expression<T1> & tri_a,
    const zjucad::matrix::matrix_expression<T1> & tri_b,
    const std::vector<zjucad::matrix::matrix<size_t> > & tri_a_connect,
    const std::vector<zjucad::matrix::matrix<size_t> > & tri_b_connect,
    const double total_area,
    const double * x,
    const size_t dim, const size_t node_num,
    const node_mapping * NM = 0)
{
  using namespace zjucad::matrix;

  matrix<double> tri_node;
  double area_a, area_b;

  if(!NM){
      itr_matrix<const double *> x0(dim, node_num, x);
      tri_node = x0(colon(), tri_a);
    }else{
      get_cell_node(x, tri_a, *NM, tri_node);
    }

  calc_tri_area_(&area_a, &tri_node[0]);

  if(!NM){
      itr_matrix<const double *> x0(dim, node_num, x);
      tri_node = x0(colon(), tri_b);
    }else{
      get_cell_node(x, tri_b, *NM, tri_node);
    }

  calc_tri_area_(&area_b, &tri_node[0]);

  double total_area_a = 0, total_area_b = 0, temp;
  for(size_t i = 0; i < tri_a_connect.size(); ++i){
      if(!NM){
          itr_matrix<const double *> x0(dim, node_num, x);
          tri_node = x0(colon(), tri_a_connect[i]);
        }else{
          get_cell_node(x, tri_a_connect[i], *NM, tri_node);
        }
      calc_tri_area_(&temp, &tri_node[0]);
      total_area_a += temp;
    }

  for(size_t i = 0; i < tri_b_connect.size(); ++i){
      if(!NM){
          itr_matrix<const double *> x0(dim, node_num, x);
          tri_node = x0(colon(), tri_b_connect[i]);
        }else{
          get_cell_node(x, tri_b_connect[i], *NM, tri_node);
        }

      calc_tri_area_(&temp, &tri_node[0]);
      total_area_b += temp;
    }

  return (area_a / total_area_b * area_b + area_b / total_area_a * area_a)/total_area;
}


class adj_normal_subtract_new_weight : public hj_func
{
public:
  adj_normal_subtract_new_weight(
      const zjucad::matrix::matrix<double> &node,
      const zjucad::matrix::matrix<size_t> &tri_a,
      const zjucad::matrix::matrix<size_t> &tri_b,
      const std::vector<zjucad::matrix::matrix<size_t> > & tri_a_connect,
      const std::vector<zjucad::matrix::matrix<size_t> > & tri_b_connect,
      const double total_area,
      const double &weight,
      const node_mapping * NM = 0)
    :node_num_(node.size(2)), weight_(weight), tri_a_(tri_a), tri_b_(tri_b),
      tri_a_connect_(tri_a_connect),tri_b_connect_(tri_b_connect),total_area_(total_area),
      node_mapping_(NM){
    std::vector<size_t> sort_points;
    for(size_t i = 0; i < tri_a.size(); ++i)
      for(size_t j = 0; j < tri_b.size(); ++j){
          if(tri_a[i] == tri_b[j]) {
              sort_points.push_back(tri_a[i]);
              break;
            }
        }
    assert(sort_points.size() == 2);
    const size_t shared_points =
        std::accumulate(sort_points.begin(), sort_points.end(), static_cast<size_t>(0));
    sort_points.push_back(
          std::accumulate(tri_a.begin(), tri_a.end(), static_cast<size_t>(0)) - shared_points);
    sort_points.push_back(
          std::accumulate(tri_b.begin(), tri_b.end(), static_cast<size_t>(0)) - shared_points);

    tp_.resize(4,1);
    std::copy(sort_points.begin(), sort_points.end(), tp_.begin());

    assemble();
  }
  virtual size_t dim_of_x(void) const {
    if(!node_mapping_)
      return node_num_*3;
    else{
        return node_mapping_->ZT.size(1);
      }
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }

  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {

    zjucad::matrix::matrix<double> tp_node;
    if(!node_mapping_){
        const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
        tp_node = T(zjucad::matrix::colon(), tp_);
      }else{
        tp_node.resize(3, tp_.size());
        get_cell_node(x, tp_, *node_mapping_, tp_node);
      }

    calc_adj_norm_sub_(f, &tp_node[0]);

    const double w = weight_ * sqrt(
          edge_area(tri_a_, tri_b_, tri_a_connect_,
                    tri_b_connect_, total_area_, &x[0], 3, node_num_, node_mapping_));
    for(int i = 0; i < 3; ++i) {
        f[i] *= w;
        jtf::math::erase_nan_inf(f[i]);
      }
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    zjucad::matrix::matrix<double> tp_node;

    if(!node_mapping_){
        const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
        tp_node = T(zjucad::matrix::colon(), tp_);
      }else{
        tp_node.resize(3, tp_.size());
        get_cell_node(x, tp_, *node_mapping_, tp_node);
      }

    zjucad::matrix::matrix<double> jac_tp(3,12);

    calc_adj_norm_sub_jac_(&jac_tp[0], &tp_node[0]);

    const double w = weight_ * sqrt(
          edge_area(tri_a_, tri_b_, tri_a_connect_, tri_b_connect_,
                    total_area_, x, 3, node_num_, node_mapping_));

    for(int fi = 0; fi < 3; ++fi) {
        if(!node_mapping_){
            ptr[fi+1] = ptr[fi] + 12;
            for(size_t pi = 0; pi < 4; ++pi) {
                for(size_t di = 0; di < 3; ++di){
                    idx[ptr[fi] + pi * 3 + di] = tp_[pi]*3+di;
                    val[ptr[fi] + pi * 3 + di] = jac_tp(fi, pi * 3 + di)*w;
                    jtf::math::erase_nan_inf(val[ptr[fi] + pi * 3 + di]);
                  }
              }
          }else{
            ptr[fi+1] = ptr[fi] + real_v_to_orig_v_coeff_.size();
            size_t vi = 0;
            for(const auto & one_real_v : real_v_to_orig_v_coeff_){
                idx[ptr[fi] + vi] = one_real_v.first;
                const std::vector<std::pair<size_t,double> > & to_orig_v = one_real_v.second;
                val[ptr[fi] + vi] = 0;
                for(const auto & one_exp : to_orig_v){
                    val[ptr[fi] + vi] += jac_tp(fi, one_exp.first)*w * one_exp.second;
                  }
                jtf::math::erase_nan_inf(val[ptr[fi]+vi]);
                ++vi;
              }
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    if(!node_mapping_)
      return 36 * dim_of_f();
    else
      return dim_of_f() * real_v_to_orig_v_coeff_.size(); //TODO: here jac nnz shoud be recalculated
  }

  void assemble(){
    if(node_mapping_){
        get_node_mapping(tp_, node_mapping_->ZT, node_mapping_of_each_variable_);
        for(size_t ni = 0; ni < node_mapping_of_each_variable_.size(); ++ni){
            const std::vector<std::pair<size_t,double> > & one_eqn =
                node_mapping_of_each_variable_[ni];
            for(const auto & one_exp : one_eqn){
                real_v_to_orig_v_coeff_[one_exp.first].push_back(
                      std::make_pair(ni, one_exp.second));
              }
          }
      }
  }
private:
  const size_t node_num_;
  zjucad::matrix::matrix<size_t> tp_;
  const double weight_;
  const zjucad::matrix::matrix<size_t> tri_a_;
  const zjucad::matrix::matrix<size_t> tri_b_;
  const std::vector<zjucad::matrix::matrix<size_t> >  tri_a_connect_;
  const std::vector<zjucad::matrix::matrix<size_t> >  tri_b_connect_;
  const double total_area_;
  const node_mapping * node_mapping_;
  std::vector<std::vector<std::pair<size_t,double> > > node_mapping_of_each_variable_;
  std::map<size_t, std::vector<std::pair<size_t,double> > > real_v_to_orig_v_coeff_;
};


template <typename T1, typename T2>
const hj_func*
build_normal_diff_constraint(
    const zjucad::matrix::matrix_expression<T1> & faces,
    const zjucad::matrix::matrix_expression<T2> & node,
    const double w,
    const node_mapping * node_mapping = 0)
{
  using namespace std;
  using namespace zjucad::matrix;

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(faces()));
  if(!ea.get()){
      std::cerr << "# [error] can not build edge2cell_adjacent." << std::endl;
      return nullptr;
    }

  vector<vector<matrix<size_t> > > face_connect(faces().size(2));

  for(size_t fi = 0; fi < faces().size(2); ++fi){
      for(size_t pi = 0; pi < faces().size(1); ++pi){
          const size_t edge_idx = ea->get_edge_idx(faces()(pi,fi), faces()((pi+1)%faces().size(1), fi));
          if(edge_idx == -1){
              cerr << "# [error] strange surface triangle contains non-surface edges." << endl;
              return nullptr;
            }
          const pair<size_t,size_t> & face_pair = ea->edge2cell_[edge_idx];
          if(ea->is_boundary_edge(face_pair)) continue;
          face_connect[fi].push_back(faces()(colon(), (face_pair.first != fi?face_pair.first:face_pair.second)));
        }
    }

  double total_area = 0.0;

  for(size_t fi = 0; fi < faces().size(2); ++fi){
      total_area += jtf::mesh::cal_face_area(faces()(colon(),fi), node);
    }

  shared_ptr<vector<hj_func_cons_ptr> > func_vec(new vector<hj_func_cons_ptr>);

  vector<matrix<size_t> > two_faces(2);
  for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
      const pair<size_t,size_t> & two_cell = ea->edge2cell_[ei];
      if(ea->is_boundary_edge(two_cell)) continue;
      const pair<size_t,size_t> & one_edge = ea->edges_[ei];

      for(size_t fi = 0; fi < 2; ++fi){
          const size_t face_idx = (fi==0?two_cell.first:two_cell.second);
          two_faces[fi].resize(3,1);
          two_faces[fi][0] = one_edge.first;
          two_faces[fi][1] = one_edge.second;
          two_faces[fi][2] =
              std::accumulate(faces()(colon(),face_idx).begin(),
                              faces()(colon(), face_idx).end(),
                              static_cast<size_t>(0))
              - one_edge.first - one_edge.second;
        }
      func_vec->push_back(
            hj_func_cons_ptr(
              new adj_normal_subtract_new_weight(
                node, two_faces.front(), two_faces.back(),
                face_connect[two_cell.first],
              face_connect[two_cell.second], total_area, w , node_mapping)));
    }

  return hj::function::new_catenated_function<double,int32_t>(func_vec);
}

#endif // NORMAL_DIFF_H
