#ifndef NORMALANTIFLIP_H
#define NORMALANTIFLIP_H

#include "../def.h"
#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/mesh.h>
#include "../../../mesh_func/tri-normal.h"
#include <hjlib/sparse/sparse.h>
#include "../../common/util.h"

class adj_normal_func : public hj_func
{
public:
  adj_normal_func(const zjucad::matrix::matrix<double> &node,
                  const zjucad::matrix::matrix<size_t> &tri_a,
                  const zjucad::matrix::matrix<size_t> &tri_b,
                  const node_mapping * node_mapping_ = 0)
    :node_num_(node.size(2)), node_mapping_(node_mapping_) {
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

    double f = 0;
    zjucad::matrix::matrix<double> tp_node = node(zjucad::matrix::colon(), tp_);
    calc_tri_normal_dot_(&f, &tp_node[0]);
    min_dot_ = -1;

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
  virtual size_t dim_of_x(void) const {
    if(!node_mapping_)
      return node_num_*3;
    else
      return node_mapping_->ZT.size(1);
  }
  virtual size_t dim_of_f(void) const {
    return 1;
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
    calc_tri_normal_dot_(f, &tp_node[0]);
    double dot_val = *f;
    *f = *f-min_dot_;

    jtf::math::erase_nan_inf(*f);
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
    zjucad::matrix::matrix<double> jac_tp(1,12);

    calc_tri_normal_dot_jac_(&jac_tp[0], &tp_node[0]);

    if(!node_mapping_){
        ptr[1] = ptr[0] + 12;
        for(size_t pi = 0; pi < 4; ++pi)
          for(size_t di = 0; di < 3; ++di){
              idx[ptr[0] + pi * 3 + di] = tp_[pi] * 3 + di;
              val[ptr[0] + pi * 3 + di] = jac_tp[pi * 3 + di];
              jtf::math::erase_nan_inf(val[ptr[0] + pi * 3 + di]);
            }
      }else{
        ptr[1] = ptr[0] + real_v_to_orig_v_coeff_.size();
        size_t vi = 0;
        for(const auto & one_real_v : real_v_to_orig_v_coeff_){
            idx[ptr[0] + vi] = one_real_v.first;
            val[ptr[0] + vi] = 0;
            const std::vector<std::pair<size_t,double> > & to_orig_v = one_real_v.second;
            for(const auto & one_exp : to_orig_v){
                val[ptr[0] + vi] += jac_tp[one_exp.first] * one_exp.second;
              }
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    if(!node_mapping_)
      return 12 * dim_of_f();
    else{
        return dim_of_f() * real_v_to_orig_v_coeff_.size();
      }
  }
private:
  const size_t node_num_;
  zjucad::matrix::matrix<size_t> tp_;
  double  min_dot_;
  const node_mapping * node_mapping_;
  std::vector<std::vector<std::pair<size_t,double> > > node_mapping_of_each_variable_;
  std::map<size_t, std::vector<std::pair<size_t,double> > > real_v_to_orig_v_coeff_;
};

/// @brief: split the triangle area according to its neighboring areas
template <typename T1, typename T2>
void calc_edge_area(const zjucad::matrix::matrix<T1> &faces,
                    const jtf::mesh::edge2cell_adjacent &ea,
                    const zjucad::matrix::matrix<T2> &face_area,
                    zjucad::matrix::matrix<T2> &edge_area)
{
  using namespace zjucad::matrix;

  edge_area = zeros<double>(ea.edges_.size(), 1);
  for(size_t fi = 0; fi < faces.size(2); ++fi) {
      matrix<size_t> nb_f_idx(3), nb_e_idx(3);
      double total_nb_area = 0;
      for(size_t ei = 0; ei < 3; ++ei) {
          nb_e_idx[ei] = ea.get_edge_idx(faces(ei, fi), faces((ei+1)%3, fi));
          const std::pair<size_t, size_t> &fp = ea.edge2cell_[nb_e_idx[ei]];
          nb_f_idx[ei] = fp.first+fp.second-fi;
          if(nb_f_idx[ei] >= faces.size(2)) {//boundary edge
              continue;
            }
          total_nb_area += face_area[nb_f_idx[ei]];
        }
      for(size_t ei = 0; ei < 3; ++ei)
        if(nb_f_idx[ei] < faces.size(2))
          edge_area[nb_e_idx[ei]] += face_area[fi]*face_area[nb_f_idx[ei]]/total_nb_area;
    }
}



template <typename T1, typename T2 >
hj::function::function_t<double, int32_t> *
build_adj_normal_func(const zjucad::matrix::matrix<T1> & node,
                      const zjucad::matrix::matrix<T2> & faces,
                      zjucad::matrix::matrix<T1> &weight,
                      const node_mapping * node_mapping_ = 0)
{
  using namespace std;
  using namespace zjucad::matrix;

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(faces));
  if(!ea.get()){
      cerr  << "# [error] can not build edge2triangle adjacent." << endl;
      return nullptr;
    }

  size_t e_num_total = 0, e_num;
  for(size_t ei = 0; ei < ea->edge2cell_.size(); ++ei) {
      const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
      if(ea->is_boundary_edge(tri_pair)) continue;
      ++e_num_total;
    }

  matrix<T1> face_area = zeros<T1>(faces.size(2),1);
  for(size_t fi = 0; fi < face_area.size(); ++fi){
      face_area[fi] = jtf::mesh::cal_face_area(faces(colon(),fi), node);
    }

  matrix<double> edge_area;
  calc_edge_area(faces, *ea, face_area, edge_area);

  const double total_area =
      std::accumulate(face_area.begin(), face_area.end(), 0.0);
  edge_area /= total_area;

  std::shared_ptr<vector<std::shared_ptr<hj::function::function_t<double, int32_t> > > >
      funcs(new vector<std::shared_ptr<hj::function::function_t<double, int32_t> > >(e_num_total));


  //cout << "adj_normal_w: " << adj_normal_w << endl;

  weight.resize(e_num_total);
  size_t edge_num = 0;
  for(size_t ei = 0; ei < ea->edge2cell_.size(); ++ei) {
      const pair<size_t,size_t> & one_edge = ea->edges_[ei];
      const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
      if(ea->is_boundary_edge(tri_pair)) continue;


      double diffusd_vale = 1.0;
      std::shared_ptr<hj::function::function_t<double, int32_t> > adj_normal
          (new adj_normal_func
           (node, faces(colon(),tri_pair.first), faces(colon(), tri_pair.second),
            node_mapping_));

      (*funcs)[edge_num] = adj_normal;
      //weight[edge_num] = (face_area[tri_pair.first] + face_area[tri_pair.second])/(3*total_area);
      weight[edge_num] = edge_area[ei];
      ++edge_num;
    }
  return hj::function::new_catenated_function<double, int32_t>(funcs);
}

#endif // NORMALANTIFLIP_H
