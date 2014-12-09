#include "polycube_surface_func.h"
#include <zjucad/matrix/io.h>
#include <fstream>

#include "util.h"
#include "../common/vtk.h"
#include "../tetmesh/tetmesh.h"
#include "../common/util.h"
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/mesh.h>
#include "../mesh_func/tri-area-normal-func.h"

#include "tet_func.h"

using namespace std;
using namespace zjucad::matrix;
using namespace hj::function;

extern double L1_sqrt_eps;
hj::function::function_t<double, int32_t> *
build_fix_zero_node_func(const matrix<double> &node,
                         const matrix<double> &zero_pos,
                         const double weight)
{
  return new fix_zero_node_func(node, zero_pos, weight);
}

hj::function::function_t<double, int32_t> *
build_surface_arap_func(const matrix<double> &node,
                        const matrix<size_t> &faces)
{
  boost::shared_ptr<vector<boost::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double, int32_t> > >(faces.size(2)));

  matrix<double> face_area = zeros<double>(faces.size(2),1);
  matrix<double> face_node;
  for(size_t fi = 0; fi < faces.size(2); ++fi){
      face_node = node(colon(), faces(colon(),fi));
      calc_tri_area_(&face_area[fi], &face_node[0]);
    }
  face_area /= std::accumulate(face_area.begin(), face_area.end(),0.0);

  for(size_t fi = 0; fi < faces.size(2); ++fi){
      (*funcs)[fi].reset(new surface_arap_func(node, faces(colon(),fi), face_area[fi]));
    }

  return new_catenated_function<double,int32_t>(funcs);
}

hj::function::function_t<double, int32_t> *
build_feature_line_func(
    const matrix<double> &node,
    const matrix<size_t> &outside_face,
    const jtf::mesh::edge2cell_adjacent &ea,
    const boost::unordered_set<pair<size_t,size_t> > &feature_edges,
    const double feature_line_w)
{
  if(feature_edges.empty()) return 0;

  boost::shared_ptr<vector<boost::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double, int32_t> > >);

  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator scit =
      feature_edges.begin(); scit != feature_edges.end(); ++scit){
      const size_t & edge_idx = ea.get_edge_idx(scit->first, scit->second);
      if(edge_idx == -1) {
          cerr << "# [error] can not find edge " << scit->first << " "
               << scit->second << endl;
          continue;
        }
      const pair<size_t,size_t> & tri_pair = ea.edge2cell_[edge_idx];
      boost::shared_ptr<normal_dot_func> ndf(
            new normal_dot_func(node, outside_face(colon(), tri_pair.first),
                                outside_face(colon(), tri_pair.second),
                                sqrt(feature_line_w)));
      funcs->push_back(ndf);
    }

  return new_catenated_function<double, int32_t>(funcs);
}

hj::function::function_t<double, int32_t> *
build_group_equation_func(
    const matrix<double> &node,
    const vector<boost::unordered_set<size_t> > &fnode_group,
    const double group_weight)
{
  boost::shared_ptr<vector<boost::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double, int32_t> > >());

  for(size_t gi = 0; gi < fnode_group.size(); ++gi){
      const boost::unordered_set<size_t> & one_group = fnode_group[gi];

      for(boost::unordered_set<size_t>::const_iterator cit = one_group.begin();
          cit != one_group.end(); ++cit){
          boost::unordered_set<size_t>::const_iterator ccit = cit;
          ++ccit;
          if(ccit == one_group.end()) continue;
          funcs->push_back(boost::shared_ptr<function_t<double, int32_t> >
                           (new equation_func(node, *cit, *ccit, sqrt(group_weight))));
        }
    }
  if(funcs->empty()) return 0;
  return new_catenated_function<double, int32_t>(funcs);
}

hj::function::function_t<double, int32_t> *
build_patch_normal_fix_func(const matrix<double> &node,
                            const vector<matrix<double> > & patch_normal,
                            const vector<matrix<size_t> > &surface_patches,
                            const double normal_fix_w,
                            const vector<vector<size_t> > * select_ptr)
{
  if(patch_normal.empty()) return 0;

  double total_area = 0.0;
  vector<matrix<double> > surface_patch_area(surface_patches.size());
  matrix<double> tri_node;
  for(size_t pi = 0; pi < surface_patches.size(); ++pi){
      const matrix<size_t> & one_patch = surface_patches[pi];
      surface_patch_area[pi].resize(one_patch.size(2));
      for(size_t fi = 0; fi < one_patch.size(2); ++fi){
          tri_node = node(colon(), one_patch(colon(),fi));
          calc_tri_area_(&surface_patch_area[pi][fi], &tri_node[0]);
          total_area += surface_patch_area[pi][fi];
        }
    }

  boost::shared_ptr<vector<boost::shared_ptr<function_t<double,int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double,int32_t> > >);
  for(size_t pi = 0; pi < surface_patches.size(); ++pi){
      const matrix<size_t> & one_patch = surface_patches[pi];
      for(size_t fi = 0; fi < one_patch.size(2); ++fi){

          boost::shared_ptr<normal_fix_func> nff(
                new normal_fix_func(
                  node, patch_normal[pi], one_patch(colon(),fi),
                  sqrt(normal_fix_w * surface_patch_area[pi][fi]/total_area)));
#if 0
          {
            const double err = jac_err(*nff, &node[0]);
            if(err > 1e-5) {
                cerr << "large error in jac: " << 1 << " " << err << endl;
              }
          }
#endif
          funcs->push_back(nff);
        }
    }

  return new_catenated_function<double, int32_t>(funcs);
}

class adj_normal_subtract : public hj::function::function_t<double,int32_t>
{
public:
  adj_normal_subtract(const zjucad::matrix::matrix<double> &node,
                      const zjucad::matrix::matrix<size_t> &tri_a,
                      const zjucad::matrix::matrix<size_t> &tri_b,
                      const double &weight)
    :node_num_(node.size(2)), weight_(weight) {
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
  }
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tp_node = T(zjucad::matrix::colon(), tp_);
    calc_adj_norm_sub_(f, &tp_node[0]);
    for(int i = 0; i < 3; ++i) {
        f[i] *= weight_;
        jtf::math::erase_nan_inf(f[i]);
      }
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tp_node = T(zjucad::matrix::colon(), tp_);
    zjucad::matrix::matrix<double> jac_tp(3,12);

    calc_adj_norm_sub_jac_(&jac_tp[0], &tp_node[0]);

    for(int fi = 0; fi < 3; ++fi) {
        ptr[fi+1] = ptr[fi] + 12;
        for(size_t pi = 0; pi < 4; ++pi) {
            for(size_t di = 0; di < 3; ++di){
                idx[ptr[fi] + pi * 3 + di] = tp_[pi] * 3 + di;
                val[ptr[fi] + pi * 3 + di] = jac_tp(fi, pi * 3 + di)*weight_;
                jtf::math::erase_nan_inf(val[ptr[fi] + pi * 3 + di]);
              }
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 36 * dim_of_f();
  }
private:
  const size_t node_num_;
  zjucad::matrix::matrix<size_t> tp_;
  const double weight_;
};

static double edge_area(const matrix<size_t> & tri_a,
                        const matrix<size_t> & tri_b,
                        const vector<matrix<size_t> > & tri_a_connect,
                        const vector<matrix<size_t> > & tri_b_connect,
                        const double total_area,
                        const double * x,
                        const size_t dim, const size_t node_num)
{    
  itr_matrix<const double *> x0(dim, node_num, x);
  matrix<double> tri_node;
  double area_a, area_b;
  tri_node = x0(colon(), tri_a);
  calc_tri_area_(&area_a, &tri_node[0]);

  tri_node = x0(colon(), tri_b);
  calc_tri_area_(&area_b, &tri_node[0]);

  double total_area_a = 0, total_area_b = 0, temp;
  for(size_t i = 0; i < tri_a_connect.size(); ++i){
      tri_node = x0(colon(), tri_a_connect[i]);
      calc_tri_area_(&temp, &tri_node[0]);
      total_area_a += temp;
    }

  for(size_t i = 0; i < tri_b_connect.size(); ++i){
      tri_node = x0(colon(), tri_b_connect[i]);
      calc_tri_area_(&temp, &tri_node[0]);
      total_area_b += temp;
    }

  return (area_a / total_area_b * area_b + area_b / total_area_a * area_a)/total_area;
}

class adj_normal_subtract_new_weight : public hj::function::function_t<double,int32_t>
{
public:
  adj_normal_subtract_new_weight(
      const zjucad::matrix::matrix<double> &node,
      const zjucad::matrix::matrix<size_t> &tri_a,
      const zjucad::matrix::matrix<size_t> &tri_b,
      const std::vector<zjucad::matrix::matrix<size_t> > & tri_a_connect,
      const std::vector<zjucad::matrix::matrix<size_t> > & tri_b_connect,
      const zjucad::matrix::matrix<size_t> & faces,
      const double total_area,
      const double &weight)
    :node_num_(node.size(2)), weight_(weight), tri_a_(tri_a), tri_b_(tri_b),
      tri_a_connect_(tri_a_connect),tri_b_connect_(tri_b_connect), faces_(faces),total_area_(total_area) {
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

  }
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }

  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {

    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tp_node = T(zjucad::matrix::colon(), tp_);
    calc_adj_norm_sub_(f, &tp_node[0]);

    const double w = weight_ * sqrt(edge_area(tri_a_, tri_b_, tri_a_connect_, tri_b_connect_, total_area_, x, 3, node_num_));
    for(int i = 0; i < 3; ++i) {
        f[i] *= w;
        jtf::math::erase_nan_inf(f[i]);
      }
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tp_node = T(zjucad::matrix::colon(), tp_);
    zjucad::matrix::matrix<double> jac_tp(3,12);

    //const double area = ac_->get_area(x, 3, node_num_, faces_);

    calc_adj_norm_sub_jac_(&jac_tp[0], &tp_node[0]);

    const double w = weight_ * sqrt(edge_area(tri_a_, tri_b_, tri_a_connect_, tri_b_connect_, total_area_, x, 3, node_num_));
    for(int fi = 0; fi < 3; ++fi) {
        ptr[fi+1] = ptr[fi] + 12;
        for(size_t pi = 0; pi < 4; ++pi) {
            for(size_t di = 0; di < 3; ++di){
                idx[ptr[fi] + pi * 3 + di] = tp_[pi] * 3 + di;
                val[ptr[fi] + pi * 3 + di] = jac_tp(fi, pi * 3 + di)*w;
                jtf::math::erase_nan_inf(val[ptr[fi] + pi * 3 + di]);
              }
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 36 * dim_of_f();
  }
private:
  const size_t node_num_;
  zjucad::matrix::matrix<size_t> tp_;
  const double weight_;
  const zjucad::matrix::matrix<size_t> tri_a_;
  const zjucad::matrix::matrix<size_t> tri_b_;
  const std::vector<zjucad::matrix::matrix<size_t> >  tri_a_connect_;
  const std::vector<zjucad::matrix::matrix<size_t> >  tri_b_connect_;
  const zjucad::matrix::matrix<size_t> & faces_;
  const double total_area_;
  //std::unique_ptr<area> ac_;
};


class adj_normal_func : public hj::function::function_t<double,int32_t>
{
public:
  adj_normal_func(const zjucad::matrix::matrix<double> &node,
                  const zjucad::matrix::matrix<size_t> &tri_a,
                  const zjucad::matrix::matrix<size_t> &tri_b)
    :node_num_(node.size(2)) {
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
    //    if(f > -0.8)
    //      min_dot_ = -0.8;
  }
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 1;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tp_node = T(zjucad::matrix::colon(), tp_);

    calc_tri_normal_dot_(f, &tp_node[0]);
    double dot_val = *f;
    // if(*f < -1) { // at nearly degenerated case, it may be true.
    //   cerr << "# normal dot < -1: " << tp_node;
    // }

    *f = *f-min_dot_;

    jtf::math::erase_nan_inf(*f);
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tp_node = T(zjucad::matrix::colon(), tp_);
    zjucad::matrix::matrix<double> jac_tp(1,12);

    calc_tri_normal_dot_jac_(&jac_tp[0], &tp_node[0]);

    ptr[1] = ptr[0] + 12;
    for(size_t pi = 0; pi < 4; ++pi)
      for(size_t di = 0; di < 3; ++di){
          idx[ptr[0] + pi * 3 + di] = tp_[pi] * 3 + di;
          val[ptr[0] + pi * 3 + di] = jac_tp[pi * 3 + di];
          jtf::math::erase_nan_inf(val[ptr[0] + pi * 3 + di]);
        }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 12 * dim_of_f();
  }
private:
  const size_t node_num_;
  zjucad::matrix::matrix<size_t> tp_;
  double  min_dot_;
};

class area_preserve_func : public hj::function::function_t<double,int32_t>
{
public:
  area_preserve_func(const size_t node_num, const matrixst & faces,
                     const double total_area, const double weight,
                     const double percent = 0.9)
    :node_num_(node_num), faces_(faces), orig_area_(total_area), w_(weight),
      per_(percent){}
  virtual size_t dim_of_x()const{
    return 3 * node_num_;
  }
  virtual size_t dim_of_f()const{
    return 1;
  }
  virtual int val(const value_type *x, value_type *f, func_ctx *ctx) const {
    *f = -1.0 * per_ * orig_area_;
    itr_matrix<const value_type*> x0(3, node_num_, x);
    matrix<double> tri_node ;
    double area = 0;
    for(size_t fi = 0; fi < faces_.size(2); ++fi){
        tri_node = x0(colon(), faces_(colon(),fi));
        calc_tri_area_(&area, &tri_node[0]);
        *f += area;
      }
    *f *= w_;
    return 0;
  }
  virtual int jac(const value_type *x, value_type *val, int_type *ptr,
                  int_type *idx, func_ctx *ctx) const{
    ptr[1] = ptr[0] + 3*node_num_;
    for(size_t i = 0; i < 3 * node_num_; ++i){
        idx[ptr[0] + i] = i;
        val[ptr[0] + i] = 0;
      }
    itr_matrix<const value_type*> x0(3, node_num_, x);
    matrix<double> tri_node, jac(1,9) ;
    for(size_t fi = 0; fi < faces_.size(2); ++fi){
        tri_node = x0(colon(), faces_(colon(),fi));
        calc_tri_area_jac_(&jac[0] , &tri_node[0]);
        for(size_t pi = 0; pi < faces_.size(1); ++pi){
            for(size_t di = 0; di < 3; ++di){
                val[ptr[0] + faces_(pi,fi) * 3 + di] += jac[3*pi+di] * w_;
              }
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void)const{
    return 3 *node_num_;
  }
private:
  const size_t node_num_;
  const matrixst & faces_;
  const double orig_area_;
  const double w_;
  const double per_;
};

/// @brief: split the triangle area according to its neighboring areas
void calc_edge_area(const matrix<size_t> &faces,
                    const jtf::mesh::edge2cell_adjacent &ea,
                    const matrix<double> &face_area,
                    matrix<double> &edge_area)
{
  edge_area = zeros<double>(ea.edges_.size(), 1);
  for(size_t fi = 0; fi < faces.size(2); ++fi) {
      matrix<size_t> nb_f_idx(3), nb_e_idx(3);
      double total_nb_area = 0;
      for(size_t ei = 0; ei < 3; ++ei) {
          nb_e_idx[ei] = ea.get_edge_idx(faces(ei, fi), faces((ei+1)%3, fi));
          const pair<size_t, size_t> &fp = ea.edge2cell_[nb_e_idx[ei]];
          nb_f_idx[ei] = fp.first+fp.second-fi;
          if(nb_f_idx[ei] >= faces.size(2)) {//boundary edge
              //        cerr << "boundary edge? " << fi << " " << fp.first << " " << fp.second
              //             << " " << ei << endl;
              //exit(0);
              continue;
            }
          total_nb_area += face_area[nb_f_idx[ei]];
        }
      for(size_t ei = 0; ei < 3; ++ei)
        if(nb_f_idx[ei] < faces.size(2))
          edge_area[nb_e_idx[ei]] += face_area[fi]*face_area[nb_f_idx[ei]]/total_nb_area;
    }
}

hj::function::function_t<double, int32_t> *
build_adj_normal_sub_func(
    const matrix<double> & node,
    const matrix<size_t> & faces,
    double weight)
{
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(faces));
  if(!ea.get()){
      cerr  << "# [error] can not build edge2triangle adjacent." << endl;
      return 0;
    }

  size_t e_num_total = 0, e_num;
  for(size_t ei = 0; ei < ea->edge2cell_.size(); ++ei) {
      const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
      if(ea->is_boundary_edge(tri_pair)) continue;
      ++e_num_total;
    }

  matrixd face_area = zeros<double>(faces.size(2),1);
  for(size_t fi = 0; fi < face_area.size(); ++fi){
      face_area[fi] = jtf::mesh::cal_face_area(faces(colon(),fi), node);
    }

  //  matrix<double> edge_area;
  //  calc_edge_area(faces, *ea, face_area, edge_area);
  //  edge_area /= std::accumulate(edge_area.begin(), edge_area.end(), 0.0);

  const double total_area =
      std::accumulate(face_area.begin(), face_area.end(), 0.0);


  boost::shared_ptr<vector<boost::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double, int32_t> > >(e_num_total));

  vector<vector<matrix<size_t> > > face_connect(faces.size(2));

  for(size_t fi = 0; fi < faces.size(2); ++fi){
      for(size_t pi = 0; pi < faces.size(1); ++pi){
          const size_t edge_idx = ea->get_edge_idx(faces(pi,fi), faces((pi+1)%faces.size(1), fi));
          if(edge_idx == -1){
              cerr << "# [error] strange surface triangle contains non-surface edges." << endl;
              return 0;
            }
          const pair<size_t,size_t> & face_pair = ea->edge2cell_[edge_idx];
          if(ea->is_boundary_edge(face_pair)) continue;
          face_connect[fi].push_back(faces(colon(), (face_pair.first != fi?face_pair.first:face_pair.second)));
        }
    }

  size_t edge_num = 0;
  for(size_t ei = 0; ei < ea->edge2cell_.size(); ++ei) {
      const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
      if(ea->is_boundary_edge(tri_pair)) continue;
      // cerr << "#  " << sqrt(weight*(face_area[tri_pair.first]+face_area[tri_pair.second])/(3*total_area)) << endl;
      boost::shared_ptr<function_t<double, int32_t> > adj_normal_sub(
            new adj_normal_subtract_new_weight
            (node, faces(colon(),tri_pair.first),
             faces(colon(), tri_pair.second),
             face_connect[tri_pair.first],
            face_connect[tri_pair.second],
          faces,
          total_area,

          sqrt(weight)));

      (*funcs)[edge_num] = adj_normal_sub;
      ++edge_num;
    }
  return new_catenated_function<double, int32_t>(funcs);
}

hj::function::function_t<double,int32_t> *
build_area_preserve_func(const size_t node_num,
                         const matrix<size_t> & faces,
                         const double total_area,
                         const double percent,
                         const double weight)
{
  return new area_preserve_func(node_num, faces, total_area, weight, percent);
}

hj::function::function_t<double, int32_t> *
build_adj_normal_func(
    const matrix<double> & node,
    const matrix<size_t> & faces,
    matrix<double> &weight)
{
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(faces));
  if(!ea.get()){
      cerr  << "# [error] can not build edge2triangle adjacent." << endl;
      return 0;
    }

  size_t e_num_total = 0, e_num;
  for(size_t ei = 0; ei < ea->edge2cell_.size(); ++ei) {
      const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
      if(ea->is_boundary_edge(tri_pair)) continue;
      ++e_num_total;
    }

  matrixd face_area = zeros<double>(faces.size(2),1);
  for(size_t fi = 0; fi < face_area.size(); ++fi){
      face_area[fi] = jtf::mesh::cal_face_area(faces(colon(),fi), node);
    }

  matrix<double> edge_area;
  calc_edge_area(faces, *ea, face_area, edge_area);

  const double total_area =
      std::accumulate(face_area.begin(), face_area.end(), 0.0);
  edge_area /= total_area;

  boost::shared_ptr<vector<boost::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double, int32_t> > >(e_num_total));


  //cout << "adj_normal_w: " << adj_normal_w << endl;

  weight.resize(e_num_total);
  size_t edge_num = 0;
  for(size_t ei = 0; ei < ea->edge2cell_.size(); ++ei) {
      const pair<size_t,size_t> & one_edge = ea->edges_[ei];
      const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
      if(ea->is_boundary_edge(tri_pair)) continue;


      double diffusd_vale = 1.0;
      boost::shared_ptr<function_t<double, int32_t> > adj_normal
          (new adj_normal_func
           (node, faces(colon(),tri_pair.first), faces(colon(), tri_pair.second)));
#if 0 // validate jac
      {
        const double err = jac_err(*adj_normal, &node[0]);
        if(err > 1e-4) {
            cerr << "large error in jac: " << ei << " " << err << endl;
          }
      }
#endif
      (*funcs)[edge_num] = adj_normal;
      //weight[edge_num] = (face_area[tri_pair.first] + face_area[tri_pair.second])/(3*total_area);
      weight[edge_num] = edge_area[ei];
      ++edge_num;
    }
  return new_catenated_function<double, int32_t>(funcs);
}

static int calc_patch_area(const matrix<size_t> &faces,
                           const matrix<double> &node,
                           vector<double> &area)
{
  area.resize(faces.size(2));
  matrix<double> tri_node;
  for(size_t fi = 0; fi < faces.size(2); ++fi){
      tri_node = node(colon(), faces(colon(),fi));
      calc_tri_area_(&area[fi], &tri_node[0]);
    }
  return 0;
}

hj::function::function_t<double, int32_t> *
build_surface_patch_align_func(
    const matrix<double> &node,
    const vector<matrix<size_t> > &surface_patch,
    const double patch_w)
{
  if(surface_patch.empty()) return 0;

  boost::shared_ptr<vector<boost::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double, int32_t> > >);

  vector<double> area;
  for(size_t pi = 0; pi < surface_patch.size(); ++pi){
      const matrix<size_t> & patch = surface_patch[pi];
      calc_patch_area(patch, node, area);
      const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
      unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
            jtf::mesh::edge2cell_adjacent::create(patch));
      for(size_t ei = 0; ei < ea->edge2cell_.size(); ++ei){
          const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
          if(ea->is_boundary_edge(tri_pair)) continue;

          boost::shared_ptr<function_t<double,int32_t> > pf(
                new normal_dot_func(
                  node, patch(colon(), tri_pair.first),
                  patch(colon(), tri_pair.second),
                  sqrt(patch_w * (area[tri_pair.first] + area[tri_pair.first])/(3*total_area) ),1.0));

          funcs->push_back(pf);
        }
    }

  return new_catenated_function<double, int32_t>(funcs);
}

//! dot(R(xyz, :), un)
class rot_unit_normal : public jtf::function::functionN1_t<double,int32_t>
{
public:
  rot_unit_normal(const zjucad::matrix::matrix<double> &un,
                  size_t xyz)
    :un_(un), xyz_(xyz) {
  }
  virtual size_t dim(void) const {
    return 3*3;
  }
  virtual int val(const double *x, double &v) {
    for(size_t i = 0; i < 3; ++i)
      v += x[xyz_+i*3]*un_[i];
    return 0;
  }
  virtual int gra(const double *x, double *g) { assert(0); }
  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx) {
    if(g == 0) {
        nnz = 3;
        return 0;
      }
    if(nnz != 3) {
        cerr << "strange." << endl;
      }
    for(size_t i = 0; i < 3; ++i) {
        idx[i] = xyz_+i*3;
        g[i] = un_[i];
      }
    return 0;
  }
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha = 1) {
    format = 1;
    if(h == 0 && ptr == 0 && idx == 0) { // query nnz
        nnz = 0;
        return 0;
      }
    if(h == 0 && ptr != 0 && idx != 0) { // query pattern
        return 0;
      }
    if(h != 0 && ptr != 0 && idx != 0) { // accumulate
        return 0;
      }
    return __LINE__;
  }
  virtual int hes_block(const double *x, double *h, double alpha = 1) {
    return 0;
  }
private:
  const size_t xyz_;
  zjucad::matrix::matrix<double> un_;
};

jtf::function::functionN1_t<double,int32_t> *
build_polycube_rot_func2(const zjucad::matrix::matrix<double> &node,
                         const zjucad::matrix::matrix<size_t> &faces,
                         const zjucad::matrix::matrix<double> &areas)
{
  const matrixd &surface_area = areas;

  const double total_area =
      std::accumulate(surface_area.begin(), surface_area.end(), 0.0);
  cerr << "total_area in build_smooth_L1_normal: " << total_area << endl;

  const size_t fn = faces.size(2);
  vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > sum(fn*3+1);
  for(size_t fi = 0; fi < fn; ++fi) {
      const matrix<double> tri = node(colon(), faces(colon(), fi));
      matrix<double> unit_normal = cross(tri(colon(), 1)-tri(colon(), 0), tri(colon(), 2)-tri(colon(), 0));
      assert(fabs(norm(unit_normal)-surface_area[fi]*2) < 1e-8);
      if(surface_area[fi] > 1e-8)
        unit_normal /= surface_area[fi]*2;
      else
        unit_normal(colon()) = 0;

      for(size_t i = 0; i < 3; ++i) {
          shared_ptr<jtf::function::functionN1_t<double,int32_t> > func
              (new rot_unit_normal(unit_normal, i));
          sum[fi*3+i].reset(new smooth_L1(func, surface_area[fi]/total_area, L1_sqrt_eps));
        }
    }
  {
    static shared_ptr<so3> func(new so3());
    sum[fn*3].reset(jtf::function::least_square_warpper(func));
    //new least_square_Gauss_Newton(*func));
  }
  unique_ptr<jtf::function::functionN1_t<double,int32_t> > rtn(new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD>(sum));
  //  cerr << "gra error: " << hj_func_opt::gra_err(*rtn, &node[0]) << endl;
  //  cerr << "hes error: " << hj_func_opt::hes_err(*rtn, &node[0]) << endl;
  return rtn.release();
}

jtf::function::functionN1_t<double,int32_t> *
build_smooth_L1_area_normal(const matrix<double> &node, const matrix<size_t> &tri,
                            const matrix<double> &areas,
                            double normal_align_w)
{
  const matrixd &surface_area = areas;

  const double total_area = std::accumulate(surface_area.begin(), surface_area.end(), 0.0);
  cerr << "total_area in build_smooth_L1_normal: " << total_area << endl;

  const size_t fn = tri.size(2);
  vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > sum(fn*3);
  const double weight = normal_align_w / total_area;
  for(size_t fi = 0; fi < fn; ++fi) {
      for(size_t i = 0; i < 3; ++i) {
          shared_ptr<jtf::function::functionN1_t<double,int32_t> > normal_func(new face_area_normal(node, tri(colon(), fi), i, 1));
          sum[fi*3+i].reset(new smooth_L1(normal_func, weight, L1_sqrt_eps*surface_area[fi]));
        }
    }
  unique_ptr<jtf::function::functionN1_t<double,int32_t> > rtn(
        new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD>(sum));
  return rtn.release();
}

hj::function::function_t<double, int32_t> *
build_polycube_function(
    const matrix<double> &orig_node,
    const matrix<double> &deform_node,
    const matrix<double> &zero_pos,
    const matrix<size_t> &tet,
    const matrix<size_t> & faces,
    hj::function::function_t<double,int32_t> *&diagnose,
    const boost::unordered_set<pair<size_t,size_t> > *feature_line,
    const vector<matrix<size_t> > *surface_patch,
    const double adj_normal_w, const double fl_w, const double patch_w)
{
  assert(tet.size(1) == 4 && faces.size(1) == 3);
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(faces));
  if(!ea.get()){
      cerr  << "# [error] can not build edge2triangle adjacent." << endl;
      return 0;
    }

  //////////////////////////////////////////////////////////////////
  ////////////////  add function term  /////////////////////////////

  shared_ptr<function_t<double, int32_t> > tetmesh_arap(
        build_tetmesh_arap_func(orig_node, tet, 1.0)); // arap_w = 1.0
  diagnose = tetmesh_arap.get();

  shared_ptr<function_t<double, int32_t> > fix_pos(
        new fix_zero_node_func(orig_node, zero_pos, 1.0)); // fix_zero pos = 1.0

  shared_ptr<function_t<double, int32_t> > adj_normal(
        build_adj_normal_sub_func(deform_node, faces, adj_normal_w));

  shared_ptr<function_t<double, int32_t> > feature_line_func;
  if(feature_line != nullptr)
    feature_line_func.reset(build_feature_line_func(orig_node, faces, *ea, *feature_line, fl_w));

  shared_ptr<function_t<double, int32_t> > surface_patch_algin_func;
  if(surface_patch != nullptr)
    surface_patch_algin_func.reset(build_surface_patch_align_func(orig_node, *surface_patch, patch_w));

  boost::shared_ptr<vector<shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<shared_ptr<function_t<double, int32_t> > >);

  if(tetmesh_arap.get()){
      funcs->push_back(tetmesh_arap);
      cerr << "# [info] add jtf::mesh::mesh arap function, weight " << 1 << endl;
    }

  if(fix_pos.get()){
      funcs->push_back(fix_pos);
      cerr << "# [info] add fix first node zero function, weight " << 1 << endl;
    }

  if(adj_normal.get() && adj_normal_w > 0){
      funcs->push_back(adj_normal);
      cerr << "# [info] add adj_normal function, weight " << adj_normal_w << endl;
    }

  if(feature_line!= nullptr && feature_line_func.get() && non_zero(fl_w)){
      funcs->push_back(feature_line_func);
      cerr << "# [info] add feature_line function, weight " << fl_w << endl;
    }

  if(surface_patch != nullptr && surface_patch_algin_func.get() && non_zero(patch_w)){
      funcs->push_back(surface_patch_algin_func);
      cerr << "# [info] add surface_patch_align function, weight " << patch_w << endl;
    }

  return new_catenated_function<double, int32_t>(funcs);
}
