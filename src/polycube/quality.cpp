#include "quality.h"

#include <zjucad/matrix/itr_matrix.h>

#include "../mesh_func/tri-area-normal.h"

using namespace std;
using namespace zjucad::matrix;

double polycube_L1_area_quality(const double *x, size_t node_num, const zjucad::matrix::matrix<size_t> &faces)
{
  const zjucad::matrix::itr_matrix<const double *> T(3, node_num, x);
  double polycube_area = 0;
  matrix<double> n(3), tri(3, 3);
  double total_area = 0;
  for(size_t fi = 0; fi < faces.size(2); ++fi) {
    tri = T(colon(), faces(colon(), fi));
    calc_tri_area_normal_(&n[0], &tri[0]);
    total_area += norm(n);
    polycube_area += dot(fabs(n), ones<double>(3, 1));
  }
  return polycube_area/total_area-1;
}

int polycube_tet_quality(const double *x, size_t node_num,
                         const zjucad::matrix::matrix<size_t> &tet,
                         const zjucad::matrix::matrix<size_t> &faces,
                         std::map<std::string, double> &stat)
{
  if(x == 0) { // query
    stat["L1-area"] = -1;
    return 0;
  }
  if(stat.find("L1-area") != stat.end())
    stat["L1-area"] = polycube_L1_area_quality(x, node_num, faces);
  return 0;
}

int polycube_tet_quality(const double *x, size_t node_num,
                         const zjucad::matrix::matrix<size_t> &faces,
                         std::map<std::string, double> &stat)
{
  if(x == 0) { // query
    stat["L1-area"] = -1;
    return 0;
  }
  if(stat.find("L1-area") != stat.end())
    stat["L1-area"] = polycube_L1_area_quality(x, node_num, faces);
  return 0;
}

quality_checker::quality_checker(
    const zjucad::matrix::matrix<double> &node,
    const zjucad::matrix::matrix<size_t> &tet,
    const zjucad::matrix::matrix<size_t> &faces,
    std::vector<std::pair<jtf_func_type *, double> > &weighted_funcs)
  :node_(node), tet_(tet), faces_(faces),
    weighted_funcs_(weighted_funcs), func_v_and_g_inf_(weighted_funcs.size())
{
  polycube_tet_quality(0, node_.size(2), tet,faces, stat_);
  g_.resize(weighted_funcs_.size());
  for(size_t i = 0; i < weighted_funcs_.size(); ++i) {
    g_[i] = zeros<double>(weighted_funcs_[i].first->dim(), 1);
  }
}

quality_checker::quality_checker(
    const zjucad::matrix::matrix<double> &node,
    const zjucad::matrix::matrix<size_t> &faces,
    std::vector<std::pair<jtf_func_type *, double> > &weighted_funcs)
  :node_(node), faces_(faces),
    weighted_funcs_(weighted_funcs), func_v_and_g_inf_(weighted_funcs.size())
{
  polycube_tet_quality(0, node_.size(2), faces, stat_);
  g_.resize(weighted_funcs_.size());
  for(size_t i = 0; i < weighted_funcs_.size(); ++i) {
    g_[i] = zeros<double>(weighted_funcs_[i].first->dim(), 1);
  }
}

int quality_checker::at_point(const double *x)
{
  polycube_tet_quality(x, node_.size(2), tet_, faces_, stat_);
  size_t i;
#pragma omp parallel for private(i)
  for(i = 0; i < weighted_funcs_.size(); ++i) {
    func_v_and_g_inf_[i].first = 0;
    weighted_funcs_[i].first->val(x, func_v_and_g_inf_[i].first);
    func_v_and_g_inf_[i].first /= weighted_funcs_[i].second;

    g_[i](colon()) = 0;
    weighted_funcs_[i].first->gra(x, &g_[i][0]);
    func_v_and_g_inf_[i].second = max(fabs(g_[i]));
    func_v_and_g_inf_[i].second /= weighted_funcs_[i].second;
  }
  return 0;
}
