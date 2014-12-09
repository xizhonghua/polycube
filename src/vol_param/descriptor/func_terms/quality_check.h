#ifndef QUALITY_CHECK_H
#define QUALITY_CHECK_H

#include <jtflib/optimizer/opt.h>
#include <zjucad/matrix/matrix.h>
#include <iostream>
#include <map>
#include "../def.h"
#include "l1-normal.h"

double polycube_L1_area_quality(const double *x, size_t node_num, const zjucad::matrix::matrix<size_t> &faces);

class quality_checker : public jtf::opt::callbacks
{
public:
  quality_checker(const zjucad::matrix::matrix<double> &node,
                  const zjucad::matrix::matrix<size_t> &tet,
                  const zjucad::matrix::matrix<size_t> &faces,
                  std::vector<std::pair<jtf_func *, double> > &weighted_funcs);
  quality_checker(const zjucad::matrix::matrix<double> &node,
                  const zjucad::matrix::matrix<size_t> &faces,
                  std::vector<std::pair<jtf_func *, double> > &weighted_funcs);
  virtual int at_point(const double *x);
protected:
  const zjucad::matrix::matrix<double> &node_;
  const zjucad::matrix::matrix<size_t> tet_, &faces_;
  std::map<std::string, double> stat_;
  const std::vector<std::pair<jtf_func *, double> > &weighted_funcs_;
  std::vector<std::pair<double, double> > func_v_and_g_inf_;
  std::vector<zjucad::matrix::matrix<double> > g_;
};


class unnormalized_normal_quality_checker : public quality_checker
{
public:
  unnormalized_normal_quality_checker(
      const zjucad::matrix::matrix<double> &node,
      const zjucad::matrix::matrix<size_t> &tet,
      const zjucad::matrix::matrix<size_t> &faces,
      std::vector<std::pair<jtf_func *, double> > &weighted_funcs,
      jtf_func &constraint)
    :quality_checker(node, tet, faces, weighted_funcs),
      L1_sqrt_eps_normalize_(2*get_l1_eps_weight()+sqrt(1+get_l1_eps_weight()*get_l1_eps_weight())),
      constraint_(constraint), step_idx_(0) {
    std::cout << "#diagnose_beg: (E_d.v, E_d.g) (E_d.v, E_d.g) (E.v, E.g) L1-area" << std::endl;
    prev_quality_ = beg_quality_ = polycube_L1_area_quality(&node[0], node.size(2), faces);
    quality_gradient_.resize(5);
    fill(quality_gradient_.begin(), quality_gradient_.end(), beg_quality_);
  }
protected:
  virtual int at_point(const double *x) {
    quality_checker::at_point(x);
    // func_v_and_g_inf_[1].first -= L1_sqrt_eps_normalize_;
    // func_v_and_g_inf_[2].first -= L1_sqrt_eps_normalize_
    //         *(weighted_funcs_[1].second/(1+weighted_funcs_[1].second));
    std::cout << "#diagnose: ";
    for(size_t i = 0; i < 3; ++i)
      std::cout << "(" << func_v_and_g_inf_[i].first << ", " << func_v_and_g_inf_[i].second << ") ";
    std::cout << stat_["L1-area"];

    quality_gradient_[step_idx_%quality_gradient_.size()] = fabs(stat_["L1-area"] - prev_quality_);
    prev_quality_ = stat_["L1-area"];
    if(step_idx_ > quality_gradient_.size()) {
        const double avg_gradient = std::accumulate(quality_gradient_.begin(), quality_gradient_.end(), 0.0) / quality_gradient_.size() / beg_quality_;
        std::cout << " " << avg_gradient << std::endl;
        if(avg_gradient < 1e-2 && step_idx_ > quality_gradient_.size()) {
            std::cout << "stop in callbacks." << std::endl;
            step_idx_ = 0;
            return 1;
          }
      }
    else
      std::cout << std::endl;

    ++step_idx_;
    return 0;
  }
  const double L1_sqrt_eps_normalize_;
  jtf_func &constraint_;
  double beg_quality_;
  std::vector<double> quality_gradient_;
  size_t step_idx_;
  double prev_quality_;
};



#endif // QUALITY_CHECK_H
