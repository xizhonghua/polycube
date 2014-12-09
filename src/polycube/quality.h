#ifndef HJ_POLYCUBE_QUALITY_H_
#define HJ_POLYCUBE_QUALITY_H_

#include <map>
#include <string>
#include <vector>

#include <zjucad/matrix/matrix.h>
#include <jtflib/optimizer/opt.h>

double polycube_L1_area_quality(const double *x, size_t node_num, const zjucad::matrix::matrix<size_t> &faces);

int polycube_tet_quality(const double *x,
                         const zjucad::matrix::matrix<size_t> &tet,
                         const zjucad::matrix::matrix<size_t> &faces,
                         std::map<std::string, double> &stat);

//! for checking the goodness of the iterations
class quality_checker : public jtf::opt::callbacks
{
public:
  typedef jtf::function::functionN1_t<double,int32_t> jtf_func_type;
  quality_checker(const zjucad::matrix::matrix<double> &node, const zjucad::matrix::matrix<size_t> &tet,
                  const zjucad::matrix::matrix<size_t> &faces, std::vector<std::pair<jtf_func_type *, double> > &weighted_funcs);
  quality_checker(const zjucad::matrix::matrix<double> &node, const zjucad::matrix::matrix<size_t> &faces,
                  std::vector<std::pair<jtf_func_type *, double> > &weighted_funcs);
  virtual int at_point(const double *x);
protected:
  const zjucad::matrix::matrix<double> &node_;
  const zjucad::matrix::matrix<size_t> tet_, &faces_;
  std::map<std::string, double> stat_;
  const std::vector<std::pair<jtf_func_type *, double> > &weighted_funcs_;
  std::vector<std::pair<double, double> > func_v_and_g_inf_;
  std::vector<zjucad::matrix::matrix<double> > g_;
};

extern double L1_sqrt_eps;

class unnormalized_normal_quality_checker : public quality_checker
{
public:
  unnormalized_normal_quality_checker(
      const zjucad::matrix::matrix<double> &node, const zjucad::matrix::matrix<size_t> &tet,
      const zjucad::matrix::matrix<size_t> &faces,
      std::vector<std::pair<jtf::function::functionN1_t<double,int32_t> *, double> > &weighted_funcs,
      const jtf::function::functionN1_t<double,int32_t> &constraint)
    :quality_checker(node, tet, faces, weighted_funcs),
      L1_sqrt_eps_normalize_(2*L1_sqrt_eps+sqrt(1+L1_sqrt_eps*L1_sqrt_eps)),
      constraint_(constraint), step_idx_(0) {
    std::cout << "#diagnose_beg: (E_d.v, E_d.g) (E_d.v, E_d.g) (E.v, E.g) L1-area" << std::endl;
    prev_quality_ = beg_quality_ = polycube_L1_area_quality(&node[0], node.size(2), faces);
    quality_gradient_.resize(5);
    fill(quality_gradient_.begin(), quality_gradient_.end(), beg_quality_);
  }

  unnormalized_normal_quality_checker(
      const zjucad::matrix::matrix<double> &node,
      const zjucad::matrix::matrix<size_t> &faces,
      std::vector<std::pair<jtf::function::functionN1_t<double,int32_t> *, double> > &weighted_funcs,
      const jtf::function::functionN1_t<double,int32_t> &constraint)
    :quality_checker(node, faces, weighted_funcs),
      L1_sqrt_eps_normalize_(2*L1_sqrt_eps+sqrt(1+L1_sqrt_eps*L1_sqrt_eps)),
      constraint_(constraint), step_idx_(0) {
    std::cout << "#diagnose_beg: (E_d.v, E_d.g) (E_d.v, E_d.g) (E.v, E.g) L1-area" << std::endl;
    prev_quality_ = beg_quality_ = polycube_L1_area_quality(&node[0], node.size(2), faces);
    quality_gradient_.resize(5);
    fill(quality_gradient_.begin(), quality_gradient_.end(), beg_quality_);
  }
protected:
  virtual int at_point(const double *x) {
    using namespace std;
    quality_checker::at_point(x);
     func_v_and_g_inf_[1].first -= L1_sqrt_eps_normalize_;
     func_v_and_g_inf_[2].first -= L1_sqrt_eps_normalize_
             *(weighted_funcs_[1].second/(1+weighted_funcs_[1].second));
    cout << "#diagnose: ";
    for(size_t i = 0; i < 3; ++i)
      cout << "(" << func_v_and_g_inf_[i].first << ", " << func_v_and_g_inf_[i].second << ") ";
    cout << stat_["L1-area"];

    quality_gradient_[step_idx_%quality_gradient_.size()] = std::fabs(stat_["L1-area"] - prev_quality_);
    prev_quality_ = stat_["L1-area"];
    if(step_idx_ > quality_gradient_.size()) {
        const double avg_gradient = std::accumulate(quality_gradient_.begin(), quality_gradient_.end(), 0.0) / quality_gradient_.size() / beg_quality_;
        cout << " " << avg_gradient << endl;
        if(avg_gradient < 1e-2 && step_idx_ > quality_gradient_.size()) {
            cout << "stop in callbacks." << endl;
            return 1;
          }
      }
    else
      cout << endl;

    ++step_idx_;
    return 0;
  }
  const double L1_sqrt_eps_normalize_;
  const jtf::function::functionN1_t<double,int32_t> &constraint_;
  double beg_quality_;
  std::vector<double> quality_gradient_;
  size_t step_idx_;
  double prev_quality_;
};
#endif
