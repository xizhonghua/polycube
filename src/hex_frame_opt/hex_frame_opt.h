#ifndef HJ_SYM_FRAME_OPT_H_
#define HJ_SYM_FRAME_OPT_H_

#include <memory>
#include <tuple>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <zjucad/matrix/matrix.h>
#include <hjlib/sparse/sparse.h>
#include <hjlib/function/function.h>

#include "hex_function_term.h"
#include "../common/def.h"


class sym_frame_opt_sys;
class gradient_function_f;
class sym_frame_opt_function;
//class gradient_function_f_inner_4face_adjacent;
class fix_function;
class fix_function_inner;
class align_function;




struct constraint_align
{
  std::shared_ptr<std::vector<std::tuple<matrixd, boost::unordered_set<size_t> > > > box_vec;
  std::shared_ptr<std::vector<std::tuple<matrixd, boost::unordered_set<size_t> > > > plane_vec;
  std::shared_ptr<std::vector<std::tuple<matrixd, boost::unordered_set<size_t> > > > line_vec;
};

class cut_inner_smooth_function;
class cut_surface_align_function;
class cut_jump_smooth_function;
class align_function_LP;
class cut_inner_smooth_linear_function;
class cut_surface_align_linear_function;
class cut_jump_smooth_linear_function;
class cut_RTR_linear_function;
class cut_surface_align_linear_function_new;
class cut_jump_smooth_linear_function_new;
class cut_RTR_linear_function_new;

class non_sym_frame_opt
{
public:
  int setup_equations_new(
      const matrixst &tet,
      const matrixd &node,
      const matrixd &fixed_frame,
      const zjucad::matrix::matrix<zjucad::matrix::idx_type> &fixed_frame_idx,
      //const matrixst &surf_aligned_type,
      const boost::unordered_map<size_t,size_t> &surface_type,
      const matrixst &aligned_surface, // face
      const zjucad::matrix::matrix<zjucad::matrix::idx_type> &aligned_idx,
      const double weight[2], //! [align, fix], default 1e6
      const matrixd &stiff, // tet weighting
      //const matrixst &face_pair,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> &jump_face_type,
      const double LP_surface,
      const double LP_smooth
      );

  int setup_equations_linear_new(
      const matrixst &tet,
      const matrixd &node,
      const matrixd &fixed_frame,
      const zjucad::matrix::matrix<zjucad::matrix::idx_type> &fixed_frame_idx,
      //const matrixst &surf_aligned_type,
      const boost::unordered_map<size_t,size_t> &surface_type,
      const matrixst &aligned_surface, // face
      const zjucad::matrix::matrix<zjucad::matrix::idx_type> &aligned_idx,
      const double weight[2], //! [align, fix], default 1e6
      const matrixd &stiff, // tet weighting
      //const matrixst &face_pair,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> &jump_face_type,
      const double LP_surface,
      const double LP_smooth,
      const double LP_rtr,
      const double RTR_w
      );

  hj::function::function *get(void) { return all_f_.get(); }
  std::vector<std::shared_ptr<hj::function::function> > funcs_;
  std::unique_ptr<hj::function::function> all_f_;
 // std::unique_ptr<sym_frame_opt_sys> sys_;
  //std::unique_ptr<non_sym_frame_opt_sys> sys_;
  std::vector<std::shared_ptr<cut_inner_smooth_function> > cut_inner_smooth_function_;
  std::vector<std::shared_ptr<cut_surface_align_function> > cut_surface_align_function_;
  //std::vector<std::shared_ptr<align_function> > cut_surface_align_function2_;
  std::vector<std::shared_ptr<align_function_LP> > cut_surface_align_function2_;
  std::vector<std::shared_ptr<cut_jump_smooth_function> > cut_jump_smooth_function_;

  std::vector<std::shared_ptr<cut_inner_smooth_linear_function> > cut_inner_smooth_linear_function_;
  std::vector<std::shared_ptr<cut_surface_align_linear_function> > cut_surface_align_linear_function_;
  std::vector<std::shared_ptr<cut_jump_smooth_linear_function> > cut_jump_smooth_linear_function_;
  std::vector<std::shared_ptr<cut_RTR_linear_function> > cut_RTR_linear_function_;
  std::vector<std::shared_ptr<cut_surface_align_linear_function_new> > cut_surface_align_linear_function_new_;
  std::vector<std::shared_ptr<cut_jump_smooth_linear_function_new> > cut_jump_smooth_linear_function_new_;
  std::vector<std::shared_ptr<cut_RTR_linear_function_new> > cut_RTR_linear_function_new_;
};

#endif
