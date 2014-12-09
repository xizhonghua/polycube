#ifndef FRAME_FUNCTION_H
#define FRAME_FUNCTION_H

#include <string>
#include <tuple>

#include <boost/unordered_set.hpp>

#include <hjlib/function/function.h>
#include <jtflib/function/function.h>
#include <zjucad/ptree/ptree.h>

#include "../common/def.h"
#include "../tetmesh/tetmesh.h"

class sym_frame_opt_sys;

class control_unit{
public:
  typedef std::vector<std::tuple<matrixd, boost::unordered_set<size_t> > > UNIT_TYPE;
  std::vector<std::pair<std::string, std::shared_ptr<UNIT_TYPE> > > units_;
};
//! @param func_type : "sh2zyz", "sh"
hj::function::function_t<double, int32_t> *
build_frame_smooth_func(const matrixst & tet, const matrixd & node,
                        std::shared_ptr<jtf::mesh::face2tet_adjacent> &fa,
                        const std::string &func_type,
                        const std::string smooth_strategy);

//! @param func_type : "sh2zyz", "sh"
hj::function::function_t<double, int32_t> *
build_normal_align_func(const matrixst &tet, const matrixd &node,
                        const matrixst & aligned_face_idx,
                        const matrixd & aligned_normal,
                        std::shared_ptr<jtf::mesh::face2tet_adjacent> & fa,
                        const std::string & func_type,
                        const double weight);

hj::function::function_t<double, int32_t> *
build_normal_align_func_tri(
    const zjucad::matrix::matrix<size_t> &tri,
    const zjucad::matrix::matrix<double> &node,
    const zjucad::matrix::matrix<double> & normal,
    const std::string strategy,
    const double w);

//void build_normal_align_func_tri_no_quadratic(
//    const zjucad::matrix::matrix<size_t> &tri,
//    const zjucad::matrix::matrix<double> &node,
//    const zjucad::matrix::matrix<double> & normal,
//    const std::string strategy,
//    std::shared_ptr<std::vector<std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > > > constraints);

//! @param aligned_frame is 3 * n matrix
hj::function::function_t<double, int32_t> *
build_frame_fix_func(const matrixst & tet, const matrixd & node,
                     const matrixst & aligned_tet_idx,
                     const matrixd & aligned_frame,
                     std::shared_ptr<jtf::mesh::face2tet_adjacent> & fa,
                     const double weight);

hj::function::function_t<double, int32_t> *
build_control_function(const matrixst & tet, const matrixd &node,
                       const control_unit &cu, const double control_weight);

std::shared_ptr<jtf::function::functionN1_t<double, int32_t> >
build_frame_surface_smooth_function_jtf_tri(
    const size_t cell_num,
    const zjucad::matrix::matrix<size_t> & tri,
    const zjucad::matrix::matrix<double> & node,
    const zjucad::matrix::matrix<double> & normal,
    const jtf::mesh::edge2cell_adjacent &ea,
    const std::string strategy,
    const double w,
    boost::property_tree::ptree &pt,
    const zjucad::matrix::matrix<size_t> * tri_idx = 0,
    const jtf::mesh::face2tet_adjacent * fa = 0);


int generate_cross_field_according_to_angle_defet(
    const jtf::mesh::edge2cell_adjacent &ea,
    const zjucad::matrix::matrix<size_t> & tri,
    const zjucad::matrix::matrix<double> & node,
    const zjucad::matrix::matrix<double> & angle_defect,
    const zjucad::matrix::matrix<double> & normal,
    const char * fv_file);
#endif // FRAME_FUNCTION_H
