#ifndef OPTIMIZE_FRAME_FIELD_L1_H
#define OPTIMIZE_FRAME_FIELD_L1_H

#include <boost/property_tree/ptree.hpp>
#include <zjucad/matrix/matrix.h>
#include "../tetmesh/tetmesh.h"

void optimize_frame_field_L1(const jtf::tet_mesh &tm,
                             zjucad::matrix::matrix<double> & frame,
                             boost::property_tree::ptree &pt);

///
/// \brief optimize_frame_field_L1_euler, variables are euler angle.
/// \param tm
/// \param zyz
/// \param pt
///
void optimize_frame_field_L1_euler(const jtf::tet_mesh & tm,
                                   zjucad::matrix::matrix<double> & zyz,
                                   boost::property_tree::ptree &pt);
#endif // OPTIMIZE_FRAME_FIELD_L1_H
