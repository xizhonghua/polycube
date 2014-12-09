#ifndef OPTIMIZE_FRAME_FIELD_WITH_TYPE_H
#define OPTIMIZE_FRAME_FIELD_WITH_TYPE_H

#include "../tetmesh/tetmesh.h"
void optimize_frame_field_with_type(const jtf::tet_mesh &tm,
                                    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_type,
                                    zjucad::matrix::matrix<double> &zyz);

#endif // OPTIMIZE_FRAME_FIELD_WITH_TYPE_H
