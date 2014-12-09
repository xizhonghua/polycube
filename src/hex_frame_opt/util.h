#ifndef HEX_FRAME_OPT_UTIL_H
#define HEX_FRAME_OPT_UTIL_H
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_set.hpp>
#include "../common/def.h"


int load_box_constraint(
    const char * box_file,
    std::vector<std::tuple<matrixd, boost::unordered_set<size_t> > > & box_vec);

int load_plane_constraint(
    const char * plane_file,
    std::vector<std::tuple<matrixd, boost::unordered_set<size_t> > > & plane_vec);

int load_line_constraint(
    const char * line_file,
    std::vector<std::tuple<matrixd, boost::unordered_set<size_t> > > & line_vec);

#endif // #ifndef HEX_FRAME_OPT_UTIL_H
