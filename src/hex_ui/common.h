#ifndef HJ_HEX_UI_COMMON_H_
#define HJ_HEX_UI_COMMON_H_

#include <zjucad/matrix/matrix.h>
#include <iostream>

#include "../common/util.h"

inline double wave_length_from_argv(const char *argv, const matrixd &node)
{
	double wave_length = atof(argv);
	if(wave_length < 0) {
		wave_length = calc_max_bounding_box_edge_length(node)/(-wave_length);
		std::cerr << "# wave length is set to: " << wave_length << std::endl;
	}
	return wave_length;
}

template <typename T>
inline void id2xyz(int id, T *xyz)
{
	for(int d = 0; d < 3; ++d)
		xyz[d] = (id >> d)&0x01;
}

#endif
