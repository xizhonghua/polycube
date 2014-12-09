#ifndef HJ_COSSERAT_H_
#define HJ_COSSERAT_H_

#include <zjucad/matrix/matrix.h>

int init_first_frame(const matrixd &pos,
		     matrixd &Fe);

int position_to_Bishop_frame_at_edge(const matrixd &pos,
				     zjucad::matrix::matrix<matrixd > &Fe);

int edge_frame_to_node_frame(const zjucad::matrix::matrix<matrixd > &Fe,
			     zjucad::matrix::matrix<matrixd > &Fn);

int position_to_Bishop_frame_at_node(const matrixd &pos,
				     zjucad::matrix::matrix<matrixd > &Fn);

#endif
