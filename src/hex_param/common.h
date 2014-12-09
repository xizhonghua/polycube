#ifndef HJ_HEX_PARAM_COMMON_H_
#define HJ_HEX_PARAM_COMMON_H_

#include <vector>
#include <algorithm>
#include <zjucad/matrix/matrix_expression.h>

#include "../tetmesh/tetmesh.h"

#define VERIFY(exp, msg) if(!(exp)) cerr << msg << endl;

template <typename E>
std::vector<typename E::value_type> to_vec(const zjucad::matrix::matrix_expression<E> &e) {
	std::vector<typename E::value_type> rtn(e().size());
	copy(e().begin(), e().end(), rtn.begin());
	return rtn;
}

inline size_t get_tet_of_outside_face(const jtf::mesh::face2tet_adjacent &fa, size_t face_id)
{
	const std::pair<size_t, size_t> &tet_id = fa.face2tet_[face_id];
	if(!fa.is_outside_face(tet_id))
		return -1;
	return (tet_id.first == -1)?tet_id.second:tet_id.first;
}

inline bool is_outside_vertex(const size_t v_idx, const matrixst &outside_face)
{
    if(std::find(outside_face.begin(),outside_face.end(),v_idx) != outside_face.end())
        return true;
    else
        return false;
}
double bary2dist(const matrixd &bary, const matrixd &tet);

/**
   F*uvw=xyz
 */
void estimate_local_frame(
	const jtf::mesh::meshes &tm,
	const matrixst &cut_tet,
    const matrixd &uvw,
    std::vector<matrixd > &local_frames);

int normal_direction(const matrixd &normal,
                     const matrixd &frame);

///
/// @brief extract inner jump type between tets
/// @param tet input tetmesh
/// @param frame input frames
/// @param fa input face_adjacent
/// @param inner_face_jump_type output inner_face_jump_type
///
void extract_inner_jump_type(
    const matrixst &tet,
    const zjucad::matrix::matrix<matrixd> &frame,
    const jtf::mesh::face2tet_adjacent &fa,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type);

#endif
