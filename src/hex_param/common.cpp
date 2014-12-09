#include <iostream>
#include <jtflib/mesh/util.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include "common.h"
#include "../common/transition.h"
#include "../common/transition_type.h"
using namespace std;
using namespace zjucad::matrix;

double bary2dist(const matrixd &bary, const matrixd &tet)
{
	const size_t dim = bary.size()-1;
	const size_t opposite_node_idx = min_element(bary.begin(), bary.end())-bary.begin();
	const double total_tet_vol = jtf::mesh::cal_tet_vol(tet);
	VERIFY(total_tet_vol > 0, "# error in tet volume.");
	const double vol = total_tet_vol*bary[opposite_node_idx];
	vector<matrixd > base_tri;
	for(size_t i = 0; i < 4; ++i) {
		if(i == opposite_node_idx) continue;
		base_tri.push_back(tet(colon(), i));
	}
	const matrixd edges[2] = {
		base_tri[1]- base_tri[0],
		base_tri[2]- base_tri[0]
	};
	const double area = norm(cross(edges[0], edges[1]))/2;
	return vol*3/area;
}

void estimate_local_frame(
	const jtf::mesh::meshes &tm,
	const matrixst &cut_tet,
	const matrixd &uvw,
	vector<matrixd > &local_frames)
{
	for(size_t ti = 0; ti < cut_tet.size(2); ++ti) {
		matrixd P = ones<double>(4, 4);
		P(colon(0, 2), colon(0, 3)) = uvw(colon(), cut_tet(colon(), ti));
        if(inv(P)) {
			cerr << "# inv fail at: " << ti << endl;
		}
		local_frames[ti] = tm.node_(colon(), tm.mesh_(colon(), ti))*P(colon(), colon(0, 2));
	}
}

int normal_direction(const matrixd &normal,
					 const matrixd &frame)
{
	vector<pair<double, int> > align_err(3);
	for(int i = 0; i < 3; ++i)
		align_err[i] = make_pair(fabs(dot(frame(colon(), i), normal)), i);
	return max_element(align_err.begin(), align_err.end())->second;
}


void extract_inner_jump_type(
    const matrixst &tet,
    const matrix<matrixd > &frame_cut_tet,
    const jtf::mesh::face2tet_adjacent &fa,
    boost::unordered_map<pair<size_t,size_t>,size_t>  &inner_face_jump_type)
{
  matrixst ptr,idx;
  matrixd rot = zeros<double>(3,3);

  tet2tet_adj_matrix(tet,ptr,idx,&fa);
  assert(ptr.size() == tet.size(2) + 1);
  for(size_t t = 1; t< ptr.size(); ++t) {
    if(ptr[t] == ptr[t-1]) continue; // strange case, seams this tet has no adjacent tets

    for(size_t j = ptr[t-1]; j < ptr[t]; ++j) {
      const matrixd &tet_from = frame_cut_tet[t-1];
      const matrixd &tet_to = frame_cut_tet[idx[j]];

      get_best_alignment(&tet_from[0],&tet_to[0], &rot[0]);
      inner_face_jump_type[make_pair(t-1,idx[j])] = type_transition1(rot);
    }
  }
}
