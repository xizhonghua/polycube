#include "../common/obj_like_format_mat.h"
#include <zjucad/matrix/matrix.h>
#include "../common/def.h"

#ifdef WIN32
#include <windows.h>
#undef max
#undef min
#endif

#include <fstream>
#include <iostream>
#include <vector>

extern "C" {
#include <hjlib/ANN_c.h>
}

using namespace std;
using namespace zjucad::matrix;

int load_vg_edge(const char *path, matrix<int>& edge_)
{
	ifstream ifs(path);
	if(ifs.fail()) {
		cerr << "open " << path << " fail." << endl;
		return __LINE__;
	}

	matrixd  pos_;
	if(read_section(ifs, 3, pos_, "node"))
	   return __LINE__;

	if(read_section(ifs, 2, edge_, "edge"))
		return __LINE__;
	edge_ -= 1;

	cerr << "nn: " << pos_.size(2) << " en: " << edge_.size(2) << endl;

	if(ifs.fail()) {
		cerr << "parse " << path << " fail." << endl;
		return __LINE__;
	}

	return 0;
}
int load_tet(const char *path, matrix<int>& tet_, matrixd& pos_)
{
	ifstream ifs(path);
	if(ifs.fail()) {
		cerr << "open " << path << " fail." << endl;
		return __LINE__;
	}

	if(read_section(ifs, 3, pos_, "node"))
	   return __LINE__;

	matrix<int> face_;
	if(read_section(ifs, 3, face_, "face"))
		return __LINE__;
	face_ -= 1;

	if(read_section(ifs, 4, tet_, "tet"))
		return __LINE__;
	tet_ -= 1;

	cerr << "nn: " << pos_.size(2) << " fn: " << face_.size(2) << " tn: " << tet_.size(2) << endl;

	if(ifs.fail()) {
		cerr << "parse " << path << " fail." << endl;
		return __LINE__;
	}
	return 0;
}
int load_vf_in_tet(const char *path, matrixd& uvw_, matrixd& frame_)
{
	ifstream ifs(path);
	if(ifs.fail()) {
		cerr << "open " << path << " fail." << endl;
		return __LINE__;
	}

	if(read_section(ifs, 3, uvw_, "uvw"))
	   return __LINE__;

	if(read_section(ifs, 9, frame_, "frame"))
	   return __LINE__;

	cerr << "pn: " << uvw_.size(2) << "fn: " << frame_.size(2) << endl;

	if(ifs.fail()) {
		cerr << "parse " << path << " fail." << endl;
		return __LINE__;
	}
	return 0;
}

double bounding_box_len(matrixd& pos_)
{
	double min_max[2][3];
	double r[3];
	for(int d = 0; d < 3; ++d) {
		min_max[0][d] = min(pos_(d, colon()));
		min_max[1][d] = max(pos_(d, colon()));
		r[d] = min_max[1][d] - min_max[0][d];
	}
	double min_len = *min_element(r, r+3);
	return min_len;
}
void cal_tet_uvw_normal_center(matrixd& each_tet_uvw_, matrixd& normal_, 
						matrixd& center_)
{
	int tet_num = each_tet_uvw_.size(2) / 4;
	normal_.resize(3, each_tet_uvw_.size(2));
	center_.resize(3, tet_num);
    matrixd edge[2];
	for(int ti = 0, ni = 0; ti < tet_num; ++ti) {
		int ti_jump = ti*4;
		for(int fi = 0; fi < 4; ++fi, ++ni) {
			edge[0] = each_tet_uvw_(colon(), (fi+2)%4+ti_jump) - each_tet_uvw_(colon(), (fi+0)%4+ti_jump);
			edge[1] = each_tet_uvw_(colon(), (fi+1)%4+ti_jump) - each_tet_uvw_(colon(), (fi+0)%4+ti_jump);
			if(fi%2==0)
				normal_(colon(), ni) = cross(edge[0], edge[1]);
			else
				normal_(colon(), ni) = cross(edge[1], edge[0]);
			const double len = norm(normal_(colon(), ni));
			if(len > 1e-9)
				normal_(colon(), ni) /= len;
			else
				cerr << "degenerated uvw normal at: " << ni << endl;
		}
		center_(colon(), ti) = each_tet_uvw_(colon(), colon(ti_jump, ti_jump+3))*ones<double>(4, 1)/4;
	}
}
void cal_tet_uvw_bounding_box(matrixd& each_tet_uvw_, matrixd& tet_bx_)
{
	int tet_num = each_tet_uvw_.size(2) / 4;
	tet_bx_.resize(3, tet_num*2);

	for(int ti = 0; ti < tet_num; ++ti) {
		matrixd one_tet_uvw = each_tet_uvw_(colon(), colon(ti*4, ti*4+3));

		for(int d = 0; d < 3; ++d) {
			tet_bx_(d, ti*2) = min(one_tet_uvw(d, colon()));
			tet_bx_(d, ti*2+1) = max(one_tet_uvw(d, colon()));
		}
	}
}
matrixd trans_uvw(matrixd& one_uvw_, matrixd& one_fun_)
{
	matrixd rot_fun_ = one_fun_(colon(), colon(0, 2));
	matrixd dis_fun_ = one_fun_(colon(), 3);

	return rot_fun_*one_uvw_+dis_fun_;
}
double cal_avg_gradient(matrix<int>& tet_, matrixd& tet_pos_, 
						matrixd& each_tet_uvw_)
{
	vector<pair<int, int> > tet_edge_index_vec(6);
	tet_edge_index_vec[0] = make_pair(0, 1);
	tet_edge_index_vec[1] = make_pair(0, 2);
	tet_edge_index_vec[2] = make_pair(0, 3);
	tet_edge_index_vec[3] = make_pair(1, 2);
	tet_edge_index_vec[4] = make_pair(1, 3);
	tet_edge_index_vec[5] = make_pair(2, 3);

	//
	double tet_len = 0;
	double uvw_len = 0;
	matrixd edge[2];
	for (int ti = 0; ti < tet_.size(2); ti++)
	{
		int ti_jump = ti*4;
		for(int ei = 0; ei < 6; ++ei) {
			pair<int, int>& tet_edge_idx = tet_edge_index_vec[ei];

			edge[0] = tet_pos_(colon(), tet_(tet_edge_idx.first, ti)) - 
				      tet_pos_(colon(), tet_(tet_edge_idx.second, ti));
			edge[1] = each_tet_uvw_(colon(), (tet_edge_idx.first+ti_jump)) - 
				      each_tet_uvw_(colon(), (tet_edge_idx.second+ti_jump));

			tet_len += norm(edge[0]);
			uvw_len += norm(edge[1]);
		}
	}
	return uvw_len / tet_len;
}
bool is_in_tet_bounding_box(matrixd& pt, matrixd& one_tet_bx_)
{
	if (pt(0, 0) >= one_tet_bx_(0, 0) && 
		pt(1, 0) >= one_tet_bx_(1, 0) &&
		pt(2, 0) >= one_tet_bx_(2, 0) &&
		pt(0, 0) <= one_tet_bx_(0, 1) && 
		pt(1, 0) <= one_tet_bx_(1, 1) &&
		pt(2, 0) <= one_tet_bx_(2, 1)) {
			return true;
	}
	return false;
}
int intege_pt_in_tet_bounding_box(matrixd& one_tet_bx_, double step_, 
								  matrixd& pt)
{
	vector<int> small_i_uvw(3), big_i_uvw(3);
	for (int i = 0; i < 3; i++) {
		small_i_uvw[i] = (int) (one_tet_bx_(i, 0) / step_) - 1;
		big_i_uvw[i] = (int) (one_tet_bx_(i, 1) / step_) + 1;
	}

	//
	int pt_num = 1;
	vector<int> dim_num_vec(3);
	for (int i = 0; i < 3; i++) {
		dim_num_vec[i] = big_i_uvw[i] - small_i_uvw[i];
		pt_num *= dim_num_vec[i];
	}

	//
	if (pt_num > 0) {
		pt.resize(3, pt_num);
		for (int i = 0; i < dim_num_vec[0]; i++) {
			for (int j = 0; j < dim_num_vec[1]; j++) {
				for (int k = 0; k < dim_num_vec[2]; k++)  {
					int index_ = i*dim_num_vec[1]*dim_num_vec[2] + j*dim_num_vec[2]+k;
					   pt(0, index_) = (big_i_uvw[0]-i)*step_;
					   pt(1, index_) = (big_i_uvw[1]-j)*step_;
					   pt(2, index_) = (big_i_uvw[2]-k)*step_;
				}
			}
		}
	}

	return pt_num;
}
bool is_in_tet(matrixd& pt, matrixd& one_tet_pos_,
			   matrixd& one_tet_normal_, matrixd& bc_coord)
{
	for(int fi = 0; fi < 4; ++fi) {
		matrixd pv_ = pt - one_tet_pos_(colon(), fi); 
		const double len = norm(pv_);
		if(len > 1e-9)
			pv_ /= len;
		else
			cerr << "degenerated p-v normal at: " << fi << endl;
		double fn_test_ = dot(pv_, one_tet_normal_(colon(), fi));
		if (fn_test_ > 0 && fabs(fn_test_) > 1e-8) {
			return false;
		}
	}

	// compute barycenter coord
	bc_coord.resize(4, 1);
	matrixd v_03 = one_tet_pos_(colon(), 0) - one_tet_pos_(colon(), 3);
	matrixd v_12 = one_tet_pos_(colon(), 1) - one_tet_pos_(colon(), 2);
	
	bc_coord(0, 0) = dot(one_tet_normal_(colon(), 1), (pt - one_tet_pos_(colon(), 3))) /
		        dot(one_tet_normal_(colon(), 1), v_03);
	bc_coord(1, 0) = dot(one_tet_normal_(colon(), 2), (pt - one_tet_pos_(colon(), 2))) /
		        dot(one_tet_normal_(colon(), 2), v_12);
	bc_coord(2, 0) = dot(one_tet_normal_(colon(), 3), (pt - one_tet_pos_(colon(), 1))) /
		        dot(one_tet_normal_(colon(), 3), -v_12);
	bc_coord(3, 0) = dot(one_tet_normal_(colon(), 0), (pt - one_tet_pos_(colon(), 0))) /
		        dot(one_tet_normal_(colon(), 0), -v_03);

	// debug.
	matrixd bc_pt = one_tet_pos_*bc_coord;
	bc_pt -= pt;
	if (norm(bc_pt) > 1e-8) {
		cerr << "error bc coord" << endl;
	}
	return true;
}
void get_pt_neighbors_on_tet_frame_cube(matrixd& pt_, int res_, double tet_step, 
									   matrixd& tet_frame, 
									   matrixd& cube_neigh)
{
	int half_res_ = (int) res_ / 2;
	cube_neigh.resize(3, res_*res_*res_);
	int index_ = 0;
	
	matrixd u_ = tet_frame(colon(), 0);
	matrixd v_ = tet_frame(colon(), 1);
	matrixd w_ = tet_frame(colon(), 2);

	for (int i = -half_res_; i <= half_res_; i++) {
		for (int j = -half_res_; j <= half_res_; j++) {
			for (int k = -half_res_; k <= half_res_; k++) {
				cube_neigh(colon(), index_++) = pt_ +
					u_*(tet_step*i) + v_*(tet_step*j) + w_*(tet_step*k);
			}
		}
	}
}

int gp2pts(int argc, char *argv[])
{
	if(argc < 7) {
		cerr << "Usage: gp2pts " << endl;
		return __LINE__;
	}

	int pts_dim = atoi(argv[1]);
	int neigh_res = atoi(argv[2]);
	if (neigh_res % 2 == 0) {
		cerr << "neighbor reslution " << neigh_res << " must be odd number.\n";
	}

	// load 
	matrix<int> vg_edge_, tet_;
	matrixd tet_pos_, each_tet_uvw_vec, each_tet_frame_vec;
	vector<matrixd > trans_fun_;

	load_vg_edge(argv[3], vg_edge_);
	load_tet(argv[4], tet_, tet_pos_);
	load_vf_in_tet(argv[5], each_tet_uvw_vec, each_tet_frame_vec);

	// set sample points
	double tet_bx_min_len = bounding_box_len(tet_pos_);
	double uvw_step =  tet_bx_min_len * 
		               cal_avg_gradient(tet_, tet_pos_, each_tet_uvw_vec)/ 
		               pts_dim;  
	double tet_step = tet_bx_min_len / pts_dim;
	printf("tet distance = %f\n", tet_step);


	// cal uvw tet normal, center, bounding box.
	matrixd uvw_tet_normal_, uvw_tet_center_, uvw_tet_bx_;
	cal_tet_uvw_normal_center(each_tet_uvw_vec, uvw_tet_normal_, uvw_tet_center_);
	cal_tet_uvw_bounding_box(each_tet_uvw_vec, uvw_tet_bx_);

	// find each sample point in uvw tet.
	int find_num = 0;
	vector<matrixd > find_pt_bc_coord_vec;
	vector<int> find_pt_ti_vec;

	for (int ti = 0; ti < tet_.size(2); ti++) {
		matrixd one_tet_bx = uvw_tet_bx_(colon(), colon(ti*2, ti*2+1));
		matrixd one_tet_pos = each_tet_uvw_vec(colon(), colon(ti*4, ti*4+3));
		matrixd one_tet_normal = uvw_tet_normal_(colon(), colon(ti*4, ti*4+3));
		matrixd bc_coord;
		matrixd int_pts_;
		int pt_num = intege_pt_in_tet_bounding_box(one_tet_bx, uvw_step, int_pts_);
		for (int vi = 0; vi < pt_num; vi++) {
			matrixd pt_ = int_pts_(colon(), vi);

			if (is_in_tet(pt_, one_tet_pos, one_tet_normal, bc_coord)) {
				find_pt_bc_coord_vec.push_back(bc_coord);
				find_pt_ti_vec.push_back(ti);
				find_num++;
			}
		}
	}
	printf("find %d sample points from uvw.\n", find_num);

	// convert uvw sample point to real 3d positions.
	matrixd pts_coord_(3, find_num);
	for (int i = 0; i < find_num; i++) {
		int ti = find_pt_ti_vec[i];
		matrixd&	bc_coord = find_pt_bc_coord_vec[i];
		matrixd one_tet_pos(3, 4);	
		for(int vi = 0; vi < 4; ++vi) {
		   one_tet_pos(colon(), vi) = tet_pos_(colon(), tet_(vi, ti));
		}
		pts_coord_(colon(), i) = one_tet_pos * bc_coord;
	}

	// build ann for searching
	vector<double *> ppts(pts_coord_.size(2));
	for(int pi = 0; pi < pts_coord_.size(2); ++pi) {
		ppts[pi] = &pts_coord_(0, pi);
	}
	void *ANNkd_tree_handle = ANNkd_tree_new(&ppts[0], pts_coord_.size(2), 3);

	double xyz[3];
	int num_neighbor = 5;
	vector<int> idx(num_neighbor);
	vector<double> dist(num_neighbor);

	// first remove repeat pt.
	double pt_nearest_len = tet_step * 1e-2;
	vector<int> pts_valid_flag(pts_coord_.size(2));
	fill(pts_valid_flag.begin(), pts_valid_flag.end(), 1);
	int remove_pt_num = 0;
	do {
		remove_pt_num = 0;
		for (int pi = 0; pi < pts_coord_.size(2); ++pi) {
			if (pts_valid_flag[pi] > 0) {
				for (int k =0; k < 3; k++) {
					xyz[k]=pts_coord_(k, pi);
				}
				ANNkd_tree_search(ANNkd_tree_handle,
					xyz, num_neighbor, &idx[0], &dist[0]);

				// ni=1, ignore itself
				for (int ni = 1; ni < num_neighbor; ni++) {
					if (pts_valid_flag[idx[ni]] > 0 &&
						dist[ni] < pt_nearest_len) {
						pts_valid_flag[idx[ni]] = -1;
						remove_pt_num++;
					}
				}
			}
		}
	}while(remove_pt_num > 0);
	ANNkd_tree_delete(ANNkd_tree_handle);

	int pt_index_ = 0;
	for (size_t i = 0; i < pts_valid_flag.size(); i++) {
		pts_valid_flag[i] = pts_valid_flag[i] < 0 ? -1 : pt_index_++;
	}
	printf("remove %ld pts, left %ld pts.\n", pts_coord_.size(2)-pt_index_, pt_index_);

	// reset ANN here.
	matrixd up_pts_coord_(3, pt_index_);
	vector<int> up_find_pt_ti_vec;
	for (int i = 0; i < pts_coord_.size(2); i++) {
		if (pts_valid_flag[i] >= 0) {
			up_pts_coord_(colon(), pts_valid_flag[i]) = pts_coord_(colon(), i);
			up_find_pt_ti_vec.push_back(find_pt_ti_vec[i]);
		}
	}
	ppts.resize(up_pts_coord_.size(2));
	for(int pi = 0; pi < up_pts_coord_.size(2); ++pi) {
		ppts[pi] = &up_pts_coord_(0, pi);
	}
	ANNkd_tree_handle = ANNkd_tree_new(&ppts[0], up_pts_coord_.size(2), 3);

	// then find (neigh_res*neigh_res*neigh_res) neighbor for cube
	double pt_farest_len = tet_step * 2;
	num_neighbor = 1;
	idx.resize(num_neighbor);
	dist.resize(num_neighbor);
	int neigh_sp_num = neigh_res*neigh_res*neigh_res;
	matrix<int> cube_neigh_pid_(neigh_sp_num, up_pts_coord_.size(2));
	for (int pi = 0; pi < up_pts_coord_.size(2); ++pi) {

		int ti = up_find_pt_ti_vec[pi];
		matrixd pt_coord = up_pts_coord_(colon(), pi);
		matrixd tf_(3, 3), cube_neigh;

		for (int ri = 0; ri < 3; ri++) {
			for (int ci = 0; ci < 3; ci++) {
				tf_(ri, ci) = each_tet_frame_vec(ri*3+ci, ti); 
			}
		}

		//
		get_pt_neighbors_on_tet_frame_cube(pt_coord, neigh_res, tet_step, tf_, cube_neigh);

		for (int i = 0; i < cube_neigh.size(2); i++) {
			for (int k =0; k < 3; k++) {
				xyz[k]=cube_neigh(k, i);
			}
			ANNkd_tree_search(ANNkd_tree_handle,
				xyz, num_neighbor, &idx[0], &dist[0]);
			if (dist[0] < pt_farest_len) {
				cube_neigh_pid_(i, pi) = idx[0];
			} else {
				cube_neigh_pid_(i, pi) = -1;
			}
		}
	}
	ANNkd_tree_delete(ANNkd_tree_handle);

	//
	ofstream ofs(argv[6]);
	write_section(ofs, up_pts_coord_.size(2), 3, &up_pts_coord_[0], "gp_node", "v");
	write_section(ofs, cube_neigh_pid_.size(2), neigh_sp_num, &cube_neigh_pid_[0], "gp_neigh", "n");
	ofs.close();
	return 0;
}
