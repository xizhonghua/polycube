#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>

#include "../common/IO.h"
extern "C" {
#include <hjlib/ANN_c.h>
}

#include "../tetmesh/tetmesh.h"

#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include <zjucad/matrix/io.h>

using namespace zjucad::matrix;
using namespace std;

//! return: distance to the tet according to the minimal barycenter
//! coordinates, < 0 means out side of the tet
static inline double distance_to_tet(const matrixd &pt, const matrixd &tet)
{
  matrixd B = ones<double>(4, 4);
  B(colon(0, 2), colon()) = tet;
  if(inv(B)) {
    cerr << "# warning degenerated tet." << endl;
    return -numeric_limits<double>::infinity();
  }
  matrixd coor = B(colon(), colon(0, 2))*pt+B(colon(), 3);
  return min(coor);
}

static void tet_centers(const jtf::mesh::meshes &tm, matrixd &pts)
{
  pts.resize(3, tm.mesh_.size(2));
  for(size_t pi = 0; pi < pts.size(2); ++pi)
    pts(colon(), pi) = tm.node_(colon(), tm.mesh_(colon(), pi))*ones<double>(4, 1)/4.0;
}

static int build_nearest(const jtf::mesh::meshes &tm, const matrixd &dst_pts, matrix<uint32_t> &nearest)
{
  matrixd src_pts;
  tet_centers(tm, src_pts);

  const size_t DIM = src_pts.size(1);
  assert(dst_pts.size(1) == DIM);
  vector<const double *> ppts(src_pts.size(2));
  for(size_t pi = 0; pi < ppts.size(); ++pi)
    ppts[pi] = &src_pts(0, pi);
  void *ANNkd_tree_handle = ANNkd_tree_new(const_cast<double **>(&ppts[0]), ppts.size(), DIM);

  const size_t MAX_SEARCH_NUM = 128;
  vector<int> idx(MAX_SEARCH_NUM);
  vector<double> dist2(MAX_SEARCH_NUM);
  vector<pair<double, int> > dist_to_tet(MAX_SEARCH_NUM);
  for(size_t dpi = 0; dpi < dst_pts.size(2); ++dpi) {
    size_t nb_num = 4;
    pair<double, int> best(0,0);
    while(nb_num < MAX_SEARCH_NUM) { // in case of poor tet shape
      // get possible tet which contains pt
      fill(dist2.begin(), dist2.end(), -numeric_limits<double>::max());
      ANNkd_tree_search(ANNkd_tree_handle,
                        const_cast<double *>(&dst_pts(0, dpi)), nb_num, &idx[0], &dist2[0]);

      // sort the quality of containing, in case of some points
      // which are outside of the whole tet mesh
      dist_to_tet.clear();
      for(size_t nbi = 0; nbi < nb_num && dist2[nbi] >= 0; ++nbi) {
        dist_to_tet.push_back(make_pair(
                                distance_to_tet(dst_pts(colon(), dpi),
                                                tm.node_(colon(), tm.mesh_(colon(), idx[nbi]))),
                                nbi));
        if(dist_to_tet[nbi].first > 0)
          break;
      }
      best = *max_element(dist_to_tet.begin(), dist_to_tet.end());
      if(best.first < 0) // still out side
        nb_num *= 2;
      else
        break;
    } // end while
    //TODO: here may introduce error, if beyond the search radius, just pick tet0
    nearest[dpi] = idx[best.second];
    if(best.first < 0) {
      cerr << "# warning: outside point " << dpi << ": " << best.first << endl;
    }
  }
  ANNkd_tree_delete(ANNkd_tree_handle);

  return 0;
}

static int interp(int argc, char *argv[], size_t dim, const matrix<uint32_t> &nearest)
{
  for(int ai = 0; ai+1 < argc; ai += 2) {
    matrixd src_val;
    if(jtf::mesh::read_matrix(argv[ai], src_val))
      return __LINE__;
    if(src_val.size(1) != dim) {
      matrixd tmp = src_val;
      src_val.resize(dim, tmp.size()/dim);
      if(src_val.size() != tmp.size()) {
        cerr << "incompatible value file: " << dim << " " << tmp.size(2) << endl;
        return __LINE__;
      }
      copy(tmp.begin(), tmp.end(), src_val.begin());
    }
    matrixd dst_val = src_val(colon(), nearest);
    if(jtf::mesh::write_matrix(argv[ai+1], dst_val))
      return __LINE__;
  }
  return 0;
}

int interp_tet(int argc, char *argv[])
{
  if(argc < 3) {
    cerr << "interp src_tet dst_tet nearest [dim value output]*" << endl;
    return __LINE__;
  }

  matrix<uint32_t> nearest;
  int has_nearest = (jtf::mesh::read_matrix(argv[3], nearest) == 0);

  if(!has_nearest) {
    cerr << "# create nearest" << endl;
    jtf::mesh::meshes src_tm, dst_tm;
    if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &src_tm.node_, &src_tm.mesh_))
      return __LINE__;
    if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &dst_tm.node_, &dst_tm.mesh_))
      return __LINE__;

    matrixd dst_pts;
    tet_centers(dst_tm, dst_pts);

    nearest.resize(dst_pts.size(2));
    if(build_nearest(src_tm, dst_pts, nearest))
      return __LINE__;
    if(jtf::mesh::write_matrix(argv[3], nearest))
      return __LINE__;
    cerr << "# [info] generate nearset size "
         << nearest.size(1) << " "
         << nearest.size(2) << endl;
  }
  if(argc > 6)
    return interp(argc-5, argv+5, atoi(argv[4]), nearest);
  return 0;
}
