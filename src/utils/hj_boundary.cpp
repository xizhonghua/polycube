#include <jtflib/mesh/mesh.h>
#include "../tetmesh/tetmesh.h"
#include "../common/IO.h"
#include "../hex_param/global_alignment.h"
#include "../common/zyz.h"
#include "../common/vtk.h"
#include "../common/transition.h"
#include "../common/transition_type.h"
#include "../numeric/util.h"
#include <iostream>

#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

using namespace std;
using namespace zjucad::matrix;

inline int load_zyz(char * zyz_file,
                    matrix<double> & zyz,
                    const size_t tet_size)
{
  if(read_zyz(zyz_file, zyz)) {
      cerr << "# [error] can not read zyz file." << endl;
      return __LINE__;
    }
  if(zyz.size(2) != tet_size){
      cerr << "# [error] impatible zyz field." << endl;
      return __LINE__;
    }
  return 0;
}

void cal_Rab(const pair<size_t,size_t> & face_pair,
             const matrix<size_t> & surface,
             const zjucad::matrix::matrix<double> &node,
             const matrix<double> & na, const matrix<double> & nb,
             matrix<double> &Rab)
{
  size_t common_edge[2];
  jtf::mesh::find_common_edge(surface(colon(), face_pair.first), surface(colon(), face_pair.second), &common_edge[0]);
  const size_t other_p = std::accumulate(surface(colon(), face_pair.first).begin(),
                                         surface(colon(), face_pair.first).end(), static_cast<size_t>(0))
      - common_edge[0] - common_edge[1];
  size_t other_p_idx = std::find(surface(colon(), face_pair.first).begin(),
                                 surface(colon(), face_pair.first).end(), other_p)
      -surface(colon(), face_pair.first).begin();
  size_t next_p = surface((other_p_idx+1)%3, face_pair.first);
  if(common_edge[1] == next_p) swap(common_edge[0], common_edge[1]);

  matrix<double> dir = node(colon(), common_edge[1]) - node(colon(), common_edge[0]);
  dir /= norm(dir);

  matrix<double> r1(3,3), r2(3,3);
  r1(colon(),0) = na;
  r1(colon(),1) = dir;
  r1(colon(),2) = cross(na, dir);

  r2(colon(),0) = nb;
  r2(colon(),1) = dir;
  r2(colon(),2) = cross(nb,dir);

  inv(r1);
  Rab = r2*r1;
}

void cal_tba(zjucad::matrix::matrix<double> &tba,
             const vector<size_t> & one_ring_tets,
             const pair<size_t,size_t> & tet_ab,
             const matrix<matrix<double> > & frame){
  assert(std::find(one_ring_tets.begin(), one_ring_tets.end(), tet_ab.first) != one_ring_tets.end());
  assert(std::find(one_ring_tets.begin(), one_ring_tets.end(), tet_ab.second) != one_ring_tets.end());
  vector<size_t> one_ring;
  for(size_t i = 0;i < one_ring_tets.size(); ++i){
      if(one_ring_tets[i] != -1) one_ring.push_back(one_ring_tets[i]);
    }
  if(one_ring.front() != tet_ab.second && one_ring.back() == tet_ab.second)
    reverse(one_ring.begin(), one_ring.end());
  assert(one_ring.front() == tet_ab.second && one_ring.back() == tet_ab.first);

  matrix<double> rot = eye<double>(3);
  tba = eye<double>(3);
  for(size_t i = 0 ; i != one_ring.size()-1; ++i){
      get_best_alignment(&frame[one_ring[i]][0], &frame[one_ring[i+1]][0], &rot[0]);
      tba = temp(tba*rot);
    }
}

int hj_boundary(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [usage] hj_boundary tet zyz." << endl;
      return __LINE__;
    }

  jtf::tet_mesh tm(argv[1]);
  matrix<double> zyz;
  if(load_zyz(argv[2], zyz, tm.tetmesh_.mesh_.size(2))){
      cerr << "# [error] wrong zyz." << endl;
      return __LINE__;
    }

  matrix<double> zyz_bkp = zyz;
  hj_frame_alignemt(tm.tetmesh_.mesh_, *tm.fa_, zyz_bkp, zyz);
  matrix<matrix<double> > frame;
  zyz2frame(zyz, frame);

  matrix<double> e(3, tm.outside_face_.size(2));
  matrix<double> eye_m = eye<double>(3);

  vector<pair<double, int> > idx(10);
  for(size_t fi = 0; fi < tm.outside_face_.size(2); ++fi){
      const pair<size_t,size_t> & tet_pair = tm.fa_->face2tet_[tm.outside_face_idx_[fi]];
      assert(tet_pair.first == -1 || tet_pair.second == -1);
      const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);

      for(size_t i = 0; i < 10; ++i){
          idx[i].first = dot(tm.outside_face_normal_(colon(), fi), frame[tet_idx]* type_transition2(i) * eye_m(colon(),2));
          idx[i].second = i;
        }
      sort(idx.begin(), idx.end());
      e(colon(), fi) = type_transition2(idx.back().second) * eye_m(colon(),2);
    }

  vector<size_t> edges;
  vector<double> edge_type;
  matrix<double> tba(3,3);
  matrix<double> Rab(3,3);
  matrix<double> axis(3,1);
  double angle = 0;
  for(size_t ei = 0; ei < tm.ea_outside_->edge2cell_.size(); ++ei){
      const pair<size_t,size_t> & one_edge = tm.ea_outside_->edges_[ei];

      const pair<size_t,size_t> & face_pair = tm.ea_outside_->edge2cell_[ei];
      pair<size_t,size_t> tet_pair;
      {
        const pair<size_t,size_t> & tet_pair0 = tm.fa_->face2tet_[tm.outside_face_idx_[face_pair.first]];
        assert(tet_pair0.first == -1 || tet_pair0.second == -1);
        tet_pair.first = (tet_pair0.first==-1?tet_pair0.second:tet_pair0.first);

        const pair<size_t,size_t> & tet_pair1 = tm.fa_->face2tet_[tm.outside_face_idx_[face_pair.second]];
        assert(tet_pair1.first == -1 || tet_pair1.second == -1);
        tet_pair.second = (tet_pair1.first==-1?tet_pair1.second:tet_pair1.first);
      }
      edges.push_back(one_edge.first);
      edges.push_back(one_edge.second);
      cal_Rab(face_pair, tm.outside_face_, tm.tetmesh_.node_,
              tm.outside_face_normal_(colon(), face_pair.first),
              tm.outside_face_normal_(colon(), face_pair.second), Rab);
      auto it = tm.ortae_.e2t_.find(one_edge);
      assert(it != tm.ortae_.e2t_.end());
      cal_tba(tba, it->second, tet_pair, frame);

      convert_rotation_matrix_to_axis_angle(Rab, axis, angle);

      matrix<double> axis_nanb = cross(e(colon(), face_pair.first), e(colon(), face_pair.second));
      const double len = norm(axis_nanb);
      if(len > 1e-6)
        axis_nanb /= norm(axis_nanb);
      else{
          size_t zi = 0;
         for(zi = 0; zi < 3; ++zi){
             if(fabs(e(zi,face_pair.first)) < 1e-6) break;
           }
         axis_nanb[zi] = 1;
         axis_nanb[(zi+1)%3] = 0;
         axis_nanb[(zi+2)%3] = 0;
        }
      from_angle_to_rotation_matrix(angle, axis_nanb, Rab);

      edge_type.push_back(dot(Rab*e(colon(), face_pair.first),
                              trans(tba)*e(colon(), face_pair.second)));
    }

  ofstream ofs("line.vtk");
  line2vtk(ofs, &tm.tetmesh_.node_[0], tm.tetmesh_.node_.size(2), &edges[0], edges.size()/2);
  cell_data(ofs, &edge_type[0], edge_type.size(), "dot");
  return 0;
}
