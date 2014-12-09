#include <stack>
#include <numeric>
#include <algorithm>
#include <sys/resource.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>

#include <zjucad/optimizer/optimizer.h>
#include <zjucad/ptree/ptree.h>
#include <zjucad/matrix/io.h>
#include <hjlib/function/func_aux.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include "common.h"
#include "hex_param.h"
#include "cut_tet.h"
#include "singularity_adjustment.h"

#include "../tetmesh/hex_io.h"
#include "../common/zyz.h"
#include "../common/IO.h"
#include "../common/vtk.h"
#include "../common/util.h"
#include "../common/transition.h"
#include "../common/transition_type.h"
#include "../common/timer.h"
#include "../numeric/util.h"
#include "topology_analysis.h"
#include "io.h"
#include "../common/visualize_tool.h"
#include "find_singularities.h"
//#include "../utils/line2cylinder.h"


using namespace std;
using namespace zjucad::matrix;
int find_edge_adjacent_faces(matrixst &edge_adjacent_faces,
                             const matrixst &outside_face,
                             const pair<size_t,size_t> &edge);
int label_surface_type(matrixst &surface_type,
                       const jtf::mesh::face2tet_adjacent &fa,
                       const matrixst &outside_face,
                       const matrixst &outside_face_idx,
                       const matrixst &tet,
                       const matrixd &node,
                       const matrixd &zyz)
{
  matrix<matrixd > rot_frame_at_face(outside_face.size(2));
  matrixst outside_face_belong_to_tet_idx(outside_face.size(2));

  { // assign the tet frame to outside_face
    for(size_t fi = 0; fi < outside_face.size(2); ++fi){
        const pair<size_t,size_t> &tet_adj = fa.face2tet_[outside_face_idx[fi]];
        outside_face_belong_to_tet_idx[fi] = (tet_adj.first == -1)? tet_adj.second:tet_adj.first;
      }

    for(size_t i = 0; i < outside_face.size(2); ++i) {
        rot_frame_at_face[i].resize(3,3);
        zyz_angle_2_rotation_matrix1(&zyz(0, outside_face_belong_to_tet_idx[i]), &rot_frame_at_face[i][0]);
      }
  }

  matrixd stander_axis(3,6);
  {
    const matrixd I = eye<double>(3);
    const matrixd& x = I(colon(),0);
    const matrixd& y = I(colon(),1);
    const matrixd& z = I(colon(),2);

    stander_axis(colon(),0) = x;
    stander_axis(colon(),1) = -x;
    stander_axis(colon(),2) = y;
    stander_axis(colon(),3) = -y;
    stander_axis(colon(),4) = z;
    stander_axis(colon(),5) = -z;
  }

  matrixd normal = zeros<double>(3,3);
  cerr << "# outside face size " << outside_face.size(1) << " " << outside_face.size(2) << endl;
  cerr << rot_frame_at_face.size() << endl;


  // ofstream ofs("test_original_surface_type");

  for(size_t t = 0; t < outside_face.size(2); ++t){
      //if(t == 4167)
      //  cerr << endl;
      const pair<size_t,size_t> &tet_pair = fa.face2tet_[outside_face_idx[t]];
      const size_t t_idx = (tet_pair.first == -1)?tet_pair.second:tet_pair.first;

      matrixst face(3);
      copy(fa.faces_[outside_face_idx[t]].begin(),fa.faces_[outside_face_idx[t]].end(),face.begin());

      if(jtf::tetmesh::calculate_face_normal(tet,node,t_idx,face,normal)) {
          cerr << "# strange degenerate face." << endl;
          continue;
        }
      //cerr << "t " << t << normal << endl;
      vector<pair<double,size_t> > resi_label(6);
      for(size_t li = 0; li < 6; ++li) {
          //resi_label[li]= make_pair(dot(normal,rot_frame_at_face[t]
          //* stander_axis(colon(),li)),li);
          resi_label[li]= make_pair(dot(normal,((li%2==0)?1:-1)*rot_frame_at_face[t](colon(),li/2)),li);
        }
      sort(resi_label.begin(),resi_label.end());
      surface_type[t] = resi_label[5].second;
#if 0 // test
      ofs << t_idx << " " << surface_type[t] << endl;
      ofs << normal << endl;
#endif
    }
  return 0;
}

// label the surface type accroding the inner jump face type
int label_surface_type_new(const matrixst &outside_face,
                           const matrixst &tet,
                           const matrixd &node,
                           const jtf::mesh::face2tet_adjacent &fa,
                           const map<pair<size_t,size_t>,size_t> &jump_type_between_tets,
                           matrixst &surf_aligned_type,
                           const matrixst &outside_face_idx,
                           const matrix<matrixd > &frame_inner,
                           const size_t seed_idx)
{
  map<size_t,size_t> surface_type; // store the face idx and its surface type

  matrix<matrixd > test_frame(tet.size(2));
  for(size_t t = 0; t < tet.size(2); ++t) test_frame[t].resize(3,3);

  // first:store the tet_idx; second.first: the current frame; second.second the rotatin frome original frame to current
  stack<pair<size_t,pair<matrixd,matrixd > > > visit_tet;
  visit_tet.push(make_pair(seed_idx,make_pair(frame_inner[seed_idx],eye<double>(3))));
  vector<bool> is_tet_visited(tet.size(2),false);
  matrixst faces_of_one_tet(3,4);
  matrixd face_normal = zeros<double>(3,1);
  vector<pair<double,size_t> > resi_label(6);


  matrixd rot_test(3,3);

  typedef map<pair<size_t,size_t>,size_t >::const_iterator mci;
  while(!visit_tet.empty()){
      pair<size_t,pair<matrixd,matrixd > > top = visit_tet.top();
      const size_t &tet_idx = top.first;
      const matrixd & frame_at_this_tet = top.second.first;
      const matrixd & total_rotation = top.second.second;
      //const matrixd total_rotation = eye<double>(3);
      visit_tet.pop();
      test_frame[tet_idx] = frame_at_this_tet;
      is_tet_visited[tet_idx] = true;

#if 0 // check the aliged frame
      get_best_alignment(&frame_inner[seed_idx](0,0),
                         &frame_at_this_tet(0,0),&rot_test(0,0));
          if(norm(rot_test - eye<double>(3)) > 1e-8){
        cerr << "# error not aligned." << endl;
      }
#endif
      // get the four face of tet
      faces_of_one_tet(0,0) = tet(0,tet_idx);
      faces_of_one_tet(1,0) = tet(1,tet_idx);
      faces_of_one_tet(2,0) = tet(2,tet_idx);

      faces_of_one_tet(0,1) = tet(0,tet_idx);
      faces_of_one_tet(1,1) = tet(1,tet_idx);
      faces_of_one_tet(2,1) = tet(3,tet_idx);

      faces_of_one_tet(0,2) = tet(0,tet_idx);
      faces_of_one_tet(1,2) = tet(2,tet_idx);
      faces_of_one_tet(2,2) = tet(3,tet_idx);

      faces_of_one_tet(0,3) = tet(1,tet_idx);
      faces_of_one_tet(1,3) = tet(2,tet_idx);
      faces_of_one_tet(2,3) = tet(3,tet_idx);

      for(size_t t = 0;t < 4; ++t){ // for each face
          const size_t face_idx = fa.get_face_idx(&faces_of_one_tet(0,t));
          const pair<size_t,size_t> &tet_pair = fa.face2tet_[face_idx];
          if(tet_pair.first == -1 || tet_pair.second == -1){ // it's a surface face, need to align the surface type
              jtf::tetmesh::calculate_face_normal(tet,node,tet_idx,faces_of_one_tet(colon(),t),face_normal);

              for(size_t li = 0; li < 6; ++li) {
                  resi_label[li]= make_pair(dot(face_normal,((li%2==0)?1:-1)*frame_at_this_tet(colon(),li/2)),li);
                }
              sort(resi_label.begin(),resi_label.end());
              surface_type[face_idx] = resi_label.back().second;
            }else{
              const size_t other_tet_idx = (tet_pair.first != tet_idx)?tet_pair.first:tet_pair.second;
              mci mci_ = jump_type_between_tets.find(make_pair(other_tet_idx,tet_idx));

              if(is_tet_visited[other_tet_idx]) {// this tet already been visited,check whether the type is correct
                  matrixd test = zeros<double>(3,1);
                  if(mci_ != jump_type_between_tets.end()){
                      test = temp(frame_inner[other_tet_idx] * type_transition2(mci_->second) ) * total_rotation ;
                    }else
                    test = frame_inner[other_tet_idx] *  total_rotation;
                  if(norm(test - test_frame[other_tet_idx]) > 1e-8){
                      cerr << "# error: this frame is set to be other type." << endl;
                    }
                  continue;
                }

              if(mci_ != jump_type_between_tets.end()){ // need a jump
                  visit_tet.push(make_pair(other_tet_idx,make_pair(temp(frame_inner[other_tet_idx] * type_transition2(mci_->second)) * total_rotation,
                                                                   type_transition2(mci_->second) * total_rotation  )));
                  test_frame[other_tet_idx] = temp(frame_inner[other_tet_idx] * type_transition2(mci_->second)) * total_rotation;

                }else{ // do not need a jump
                  visit_tet.push(make_pair(other_tet_idx,make_pair(frame_inner[other_tet_idx] * total_rotation,
                                                                   total_rotation)));
                  test_frame[other_tet_idx] = frame_inner[other_tet_idx] * total_rotation;
                }
            }
        }
    }

  {// record the surface type
    if(count(is_tet_visited.begin(),is_tet_visited.end(),false) != 0) cerr << "# have not visist all tets." << endl;
    typedef map<size_t,size_t>::const_iterator mi;
    for(size_t t = 0; t < outside_face_idx.size(); ++t){
        mi mi_ = surface_type.find(outside_face_idx[t]);
        if(mi_ == surface_type.end()) {
            cerr << "# strange can not find this outside face." << endl;
            continue;
          }else{
            surf_aligned_type[t] = mi_->second;
          }
      }
  }

#if 1
  {// to vtk
    ofstream ofs("surface-aligned-frame-residual.vtk");
    vector<vector<size_t> > tet_adjacent_tet;
    tet_adjacent_tet.resize(tet.size(2));
    const size_t face_num = fa.faces_.size();
    for(size_t fi = 0; fi < face_num; ++fi)
      {
        const size_t vi = fa.face2tet_[fi].first;
        const size_t vj = fa.face2tet_[fi].second;
        if(vi != -1 && vj != -1)
          {
            tet_adjacent_tet[vi].push_back(vj);
            tet_adjacent_tet[vj].push_back(vi);
          }
      }
    matrixd residual(tet.size(2),1);

    for(size_t vi = 0; vi < tet.size(2); ++vi)
      {
        double residual_each_tet = 0.0;
        for(size_t advi = 0; advi < tet_adjacent_tet[vi].size(); ++advi)    {
            residual_each_tet += norm(test_frame[vi] - test_frame[tet_adjacent_tet[vi][advi]]);
          }
        if(tet_adjacent_tet[vi].size())
          residual[vi] = residual_each_tet / tet_adjacent_tet[vi].size();
      }

    tet2vtk(ofs,&node[0],node.size(2),&tet[0],tet.size(2));
    cell_data(ofs,&residual[0],residual.size(),"new_aligned_frame");
  }
#endif
  return 0;
}

size_t get_best_aligned_axis(const matrixd &normal,
                             const matrixd &frame)
{
  vector<pair<double,size_t> > resi_label(6);
  for(size_t li = 0; li < 6; ++li) {
      resi_label[li]= make_pair(dot(normal,((li%2==0)?1:-1)*frame(colon(),li/2)),li);
    }
  sort(resi_label.begin(),resi_label.end());
  return resi_label.back().second;
}

int find_outside_face_of_tet_with_edge(const matrixst &tet,
                                       const jtf::mesh::face2tet_adjacent &fa,
                                       matrixst &face,
                                       const pair<size_t,size_t> &edge)
{
  assert(tet.size() == 4);
  matrixst faces_of_one_tet(3,4);
  faces_of_one_tet(0,0) = tet[0];
  faces_of_one_tet(1,0) = tet[1];
  faces_of_one_tet(2,0) = tet[2];

  faces_of_one_tet(0,1) = tet[0];
  faces_of_one_tet(1,1) = tet[1];
  faces_of_one_tet(2,1) = tet[3];

  faces_of_one_tet(0,2) = tet[0];
  faces_of_one_tet(1,2) = tet[2];
  faces_of_one_tet(2,2) = tet[3];

  faces_of_one_tet(0,3) = tet[1];
  faces_of_one_tet(1,3) = tet[2];
  faces_of_one_tet(2,3) = tet[3];

  for(size_t t = 0; t < 4; ++t){
      if(fa.is_outside_face(fa.face2tet_[fa.get_face_idx(&faces_of_one_tet(0,t))])
         && find(faces_of_one_tet(colon(),t).begin(),faces_of_one_tet(colon(),t).end(),edge.first) != faces_of_one_tet(colon(),t).end()
         && find(faces_of_one_tet(colon(),t).begin(),faces_of_one_tet(colon(),t).end(),edge.second) != faces_of_one_tet(colon(),t).end()){
          face = faces_of_one_tet(colon(),t);
          return 0;
        }
    }
  return __LINE__;
}

int find_outside_face_of_tet(const matrixst &tet,
                             const jtf::mesh::face2tet_adjacent &fa,
                             matrixst &face)
{
  assert(tet.size() == 4);
  matrixst faces_of_one_tet(3,4);
  faces_of_one_tet(0,0) = tet[0];
  faces_of_one_tet(1,0) = tet[1];
  faces_of_one_tet(2,0) = tet[2];

  faces_of_one_tet(0,1) = tet[0];
  faces_of_one_tet(1,1) = tet[1];
  faces_of_one_tet(2,1) = tet[3];

  faces_of_one_tet(0,2) = tet[0];
  faces_of_one_tet(1,2) = tet[2];
  faces_of_one_tet(2,2) = tet[3];

  faces_of_one_tet(0,3) = tet[1];
  faces_of_one_tet(1,3) = tet[2];
  faces_of_one_tet(2,3) = tet[3];

  for(size_t t = 0; t < 4; ++t){
      if(fa.is_outside_face(fa.face2tet_[fa.get_face_idx(&faces_of_one_tet(0,t))])){
          face = faces_of_one_tet(colon(),t);
          return 0;
        }
    }
  return __LINE__;
}

int draw_surface_singularity_edge(
    const matrixst &outside_face,
    const matrixst &tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent &fa,
    const map<pair<size_t,size_t>,size_t> &jump_type_between_tets,
    matrixst &surf_aligned_type,
    const matrixst &outside_face_idx,
    const matrix<matrixd > &frame_inner,
    const jtf::mesh::one_ring_tet_at_edge &ortae)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oci;
  typedef map<pair<size_t,size_t>,size_t>::const_iterator mci;
  matrixd rot = zeros<double>(3,3);
  matrixd normal = zeros<double>(3,1);
  matrixst face = zeros<size_t>(3,1);
  vector<size_t> surface_singularity_edge;

  for(oci ci = ortae.e2t_.begin(); ci != ortae.e2t_.end(); ++ci){
      const vector<size_t> &loop = ci->second;
      if(loop.front() != loop.back() && loop.size() != 2){
          cerr << "# error edge. " << endl;
          continue;
        }
      if(loop.size() == 2){ // sharp edge
#if 0
          surface_singularity_edge.push_back(ci->first.first);
          surface_singularity_edge.push_back(ci->first.second);
#endif
          continue;
        }
      if(loop.front() == -1 && loop.back() == -1){ // it's a surface edge
          rot = eye<double>(3);
          if(find_outside_face_of_tet(tet(colon(),loop[1]),fa,face)){
              cerr << "# can not find outside face on this tet." << endl;
              return __LINE__;
            }
          jtf::tetmesh::calculate_face_normal(tet,node,loop[1],face,normal);
          const size_t axis_type_1 = get_best_aligned_axis(normal,frame_inner[loop[1]]);

          for(size_t t = loop.size() - 2; t != 1; --t){
              mci it = jump_type_between_tets.find(make_pair(loop[t],loop[t-1]));
              if(it == jump_type_between_tets.end()) continue;
              rot = temp(rot * type_transition2(it->second));
            }

          if(find_outside_face_of_tet(tet(colon(),loop[loop.size() - 2]),fa,face)){
              cerr << "# can not find outside face on this tet." << endl;
              return __LINE__;
            }
          jtf::tetmesh::calculate_face_normal(tet,node,loop[loop.size() - 2],face,normal);
          const size_t axis_type_2 = get_best_aligned_axis(normal,frame_inner[loop[loop.size() - 2]] * rot);
          if((axis_type_2 != axis_type_1 )
             /*  && (axis_type_2/2 == axis_type_1/2) */){
              surface_singularity_edge.push_back(ci->first.first);
              surface_singularity_edge.push_back(ci->first.second);
            }
        }
    }

  { // to vtk
    ofstream ofs("surface_singualrity_edge.vtk");
    line2vtk(ofs,&node[0],node.size(2),&surface_singularity_edge[0],surface_singularity_edge.size()/2);
  }
  return 0;

}

int draw_surface_singularity_edge_no_jump_type(
    const matrixst &outside_face,
    const matrixst &tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent &fa,
    const map<pair<size_t,size_t>,size_t> &jump_type_between_tets,
    matrixst &surf_aligned_type,
    const matrixst &outside_face_idx,
    const matrix<matrixd > &frame_inner,
    const jtf::mesh::one_ring_tet_at_edge &ortae)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oci;
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mci;
  matrixd rot = zeros<double>(3,3);
  matrixd normal = zeros<double>(3,1);
  matrixst face = zeros<size_t>(3,1);
  vector<size_t> surface_singularity_edge;

  for(oci ci = ortae.e2t_.begin(); ci != ortae.e2t_.end(); ++ci){
      const vector<size_t> &loop = ci->second;
      if(loop.front() != loop.back() && loop.size() != 2){
          cerr << "# error edge. " << endl;
          continue;
        }
      if(loop.size() == 2){ // sharp edge
#if 0
          surface_singularity_edge.push_back(ci->first.first);
          surface_singularity_edge.push_back(ci->first.second);
#endif
          continue;
        }
      if(loop.front() == -1 && loop.back() == -1){ // it's a surface edge
          rot = eye<double>(3);
          if(find_outside_face_of_tet(tet(colon(),loop[1]),fa,face)){
              cerr << "# can not find outside face on this tet." << endl;
              return __LINE__;
            }
          jtf::tetmesh::calculate_face_normal(tet,node,loop[1],face,normal);
          const size_t axis_type_1 = get_best_aligned_axis(normal,frame_inner[loop[1]]);

          matrixd rot_temp = zeros<double>(3,3);
          for(size_t t = loop.size() - 2; t != 1; --t){
              //                mci it = jump_type_between_tets.find(make_pair(loop[t],loop[t-1]));
              //                if(it == jump_type_between_tets.end()) continue;
              //                rot = temp(rot * type_transition2(it->second));
              get_best_alignment(&frame_inner[loop[t]][0],
                  &frame_inner[loop[t-1]][0],
                  &rot_temp[0]);
              rot = temp(rot * rot_temp);
            }

          if(find_outside_face_of_tet(tet(colon(),loop[loop.size() - 2]),fa,face)){
              cerr << "# can not find outside face on this tet." << endl;
              return __LINE__;
            }
          jtf::tetmesh::calculate_face_normal(tet,node,loop[loop.size() - 2],face,normal);
          const size_t axis_type_2 = get_best_aligned_axis(normal,frame_inner[loop[loop.size() - 2]] * rot);
          if((axis_type_2 != axis_type_1 )
             && (axis_type_2/2 == axis_type_1/2) ){
              surface_singularity_edge.push_back(ci->first.first);
              surface_singularity_edge.push_back(ci->first.second);
            }
        }
    }

  { // to vtk
    ofstream ofs("surface_singualrity_edge_no_jump_type.vtk");
    line2vtk(ofs,&node[0],node.size(2),&surface_singularity_edge[0],surface_singularity_edge.size()/2);
  }
  return 0;
}

int rotationFromTo(const matrixd& from,
                   const matrixd &to,
                   matrixd &rot)
{

  const matrixd from_ = from / norm(from);
  const matrixd to_ = to / norm(to);

  const double d = dot(from_,to_);

  if(d - 1 > 1e-8){
      rot = eye<double>(3);
      return 0;
    }
  const double s = sqrtf( (1+d)*2 ); // optimize inv_sqrt
  const double invs = 1.f / s;
  const matrixd c = cross(from_,to_) * invs;// v0.crossProduct(v1)*invs;

  matrixd quaternion(4,1); // x,y,z,w

  quaternion[0] = c[0];
  quaternion[1] = c[1];
  quaternion[2] = c[2];
  quaternion[3] = s * 0.5;

  const double &X = quaternion[0];
  const double &Y = quaternion[1];
  const double &Z = quaternion[2];
  const double &W = quaternion[3];

  rot(0,0) = 1.0f - 2.0f*Y*Y - 2.0f*Z*Z;
  rot(1,0) = 2.0f*X*Y + 2.0f*Z*W;
  rot(2,0) = 2.0f*X*Z - 2.0f*Y*W;

  rot(0,1) = 2.0f*X*Y - 2.0f*Z*W;;
  rot(1,1) = 1.0f - 2.0f*X*X - 2.0f*Z*Z;
  rot(2,1) = 2.0f*Z*Y + 2.0f*X*W;

  rot(0,2) = 2.0f*X*Z + 2.0f*Y*W;
  rot(1,2) = 2.0f*Z*Y - 2.0f*X*W;
  rot(2,2) = 1.0f - 2.0f*X*X - 2.0f*Y*Y;

  return 0;
}
int get_best_alignment_normal_and_tangent(const matrixd &frame_from,
                                          const matrixd &normal_from,
                                          const matrixd &frame_to,
                                          const matrixd &normal_to,
                                          matrixd &rot,
                                          const matrixd &axis)
{
  matrixd Re = eye<double>(3);

  rotationFromTo(normal_to,normal_from,Re);

  size_t normal_align_axis_type_from = get_best_aligned_axis(normal_from,frame_from);
  vector<pair<double,size_t> > align_normal_type;
  for(size_t t = 0; t < 24; ++t){
      const matrixd new_frame_to = Re * frame_to * type_transition2(t);
      if(get_best_aligned_axis(normal_to, new_frame_to) == normal_align_axis_type_from){
          align_normal_type.push_back(make_pair(norm(new_frame_to - frame_from),t));
        }
    }
  sort(align_normal_type.begin(),align_normal_type.end());
  rot = type_transition2(align_normal_type.front().second);
  return 0;
}

int draw_surface_singularity_edge_normal_and_tangent(
    const matrixst &outside_face,
    const matrixst &tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent &fa,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &jump_type_between_tets,
    matrixst &surf_aligned_type,
    const matrixst &outside_face_idx,
    const matrix<matrixd > &frame_inner,
    const jtf::mesh::one_ring_tet_at_edge &ortae)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oci;
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mci;
  matrixd rot = zeros<double>(3,3);
  //matrixd normal = zeros<double>(3,1);
  //matrixst face = zeros<size_t>(3,1);
  vector<size_t> surface_singularity_edge;

  for(oci ci = ortae.e2t_.begin(); ci != ortae.e2t_.end(); ++ci){
      const vector<size_t> &loop = ci->second;
      if(loop.front() != loop.back() && loop.size() != 2){
          cerr << "# error edge. " << endl;
          continue;
        }
      if(loop.size() == 2){ // sharp edge
#if 0
          surface_singularity_edge.push_back(ci->first.first);
          surface_singularity_edge.push_back(ci->first.second);
#endif
          continue;
        }
      if(loop.front() == -1 && loop.back() == -1){ // it's a surface edge
          rot = eye<double>(3);
          matrixst face_0 = zeros<size_t>(3,1);
          matrixst face_1 = zeros<size_t>(3,1);

          if(find_outside_face_of_tet(tet(colon(),loop[1]),fa,face_0)){
              cerr << "# can not find outside face on this tet." << endl;
              return __LINE__;
            }

          if(find_outside_face_of_tet(tet(colon(),loop[loop.size() - 2]),fa,face_1)){
              cerr << "# can not find outside face on this tet." << endl;
              return __LINE__;
            }

          matrixd face_normal_0 = zeros<double>(3,1);
          matrixd face_normal_1 = zeros<double>(3,1);

          jtf::tetmesh::calculate_face_normal(tet,node,loop[1],face_0,face_normal_0);
          jtf::tetmesh::calculate_face_normal(tet,node,loop[loop.size() - 2],face_1,face_normal_1);

          for(size_t t = loop.size() - 2; t != 1; --t){
              mci it = jump_type_between_tets.find(make_pair(loop[t],loop[t-1]));
              if(it == jump_type_between_tets.end()) continue;
              rot = temp(rot * type_transition2(it->second));
            }

          matrixd jump_rot = zeros<double>(3,1);
          matrixd current_edge = node(colon(),ci->first.second) - node(colon(),ci->first.first);
          get_best_alignment_normal_and_tangent(frame_inner[loop[loop.size() - 2]] * rot, face_normal_1,
              frame_inner[loop[1]],face_normal_0,
              jump_rot,
              current_edge);
          if(norm(  jump_rot - eye<double>(3)) > 1e-8 ){
              surface_singularity_edge.push_back(ci->first.first);
              surface_singularity_edge.push_back(ci->first.second);
            }
        }
    }

  { // to vtk
    ofstream ofs("surface_singualrity_edge_norm&tangent.vtk");
    line2vtk(ofs,&node[0],node.size(2),&surface_singularity_edge[0],surface_singularity_edge.size()/2);
  }
  return 0;

}

int label_surface_type_by_connection(
    const matrixst &outside_face,
    const matrixst &tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent &fa,
    const map<pair<size_t,size_t>,size_t> &jump_type_between_tets,
    matrixst &surf_aligned_type,
    const matrixst &outside_face_idx,
    const matrix<matrixd > &frame_inner,
    const jtf::mesh::one_ring_tet_at_edge &ortae)
{
  // store the outsideface_idx,tet_idx and the matrix aligned in face_idx tet, and rotation from the seed
  stack<pair<pair<size_t,size_t>,
      pair<matrixd,matrixd >
      > > visit_tet;
  map<size_t,size_t> surface_type;
  pair<pair<size_t,size_t>,
      pair<matrixd,matrixd >
      > seed;
  const pair<size_t,size_t> &tet_pair = fa.face2tet_[outside_face_idx[0]];
  matrix<matrixd > aligned_frame(tet.size(2));
  for(size_t t = 0; t < tet.size(2); ++t) {
      aligned_frame[t].resize(3,3);
      aligned_frame[t] = frame_inner[t];
    }

  seed.first.first = outside_face_idx[0];
  seed.first.second = (tet_pair.first == -1)?tet_pair.second:tet_pair.first;
  seed.second.first = frame_inner[seed.first.second];
  seed.second.second = eye<double>(3);
  vector<bool> is_face_visited(fa.faces_.size(),false);
  visit_tet.push(seed);


  vector<size_t> error_face;

  matrixst edges_of_one_face(2,3);
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator ocit;
  typedef map<pair<size_t,size_t>,size_t>::const_iterator mcit;
  matrixd face_nromal = zeros<double>(3,1);

  while(!visit_tet.empty()){
      const pair<pair<size_t,size_t>,pair<matrixd,matrixd > > top = visit_tet.top();
      visit_tet.pop();
      is_face_visited[top.first.first] = true;
      aligned_frame[top.first.second] = top.second.first;

      assert(fa.is_outside_face(fa.face2tet_[top.first.first]));
      const vector<size_t> &current_face = fa.faces_[top.first.first];
      for(size_t t = 0 ; t < 3; ++t){
          edges_of_one_face(0,t) =  current_face[t];
          edges_of_one_face(1,t) =  current_face[(t+1)%3];
        }

      for(size_t t = 0; t < 3; ++t){
          ocit it = ortae.e2t_.find(make_pair(edges_of_one_face(0,t),
                                              edges_of_one_face(1,t)));
          if(it == ortae.e2t_.end()) it = ortae.e2t_.find(make_pair(
                                                            edges_of_one_face(1,t),
                                                            edges_of_one_face(0,t)));
          if(it == ortae.e2t_.end()){
              cerr << "# strange can not find this edge." << endl;
              continue;
            }else{ // find this edge
              const vector<size_t> & loop = it->second;
              if(loop.front() != loop.back() && loop.size() != 2){
                  cerr << "# error edge: " << edges_of_one_face(0,t) << " " << edges_of_one_face(1,t) << endl;
                  continue;
                }
              if(loop.size() == 2){ // sharp edge, which will contain at least two outside face on one tet
                  const size_t& tet_idx = (loop[0] == -1)?loop[1]:loop[0];
                  assert(tet_idx == top.first.second);
                  size_t face_need_to_handle;
                  { // to find an adjacent face with current face
                    matrixst two_outside_face(3,2);
                    two_outside_face(0,0) = two_outside_face(0,1) = edges_of_one_face(0,t);
                    two_outside_face(1,0) = two_outside_face(1,1) = edges_of_one_face(1,t);

                    for(size_t i = 0,j = 0; i < 4; ++i){
                        if(tet(i,tet_idx) != edges_of_one_face(0,t)
                           && tet(i,tet_idx) != edges_of_one_face(1,t)){
                            two_outside_face(2,j) = tet(i,tet_idx);
                            j++;
                          }
                      }

                    const size_t face_idx0 = fa.get_face_idx(&two_outside_face(0,0));
                    const size_t face_idx1 = fa.get_face_idx(&two_outside_face(0,1));

                    //assert(top.first.first == face_idx0 || top.first.first == face_idx1);
                    if(fa.is_outside_face(fa.face2tet_[face_idx0])) face_need_to_handle = face_idx0;
                    else
                      face_need_to_handle = face_idx1;

                  }

                  if(is_face_visited[face_need_to_handle]){ // has been visited, need to check
                      jtf::tetmesh::calculate_face_normal_wrapper(tet,node,tet_idx,&fa.faces_[face_need_to_handle][0],face_nromal);
                      const size_t face_type = get_best_aligned_axis(face_nromal,aligned_frame[tet_idx]);
                      if(face_type != surface_type[face_need_to_handle]){
                          cerr << "# conflict face type assign. " << __LINE__ << endl;
#if 1
                          for(size_t j = 0; j < 3; ++j){
                              error_face.push_back(fa.faces_[face_need_to_handle][j]);
                            }
#endif
                          continue;
                        }
                    }else
                    {
                      jtf::tetmesh::calculate_face_normal_wrapper(tet,node,tet_idx,&fa.faces_[face_need_to_handle][0],face_nromal);
                      const size_t face_type = get_best_aligned_axis(face_nromal,aligned_frame[tet_idx]);
                      //is_face_visited[face_need_to_handle] = true;
                      surface_type[face_need_to_handle] = face_type;

                      visit_tet.push(make_pair(
                                       make_pair(face_need_to_handle,tet_idx),
                                       make_pair(aligned_frame[tet_idx],top.second.second))
                                     );
                    }
                } // end sharp edge
              else{ // normal outside edge, loop is : -1 t0 t1 ... tn -1
                  assert(loop.front() == -1 && loop.back() == -1);
                  matrixst surface_tet_idx(2);
                  surface_tet_idx[0] = loop[1];
                  surface_tet_idx[1] = loop[loop.size() - 2];

                  const size_t& current_tet = top.first.second;

                  if(current_tet == surface_tet_idx[0]){ // the other tet is at the end of loop
                      matrixd rot = eye<double>(3);
                      for(size_t k = loop.size()-2; k != 1; --k){
                          mcit cit = jump_type_between_tets.find(make_pair(loop[k],loop[k-1]));
                          if(cit == jump_type_between_tets.end()) continue; // this jump type is identity
                          rot = temp(rot * type_transition2(cit->second));
                        }
                      rot = temp(rot * top.second.second); // apply the rotation frome current tet to seed tet
                      matrixd new_aligned_frame = frame_inner[surface_tet_idx[1]] * rot; // align to the seed tet
                      // to find the outside face of the other tet

                      size_t face_need_to_handle;
                      { // to find an adjacent face with current face
                        matrixst two_outside_face(3,2);
                        two_outside_face(0,0) = two_outside_face(0,1) = edges_of_one_face(0,t);
                        two_outside_face(1,0) = two_outside_face(1,1) = edges_of_one_face(1,t);

                        for(size_t i = 0,j = 0; i < 4; ++i){
                            if(tet(i,surface_tet_idx[1]) != edges_of_one_face(0,t)
                               && tet(i,surface_tet_idx[1]) != edges_of_one_face(1,t)){
                                two_outside_face(2,j) = tet(i,surface_tet_idx[1]);
                                j++;
                              }
                          }
                        const size_t face_idx0 = fa.get_face_idx(&two_outside_face(0,0));
                        const size_t face_idx1 = fa.get_face_idx(&two_outside_face(0,1));

                        //assert(top.first.first == face_idx0 || top.first.first == face_idx1);
                        if(fa.is_outside_face(fa.face2tet_[face_idx0])) face_need_to_handle = face_idx0;
                        else
                          face_need_to_handle = face_idx1;
                      }

                      const vector<size_t> &face_ = fa.faces_[face_need_to_handle];
                      if(is_face_visited[face_need_to_handle]) { // need to check whether this face is set correct
                          jtf::tetmesh::calculate_face_normal_wrapper(tet,node,surface_tet_idx[1],&face_[0],face_nromal);
                          const size_t face_type = get_best_aligned_axis(face_nromal,new_aligned_frame);
                          if(face_type != surface_type[face_need_to_handle]){
                              cerr << "# conflict face type assign. " << __LINE__ << endl;
#if 1
                              for(size_t j = 0; j < 3; ++j){
                                  error_face.push_back(fa.faces_[face_need_to_handle][j]);
                                }
#endif
                              continue;
                            }
                        }else{
                          jtf::tetmesh::calculate_face_normal_wrapper(tet,node,surface_tet_idx[1],&face_[0],face_nromal);
                          const size_t face_type = get_best_aligned_axis(face_nromal,new_aligned_frame);
                          //is_face_visited[face_idx_] = true;
                          surface_type[face_need_to_handle] = face_type;
                          //aligned_frame[ surface_tet_idx[1]] = new_aligned_frame;
                          visit_tet.push(make_pair(
                                           make_pair(face_need_to_handle, surface_tet_idx[1]),
                                         make_pair(new_aligned_frame,rot)));
                        }

                    }else // current_tet == surface_tet_idx[1]
                    {
                      matrixd rot = eye<double>(3);

                      for(size_t k = 1; k != loop.size() - 2; ++k){
                          mcit cit = jump_type_between_tets.find(make_pair(loop[k],loop[k+1]));
                          if(cit == jump_type_between_tets.end()) continue; // this jump type is identity
                          rot = temp(rot * type_transition2(cit->second));
                        }
                      rot = temp(rot * top.second.second);

                      matrixd new_aligned_frame = frame_inner[surface_tet_idx[0]] * rot;
                      // to find the outside face of the other tet
                      size_t face_need_to_handle;
                      { // to find an adjacent face with current face
                        matrixst two_outside_face(3,2);
                        two_outside_face(0,0) = two_outside_face(0,1) = edges_of_one_face(0,t);
                        two_outside_face(1,0) = two_outside_face(1,1) = edges_of_one_face(1,t);

                        for(size_t i = 0,j = 0; i < 4; ++i){
                            if(tet(i,surface_tet_idx[1]) != edges_of_one_face(0,t)
                               && tet(i,surface_tet_idx[1]) != edges_of_one_face(1,t)){
                                two_outside_face(2,j) = tet(i,surface_tet_idx[1]);
                                j++;
                              }
                          }
                        const size_t face_idx0 = fa.get_face_idx(&two_outside_face(0,0));
                        const size_t face_idx1 = fa.get_face_idx(&two_outside_face(0,1));

                        //assert(top.first.first == face_idx0 || top.first.first == face_idx1);
                        if(fa.is_outside_face(fa.face2tet_[face_idx0])) face_need_to_handle = face_idx0;
                        else
                          face_need_to_handle = face_idx1;
                      }

                      const vector<size_t> &face_ = fa.faces_[face_need_to_handle];

                      if(is_face_visited[face_need_to_handle]) { // need to check whether this face is set correct
                          jtf::tetmesh::calculate_face_normal_wrapper(tet,node,surface_tet_idx[0],&face_[0],face_nromal);
                          const size_t face_type = get_best_aligned_axis(face_nromal,new_aligned_frame);
                          if(face_type != surface_type[face_need_to_handle]){
                              cerr << "# conflict face type assign. " <<  __LINE__  << endl;
#if 1
                              for(size_t j = 0; j < 3; ++j){
                                  error_face.push_back(fa.faces_[face_need_to_handle][j]);
                                }
#endif
                              continue;
                            }
                        }else{
                          jtf::tetmesh::calculate_face_normal_wrapper(tet,node,surface_tet_idx[0],&face_[0],face_nromal);
                          const size_t face_type = get_best_aligned_axis(face_nromal,new_aligned_frame);
                          //is_face_visited[face_idx_] = true;
                          surface_type[face_need_to_handle] = face_type;
                          visit_tet.push(make_pair(
                                           make_pair(face_need_to_handle,surface_tet_idx[0]),
                                         make_pair(new_aligned_frame,rot)));
                        }
                    }
                } // end normal edge
            }
        }
    }

  {//
    typedef map<size_t,size_t>::const_iterator mci;
    for(size_t t = 0; t < surf_aligned_type.size(); ++t){
        mci ci = surface_type.find(outside_face_idx[t]);
        if(ci == surface_type.end()) {
            cerr << "# strange, can not find this outside face." << endl;
            continue;
          }
        surf_aligned_type[t] = ci->second;
      }
  }

  {// error face
    ofstream ofs("error_face.vtk");
    tri2vtk(ofs,&node[0],node.size(2),&error_face[0],error_face.size()/3);
  }

  return 0;
}

int get_jump_type(
    matrixst &common_face,
    const matrixst &tet,
    const matrixd &node,
    const boost::unordered_map<std::pair<size_t,size_t>, std::vector<pair<size_t,size_t> > > &edge_loop_pair,
    const vector<deque<size_t> > &singularity_chain,
    const vector<deque<size_t> > &singularity_type,
    const size_t t_i,
    const size_t t_j)
{
  typedef deque<size_t>::const_iterator dci;
  typedef deque<pair<size_t,size_t> >::iterator dcip;
  typedef vector<size_t>::const_iterator vci;
  typedef boost::unordered_map<std::pair<size_t,size_t>, std::vector<pair<size_t,size_t> > > e2tet_type_p;
  typedef e2tet_type_p::const_iterator ecip;

  size_t other_vertex_in_ti = -1;
  for(size_t t = 0; t < 4; ++t){
      if(tet(t,t_i) != common_face[0]
         && tet(t,t_i) != common_face[1]
         && tet(t,t_i) != common_face[2])
        {
          other_vertex_in_ti = tet(t,t_i);
          break;
        }
    }


  if(other_vertex_in_ti == -1){
      cerr << "# error tet." << endl;
    }

  matrixd edge20 = node(colon(),common_face[0]) - node(colon(),common_face[2]);
  matrixd edge01 = node(colon(),common_face[1]) - node(colon(),common_face[0]);
  matrixd edge_other = node(colon(),common_face[0]) - node(colon(),other_vertex_in_ti);

  if(dot(edge_other, cross(edge20,edge01)) < -1e-8) swap(common_face[1],common_face[2]); // reverse the face to make the direction di->dj is right hand

  vector<pair<size_t,size_t> > edges;
  edges.push_back(make_pair(common_face[0],common_face[1]));
  edges.push_back(make_pair(common_face[1],common_face[2]));
  edges.push_back(make_pair(common_face[2],common_face[0]));


  //size_t type;
  for(size_t t = 0; t < singularity_chain.size(); ++t){
      const deque<size_t> &chain = singularity_chain[t];
      // TODO: need to speed up
      deque<pair<size_t,size_t> > chain_pair;

      for(size_t i = 0; i < chain.size(); i+=2) {
          chain_pair.push_back(make_pair(chain[i],chain[i+1]));
        }

      bool is_edge_up_down = false;
      for(size_t i = 0; i < 3; ++i){ // for each edge on face
          is_edge_up_down = false;
          const pair<size_t,size_t> &edge = edges[i];

          dcip dit_ = find(chain_pair.begin(),chain_pair.end(),edge);

          if(dit_ == chain_pair.end()) dit_ = find(chain_pair.begin(),chain_pair.end(),make_pair(edge.second,edge.first));
          if(dit_ == chain_pair.end()) // it's not singularity edge
            continue;

          ecip ecip_ = edge_loop_pair.find(edge); // to find the around tet of edge_new

          if(ecip_ == edge_loop_pair.end()) {
              ecip_ = edge_loop_pair.find(make_pair(edge.second,edge.first));
              is_edge_up_down = true;
            }
          if(ecip_ != edge_loop_pair.end()){ // find the edge

              const vector<pair<size_t,size_t> > &v_loop_p = ecip_->second;

              if(find(v_loop_p.begin(),v_loop_p.end(),make_pair(t_i,t_j)) != v_loop_p.end())
                {
                  size_t type_ = singularity_type[t][static_cast<size_t>(dit_ - chain_pair.begin())];
                  if(is_edge_up_down)
                    return type_;
                  else
                    {
                      if(type_ - type_/3 * 3 == 0)
                        type_ += 2;
                      else if(type_ - type_/3 * 3 == 2)
                        type_ -= 2;
                      return type_;
                    }
                }
              if(find(v_loop_p.begin(),v_loop_p.end(),make_pair(t_j,t_i)) != v_loop_p.end()) // need to revers type 0-->2 3-->5 6-->8
                {
                  size_t type_ = singularity_type[t][static_cast<size_t>(dit_ - chain_pair.begin())];
                  if(!is_edge_up_down)
                    return type_;
                  else
                    {
                      if(type_ - type_/3 * 3 == 0) type_ += 2;
                      else if(type_ - type_/3 * 3 == 2) type_ -= 2;
                      return type_;
                    }
                }
            }else
            {
              cerr << "# strange, it's an edge, but can not find it in map." <<endl;
            }
        }
    }
  return 10; // means no need jump
}

int find_adjacent_face(const matrixst &outside_face_cut,
                       const vector<bool> &is_outside_face_cut_visited,
                       const vector<size_t> &jump_face,
                       const size_t face_idx,
                       vector<size_t> &tmp_adjacent_face,
                       const vector<pair<size_t,size_t> > &singularity_edges,
                       const matrixst &cut_tet2tet) // notice: the vertex index store in singularity is the original tet idx
{
  tmp_adjacent_face.clear();
  vector<pair<size_t,size_t> > need_to_search_edges;
  const matrixst &face = outside_face_cut(colon(),face_idx);
  matrixst face_idx_original(3,1);

  for(size_t t = 0;t < 3; ++t){
      face_idx_original[t] = cut_tet2tet[face[t]];
    }

  for(size_t i = 0; i < 3; ++i) // check whether the face is along the singularity, the singularity edge should be boundary
    if(find(singularity_edges.begin(),singularity_edges.end(),make_pair(face_idx_original[i],face_idx_original[(i+1)%3])) == singularity_edges.end()
       && find(singularity_edges.begin(),singularity_edges.end(),make_pair(face_idx_original[(i+1)%3],face_idx_original[i])) == singularity_edges.end())
      need_to_search_edges.push_back(make_pair(face[i],face[(i+1)%3]));

  if(need_to_search_edges.empty()) return 0;

  for(size_t i = 0; i < need_to_search_edges.size(); ++i){
      const pair<size_t,size_t> &edge_ = need_to_search_edges[i];

      for(size_t t = 0; t < jump_face.size(); ++t){
          if(jump_face[t] != face_idx && !is_outside_face_cut_visited[jump_face[t]] ){
              const matrixst& each_face = outside_face_cut(colon(),jump_face[t]);
              if(find(each_face.begin(),each_face.end(),edge_.first) != each_face.end()
                 && find(each_face.begin(),each_face.end(),edge_.second) != each_face.end())
                {
                  tmp_adjacent_face.push_back(jump_face[t]);
                  break;
                }
            }
        }
    }
  return 0;
}

int type_matrix_map_test(const size_t type, matrixd &rot);
int label_jump_type_for_each_patch(
    const vector<size_t> &face_of_patch, // do not contain the pair faces in one patch
    const matrixst &cut_face_pair,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const matrixst &tet,
    const matrixd &node,
    const matrixd &cut_node,
    const matrixst &cut_tet2tet,
    matrixst &cut_jump_face_type,
    const matrixst &outside_face_cut,
    const matrixst &outside_face_cut_idx,
    const vector<deque<size_t> > &singularity_chain,
    const vector<deque<size_t> > &singularity_type,
    const boost::unordered_map<pair<size_t,size_t>, vector<pair<size_t,size_t> > > &edge_loop_pair,
    const matrix<matrixd > &frame_cut_tet)
{
  set<pair<size_t,size_t> > edges;
  size_t type;
  for(size_t t = 0; t < face_of_patch.size(); ++t){

      if(cut_face_pair[face_of_patch[t]] == -1) cerr << "# strange, can not be non jump surface." << endl;

      const pair<size_t,size_t> &cut_tet_pair_i = fa_cut.face2tet_[outside_face_cut_idx[face_of_patch[t]]];
      const pair<size_t,size_t> &cut_tet_pair_j = fa_cut.face2tet_[outside_face_cut_idx[cut_face_pair[face_of_patch[t]]]];

      assert(cut_tet_pair_i.first == -1 || cut_tet_pair_i.second == -1);
      assert(cut_tet_pair_j.first == -1 || cut_tet_pair_j.second == -1);

      const size_t t_i = (cut_tet_pair_i.first == -1)?cut_tet_pair_i.second:cut_tet_pair_i.first;
      const size_t t_j = (cut_tet_pair_j.first == -1)?cut_tet_pair_j.second:cut_tet_pair_j.first;

      matrixst common_face(3);
      common_face = outside_face_cut(colon(),face_of_patch[t]);
      matrixst face_in_origin_idx(3);
      for(size_t t = 0; t < 3; ++t)
        face_in_origin_idx[t] = cut_tet2tet[common_face[t]];

      //        if(find_common_face(tet(colon(),t_i),tet(colon(),t_j),common_face))
      //        {
      //            cerr << "# strange can not find common_face " << endl;
      //            cerr << "# tet " << t_i  << tet(colon(),t_i) << endl;
      //            cerr << "# tet " << t_j << tet(colon(),t_j) << endl;
      //            continue;
      //        }

      //        if(norm(face_in_origin_idx - common_face) != 0)
      //        {
      //            cerr << "face error." << endl;
      //            cerr << common_face << endl;
      //            cerr << face_in_origin_idx << endl;
      //        }
      type = get_jump_type(face_in_origin_idx, tet,node,edge_loop_pair, singularity_chain, singularity_type,t_i,t_j);
      assert(type != 9 && type != -1);
      if(type == 10) continue;
      else
        break;
    }
  if(type == 10) cerr << "# strange: this patch doesn't have singularities." << endl;

#if 0 // test set the type constant to be 0
  //static size_t tmp = 0;
  type = 5;
#endif
  for(size_t t = 0; t < face_of_patch.size(); ++t)
    cut_jump_face_type[face_of_patch[t]] = type;
#if 1 //check the patch is right
  {
    static size_t path_NO = 0;
    size_t patch_idx = 0;
    matrixd rot0(3,3);
    matrixd rot1(3,3);
    matrixd type_m(3,3);
    for(size_t t = 0; t < face_of_patch.size(); ++t){
        const size_t& jump_type = cut_jump_face_type[face_of_patch[t]];
        const pair<size_t,size_t> & tet_pair_0 = fa_cut.face2tet_[outside_face_cut_idx[face_of_patch[t]]];
        const pair<size_t,size_t> & tet_pair_1 = fa_cut.face2tet_[outside_face_cut_idx[cut_face_pair[face_of_patch[t]]]];

        assert(tet_pair_0.first == -1 || tet_pair_0.second == -1);
        assert(tet_pair_1.first == -1 || tet_pair_1.second == -1);

        const size_t t0 = (tet_pair_0.first == -1)?tet_pair_0.second:tet_pair_0.first;
        const size_t t1 = (tet_pair_1.first == -1)?tet_pair_1.second:tet_pair_1.first;

        get_best_alignment(&frame_cut_tet[t0](0,0),&frame_cut_tet[t1](0,0),&rot0[0]);
            get_best_alignment(&frame_cut_tet[t1](0,0),&frame_cut_tet[t0](0,0),&rot1[0]);
            type_matrix_map_test(jump_type,type_m);
        //            cerr << "# t " << t << endl;
        //            cerr << "# Before adjust : " << endl;
        //            cerr << "# jump_face difference t0-->t1 " << norm(rot0 - type_m) << endl;
        //            cerr << "# jump_face difference t1-->t0 " << norm(rot1 - type_m) << endl;
        size_t type_ = jump_type;
        if(norm(rot0 - type_m) > norm(rot1 - type_m)){
            ++patch_idx;
            if(type_ - type_/3 * 3 == 2)
              type_ -= 2;
            else if(type_ - type/3 * 3 == 0) type_ += 2;
            cut_jump_face_type[face_of_patch[t]] = type_;
          }
        //type_matrix_map_test(type_,type_m);
        //            cerr << "# After adjust : " << endl;
        //            cerr << "# jump_face difference t0-->t1 " << norm(rot0 - type_m) << endl;
        //            cerr << "# jump_face difference t1-->t0 " << norm(rot1 - type_m) << endl;
      }
    cerr << "# patch " << path_NO++ << " contains " << face_of_patch.size() << " faces." << endl;
    cerr << "# modified " << patch_idx << "face in patch." << endl;
  }
#endif

#if 1 // check dump out the face patch
  {
    static size_t path_NO_ = 0;
    std::stringstream ss;
    ss << path_NO_++;
    string vtk_file = "patch_";
    vtk_file += ss.str();
    vtk_file += ".vtk";
    matrixst tri_face(3,face_of_patch.size());
    for(size_t t = 0; t < face_of_patch.size(); ++t){
        tri_face(colon(),t) = outside_face_cut(colon(),face_of_patch[t]);
      }
    ofstream ofs(vtk_file.c_str());
    tri2vtk(ofs,&cut_node[0],cut_node.size(2),&tri_face[0],tri_face.size(2));
  }

#endif
  return 0;

}

int set_jump_face_type_new(const vector<deque<size_t> > &singularity_chain,
                           const vector<deque<size_t> > &singularity_type,
                           const jtf::mesh::face2tet_adjacent &fa,
                           const jtf::mesh::face2tet_adjacent &fa_cut,
                           const matrixst &tet,
                           const matrixd &node,
                           const matrixd &cut_node,
                           const matrixst &outside_face_cut,
                           const matrixst &outside_face_cut_idx,
                           matrixst &cut_jump_face_type,
                           matrixst &cut_face_pair,
                           const matrix<matrixd > &frame_cut_tet,
                           const matrixst &cut_tet2tet)
{

  vector<pair<size_t,size_t> >  singularity_edges;
  {//convert the singularity chain
    for(size_t t = 0; t < singularity_chain.size(); ++t) {
        for(size_t i = 0; i < singularity_chain[t].size(); i += 2){
            singularity_edges.push_back(make_pair(singularity_chain[t][i],singularity_chain[t][i+1]));
          }
      }
  }

  vector<bool> is_outside_face_cut_visited(cut_face_pair.size(),false); // is this patch face visited
  for(size_t t = 0; t < cut_face_pair.size(); ++t){
      if(cut_face_pair[t] == -1) is_outside_face_cut_visited[t] = true; // remove all non-jump face
    }

  vector<size_t> jump_face;
  //    cerr  << is_outside_face_cut_visited.size() << endl;
  for(size_t t = 0; t < is_outside_face_cut_visited.size(); ++t){
      //        if(t == 2365)
      //        cerr << t << endl;
      if(is_outside_face_cut_visited[t]) continue;
      //cut_face_pair_new[cut_face_pair_new[t]] = -1; // do not need the another jump face, i only neea one layer faces
      jump_face.push_back(t); // insert jump face Num of outside_face_cut_idx, for now, the jump_face contains all jump face_pair, which i only need one of each pair
    }
  //size_t path_idx = 0;
  vector<vector<size_t> > patch_faces;
  while(!jump_face.empty()){
      stack<size_t> face_adj_;
      face_adj_.push(jump_face.front()); // push a seed
      set<size_t> one_patch;

      if(cut_face_pair[jump_face.front()] != -1){
          is_outside_face_cut_visited[face_adj_.top()] = true;
          is_outside_face_cut_visited[cut_face_pair[face_adj_.top()]] = true; // the pair face of current face is not needed
          cut_face_pair[cut_face_pair[face_adj_.top()]] = -1; // no more need the pair face
        }

      //one_patch.insert(jump_face.front());
      while(!face_adj_.empty()){ // the patch will grow from the seed, end while reach the singularity edge

          //outside_face_cut_patch[face_adj_.top()] = path_idx; // label the path idx
          if(cut_face_pair[face_adj_.top()] == -1) {
              cerr << "# strange, this situation execpt to be not happer." << endl;
              face_adj_.pop();
              continue;
            }
          one_patch.insert(face_adj_.top());
          //            const size_t patch_num = one_patch.size();
          const size_t face_idx = face_adj_.top();

          face_adj_.pop();
#if 1 // check
          if(cut_face_pair[cut_face_pair[face_idx]] != -1
             || !is_outside_face_cut_visited[face_idx]
             || !is_outside_face_cut_visited[cut_face_pair[face_idx]])
            cerr << "# error, this face should be visited." << endl;
#endif
          vector<size_t> tmp_adjacent_face;
          find_adjacent_face(outside_face_cut,is_outside_face_cut_visited,jump_face,face_idx,tmp_adjacent_face,singularity_edges,cut_tet2tet);
          for(size_t t = 0; t < tmp_adjacent_face.size(); ++t){
              const size_t pre_size = one_patch.size();
              one_patch.insert(tmp_adjacent_face[t]);
              if(one_patch.size() != pre_size)
                {
                  face_adj_.push(tmp_adjacent_face[t]);
                  is_outside_face_cut_visited[tmp_adjacent_face[t]] = true; // tag this face visited and it's pair face also visited
                  is_outside_face_cut_visited[cut_face_pair[tmp_adjacent_face[t]]] = true;
                  cut_face_pair[cut_face_pair[tmp_adjacent_face[t]]] = -1; // break the pair face link, we do not the pair face's back link to itself
                }
            }
          //            if(patch_num == one_patch.size()) { // can not insert any more
          //                break;
          //            }
        }

      vector<size_t> jump_face_left; // find the left jump face
      for(size_t t = 0; t < jump_face.size(); ++t){
          if(!is_outside_face_cut_visited[jump_face[t]] && cut_face_pair[jump_face[t]] != -1)
            jump_face_left.push_back(jump_face[t]);
        }
      swap(jump_face_left,jump_face);
      //++path_idx;
      if(!one_patch.empty())
        {
          vector<size_t> tmp(one_patch.size());
          copy(one_patch.begin(),one_patch.end(),tmp.begin());
          patch_faces.push_back(tmp);
        }
    }

#if 0 // check the patch is correct
  size_t error_face_num = 0;
  for(size_t t = 0; t < patch_faces.size(); ++t){
      for(size_t i = 0; i < patch_faces[t].size(); ++i)
        if(cut_face_pair[cut_face_pair[patch_faces[t][i]]] != -1){
            cerr << "# incorrect face pair, not modified complete. patch " << t << " " << patch_faces[t][i] << "," << cut_face_pair[patch_faces[t][i]] << endl;
            ++error_face_num;
          }
    }
  cerr << "# error face pair " << error_face_num << endl;
#endif

  cerr << "# find " << patch_faces.size() << " patches." << endl;

  jtf::mesh::one_ring_tet_at_edge ortae;
  for(size_t ti = 0; ti < tet.size(2); ++ti)
    ortae.add_tets(tet(colon(), ti), fa);

  if(ortae.sort_into_loop(tet,node)) {
      cerr << "# sort error." << endl;
      return __LINE__;
    }

  typedef boost::unordered_map<std::pair<size_t,size_t>, std::vector<size_t> > e2tet_type;
  typedef e2tet_type::const_iterator eci;
  typedef boost::unordered_map<std::pair<size_t,size_t>, std::vector<pair<size_t,size_t> > > e2tet_type_p;
  e2tet_type_p edge_loop_pair; // initial the edge_loop_pair
  {// TODO: need to speed up
    for(eci eci_ = ortae.e2t_.begin(); eci_ != ortae.e2t_.end(); ++eci_){
        for(size_t t = 0; t < eci_->second.size() - 1; ++t)
          edge_loop_pair[eci_->first].push_back(make_pair(eci_->second[t],eci_->second[t+1]));
      }
  }

  //cut_jump_face_type.resize(cut_face_pair.size());
  cut_jump_face_type = -1 * ones<size_t>(cut_face_pair.size());


  for(size_t t = 0; t < patch_faces.size(); ++t)
    label_jump_type_for_each_patch(patch_faces[t],cut_face_pair,fa_cut,tet,node,cut_node,
                                   cut_tet2tet,cut_jump_face_type,outside_face_cut,outside_face_cut_idx,
                                   singularity_chain, singularity_type,edge_loop_pair,frame_cut_tet);

  return 0;
}


int set_jump_face_type(
    const vector<deque<size_t> > &singularity_chain,
    const vector<deque<size_t> > &singularity_type,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face_cut_idx,
    matrixst &cut_jump_face_type,
    const matrixst &cut_face_pair)
{
  jtf::mesh::one_ring_tet_at_edge ortae;
  for(size_t ti = 0; ti < tet.size(2); ++ti)
    ortae.add_tets(tet(colon(), ti), fa);

  if(ortae.sort_into_loop(tet,node)) {
      cerr << "# sort error." << endl;
      return __LINE__;
    }

  typedef boost::unordered_map<std::pair<size_t,size_t>, std::vector<size_t> > e2tet_type;
  typedef e2tet_type::const_iterator eci;
  typedef boost::unordered_map<std::pair<size_t,size_t>, std::vector<pair<size_t,size_t> > > e2tet_type_p;
  e2tet_type_p edge_loop_pair; // initial the edge_loop_pair
  {// TODO: need to speed up
    for(eci eci_ = ortae.e2t_.begin(); eci_ != ortae.e2t_.end(); ++eci_){
        for(size_t t = 0; t < eci_->second.size() - 1; ++t)
          edge_loop_pair[eci_->first].push_back(make_pair(eci_->second[t],eci_->second[t+1]));
      }
  }

  //cut_jump_face_type.resize(cut_face_pair.size());
  cut_jump_face_type = -1 * ones<size_t>(cut_face_pair.size());

  for(size_t t = 0; t < cut_face_pair.size(); ++t)
    {
      if(cut_face_pair[t] == -1) continue; // means original surface
      const pair<size_t,size_t> &cut_tet_pair_i = fa_cut.face2tet_[outside_face_cut_idx[t]];
      const pair<size_t,size_t> &cut_tet_pair_j = fa_cut.face2tet_[outside_face_cut_idx[cut_face_pair[t]]];

      assert(cut_tet_pair_i.first == -1 || cut_tet_pair_i.second == -1);
      assert(cut_tet_pair_j.first == -1 || cut_tet_pair_j.second == -1);

      const size_t t_i = (cut_tet_pair_i.first == -1)?cut_tet_pair_i.second:cut_tet_pair_i.first;
      const size_t t_j = (cut_tet_pair_j.first == -1)?cut_tet_pair_j.second:cut_tet_pair_j.first;

      matrixst common_face(3);

      if(jtf::mesh::find_common_face(tet(colon(),t_i),tet(colon(),t_j),common_face))
        {
          cerr << "# strange can not find common_face " << endl;
          cerr << "# tet " << t_i  << tet(colon(),t_i) << endl;
          cerr << "# tet " << t_j << tet(colon(),t_j) << endl;
          continue;
        }
      //        if(t == 32 && t_i == 22 && t_j == 3253)
      //            cerr << "# error " << t << endl;
      const size_t type = get_jump_type(common_face, tet,node,edge_loop_pair, singularity_chain, singularity_type,t_i,t_j);
      //const size_t type = 3;

      assert(type < 9 && type != -1);
      cut_jump_face_type[t] = type;
    }
  return 0;
}

int type_matrix_map_test(const size_t type, matrixd &rot)
{
  if(type == 9 || type == -1) return __LINE__;
  if(type == 10) // no need jump
    {
      rot = eye<double>(3);
      return 0;
    }
  static const double stander_rot[9][3][3] = {
    {{1,0,0},{0,0,1},{0,-1,0}},  // u1
    {{1,0,0},{0,-1,0},{0,0,-1}}, // u2
    {{1,0,0},{0,0,-1},{0,1,0}},  // u3

    {{0,0,-1},{0,1,0},{1,0,0}},  // v1
    {{-1,0,0},{0,1,0},{0,0,-1}}, // v2
    {{0,0,1},{0,1,0},{-1,0,0}},  // v3

    {{0,1,0},{-1,0,0},{0,0,1}},  // w1
    {{-1,0,0},{0,-1,0},{0,0,1}}, // w2
    {{0,-1,0},{1,0,0},{0,0,1}}}; // w3
  itr_matrix<const double *> rot_(3,3,&stander_rot[type][0][0]);
  rot = rot_;
  return 0;
}

double get_current_time(){
  timeval tv;
  gettimeofday(&tv,0);
  return tv.tv_sec+1e-6*tv.tv_usec;
}

int init_zyz_as_linear_optimization_for_refine_process(
    const matrixst &tet,
    const matrixd &node,
    const matrixd &fixed_frame,
    const matrixst &fixed_frame_idx,
    const boost::unordered_map<size_t,size_t> &surface_type,
    const matrixst &outside_face,
    const matrixst &outside_face_idx,
    const double weight[2],
const matrixd &stiff,
const boost::unordered_map<pair<size_t,size_t>,size_t > &inner_face_jump_type,
const double LP_surface,
const double LP_smooth,
const double LP_rtr,
const double rtr_w,
const matrix<matrixd > &frame_cut_tet,
matrixd &new_frame, // 9 * n
boost::property_tree::ptree &pt)
{
  non_sym_frame_opt func_linear;
  func_linear.setup_equations_linear_new(tet,node,fixed_frame,fixed_frame_idx, surface_type, outside_face, outside_face_idx,
                                         weight, stiff,inner_face_jump_type, LP_surface, LP_smooth, LP_rtr,rtr_w);

  new_frame.resize(9,frame_cut_tet.size());
  for(size_t t = 0;t < frame_cut_tet.size(); ++t){
      new_frame(colon(),t) = frame_cut_tet[t](colon());
    }

#if 1
  matrixd residual(func_linear.get()->dim_of_f());
  cerr << "#### frame func num = " << func_linear.get()->dim_of_f() << endl;
  zjucad::optimize(*func_linear.get(), new_frame, residual, pt);
#endif
  return 0;
}

int dump_out_angle_of_frame_field(
    boost::property_tree::ptree &pt,
    const jtf::mesh::face2tet_adjacent &fa,
    const matrixst &tet,
    const matrixd &node,
    const matrix<matrixd > &frame,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const char* tag)
{
  matrixd jump_error_angle_each_face = zeros<double>(fa.face2tet_.size(),1);

  matrixd jump_error_angle = zeros<double>(3,1);

  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mci;

  for(size_t t = 0; t < fa.face2tet_.size(); ++t){
      const pair<size_t,size_t> &tet_pair = fa.face2tet_[t];
      if(tet_pair.first == -1 || tet_pair.second == -1) {
          jump_error_angle_each_face[t] = 0;
        }else{
          mci ci = inner_face_jump_type.find(tet_pair);
          if(ci == inner_face_jump_type.end()){ // the jump type may be identity
              for(size_t i =0; i < 3; ++i){
                  jump_error_angle[i] =
                      fabs(acos(std::min(1.0,std::max(-1.0,dot(frame[tet_pair.first](colon(),i), frame[tet_pair.second])
                                / (norm(frame[tet_pair.first](colon(),i)) * norm(frame[tet_pair.second]))
                      ))) * 180.0 / My_PI()) ;
                }
              jump_error_angle_each_face[t] = *max_element(jump_error_angle.begin(),jump_error_angle.end());
            }
          else{ // jump type not equals identity
              const size_t type = ci->second;
              const matrixd type_trans = type_transition2(type);

              const matrixd first = frame[tet_pair.first] * type_trans;
              for(size_t i =0; i < 3; ++i){
                  jump_error_angle[i] =  fabs(acos(min(1.0,max(-1.0,dot(first(colon(),i), frame[tet_pair.second](colon(),i))
                                                       /(norm(first(colon(),i)) * norm(frame[tet_pair.second](colon(),i))))
                                              )) * 180.0 / My_PI()) ;
                }

              jump_error_angle_each_face[t] = *max_element(jump_error_angle.begin(),jump_error_angle.end());
            }
        }
    }

  const double max_err_angle = *max_element(jump_error_angle_each_face.begin(),jump_error_angle_each_face.end());

  assert(jump_error_angle_each_face.size());

  const double average_err_angle = accumulate(jump_error_angle_each_face.begin(),jump_error_angle_each_face.end(),0.0)
      / jump_error_angle_each_face.size();

  string jump_err_angle_perface = pt.get<string>("output_zyz.value");
  jump_err_angle_perface += ".jump_err_angle_perface_";
  jump_err_angle_perface += tag;
  ofstream ofs_jump_err_angle_perface(jump_err_angle_perface.c_str());

  for(size_t t = 0; t < jump_error_angle_each_face.size(); ++t){
      ofs_jump_err_angle_perface << jump_error_angle_each_face[t] << endl;
    }
  ofs_jump_err_angle_perface << "max angle error:" << max_err_angle << endl;
  ofs_jump_err_angle_perface << "average angle error:" << average_err_angle << endl;


  cerr << "# secess calculate the jump face error for original frame under new type and new frame under new type" << endl;

  return 0;
}

int dump_out_angle_of_frame_field_zyz(
    boost::property_tree::ptree &pt,
    const jtf::mesh::face2tet_adjacent &fa,
    const matrixst &tet,
    const matrixd &node,
    const matrixd &new_zyz,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const char * tag)
{
  matrix<matrixd > frame(tet.size(2));
  for(size_t t = 0; t < tet.size(2); ++t){
      frame[t].resize(3,3);
      zyz_angle_2_rotation_matrix1(&new_zyz(0,t),&frame[t][0]);
    }
  dump_out_angle_of_frame_field(pt,fa,tet,node,frame,inner_face_jump_type,tag);
  return 0;
}

int refine_frame_field_after_aligned_new(
    boost::property_tree::ptree &pt, jtf::tet_mesh &tm,
    const matrixd &zyz, matrixd & new_zyz)
{
  matrix<matrixd > frame_cut_tet;

  tetmesh_cutter tmc(tm);
  new_zyz = zyz;
  tmc.cut(new_zyz);

  zyz2frame(new_zyz, frame_cut_tet);

  pt.put("new_tet.desc", "new_tet desc");
  pt.put("temp_zyz.desc", "temp_zyz which is corresponding to zyz");
  pt.put("dump_out_inner_face_jump_type.desc", "output inner face jump type");
  pt.put("dump_out_surface_type.desc", "dump out surface align type");

  // initial process
  matrixst fixed_frame_idx,aligned_idx;
  matrixd fixed_frame,aligned;

  const double weight[2] = {
    pt.get<double>("surface_align_w.value"),
    pt.get<double>("jump_smooth_w.value")
  };

  matrixd stiff = ones<double>(tmc.cut_tm_.mesh_.size(2), 1); // do not know how to use the stiff

  if(zjucad::has("stiff.value",pt)){
      ifstream ifs(pt.get<string>("stiff.value").c_str(), ifstream::binary);
      if(ifs.fail()) {
          cerr << "# open stiff fail." << endl;
          return __LINE__;
        }
      cerr << "# use stiff" << endl;
      jtf::mesh::read_matrix(ifs, stiff);
    }

  if(load_from_tet_inner(tm.tetmesh_.node_, tm.tetmesh_.mesh_,
                         fixed_frame, fixed_frame_idx,
                         aligned, aligned_idx,0,0,0)) {
      cerr << "load fail." << endl;
      return __LINE__;
    }

  boost::unordered_map<pair<size_t,size_t>,size_t> inner_face_jump_type;
  extract_inner_jump_type(tm.tetmesh_.mesh_,frame_cut_tet,*tm.fa_,inner_face_jump_type);

  boost::unordered_map<pair<size_t,size_t>,size_t> init_inner_face_jump_type =
      inner_face_jump_type;

  if(zjucad::has("inner_type.value", pt)){
      if(load_inner_face_jump_type(pt.get<string>("inner_type.value").c_str(),
                                   init_inner_face_jump_type))
        return __LINE__;
      cerr << "# load inner type." << endl;
    }

  ///////////////////////////////////////////////////////////////////

  vector<deque<pair<size_t,size_t> > > chain_list;
  vector<deque<size_t>  > singularities_type;

  singularity_extractor se(tm);
  se.extract(init_inner_face_jump_type, chain_list, singularities_type);

  dump_singularity_chain_to_vtk_2("before_modify_inner_singularity.vtk",
                                  tm.tetmesh_.node_,chain_list,singularities_type);

  dump_singularity_to_cylinder("original_inner_singularity.obj",
                               tm.tetmesh_.node_,chain_list,0.002);

  relabel_singularity_chain_by_modify_face_jump_type(
        tm.tetmesh_.mesh_,tm.tetmesh_.node_,tm.outside_face_,tm.outside_face_idx_,frame_cut_tet,*tm.fa_,tm.ortae_,
        chain_list,singularities_type, inner_face_jump_type);

  //matrixst surf_aligned_type;
  boost::unordered_map<size_t,size_t> surface_type;
  //surf_aligned_type.resize(outside_face_idx.size());

  remove_surface_zigzag(tm.outside_face_,tm.tetmesh_.mesh_,tm.tetmesh_.node_,*tm.fa_,inner_face_jump_type,
                        tm.outside_face_idx_,frame_cut_tet, tm.ortae_,pt,surface_type);

  // remove surface zigzag once more
  remove_surface_zigzag(tm.outside_face_,tm.tetmesh_.mesh_,tm.tetmesh_.node_,*tm.fa_,inner_face_jump_type,
                        tm.outside_face_idx_,frame_cut_tet, tm.ortae_,pt,surface_type);

  dump_inner_face_jump_type(
        pt.get<string>("dump_out_inner_face_jump_type.value").c_str(),
        inner_face_jump_type);

  dump_surface_normal_align_type_map(
        pt.get<string>("dump_out_surface_type.value").c_str(),
        tm.outside_face_,tm.outside_face_idx_, surface_type);

  jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("new_tet.value").c_str(),
                                      &tm.tetmesh_.node_, &tm.tetmesh_.mesh_);
  {
    assert(frame_cut_tet.size() == tm.tetmesh_.mesh_.size(2));
    matrixd temp_zyz;
    frame2zyz(frame_cut_tet, temp_zyz);
    write_zyz(pt.get<string>("temp_zyz.value").c_str(), temp_zyz);
  }

  const double LP_surface = pt.get<double>("LP_surface.value");
  const double LP_smooth = pt.get<double>("LP_smooth.value");
  const double LP_rtr = pt.get<double>("LP_rtr.value");
  const double rtr_w = pt.get<double>("rtr_w.value");

  pt.put("LP_surface.desc","the surface alignment LP num.");
  pt.put("LP_smooth.desc","the inner smoothing LP num.");
  pt.put("LP_rtr.desc","the RTR LP num");
  pt.put("rtr_w.desc","weighting of rtr");

  cerr << "# optimization : LP_surface " << LP_surface << endl;
  cerr << "# optimization : LP_smooth " << LP_smooth << endl;
  cerr << "# optimization : LP_rtr " << LP_rtr << endl;
  cerr << "# optimization : rtr_w " << rtr_w << endl;

  if(zjucad::has("init_zyz.value",pt)){
      ifstream ifs(pt.get<string>("init_zyz.value").c_str(),ifstream::binary);
      if(ifs.fail()){
          cerr << "# can not open init zyz." << endl;
        }else{
          jtf::mesh::read_matrix(ifs,new_zyz);
          if(!new_zyz.size()){
              cerr << "# read zyz fail." << endl;
              return __LINE__;
            }else if(new_zyz.size(2) != tm.tetmesh_.mesh_.size(2)){
              cerr << "# error zyz format." << endl;
              return __LINE__;
            }
        }
    }else{
      matrixd new_frame(9,tm.tetmesh_.mesh_.size(2));

      // there need two iterations: one use rtr_w as 0
      // second use rtr_w as 100
      cerr << "# [info] start the first iterator using rtr_w as 0" << endl;
      init_zyz_as_linear_optimization_for_refine_process(
            tm.tetmesh_.mesh_,tm.tetmesh_.node_,fixed_frame,fixed_frame_idx, surface_type,
            tm.outside_face_, tm.outside_face_idx_, weight,stiff, inner_face_jump_type,LP_surface,
            LP_smooth, LP_rtr, 0, frame_cut_tet,new_frame,pt);

      cerr << "# [info] start the first iterator using rtr_w as " << rtr_w << endl;
      init_zyz_as_linear_optimization_for_refine_process(
            tm.tetmesh_.mesh_, tm.tetmesh_.node_,fixed_frame,fixed_frame_idx, surface_type,
            tm.outside_face_, tm.outside_face_idx_, weight,stiff, inner_face_jump_type,LP_surface,
            LP_smooth, LP_rtr, rtr_w, frame_cut_tet,new_frame,pt);
      for(size_t ti = 0; ti < new_frame.size(2); ++ti){
          rotation_matrix_2_zyz_angle(&new_frame(0,ti), &new_zyz(0,ti), 0);
        }
    }

  {// resolve the equations
    non_sym_frame_opt  func;
    func.cut_inner_smooth_function_.clear();
    func.cut_jump_smooth_function_.clear();
    func.cut_surface_align_function_.clear();
    func.funcs_.clear();
    func.setup_equations_new(tm.tetmesh_.mesh_,tm.tetmesh_.node_,fixed_frame,fixed_frame_idx, surface_type,
                             tm.outside_face_, tm.outside_face_idx_,
                             weight, stiff,inner_face_jump_type, LP_surface, LP_smooth);

    matrixd residual(func.get()->dim_of_f());
    cerr << "#### frame func num = " << func.get()->dim_of_f() << endl;
    zjucad::optimize(*func.get(), new_zyz, residual, pt);
  }

  matrixd ori_zyz;
  frame2zyz(frame_cut_tet, ori_zyz);
  dump_out_angle_of_frame_field_zyz(pt,*tm.fa_,tm.tetmesh_.mesh_,tm.tetmesh_.node_,ori_zyz,inner_face_jump_type,
                                    "old_frame_under_new_type");

  dump_out_angle_of_frame_field_zyz(pt,*tm.fa_,tm.tetmesh_.mesh_,tm.tetmesh_.node_,new_zyz,inner_face_jump_type,
                                    "new_frame_under_new_type");
  return 0;
}

//! @brief this function is used to remove black lines and zigzag by splitting tets
// this method can make sure to clean all black lines and zigzag,
// but the final singularity distribution is not good enough
int refine_frame_field_after_aligned_split(
    boost::property_tree::ptree &pt,
    jtf::tet_mesh &tm,
    const matrixd &zyz,
    matrixd & new_zyz)
{
  matrix<matrixd > frame(zyz.size(2));
  boost::unordered_map<pair<size_t,size_t>,size_t> inner_face_jump_type;

  vector<deque<pair<size_t,size_t> > > chain_list;
  vector<deque<size_t> > singularities_type_;

  { // preprocess
    pt.put("use_mst.desc","if use minimal spanning tree to adjust frame,[y/n]");
    string use_mst = pt.get<string>("use_mst.value","n");
    if(use_mst == "y" || use_mst == "Y"
       || use_mst == "yes" || use_mst == "Yes"
       || use_mst == "YES"){
        tetmesh_cutter tmc(tm);
        new_zyz = zyz;
        tmc.cut(new_zyz);
        zyz2frame(new_zyz, frame);
        cerr << "# [info] used MST" << endl;
      }else{
        zyz2frame(zyz, frame);
      }

    singularity_extractor se(tm);
    se.extract(frame, chain_list, singularities_type_);
  }

  {// remove black line and zigzag and near surface
    pt.put("inner_face_jump_type.desc","init inner face jump type file");
    if( zjucad::has("inner_face_jump_type.value",pt)){
        if(load_inner_face_jump_type(
             pt.get<string>("inner_face_jump_type.value").c_str(),
             inner_face_jump_type))
          return __LINE__;

      }else{
        extract_inner_jump_type(tm.tetmesh_.mesh_,frame,*tm.fa_,inner_face_jump_type);
#if 1 // check
        {
          vector<deque<pair<size_t,size_t> > > chain_;
          vector<deque<size_t> > type_;
          singularity_extractor se(tm);
          se.extract(inner_face_jump_type, chain_list, singularities_type_);
          dump_singularity_chain_to_vtk_2(
                "original_singularity_edge_use_type.vtk",tm.tetmesh_.node_,chain_,type_);
        }
#endif
        pt.put("dump_tet_rot.desc","filename: global rotation defined in tet to root");
        pt.put("tet_rot.desc","filename of tet rot type");
        vector<size_t> tet_rot_to_root;
        if(zjucad::has("tet_rot.value",pt)){
            load_tet_rotation_array(pt.get<string>("tet_rot.value").c_str(),tet_rot_to_root);
          }else{
            const string dump_tet_rot_str = pt.get<string>("dump_tet_rot.value");
            calc_tet_rotation(tm,inner_face_jump_type,tet_rot_to_root);
            dump_tet_rotation_from_array(dump_tet_rot_str.c_str(),tet_rot_to_root);
          }
        relabel_singularity_chain_by_splitting(
              tm,frame,chain_list, singularities_type_, inner_face_jump_type,pt,tet_rot_to_root);
        dump_out_singularity_chain_3("singularity_after_relabel.sc", chain_list,singularities_type_);
      }
  }

  {// handle the synchronizing
    if(zjucad::has("inner_face_jump_type.value",pt) ||
       zjucad::has("surface_type.value",pt)){
        if(!(zjucad::has("new_tet.value",pt)) ||
           !(zjucad::has("new_zyz.value",pt))){
            cerr << "# [error] if inner face jump type and surface type are read from file,"
                    " the new_tet and new_zyz needs to be given as well." << endl;
            return __LINE__;
          }
        tm.load(pt.get<string>("new_tet.value").c_str());
        matrixd new_zyz;
        if(read_zyz(pt.get<string>("new_zyz.value").c_str(), new_zyz)){
            cerr << "# [error] can not load new_zyz." << endl;
            return __LINE__;
          }

        if(new_zyz.size(2) != tm.tetmesh_.mesh_.size(2)){
            cerr << "# [error] wrong zyz file." << endl;
            return __LINE__;
          }
        frame.resize(new_zyz.size(2));
        zyz2frame(new_zyz, frame);
      }
  }

  pt.put("relabel_iter.desc","relabel singualrity chain iter num");
  const size_t relabel_iter_num = pt.get<size_t>("relable_iter.value");
  {
    vector<size_t> tet_rot_to_root;
    singularity_extractor se(tm);

    for(size_t i = 0; i < relabel_iter_num; ++i){
        vector<deque<pair<size_t,size_t> > > chain_;
        vector<deque<size_t> > type_;
        se.extract(inner_face_jump_type, chain_, type_);
        calc_tet_rotation(tm,inner_face_jump_type,tet_rot_to_root);
        relabel_singularity_chain_by_splitting(tm,frame,chain_list,singularities_type_,inner_face_jump_type,pt,tet_rot_to_root);

      }
    dump_out_singularity_chain_3("singularity_after_relabel.sc",chain_list,singularities_type_);
  }

  boost::unordered_map<size_t,size_t> surface_type;

  //surf_aligned_type.resize(outside_face_idx.size());
  { // relabel surface
    pt.put("surface_type.desc","init surface type file");
    if(zjucad::has("surface_type.value",pt) && relabel_iter_num == 0){
        ifstream ifs(pt.get<string>("surface_type.value").c_str());
        if(ifs.fail()){
            cerr << "# [error] can not find surface_type file." << endl;
            return __LINE__;
          }
        size_t triface[3],type;
        while(!ifs.eof()){
            ifs >> triface[0] >> triface[1] >> triface[2] >> type;
            size_t face_idx = tm.fa_->get_face_idx(triface);
            surface_type[face_idx] = type;
          }
      }else{
        //! there is no need to remove surface zigzag
        remove_surface_zigzag(tm.outside_face_,tm.tetmesh_.mesh_,tm.tetmesh_.node_,*tm.fa_,inner_face_jump_type
                              ,tm.outside_face_idx_,frame, tm.ortae_,pt,surface_type);
      }
  }

  {// data dump out
    pt.put("dump_inner_face_jump_type.desc","dump out inner face jump type");
    if(zjucad::has("dump_inner_face_jump_type.value",pt) ||
       relabel_iter_num != 0){
        ofstream ofs(pt.get<string>("dump_inner_face_jump_type.value").c_str());
        if(ofs.fail()){
            cerr << "# [error] can not open dump_inner_face_jump_type file." << endl;
            return __LINE__;
          }
        for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mcit = inner_face_jump_type.begin();
            mcit != inner_face_jump_type.end(); ++mcit){
            if(is_trivial_type(mcit->second)) continue;
            ofs << mcit->first.first << " " << mcit->first.second << " " << mcit->second << endl;
          }
        cerr << "# [info] success dump out file " << pt.get<string>("dump_inner_face_jump_type.value") << endl;
      }
    pt.put("dump_surface_type.desc","dump out surface type");
    if(zjucad::has("dump_surface_type.value",pt) ||
       relabel_iter_num != 0){
        ofstream ofs(pt.get<string>("dump_surface_type.value").c_str());
        if(ofs.fail()){
            cerr << "# [error] can not open dump_surface_type file." << endl;
            return __LINE__;
          }
        for(boost::unordered_map<size_t,size_t>::const_iterator mcit = surface_type.begin();
            mcit != surface_type.end(); ++mcit){
            const vector<size_t> &tri = tm.fa_->faces_[mcit->first];
            ofs << tri[0] << " " << tri[1] << " " << tri[2] << " " << mcit->second << endl;
          }
        cerr << "# [info] success dump out file " << pt.get<string>("dump_surface_type.value") << endl;
      }
    pt.put("dump_new_tet.desc","dump out new tet file");
    pt.put("dump_new_zyz.desc","dump out new zyz file, this zyz is synchronized with splitted tet");
    if(zjucad::has("dump_new_tet.value",pt) ||
       relabel_iter_num !=0 ){
        jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("dump_new_tet.value").c_str(),&tm.tetmesh_.node_,&tm.tetmesh_.mesh_);
      }
    if(zjucad::has("dump_new_zyz.value",pt) ||
       relabel_iter_num !=0){
        matrixd dump_zyz(3,frame.size());
        frame2zyz(frame, dump_zyz);
        write_zyz(pt.get<string>("dump_new_zyz.value").c_str(), dump_zyz);
        cerr << "# [info] success dump out file " << pt.get<string>("dump_new_zyz.value") << endl;
      }
  }

  { // optimization
    // initialization
    matrixst fixed_frame_idx,aligned_idx;
    matrixd fixed_frame,aligned;
    pt.put("surface_align_w.desc","weight of surface alignment");
    pt.put("jump_smooth_w.desc","weight of jump smoothing");
    const double weight[2] = {
      pt.get<double>("surface_align_w.value"),
      pt.get<double>("jump_smooth_w.value")
    };

    if(load_from_tet_inner(tm.tetmesh_.node_, tm.tetmesh_.mesh_,
                           fixed_frame, fixed_frame_idx,
                           aligned, aligned_idx,0,0,0)) {
        cerr << "load fail." << endl;
        return __LINE__;
      }

    matrixd stiff = ones<double>(tm.tetmesh_.mesh_.size(2), 1);

    const double LP_surface = pt.get<double>("LP_surface.value");
    const double LP_smooth = pt.get<double>("LP_smooth.value");
    const double LP_rtr = pt.get<double>("LP_rtr.value");
    const double rtr_w = pt.get<double>("rtr_w.value");

    pt.put("LP_surface.desc","the surface alignment LP num.");
    pt.put("LP_smooth.desc","the inner smoothing LP num.");
    pt.put("LP_rtr.desc","the RTR LP num");
    pt.put("rtr_w.desc","weighting of rtr");

    cerr << "# optimization : LP_surface " << LP_surface << endl;
    cerr << "# optimization : LP_smooth " << LP_smooth << endl;
    cerr << "# optimization : LP_rtr " << LP_rtr << endl;
    cerr << "# optimization : rtr_w " << rtr_w << endl;

    pt.put("init_zyz.desc","init zyz for non_linear optimization");
    if(zjucad::has("init_zyz.value",pt)){
        ifstream ifs(pt.get<string>("init_zyz.value").c_str(),ifstream::binary);
        if(ifs.fail()){
            cerr << "# can not open init zyz." << endl;
          }else{
            jtf::mesh::read_matrix(ifs,new_zyz);
            if(!new_zyz.size()){
                cerr << "# read zyz fail." << endl;
                return __LINE__;
              }else if(new_zyz.size(2) != tm.tetmesh_.mesh_.size(2)){
                cerr << "# error zyz format." << endl;
                return __LINE__;
              }
          }
      }else{
        matrixd new_frame(9,tm.tetmesh_.mesh_.size(2));
        init_zyz_as_linear_optimization_for_refine_process(
              tm.tetmesh_.mesh_,
              tm.tetmesh_.node_,fixed_frame,fixed_frame_idx, surface_type, tm.outside_face_,
              tm.outside_face_idx_,
              weight,stiff, inner_face_jump_type,LP_surface,LP_smooth, LP_rtr, rtr_w, frame,new_frame,pt);
        new_zyz.resize(3, tm.tetmesh_.mesh_.size(2));
        for(size_t ti = 0; ti < new_frame.size(2); ++ti){
            rotation_matrix_2_zyz_angle(&new_frame(0,ti), &new_zyz(0,ti), 0);
          }
      }

#if 1
    {// resolve the equations
      non_sym_frame_opt  func;
      func.cut_inner_smooth_function_.clear();
      func.cut_jump_smooth_function_.clear();
      func.cut_surface_align_function_.clear();
      func.funcs_.clear();
      func.setup_equations_new(tm.tetmesh_.mesh_,tm.tetmesh_.node_,fixed_frame,fixed_frame_idx, surface_type,
                               tm.outside_face_, tm.outside_face_idx_,
                               weight, stiff,inner_face_jump_type, LP_surface, LP_smooth);

      matrixd residual(func.get()->dim_of_f());
      cerr << "#### frame func num = " << func.get()->dim_of_f() << endl;
      zjucad::optimize(*func.get(), new_zyz, residual, pt);
    }
#endif

    matrixd ori_zyz;
    frame2zyz(frame, ori_zyz);

    dump_out_angle_of_frame_field_zyz(pt,*tm.fa_,tm.tetmesh_.mesh_,tm.tetmesh_.node_,ori_zyz,inner_face_jump_type,"old_frame_under_new_type");

    dump_out_angle_of_frame_field_zyz(pt,*tm.fa_,tm.tetmesh_.mesh_, tm.tetmesh_.node_,new_zyz,inner_face_jump_type,"new_frame_under_new_type");
  }

  return 0;
}

bool is_constructing_chain(const deque<pair<size_t,size_t> > &chain,
                           map<size_t,size_t> &vertex_counts)
{
  if(chain.empty()) return true;
  if(vertex_counts[chain.front().first] !=  2 // can not use const ?
     && vertex_counts[chain.back().second] !=  2)
    return false;
  if(chain.front().first == chain.back().second)
    return false;
  return true;
}


int cut_singularity_edge_into_chain(const vector<pair<pair<size_t,size_t>,size_t> > &surface_singularity_edge,
                                    vector<deque<pair<size_t,size_t> > > &singularity_chain,
                                    vector<deque<size_t> > &singularity_type)
{
  singularity_chain.clear();
  singularity_type.clear();
  map<pair<size_t,size_t>,size_t> surface_singularity_edge_map;
  map<size_t,size_t> vertex_counts;
  typedef map<size_t,size_t>::iterator mit;
  typedef map<pair<size_t,size_t>,size_t>::const_iterator mci;
  typedef set<pair<size_t,size_t> >::iterator sit;
  set<pair<size_t,size_t> > edge_segments;
  for(size_t t = 0; t < surface_singularity_edge.size(); ++t){
      surface_singularity_edge_map[surface_singularity_edge[t].first] = surface_singularity_edge[t].second;
      edge_segments.insert(surface_singularity_edge[t].first);
      mit it = vertex_counts.find(surface_singularity_edge[t].first.first);
      if(it == vertex_counts.end())
        vertex_counts[surface_singularity_edge[t].first.first] = 1;
      else
        it->second += 1;
      it = vertex_counts.find(surface_singularity_edge[t].first.second);
      if(it == vertex_counts.end()) vertex_counts[surface_singularity_edge[t].first.second] = 1;
      else
        it->second += 1;
    }


  while(!edge_segments.empty())
    {
      deque<pair<size_t,size_t> > chain;
      while(is_constructing_chain(chain,vertex_counts))
        {
          for(sit it = edge_segments.begin(); it != edge_segments.end(); ++it){
              if(chain.empty()){
                  chain.push_back(*it);
                  edge_segments.erase(it);
                  continue;
                }
              if(vertex_counts[chain.back().second] == 2){
                  if(chain.back().second == it->first){
                      chain.push_back(*it);
                      edge_segments.erase(*it);
                      continue;
                    }
                  if(chain.back().second == it->second){
                      chain.push_back(make_pair(it->second,it->first));
                      edge_segments.erase(*it);
                      continue;
                    }
                }
              if(vertex_counts[chain.front().first] == 2){
                  if(chain.front().first == it->second){
                      chain.push_front(*it);
                      edge_segments.erase(*it);
                      continue;
                    }
                  if(chain.front().first == it->first){
                      chain.push_front(make_pair(it->second,it->first));
                      edge_segments.erase(*it);
                      continue;
                    }
                }
              if(vertex_counts[chain.front().first] != 2 && vertex_counts[chain.back().second] != 2){
                  break;
                }
            }
        }
      singularity_chain.push_back(chain);
    }

  singularity_type.resize(singularity_chain.size());
  for(size_t t = 0; t < singularity_chain.size(); ++t){
      for(size_t j = 0; j < singularity_chain[t].size(); ++j){
          mci ci = surface_singularity_edge_map.find(singularity_chain[t][j]);
          if(ci != surface_singularity_edge_map.end()){
              singularity_type[t].push_back(ci->second);
            }else{
              ci = surface_singularity_edge_map.find(make_pair(singularity_chain[t][j].second,
                                                               singularity_chain[t][j].first));
              if(ci == surface_singularity_edge_map.end()){
                  cerr << "# strange: can not find this singularity edge." << endl;
                  continue;
                }
              matrixd rot = type_transition2(ci->second);
              singularity_type[t].push_back(type_transition1(trans(rot)));
            }
        }
    }
  return 0;
}

int reorder_edges(const pair<size_t,size_t> &edge0,
                  const pair<size_t,size_t> &edge1,
                  pair<size_t,size_t> &edge0_new,
                  pair<size_t,size_t> &edge1_new,
                  const matrixst &tet,
                  const jtf::mesh::face2tet_adjacent &fa,
                  const matrixd &node)
{
  matrixst face_a(3,1);
  face_a[0] = edge0.first;
  face_a[1] = edge0.second;
  face_a[2] = edge1.second;

  matrixd edge0_ = node(colon(),edge0.second) - node(colon(),edge0.first);
  matrixd edge1_ = node(colon(),edge1.second) - node(colon(),edge1.first);

  const pair<size_t,size_t> & tet_pair = fa.face2tet_[fa.get_face_idx(&face_a[0])];
  const size_t tet_idx = (tet_pair.first == -1)?tet_pair.second:tet_pair.first;
  matrixd face_normal_a(3,1);
  jtf::tetmesh::calculate_face_normal(tet,node,tet_idx,face_a,face_normal_a);

  if(dot(face_normal_a,cross(edge0_,edge1_)) > 0.0){
      edge0_new = edge1;
      edge1_new = edge0;
    }else{
      edge0_new = edge0;
      edge1_new = edge1;
    }
  return 0;
}

size_t get_edge_type_at_face_a2(
    const pair<size_t,size_t> &edge0_new,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::edge2cell_adjacent &ea,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &face_pair_tangent_jump,
    const boost::unordered_map<size_t,size_t> &face_jump,
    const matrixst &outside_face_idx,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &jump_type_between_tets)
{
  vector<size_t> loop;

  const size_t edge_idx = ea.get_edge_idx(edge0_new.first, edge0_new.second);
  assert(edge_idx != -1);

  auto it = ortae.e2t_.find(edge0_new);
  if(it != ortae.e2t_.end()) loop = it->second;
  else {
      it = ortae.e2t_.find(make_pair(edge0_new.second, edge0_new.first));
      if(it == ortae.e2t_.end()) {
          cerr << "# [error] can not find edge in ortae." << endl;
          return -1;
        }
      loop = it->second;
      reverse(loop.begin(), loop.end());
    }
  assert(ortae.is_inner_edge(loop) == false);
  if(loop.front() == -1) loop.erase(loop.begin(), loop.begin()+1);
  if(loop.back() == -1) loop.pop_back();

  pair<size_t,size_t> edge_cell = ea.edge2cell_[edge_idx];
  edge_cell.first = outside_face_idx[edge_cell.first];
  edge_cell.second = outside_face_idx[edge_cell.second];
  pair<size_t,size_t> tet_pair = jtf::tetmesh::face_pair2tet_pair(edge_cell, fa);
  assert(tet_pair.first == loop.front() || tet_pair.second == loop.front());
  assert(tet_pair.first == loop.back() || tet_pair.second == loop.back());

  if(tet_pair.second == loop.front()) {
      reverse(loop.begin(), loop.end());
    }

  matrix<double> rot = eye<double>(3);
  for(size_t ti = 0; ti != loop.size()-1; ++ti){
      auto it = jump_type_between_tets.find(make_pair(loop[ti], loop[ti+1]));
      if(it != jump_type_between_tets.end())
        rot = temp(rot * type_transition2(it->second));
    }

  auto ita = face_jump.find(edge_cell.first);
  auto itb = face_jump.find(edge_cell.second);
  assert(ita != face_jump.end() || itb != face_jump.end());

  rot = temp(trans(type_transition2(ita->second)) * rot * type_transition2(itb->second));

  auto itba = face_pair_tangent_jump.find(make_pair(edge_cell.second, edge_cell.first));
  rot = temp(rot * type_transition2(itba->second));
  return type_transition1(rot);
}

bool is_surface_singularity_zigzag_new(const pair<size_t,size_t> &edge0,
                                       const pair<size_t,size_t> &edge1,
                                       const jtf::mesh::face2tet_adjacent &fa,
                                       const matrixst &outside_face)
{
  set<size_t> face;
  face.insert(edge0.first);
  face.insert(edge0.second);
  face.insert(edge1.first);
  face.insert(edge1.second);

  if(face.size() != 3)
    return false;
  vector<size_t> face_(3);
  copy(face.begin(),face.end(),face_.begin());
  size_t face_idx = fa.get_face_idx(&face_[0]);
  if(face_idx == -1) return false;

  const pair<size_t,size_t> &two_tets = fa.face2tet_[face_idx];
  if(fa.is_outside_face(two_tets))
    return true;
  else
    return false;
}

int apply_rotation_to_face_to_remove_zigzag_new(
    const pair<size_t,size_t> &edge0,
    const pair<size_t,size_t> &edge1,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::face2tet_adjacent &fa,
    const matrixst & outside_face,
    const boost::unordered_map<pair<size_t,size_t>,size_t > &jump_type_between_tets,
    boost::unordered_map<size_t,size_t> &face_jump,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &face_pair_tangent_jump)
{
  assert(edge0.second == edge1.first);
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mpcit;
  typedef boost::unordered_map<size_t,size_t>::const_iterator mcit;

  size_t face_a_idx = fa.get_face_idx(edge0.first,edge0.second,edge1.second);
  if(face_a_idx == -1) {
      cerr << "# this two edges do not belong to one face." << endl;
      cerr << "# [error] edge0 <" << edge0.first << "," << edge0.second << ">." << endl;
      cerr << "# [error] edge1 <" << edge1.first << "," << edge1.second << ">." << endl;
      return __LINE__;
    }

  oecit it = ortae.e2t_.find(edge0);
  vector<size_t> around_tets;
  if(it != ortae.e2t_.end()){
      around_tets = it->second;
    }else {
      it = ortae.e2t_.find(make_pair(edge0.second,edge0.first));
      around_tets = it->second;
      reverse(around_tets.begin(),around_tets.end());
    }

  matrixst edge0_adjacent_face(3,2);
  vector<size_t> edge0_adjacent_face_idx(2);
  size_t face_b_idx;

  find_edge_adjacent_faces(edge0_adjacent_face,outside_face,edge0);
  for(size_t t = 0; t < 2; ++t){
      edge0_adjacent_face_idx[t] = fa.get_face_idx(&edge0_adjacent_face(0,t));
    }
#if 1 // check
  if(find(edge0_adjacent_face_idx.begin(),edge0_adjacent_face_idx.end(),face_a_idx) == edge0_adjacent_face_idx.end())
    {
      cerr << "# [error] strange this two adjacent faces doesn't include face_a." << endl;
      cerr << "# [error] two faces: " << edge0_adjacent_face_idx[0] << " " <<  edge0_adjacent_face_idx[1] << endl;
      cerr << "# [error] face_a " << face_a_idx << endl;
    }
#endif
  face_b_idx = (edge0_adjacent_face_idx[0] == face_a_idx)?edge0_adjacent_face_idx[1]:edge0_adjacent_face_idx[0];

  const pair<size_t,size_t> &a_tets = fa.face2tet_[face_a_idx];
  const pair<size_t,size_t> &b_tets = fa.face2tet_[face_b_idx];
  const size_t t0 = (a_tets.first == -1)?a_tets.second:a_tets.first;
  const size_t tn = (b_tets.first == -1)?b_tets.second:b_tets.first;

  matrixd rot_t02tn = eye<double>(3);
  matrixd b2a = eye<double>(3);
  matrixd t02a = eye<double>(3);
  matrixd tn2b = eye<double>(3);

  mpcit pit = face_pair_tangent_jump.find(make_pair(face_b_idx,face_a_idx));
  if(pit != face_pair_tangent_jump.end()){
      b2a = type_transition2(pit->second);
    }

  mcit ita = face_jump.find(face_a_idx);
  mcit itb = face_jump.find(face_b_idx);
  if(ita != face_jump.end()){
      t02a = type_transition2(ita->second);
    }
  if(itb != face_jump.end()){
      tn2b = type_transition2(itb->second);
    }

  if(t0 != tn) // not sharp edge, need to calculate the inner tets rot
    {
      vector<size_t> inner_tets(around_tets.size() - 2); // remove the two -1 at beginning and ending
      copy(around_tets.begin() + 1, around_tets.end() -1, inner_tets.begin());
      if(inner_tets.front() != t0) reverse(inner_tets.begin(),inner_tets.end());
      for(size_t i = 0; i < inner_tets.size() - 1; ++i) {
          mpcit it_ = jump_type_between_tets.find(make_pair(inner_tets[i],inner_tets[i+1]));
          if(it_ == jump_type_between_tets.end()) continue;
          rot_t02tn = temp(rot_t02tn * type_transition2(it_->second));
        }
    }

  matrixd type_matrix = temp(b2a *tn2b) * temp(rot_t02tn * trans(t02a));

  face_jump[face_a_idx] = type_transition1(type_matrix * t02a); // modify the t02a

  return 0;
}

int solve_face_jump_new(
    boost::unordered_map<size_t,size_t> &face_jump,
    vector<deque<pair<size_t,size_t> > > &singularity_chain,
    const vector<deque<size_t> > &singularity_type,
    const matrixst &tet,
    const matrixd &node,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    const jtf::mesh::face2tet_adjacent &fa,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &jump_type_between_tets,
    const matrix<matrixd > &frame_inner,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &face_pair_tangent_jump,
    const matrixst &outside_face)
{
  vector<pair<pair<size_t,size_t>,size_t > > error_edge;
  size_t zigzag_num = 0;
  for(size_t t = 0; t < singularity_chain.size(); ++t){
      deque<pair<size_t,size_t> > &single_chain = singularity_chain[t];
      deque<pair<size_t,size_t> > new_chain;
      new_chain.push_back(single_chain[0]);

      for(size_t i = 0; i < single_chain.size() - 1 ; ++i){
          if(is_surface_singularity_zigzag_new(new_chain.back(),single_chain[i+1],
                                               fa,outside_face))
            {
              apply_rotation_to_face_to_remove_zigzag_new(new_chain.back(),
                                                          single_chain[i+1],
                  ortae,fa,outside_face,jump_type_between_tets, face_jump,
                  face_pair_tangent_jump);
              pair<size_t,size_t> back_ = new_chain.back();
              new_chain.pop_back();
              new_chain.push_back(make_pair(back_.first,single_chain[i+1].second));
              ++zigzag_num;
            }else{
              new_chain.push_back(single_chain[i+1]);
            }

        }
      if(new_chain.size() != single_chain.size()){ // detect zigzag, need to replace the chain
          swap(single_chain,new_chain);
        }
    }

  cerr << "# [surface zigzag] surface_zigzag_num: " << zigzag_num << endl;
  return 0;
}

int extract_singularity_edge(
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    const boost::unordered_map<size_t,size_t> &face_jump,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::edge2cell_adjacent &ea,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &jump_type_between_tets,
    vector<deque<pair<size_t,size_t> > > &singularity_chain,
    vector<deque<size_t> > &singularity_type,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &face_pair_tangent_jump,
    const matrixst &outside_face_idx)
{
  typedef boost::unordered_map<size_t,size_t>::const_iterator mcit;
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mci;
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oci;
  vector<pair<pair<size_t,size_t>,size_t> > surface_singularity_edge_modificatoin;

  for(oci ci = ortae.e2t_.begin(); ci != ortae.e2t_.end(); ++ci){
      const vector<size_t> &loop = ci->second;
      if(loop.front() != loop.back() && loop.size() != 2){
          cerr << "# error edge. " << endl;
          continue;
        }
      if(loop.front() != -1 && loop.back() != -1) continue; // inner edge
      const size_t type =
          get_edge_type_at_face_a2(ci->first, ortae, fa, ea, face_pair_tangent_jump,
                                   face_jump, outside_face_idx, jump_type_between_tets);
      if(type != TRIVIAL_TYPE){ // not identity
          surface_singularity_edge_modificatoin.push_back(
                make_pair(make_pair(ci->first.first, ci->first.second),type));
        }
    }
  cut_singularity_edge_into_chain(surface_singularity_edge_modificatoin,singularity_chain,singularity_type);
  return 0;
}

int apply_global_U_rotation_with_normal(
    const matrixst &tet,
    const matrixd &node,
    boost::unordered_map<size_t,size_t> &face_jump,
    const matrix<matrixd > &frame_inner,
    const jtf::mesh::face2tet_adjacent & fa,
    const matrixst &outside_face_idx,
    const matrixd & outside_face_normal)
{
  assert(outside_face_normal.size(1) == 3);
  assert(outside_face_normal.size(2) == outside_face_idx.size());

  matrixd temp(3,3);
  vector<pair<double,size_t> > resi(24);

  for(size_t t = 0; t < outside_face_idx.size(); ++t){
      const pair<size_t,size_t> &tet_pair = fa.face2tet_[outside_face_idx[t]];
      assert(tet_pair.second == -1 || tet_pair.first == -1);

      const size_t &tet_idx = (tet_pair.first ==
                               -1)?tet_pair.second:tet_pair.first;

      for(size_t i = 0; i < 24; ++i){
          temp = frame_inner[tet_idx] * type_transition2(i);
          const matrixd U_vec = temp(colon(),0);
          resi[i] = make_pair(dot(U_vec,outside_face_normal(colon(),t)),i);
        }
      sort(resi.begin(),resi.end());

      for(size_t i = 23; i != -1; --i){
          if(resi[i].second < 10){ // choose only one axis rotaion
              face_jump[outside_face_idx[t]] = resi[i].second;
              break;
            }
        }
    }
  return 0;
}

int apply_global_U_rotation(
    const matrixst &tet,
    const matrixd &node,
    boost::unordered_map<size_t,size_t> &face_jump,
    const matrix<matrixd > &frame_inner,
    const jtf::mesh::face2tet_adjacent & fa,
    const matrixst &outside_face_idx)
{
  matrixd face_normal = zeros<double>(3,1);

  matrixd temp(3,3);
  vector<pair<double,size_t> > resi(24);

  for(size_t t = 0; t < outside_face_idx.size(); ++t){
      const pair<size_t,size_t> &tet_pair = fa.face2tet_[outside_face_idx[t]];
      assert(tet_pair.second == -1 || tet_pair.first == -1);
      const size_t &tet_idx = (tet_pair.first == -1)?tet_pair.second:tet_pair.first;
      const vector<size_t> &face = fa.faces_[outside_face_idx[t]];
      jtf::tetmesh::calculate_face_normal_wrapper(tet,node,tet_idx,&face[0],face_normal);

      for(size_t i = 0; i < 24; ++i){
          temp = frame_inner[tet_idx]*type_transition2(i);
          resi[i] = make_pair(dot(face_normal,temp(colon(),0)),i);
        }
      sort(resi.begin(),resi.end());

      for(size_t i = 23; i != -1; --i){
          if(resi[i].second < 10){ // choose one axis rotaion
              face_jump[outside_face_idx[t]] = resi[i].second;
              break;
            }
        }
    }

  return 0;
}

int find_edge_adjacent_faces(matrixst &edge_adjacent_faces,
                             const matrixst &outside_face,
                             const pair<size_t,size_t> &edge)
{
  assert(edge_adjacent_faces.size(1) == 3);
  assert(edge_adjacent_faces.size(2) == 2);
  for(size_t t = 0,j = 0; t < outside_face.size(2); ++t){
      const matrixst &face= outside_face(colon(),t);
      if(find(face.begin(),face.end(),edge.first) != face.end()
         && find(face.begin(),face.end(),edge.second) != face.end())
        {
          edge_adjacent_faces(colon(),j) = face;
          ++j;
          if(j == 2) return  0;
        }
    }
  cerr << "# error can not find two adjacent faces for edge: " << edge.first << "," << edge.second << endl;
  return __LINE__;
}

int calculate_tangent_rotation_for_each_edge(
    boost::unordered_map<pair<size_t,size_t>,size_t> &face_pair_tangent_jump,
    const jtf::mesh::face2tet_adjacent& fa,
    const jtf::mesh::edge2cell_adjacent &ea,
    const matrixst &outside_face,
    const matrixst &outside_face_idx,
    const boost::unordered_map<size_t,size_t> &face_jump,
    const matrix<matrixd > &frame_inner,
    const matrixd &node)
{
  matrix<double> dir(3,1),n1(3,1),n2(3,1);
  matrix<double> rot(3,3),rot2(3,3), RFP(3,3), FP(3,3);
  vector<pair<double,size_t> > candidate;
  for(size_t ei = 0; ei < ea.edges_.size(); ++ei){
      const pair<size_t,size_t>& one_edge_pair = ea.edges_[ei];
      const pair<size_t,size_t>& cell_pair = ea.edge2cell_[ei];
      dir = node(colon(), one_edge_pair.second) - node(colon(), one_edge_pair.first);
      dir /= norm(dir);
      jtf::mesh::cal_face_normal(outside_face(colon(), cell_pair.first), node,  n1);
      jtf::mesh::cal_face_normal(outside_face(colon(), cell_pair.second), node,  n2);
      const double dihedral_angle_degree = calculate_dihedral_angle_degree(n1,n2);
      from_angle_to_rotation_matrix(jtf::math::My_PI()*dihedral_angle_degree/180.,dir, rot);

      const pair<size_t,size_t> & tet_pair0 = fa.face2tet_[outside_face_idx[cell_pair.first]];
      const pair<size_t,size_t> & tet_pairk = fa.face2tet_[outside_face_idx[cell_pair.second]];

      const size_t t0 = (tet_pair0.first==-1?tet_pair0.second:tet_pair0.first);
      const size_t tk = (tet_pairk.first==-1?tet_pairk.second:tet_pairk.first);

      auto it0 = face_jump.find(outside_face_idx[cell_pair.first]);
      auto itk = face_jump.find(outside_face_idx[cell_pair.second]);
      assert(it0 != face_jump.end());
      assert(itk != face_jump.end());

      RFP = rot * frame_inner[t0] * type_transition2(it0->second);
      FP = frame_inner[tk] * type_transition2(itk->second);

      candidate.clear();
      for(size_t di = 0; di < 24; ++di){
          rot2 = type_transition2(di);
          if(fabs(rot2[0]-1) > 1e-6) continue;
          candidate.push_back(make_pair(norm(RFP*rot2-FP),di));
        }
      auto min_e = min_element(candidate.begin(), candidate.end());
      face_pair_tangent_jump[make_pair(outside_face_idx[cell_pair.first],
          outside_face_idx[cell_pair.second])] = min_e->second;
      face_pair_tangent_jump[make_pair(outside_face_idx[cell_pair.second],
          outside_face_idx[cell_pair.first])] = get_trans_type(min_e->second);
    }
  return 0;
}

void calculate_singularity_point_on_surface(
    const matrix<size_t> &outside_face,
    const matrix<double> &node,
    const matrix<size_t> &outside_face_idx,
    const boost::unordered_map<pair<size_t,size_t>,size_t>  &face_pair_tangent_jump)
{
  jtf::mesh::one_ring_face_at_point orfap;
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(outside_face));

  orfap.add_all_faces(outside_face, *ea);
  orfap.sort_int_loop(outside_face, node);

  matrix<double> rot = eye<double>(3);
  matrix<double> point_type(node.size(2),1);
  for(const auto & one_p : orfap.p2f_){
      rot = eye<double>(3);
      const vector<size_t> & face_loop = one_p.second;
      assert(face_loop.front() == face_loop.back());
      for(size_t fi = 0; fi < face_loop.size()-1; ++fi){
          auto it = face_pair_tangent_jump.find(
                make_pair(outside_face_idx[face_loop[fi]],
                outside_face_idx[face_loop[fi+1]]));
          assert(it != face_pair_tangent_jump.end());
          rot = temp(rot * type_transition2(it->second));
        }
      const size_type p_type = type_transition1(rot);
      if(p_type != TRIVIAL_TYPE)
        point_type[one_p.first] = p_type;
      else point_type[one_p.first] = -1.;
    }
  ofstream ofs("surface_singularity_point.vtk");
  tri2vtk(ofs, &node[0], node.size(2), &outside_face[0], outside_face.size(2));
  point_data(ofs, &point_type[0], point_type.size(), "singularity");
}

int remove_surface_zigzag(
    const matrixst &outside_face,
    const matrixst &tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent &fa,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &jump_type_between_tets,
    const matrixst &outside_face_idx,
    const matrix<matrixd > &frame_inner,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    const boost::property_tree::ptree &pt,
    boost::unordered_map<size_t,size_t> &face_jump)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oci;
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mci;

  //map<size_t,size_t> face_jump; // store the face idx, and it's jump accroding to remove surface zigzag
  face_jump.clear();
  boost::unordered_map<pair<size_t,size_t>,size_t> face_pair_tangent_jump;
  vector<deque<pair<size_t,size_t> > > singularity_chain;
  vector<deque<size_t> > singularity_type;

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()) {
      cerr << "# [error] can not create edge2cell_adjacent." << endl;
      return __LINE__;
    }

  // adjust the u vecto align the normal
  apply_global_U_rotation(tet,node,face_jump,frame_inner,fa,outside_face_idx);

  // calculate \Pi_{b,a} for each edge
  calculate_tangent_rotation_for_each_edge(
        face_pair_tangent_jump,fa,*ea, outside_face,outside_face_idx,
        face_jump,frame_inner,node);

  {
    calculate_singularity_point_on_surface(outside_face, node, outside_face_idx, face_pair_tangent_jump);
  }

  extract_singularity_edge(ortae,face_jump, fa,*ea,jump_type_between_tets,singularity_chain,
                           singularity_type, face_pair_tangent_jump,outside_face_idx);

  dump_singularity_chain_to_vtk_2("surface_singularity_edge_before_modify.vtk", node, singularity_chain, singularity_type);

  dump_singularity_to_cylinder("surface_singularity_edge_before_modify.obj", node,singularity_chain,0.002);

  timer tr;
  tr.start();
  solve_face_jump_new(face_jump,singularity_chain,singularity_type,tet,node,ortae,fa,jump_type_between_tets,frame_inner,face_pair_tangent_jump,outside_face);
  tr.finish();
  cerr << "# [surface zigzag time] cost " << tr.result()
       << "ms to remove surface zigzag " << endl;
  //dump_singularity_to_vtk("surface_singularity_after_remove_zigzag_line.vtk",node,singularity_chain);

  extract_singularity_edge(ortae,face_jump, fa, *ea, jump_type_between_tets,singularity_chain,
                           singularity_type, face_pair_tangent_jump,outside_face_idx);

  dump_singularity_chain_to_vtk_2(
        "surface_singularity_after_remove_zigzag_line.vtk",node,
        singularity_chain,singularity_type);

  dump_singularity_to_cylinder("surface_singularity_after_remove_zigzag_line.obj",
                               node,singularity_chain,0.002);

  return 0;

}
