#include <memory>
#include <vector>
#include <iostream>
#include <numeric>
#include <set>

#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/util.h>

#include "hex_process.h"
#include <jtflib/mesh/mesh.h>
#include "../hexmesh/util.h"


using namespace std;
using namespace zjucad::matrix;

//! @brief this function will check each edge of this hex
//  if the hex loop around this edge will be separated and introduced more gaps, then the removing process will
//  effect the topology
static bool is_hex_del_effect_topology(
    const matrixst &hex,
    const matrixd &ori_hex_node,
    const jtf::mesh::face2hex_adjacent &fa,
    const jtf::mesh::one_ring_hex_at_edge & ornae,
    const size_t &hex_id,
    const matrixst &hex_tri_hit)
{
  using namespace jtf::hexmesh;
  matrixst faces;
  jtf::mesh::get_faces_for_one_hex(hex(colon(),hex_id),faces);

  // if this hex is inside, then topology must be changed when this hex is deleted
  matrixst face_state(6);
  for(size_t t = 0; t < faces.size(2); ++t){
    const pair<size_t,size_t> both_hex = fa.query(&faces(0,t));
    if(fa.is_outside_face(both_hex)) face_state[t] = 0; // outside face
    else{
      const size_t other_hex_id = both_hex.first + both_hex.second - hex_id;
      if(hex_tri_hit[other_hex_id] == 0) face_state[t] = 0;
      else
        face_state[t] = 1; //inside face
    }
  }

  const size_t total_state = accumulate(face_state.begin(),face_state.end(),0);
  // all faces are outside, means this hex is an isolated one, this case will happen when it's arounding is marked as out, then removing becomes acceptable
  if (total_state == 0) return false;

  if( total_state == face_state.size()) return true; // all faces are inside, then delete this hex must introduce a hole

  vector<size_t> potential_problem_face;
  for(size_t t = 0; t < 6; t += 2){
    if(face_state[t] + face_state[t+1] == 1) {
      potential_problem_face.push_back((face_state[t] == 1)?face_state[t]:face_state[t+1]);
    }
  }
  if(potential_problem_face.empty()) return true; // means this hex connect with others only throung pair face, removing it must introduce a hole

  // each non-pair inside face, check the edges

  set<pair<size_t,size_t> > edges;
  pair<size_t,size_t> edge;

  // gather all "may introduce problems" edges
  for(size_t t = 0; t < potential_problem_face.size(); ++t){
    for(size_t i = 0; i < faces.size(1); ++i){
      edge.first = faces(i,potential_problem_face[t]);
      edge.second = faces((i+1)%4,potential_problem_face[t]);
      if(edge.first > edge.second) swap(edge.first,edge.second);
      edges.insert(edge);
    }
  }

  typedef jtf::mesh::one_ring_hex_at_edge::e2hex_type::const_iterator oecit;
  // for each edge, we need to determin:
  // the num of how much segments will be separated around each edge after this hex is deleted
  for(set<pair<size_t,size_t> >::const_iterator scit = edges.begin();
      scit != edges.end(); ++scit){
    oecit it = ornae.e2h_.find(*scit);
    if(it == ornae.e2h_.end()) it = ornae.e2h_.find(make_pair(scit->second,scit->first));
    if(it == ornae.e2h_.end()) {
      cerr << "# [error] can not find edge <" << scit->first << ","<< scit->second << "> in one_ring_neigubour_at_edge." << endl;
      continue;
    }
    const vector<size_t> & loop = it->second;
    if(loop.front() != -1 && loop.back() != -1) continue; // if this loop is closed , delete this outside hex may not effect the topology
    vector<size_t>::const_iterator idx_it = find(loop.begin(),loop.end(),hex_id);
    if(idx_it == loop.end()) {
      cerr << "# [error] edge <" << scit->first << ","<< scit->second << "> is not one edge of hex "<< hex_id << "." << endl;
      continue;
    }

    if(idx_it == loop.begin() + 1 || idx_it == loop.end() -1) continue; // will not effect topology
    return true; // will introduce new gap in this loop
  }
  return false;
}


static int find_outside_hex_idx(const jtf::mesh::face2hex_adjacent & fa,
                                const matrixst & outside_face_idx,
                                matrixst &outside_hex_idx)
{
  set<size_t> outside_hex;
  for(size_t t = 0; t < outside_face_idx.size(); ++t){
    const pair<size_t,size_t> & two_hex = fa.face2hex_[outside_face_idx[t]];
    assert((two_hex.first == -1) ^ (two_hex.second == -1));
    outside_hex.insert( (two_hex.first == -1)?two_hex.second:two_hex.first);
  }

  outside_hex_idx.resize(outside_hex.size());
  copy(outside_hex.begin(),outside_hex.end(),outside_hex_idx.begin());
  return 0;
}

int pare_hexmesh_from_surface(const matrixst &hex,
                              const matrixd &node,
                              matrixst &hex_new)
//matrixd &node_new)
{
  using namespace jtf::hexmesh;
  unique_ptr<jtf::mesh::face2hex_adjacent> fa(jtf::mesh::face2hex_adjacent::create(hex));
  if(fa.get() == 0) return __LINE__;
  matrixst outside_face_idx;
  jtf::mesh::get_outside_face_idx(*fa,outside_face_idx);
  matrixst outside_hex_idx;
  find_outside_hex_idx(*fa,outside_face_idx,outside_hex_idx);

  cerr << "# [info] total hex: " << hex.size(2) << endl;
  cerr << "# [info] outside hex: " << outside_hex_idx.size() << endl;


  vector<bool>  remove_flag(hex.size(2),false);

  //cerr << "# [info] remove_flag num " << remove_flag.size() << endl;
  for(size_t t = 0; t < outside_hex_idx.size(); ++t) remove_flag[outside_hex_idx[t]] = true;

  hex_new.resize(hex.size(1),hex.size(2) - outside_hex_idx.size());
  for(size_t t = 0,j = 0; t < remove_flag.size(); ++t){
    if(!remove_flag[t]){ // do not need to remove
      hex_new(colon(),j) = hex(colon(),t);
      ++j;
    }
  }
  return 0;
}

#if 0 // no used
//! @brief as the following fig shows, for each point of hex, we find a direction from triangle bary center to this point,
//  then if dot(direction,face_normal) > 0, this point is outside along the normal
//       if dot(direction,face_normal) < 0, this point is outside along the anti-normal
//  if all points are with the same sign, this hex is outside the triangle,
//  else, this hex intersect with the triangle
//    e-----f
//   /|    /|
//  a-|---b |
//  | g---|-h
//  |/    |/
//  c-----d
//  \  A  /
//   \/ \/
//   /\ /\
//  /  I  \
// B-------C
//
size_t cal_adjacent_state(const matrixst &each_hex,
                          const matrixd &hex_node,
                          const matrixst &face,
                          const matrixd &face_normal,
                          const matrixd &tet_node)
{
  matrixd tri_bary = zeros<double>(3,1);

  for(size_t t = 0; t < face.size(); ++t)
    tri_bary  += tet_node(colon(),face[t]);
  tri_bary /= face.size();

  vector<size_t> dir_state(each_hex.size());

  matrixd dir = zeros<double>(3,1);
  for(size_t t = 0; t < each_hex.size(); ++t){
    dir = hex_node(colon(),each_hex[t]) - tri_bary;
    if(dot(dir,face_normal) > 0) dir_state[t] = 0;
    else
      dir_state[t] = 1;
  }

  size_t state = accumulate(dir_state.begin(),dir_state.end(),0);
  if(state == 0) // 0: this hex outside along the triangle normal
    return 0;
  else if(state == each_hex.size()) // 8: this hex is outside along the anti-normal
    return 1;
  else return 2; // intersect
}
#endif

int pare_hex_outside_tet_surface(const matrixst &tet,
                                 const matrixd &tet_node,
                                 const matrixst &ori_hex,
                                 const matrixd &ori_hex_node,
                                 hit_func &func, // func is used to varify whether a hex around surface is in or out
                                 matrixst &new_hex,
                                 matrixd &new_hex_node)
{
  using namespace jtf::hexmesh;
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  unique_ptr<jtf::mesh::face2hex_adjacent> fa_hex(jtf::mesh::face2hex_adjacent::create(ori_hex));
  matrixst outside_face,outside_face_idx;
  get_outside_face(*fa,outside_face);
  get_outside_face_idx(*fa,outside_face_idx);

  matrixd face_normal = zeros<double>(3,outside_face.size(2));
  jtf::mesh::cal_face_normal(outside_face,tet_node,face_normal);
  jtf::tetmesh::orient_face_normal_outside_tetmesh(tet,tet_node,outside_face,outside_face_idx,*fa,face_normal);

  // 0: outside along the triangle normal
  // 1: inside the tet model
  // 2: intersect
  matrixst hex_tri_hit(ori_hex.size(2));


  func.hex_tri_hit_func(outside_face,tet_node,ori_hex,ori_hex_node,face_normal,hex_tri_hit);

  const size_t outside_hex_num = count(hex_tri_hit.begin(),hex_tri_hit.end(),0);
  const size_t left_hex_num = ori_hex.size(2) - outside_hex_num;

#if 1
  const size_t inside_hex_num = count(hex_tri_hit.begin(),hex_tri_hit.end(),1);
  const size_t surface_hex_num = count(hex_tri_hit.begin(),hex_tri_hit.end(),2);
  cerr << "# [info] outside hex num = " << outside_hex_num << endl;
  cerr << "# [info] inside hex num = " << inside_hex_num << endl;
  cerr << "# [info] surface hex num = " << surface_hex_num << endl;
#endif

  jtf::mesh::one_ring_hex_at_edge ornae;
  for(size_t t = 0; t < ori_hex.size(2); ++t){
    ornae.add_hex(ori_hex(colon(),t),ori_hex_node,*fa_hex);
  }
  ornae.sort_into_loop(ori_hex,ori_hex_node,*fa_hex);
  // to determin the surface hex should be removed or not
  for(size_t t = 0; t < ori_hex.size(2); ++t)
    if(hex_tri_hit[t] == 2){
      if(is_hex_del_effect_topology(ori_hex,ori_hex_node,*fa_hex,ornae,t,hex_tri_hit)) hex_tri_hit[t] = 1;
      else hex_tri_hit[t] = 0;
    }

  new_hex.resize(8,left_hex_num);
  for(size_t t = 0,i = 0; t < ori_hex.size(2); ++t){
    if(hex_tri_hit[t] != 0)
      new_hex(colon(),i++) = ori_hex(colon(),t);
  }

  new_hex_node = ori_hex_node;

  return 0;
}


size_t trivial_hit_func::ray_casting_determin(const matrixd &point,
                                              const matrixd &direction,
                                              const matrixst &outside_tet_face,
                                              const matrixd &face_normal,
                                              const matrixd &tet_node)
{
  // bounding_box: p_min, p_max
  if(bounding_box.size() != 6) // means this bounding box need to be recalculed
  {
    bounding_box.resize(3,2);
    bounding_box(colon(),0) = tet_node(colon(),outside_tet_face[0]);
    bounding_box(colon(),1) = tet_node(colon(),outside_tet_face[0]);

    set<size_t> outside_point(outside_tet_face.begin(),outside_tet_face.end());
    for(set<size_t>::const_iterator sit = outside_point.begin();
        sit != outside_point.end(); ++sit){
      for(size_t i = 0; i < 3; ++i){
        if(tet_node(i,*sit) < bounding_box(i,0)) bounding_box(i,0) = tet_node(i,*sit);
        if(tet_node(i,*sit) > bounding_box(i,1)) bounding_box(i,1) = tet_node(i,*sit);
      }
    }
  }

  // outside the bounding box
  if(point[0] > bounding_box(0,1) || point[1] > bounding_box(1,1) || point[2] > bounding_box(2,1)) return 0;
  if(point[0] < bounding_box(0,0) || point[1] < bounding_box(1,0) || point[2] < bounding_box(2,0)) return 0;

  // calculate the intersect point between the line and triangle
  // line: p = point + direction * t; t\in [-inf,+inf]
  // tri: n * (p - p0) = 0; n is the surface normal and p0 is one point of triangle
  // the interected point is: t = dot(n,p0-point)/dot(n,direction);
  size_t count_intersect = 0;
  matrixd inter_point = zeros<double>(3,1);
  for(size_t t = 0; t < outside_tet_face.size(2); ++t){
    const double nd = dot(face_normal(colon(),t),direction);
    if(fabs(nd) < 1e-8) continue; // this means t = inf, and this line is parallel with the triangle
    const double dt = dot(face_normal(colon(),t),
                          tet_node(colon(),outside_tet_face(0,t)) - point) / nd;
    if(dt > 0) continue;
    inter_point = point + direction * dt;
    if(is_point_in_triangle(inter_point,outside_tet_face(colon(),t),tet_node)) ++count_intersect;
  }

  if(count_intersect % 2 == 1) // inside the tet
    return 1;
  else
    return 0;
}


bool trivial_hit_func::is_point_in_triangle(const matrixd & point_,
                                            const matrixst &triangle,
                                            const matrixd &node)const
{
  const matrixd v0 = node(colon(),triangle[2]) - node(colon(),triangle[0]);
  const matrixd v1 = node(colon(),triangle[1]) - node(colon(),triangle[0]);
  const matrixd v2 = point_ - node(colon(),triangle[0]);

  const double dot00 = dot(v0, v0);
  const double dot01 = dot(v0, v1);
  const double dot02 = dot(v0, v2);
  const double dot11 = dot(v1, v1);
  const double dot12 = dot(v1, v2);

  const double Denom = (dot00 * dot11 - dot01 * dot01);
  if(fabs(Denom) < 1e-8) {
    cerr << "# [error] this triangle is degenerated." << endl;
    return false;
  }
  const double u = (dot11 * dot02 - dot01 * dot12) /Denom;
  const double v = (dot00 * dot12 - dot01 * dot02) /Denom;

  return (u > 0 || fabs(u) < 1e-8) && (v > 0 || fabs(v) < 1e-8) && (u + v < 1);
}

int trivial_hit_func:: hex_tri_hit_func(const matrixst &outside_tet_face,
                                        const matrixd &tet_node,
                                        const matrixst &hex,
                                        const matrixd &hex_node,
                                        const matrixd &face_normal,
                                        matrixst &hex_hit)
{
  // node num may be redundant
  set<size_t> hex_point(hex.begin(),hex.end());
  assert(hex_point.size() != hex.size());

  map<size_t,size_t> point_hit_map; // store the point is inside the tet or not, 0: out, 1: inside
  const double direction[] = {1,0,0}; // along the x axis
  itr_matrix<const double*> dir_mat(3,1,direction);
  for(set<size_t>::const_iterator scit = hex_point.begin();
      scit != hex_point.end(); ++scit){
    point_hit_map[*scit] = ray_casting_determin(hex_node(colon(),*scit),dir_mat,outside_tet_face,face_normal,tet_node);
  }

  hex_hit.resize(hex.size(2));
  typedef map<size_t,size_t>::const_iterator mcit;

  vector<size_t> each_hex_flag(hex.size(1),0);
  for(size_t t = 0; t < hex.size(2); ++t){
    for(size_t i = 0; i < hex.size(1); ++i){
      mcit it = point_hit_map.find(hex(i,t));
      if(it == point_hit_map.end()){
        cerr << "# [error] can not find the point in the map." << endl;
        return __LINE__;
      }else{
        each_hex_flag[i] = it->second;
      }
    }// end for each hex

    const size_t state = accumulate(each_hex_flag.begin(),each_hex_flag.end(),0);
    if(state == 0) {// this hex is outside
      hex_hit[t] = 0;
    }else if(state == each_hex_flag.size()) {
      hex_hit[t] = 1; // inside the model
    }else  hex_hit[t] = 2; // interect with tet surface
  }
  return 0;
}
