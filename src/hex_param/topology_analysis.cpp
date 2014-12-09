#include <set>
#include <stack>
#include <fstream>
#include <numeric>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/filesystem.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>


#include "topology_analysis.h"
#include "topology_operation.h"
#include "remove_surface_wedge.h"
#include "common.h"
#include "wyz_format.h"
#include "io.h"
#include "hex_param.h"
#include "find_singularities.h"
#include "singularity_adjustment.h"
#include "cut_tet.h"
#include "solution_searching.h"
#include "../common/vtk.h"
#include "../common/util.h"
#include "../common/zyz.h"
#include "../common/IO.h"
#include "../common/extract_loop_from_undirected_edges.h"
#include "../common/transition_type.h"
#include "../common/transition.h"
#include <jtflib/util/container_operation.h>

#include "../tetmesh/tetmesh.h"
#include "../numeric/util.h"
#include "../common/visualize_tool.h"
#include "../tetmesh_refine/tetmesh_refine.h"
#include "../tetmesh/util.h"
#include "../numeric/util.h"
#include "../equation_graph/equation_graph.h"

using namespace std;
using namespace zjucad::matrix;
using namespace boost;


#define debug
#ifdef debug
static int num = 0;
//vector<deque<pair<size_t,size_t> > > singularity_edges_of_cut_tet(1);
//vector<deque<size_t> > singularity_type_of_cut_tet(1);

//vector<deque<pair<size_t,size_t> > > compound_edges(1);
//vector<deque<pair<size_t,size_t> > > loop_edges;
//vector<size_t> jump_faces;

#endif

extern int boundary_issue_num;
extern int degree_issue_num;

int cal_singularity_type_at_given_tet(
    const vector<size_t> & tet_loop,
    const size_t &begin_tet,
    const vector<size_t> &rot_type,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & rot_type_map_idx,
    std::vector<size_t> * reordered_tet_loop_ptr = 0);

template <typename T1, typename T2>
int map_insert(map<T1, T2> & map_, const  T1 & pair_, const T2 & type_)
{
  typename map<T1,T2>::const_iterator mcit = map_.find(pair_);
  if(mcit == map_.end())
    map_.insert(make_pair(pair_,type_));
  else
    return __LINE__;
  return 0;
}

template <typename T1, typename T2>
int map_insert(boost::unordered_map<T1, T2> & map_, const  T1 & pair_, const T2 & type_)
{
  typename boost::unordered_map<T1,T2>::const_iterator mcit = map_.find(pair_);
  if(mcit == map_.end())
    map_.insert(make_pair(pair_,type_));
  else
    return __LINE__;
  return 0;
}

template <typename T1, typename T2>
int map_erase(map<T1, T2> & map_, const T1 & pair_, const T2 & type_)
{
  typename map<T1,T2>::iterator mcit = map_.find(pair_);
  if(mcit == map_.end())
    return __LINE__;
  else {
      if(mcit->second !=  type_)
        return __LINE__;
      map_.erase(mcit);
    }
  return 0;
}

template <typename T1, typename T2>
int map_erase(boost::unordered_map<T1, T2> & map_, const T1 & pair_, const T2 & type_)
{
  typename boost::unordered_map<T1,T2>::iterator mcit = map_.find(pair_);
  if(mcit == map_.end())
    return __LINE__;
  else {
      if(mcit->second !=  type_)
        return __LINE__;
      map_.erase(mcit);
    }
  return 0;
}

struct state{
  pair<size_t,size_t> tet_pair_;
  vector<size_t> possible_rots_;
  size_t current_rot_idx_;
  //size_t edge_num_;
  deque<pair<size_t,size_t> >  edges_need_to_handle_;
  vector<double> rot0_,rot1_;
  string undo_process_type_;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & tet_pair_;
    ar & possible_rots_;
    ar & current_rot_idx_;
    ar & edges_need_to_handle_;
    ar & rot0_;
    ar & rot1_;
    ar & undo_process_type_;
  }

  //  matrix_expression<double> & rot0()
  //  {
  //    itr_matrix<double*> rot_itr0(3,3,&rot0_[0]);
  //    return rot_itr0;
  //  }

  //  matrix_expression<double> & rot1()
  //  {
  //    itr_matrix<double*>rot_itr1(3,3,&rot1_[0]);
  //    return rot_itr1;
  //  }
  //  //! WARNING !!! this function must be bind with undo_process, since this
  //  // functin can only be called while searching back
  //  boost::function< int (const pair<size_t,size_t>&, const size_t& ) > undo_process_;
};



static int recover_cut_tet_mesh(
    matrixst & cut_mesh,
    const vector<set<matrixst::value_type *> > acc_table)
{
  cut_mesh(colon()) = colon(0, cut_mesh.size()-1);
  for(size_t t = 0; t < acc_table.size(); ++t){
      if(acc_table[t].empty()) continue;
      const set<matrixst::value_type *> & stv = acc_table[t];
      for(set<matrixst::value_type*>::const_iterator scit = stv.begin();
          scit != stv.end(); ++scit){
          *(*scit) = t;
        }
    }
  return 0;
}


int cal_frame_geometry_priority_queue(
    const pair<size_t,size_t>& tet_pair,
    const zjucad::matrix::matrix<matrixd > & frame,
    vector<pair<double,size_t> > & priority_queue)
{
  priority_queue.resize(24); // 24 rotations
  for(size_t t = 0; t < 24; ++t){
      priority_queue[t] =
          make_pair(norm(frame[tet_pair.first] * type_transition2(t)
                    - frame[tet_pair.second]), t);
    }
  sort(priority_queue.begin(), priority_queue.end());
  return 0;
}

static  int convert_face_normal_to_cut_face_normal(
    const matrixd & face_normal,
    const matrixst & outside_face,
    const matrixst & outside_face_cut,
    const matrixst & cut_tet2tet,
    const matrixst & outside_face_idx,
    const matrixst & outside_face_idx_in_cut,
    boost::unordered_map<size_t, matrixd > &outside_face_normal_cut,
    boost::unordered_map<size_t,size_t> & orig_face_idx2cut_face)
{
  assert(face_normal.size(1) == 3);
  assert(face_normal.size(2) == outside_face.size(2));
  assert(outside_face_cut.size(2) == outside_face_idx_in_cut.size());
  matrixst sort_outside_face = outside_face;
  for(size_t fi = 0; fi < sort_outside_face.size(2); ++fi){
      sort(sort_outside_face(colon(),fi).begin(),
           sort_outside_face(colon(),fi).end());
    }

  for(size_t fi = 0; fi < outside_face_cut.size(2); ++fi){
      matrixst orig_face = cut_tet2tet(outside_face_cut(colon(),fi));
      sort(orig_face.begin(), orig_face.end());


      for(size_t ofi = 0; ofi < sort_outside_face.size(2); ++ofi){

          if(norm( orig_face - sort_outside_face(colon(), ofi)) < 1e-6){
              outside_face_normal_cut[outside_face_idx_in_cut[fi]] =
                  face_normal(colon(),ofi);

              orig_face_idx2cut_face[outside_face_idx[ofi]]
                  = outside_face_idx_in_cut[fi];
              break;
            }
        }
      //cerr << fi << endl;
    }
  return 0;
}



int aggressive_assign_transition_with_origin_type(
    const matrixst & orig_tet,
    const matrixd & node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type,
    boost::property_tree::ptree & pt,
    const bool no_surface,
    matrixd * zyz_ptr)
{
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(orig_tet));
  pt.put("surface_type.desc","surface normal align type");

  bool is_restrcited_type = false;
  boost::unordered_map<size_t,size_t> surface_type;
  if(load_surface_type(pt.get<string>("surface_type.value").c_str(),
                       surface_type)){
      load_surface_type(pt.get<string>("surface_type.value").c_str(),
                        surface_type,fa.get());
    }
  else{
      is_restrcited_type = true;
    }

  //  {// hack
  //    surface_type.clear();
  //    for(size_t fi = 0; fi < fa->face2tet_.size(); ++fi){
  //      if(fa->is_outside_face(fa->face2tet_[fi])){
  //        surface_type[fi] = 0;
  //      }
  //    }
  //  }

  // build mst and orient the input types
  matrixst cut_tet;
  matrix<matrixd> frame;

  if(zyz_ptr){
      frame.resize(zyz_ptr->size(2));
      for(size_t ti = 0; ti < zyz_ptr->size(2); ++ti){
          frame[ti].resize(3,3);
          zyz_angle_2_rotation_matrix1(&(*zyz_ptr)(0,ti), &frame[ti][0]);
        }
    }

  cut_tetmesh(
        orig_tet, node, *fa, inner_face_type, surface_type, is_restrcited_type,
        cut_tet, (zyz_ptr == 0?0:&frame));

#if 1
  {
    if(frame.size() != 0){
        matrixd zyz(3, frame.size());
        for(size_t fi = 0; fi < frame.size(); ++fi)
          rotation_matrix_2_zyz_angle(&frame[fi][0], &zyz(0,fi),0);

        jtf::mesh::write_matrix("frame_after_aligned.zyz",zyz);
      }
  }
#endif


  /// WARNING: NOTICE that the surface used in aggressive_assign_transition with
  // type is different with the input, in the following function , surface type
  // is 0-5 which stands for +u,-u,+v,-v,+w,-w
  // so here needs to convert the surface type which is represented by rotation
  // u to normal to real axis
  if(!is_restrcited_type){
      for(boost::unordered_map<size_t,size_t>::iterator bumit = surface_type.begin();
          bumit != surface_type.end(); ++bumit)
        bumit->second = convert_surface_rotation_to_axis_type(bumit->second)/2;
    }

  aggressive_assign_transition_with_type(
        orig_tet, cut_tet, node, inner_face_type, surface_type,
        no_surface, (frame.size() == 0)?0:&frame);
  return 0;
}

int aggressive_assign_transition(
    const matrixst & orig_tet,
    const matrixd & node,
    const matrixd & zyz,
    const bool no_surface,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_type)
{
  jtf::tet_mesh tm(orig_tet, node);
  tetmesh_cutter tmc(tm);
  matrixd new_zyz = zyz;
  tmc.cut(new_zyz);

  jtf::mesh::write_matrix("frame_after_aligned.zyz",new_zyz);

  matrix<matrix<double> > frame_inner;
  zyz2frame(new_zyz, frame_inner);

  matrixd rot = eye<double>(3);
  for(size_t fi = 0; fi < tm.fa_->face2tet_.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = tm.fa_->face2tet_[fi];
      if(tm.fa_->is_outside_face(tet_pair)) continue;
      assert(tet_pair.first != -1 && tet_pair.second != -1);
      get_best_alignment(&frame_inner[tet_pair.first](0,0),
          &frame_inner[tet_pair.second](0,0), &rot[0]);
      face_type[tet_pair] = type_transition1(rot);
      face_type[make_pair(tet_pair.second, tet_pair.first)]
          = get_trans_type(type_transition1(rot));
    }

  boost::unordered_map<size_t,size_t> surface_type;
  {
    vector<matrixd> frame_vec(frame_inner.begin(), frame_inner.end());
    vector<size_t> surface_normal_align_type;
    get_surface_normal_align_type(surface_normal_align_type, orig_tet,
                                  node,tm.outside_face_, tm.outside_face_idx_, *tm.fa_,
                                  frame_vec );

    assert(surface_normal_align_type.size() == tm.outside_face_.size(2));

    for(size_t fi = 0; fi < surface_normal_align_type.size(); ++fi){
        surface_type[tm.outside_face_idx_[fi]] = surface_normal_align_type[fi];
      }
  }

  aggressive_assign_transition_with_type(
        orig_tet, tmc.cut_tm_.mesh_, node, face_type, surface_type,
        no_surface, &frame_inner);

  return 0;
}

/**
 * @brief
 *
 * @param filename
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param integer_constraints
 * @return int
 */
int load_integer_constraint(
    const char* filename,
    std::map<std::pair<size_t,size_t>,size_t> & integer_constraints)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open integer constraint file." << endl;
      return __LINE__;
    }
  integer_constraints.clear();
  size_t point_idx = -1,fix_uvw = -1;
  string trash;
  size_t linked_points_num = 0;
  size_t link_p = -1;
  size_t group_idx = 0;
  while(!ifs.eof()){
      ifs >> trash >> point_idx >> trash >> fix_uvw; // "point 0 adis u"
      ifs >> linked_points_num;
      if(point_idx == -1 || fix_uvw == -1){
          cerr << "# [info] file seems empty." << endl;
          return 0;
        }
      set<size_t> linked_points;
      for(size_t t = 0 ; t < linked_points_num; ++t){
          ifs >> link_p;
          linked_points.insert(link_p);
        }
      if(!linked_points.empty()){
          // inner_integer_constraints[make_pair(point_idx,fix_uvw)] = linked_points;
          for(set<size_t>::const_iterator sit = linked_points.begin();
              sit != linked_points.end(); ++sit){
              integer_constraints[make_pair(*sit,fix_uvw)] = group_idx;
            }
          ++group_idx;
        }
    }
  return 0;
}


int dump_singularity_degenerated_points(
    const char * filename,
    const std::vector<std::vector<size_t> > & degenerated_points)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      cerr << "# [error] can not open singularity degeneration points file." << endl;
      return __LINE__;
    }
  ofs << degenerated_points.size() << endl;
  for(size_t t = 0; t < degenerated_points.size(); ++t){
      const vector<size_t> & one_group = degenerated_points[t];
      ofs << one_group.size() << endl;
      copy(one_group.begin(),one_group.end(),ostream_iterator<size_t>(ofs," "));
      ofs << endl;
    }
  return 0;
}

int load_singularity_degenerated_points(
    const char * filename,
    std::vector<std::vector<size_t> > & degenerated_points)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open singularity degeneration points file." << endl;
      return __LINE__;
    }
  size_t points_vec_size = 0;
  size_t one_group_size = 0;

  ifs >> points_vec_size;

  degenerated_points.resize(points_vec_size);

  for(size_t t = 0; t < degenerated_points.size(); ++t){
      vector<size_t> & one_group = degenerated_points[t];
      ifs >> one_group_size;
      one_group.resize(one_group_size);
      for(size_t i = 0; i < one_group_size; ++i){
          ifs >> one_group[i];
        }
    }
  return 0;
}

/**
 * @brief
 *
 * @param filename
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param integer_constraints
 * @param group
 * @return int
 */
int load_integer_constraint_group_info(
    const char* filename,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & integer_constraints,
    std::vector<std::set<size_t> > & group)
{
  integer_constraints.clear();
  group.clear();

  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open integer constraint file." << endl;
      return __LINE__;
    }
  integer_constraints.clear();
  size_t point_idx = -1,fix_uvw = -1;
  string trash;
  size_t linked_points_num = 0;
  size_t link_p = -1;
  size_t group_idx = 0;
  while(!ifs.eof()){
      ifs >> trash >> point_idx >> trash >> fix_uvw; // "point 0 adis u"
      ifs >> linked_points_num;
      if(point_idx == -1 || fix_uvw == -1){
          cerr << "# [info] file seems empty." << endl;
          return 0;
        }
      set<size_t> linked_points;
      for(size_t t = 0 ; t < linked_points_num; ++t){
          ifs >> link_p;
          linked_points.insert(link_p);
        }
      group.push_back(linked_points);
      if(!linked_points.empty()){
          // inner_integer_constraints[make_pair(point_idx,fix_uvw)] = linked_points;
          for(set<size_t>::const_iterator sit = linked_points.begin();
              sit != linked_points.end(); ++sit){
              integer_constraints[make_pair(*sit,fix_uvw)] = group_idx;
            }
          ++group_idx;
        }
    }

  return 0;
}

/**
 * @brief
 *
 * @param filename
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param integer_cons
 * @return int
 */
int load_integer_constraint_linking_info(
    const char* filename,
    std::map<std::pair<size_t,size_t>,std::set<size_t> > & integer_cons)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open linking file." << endl;
      return __LINE__;
    }
  string trash;
  size_t point_idx,axis;
  size_t linked_num = 0;
  size_t linked_point = -1;
  while(!ifs.eof())
    {
      ifs >> trash >> point_idx >> trash >> axis; // point n axis u
      ifs >> linked_num;
      for(size_t i = 0; i < linked_num; ++i){
          ifs >> linked_point;
          integer_cons[make_pair(point_idx,axis)].insert(linked_point);
        }
    }
  return 0;
}

template <typename E>
bool is_mat_equal(const zjucad::matrix::matrix_expression<E> &a,
                  const zjucad::matrix::matrix_expression<E> &b)
{
  if(a().size(1) != b().size(1) ||
     a().size(2) != b().size(2))
    return false;
  for(size_t t = 0; t < a().size(); ++t){
      if(a()[t] > b()[t] || a()[t] < b()[t])
        return false;
    }
  return true;
}

/**
 * @brief
 *
 * @param faces
 * @param outside_face_type
 * @return bool
 */
bool is_patch_with_same_axis_fix(const std::vector<size_t> & faces,
                                 const matrixst &outside_face_type)
{
  set<size_t> fix_axes;
  for(size_t t = 0; t < faces.size(); ++t){
      fix_axes.insert(get_surface_fix_axis(outside_face_type[faces[t]]));
    }
  if(fix_axes.size() == 1) return true;
  return false;
}

/**
 * @brief
 *
 * @param degenerated_points
 * @param linked_near_surface_chain
 * @return bool
 */
bool is_near_surface_degenerated_case(
    const std::vector<size_t> & degenerated_points,
    std::vector<size_t> & linked_near_surface_chain)
{
  cerr << "# [error] is_near_surface_degenerated_case not finished." << endl;
  return false;
}

/**
 * @brief
 *
 * @param faces
 * @return bool
 */
bool is_simple_connected_patch(const matrixst & faces)
{
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(faces));
  set<size_t> points(faces.begin(),faces.end());
  const size_t vertex_num = points.size();
  const size_t edge_num = ea->edges_.size();
  const size_t face_num = faces.size(2);
  if(vertex_num + face_num - edge_num != 1) // V+F-E != 1 is not a simple connected  patch
    return false;
  else
    return true;
}

/**
 * @brief
 *
 * @param node
 * @param boundary0
 * @param boundary1
 * @param ea
 * @param outside_face
 * @param smaller_face_patch_with_boundary
 * @return int
 */
int find_face_patch_with_two_closed_boundary(
    const matrixd & node,
    const std::vector<size_t> & boundary0,
    const std::vector<size_t> & boundary1,
    const jtf::mesh::edge2cell_adjacent &ea,
    const matrixst &outside_face,
    std::vector<size_t> &smaller_face_patch_with_boundary)
{
  if(boundary0.front() != boundary1.front() &&
     boundary0.front() != boundary1.back()){
      cerr << "# [error] invalid open boundary." << endl;
      return __LINE__;
    }

  if(boundary0.back() != boundary1.front() &&
     boundary0.back() != boundary1.back()){
      cerr << "# [error] invalid open boundary." << endl;
      return __LINE__;
    }

  set<pair<size_t,size_t> > boundary_edges;
  for(size_t t = 0; t < boundary0.size()-1; ++t){
      if(boundary0[t] > boundary0[t+1])
        boundary_edges.insert(make_pair(boundary0[t+1],boundary0[t]));
      else
        boundary_edges.insert(make_pair(boundary0[t],boundary0[t+1]));
    }
  for(size_t t = 0; t < boundary1.size()-1; ++t){
      if(boundary1[t] > boundary1[t+1])
        boundary_edges.insert(make_pair(boundary1[t+1],boundary1[t]));
      else
        boundary_edges.insert(make_pair(boundary1[t],boundary1[t+1]));
    }

  const pair<size_t,size_t> start_face_pair =
      ea.query(boundary_edges.begin()->first,
               boundary_edges.begin()->second);

  if(start_face_pair.first == -1 || start_face_pair.second == -1)
    {
      cerr << "# [error] meet boundary edges." << endl;
      return __LINE__;
    }

  size_t other_points [2] ;

  other_points[0] = outside_face(0,start_face_pair.first) +
      outside_face(1,start_face_pair.first) +
      outside_face(2,start_face_pair.first) -
      boundary_edges.begin()->first -
      boundary_edges.begin()->second;

  other_points[1] = outside_face(0,start_face_pair.second) +
      outside_face(1,start_face_pair.second) +
      outside_face(2,start_face_pair.second) -
      boundary_edges.begin()->first -
      boundary_edges.begin()->second;

  set<size_t> boundary_points(boundary0.begin(),boundary0.end());
  boundary_points.insert(boundary1.begin(),boundary1.end());
  matrixd average_node = zeros<double>(3,1);
  for(set<size_t>::const_iterator scit =boundary_points.begin();
      scit != boundary_points.end(); ++scit){
      average_node +=  node(colon(),*scit);
    }
  average_node /= boundary_points.size();

  size_t start_face_idx = -1;
  if(norm(node(colon(),other_points[0]) - average_node)
     < norm(node(colon(),other_points[1]) - average_node))
    start_face_idx = start_face_pair.first;
  else
    start_face_idx = start_face_pair.second;

  vector<bool> face_visited_flag(outside_face.size(2),false);

  stack<size_t> grow_faces;
  grow_faces.push(start_face_idx);
  while(!grow_faces.empty()){
      size_t current_face_idx = grow_faces.top();
      grow_faces.pop();
      //assert(face_visited_flag[current_face_idx] == false);
      if(face_visited_flag[current_face_idx]) continue;
      face_visited_flag[current_face_idx] = true;
      for(size_t t = 0; t < 3; ++t){
          pair<size_t,size_t> edge(
                outside_face(t,current_face_idx),
                outside_face((t+1)%outside_face.size(1),current_face_idx));
          const size_t edge_idx = ea.get_edge_idx(edge.first,edge.second);
          if(edge_idx == -1){
              cerr << "# [error] strange can not find such edge <"
                   << outside_face(t,current_face_idx) << ","
                   << outside_face((t+1)%outside_face.size(1),current_face_idx)
                   << ">" << endl;
              return __LINE__;
            }
          if(edge.first > edge.second)
            swap(edge.first,edge.second);
          if(find(boundary_edges.begin(),boundary_edges.end(),edge) !=
             boundary_edges.end()) continue; // meet the boundary
          const pair<size_t,size_t> & face_pair = ea.edge2cell_[edge_idx];
          assert(face_pair.first != -1 && face_pair.second != -1);
          assert(face_pair.first == current_face_idx || face_pair.second == current_face_idx);
          const size_t other_face_idx = face_pair.first + face_pair.second - current_face_idx;

          if(face_visited_flag[other_face_idx] == true) continue;
          grow_faces.push(other_face_idx);
        }
    }


  {// check wether the visited face can make a simple connected patch,
    // if original mesh is not a closed mesh, it must be divied into simple connected one and one with a hole
    // if original mesh is a closed mesh, the two patches are both simple connected, here we take an assumption
    // the one we need boundaried with two chains must be with fewer faces

    //size_t visited_faces_num = count(face_visited_flag.begin(),face_visited_flag.end(),true);
    vector<size_t> visited_faces,non_visited_faces;
    for(size_t t  = 0; t < face_visited_flag.size(); ++t) {
        if(face_visited_flag[t])
          visited_faces.push_back(t);
        else
          non_visited_faces.push_back(t);
      }
    matrixst patch_visited,  patch_non_visited;
    itr_matrix<size_t*> visited_faces_(1,visited_faces.size(),&visited_faces[0]);
    itr_matrix<size_t*> non_visited_faces_(1,non_visited_faces.size(),&non_visited_faces[0]);

    patch_visited = outside_face(colon(),visited_faces_);
    patch_non_visited = outside_face(colon(),non_visited_faces_);

    bool is_visited_patch_simgple =  is_simple_connected_patch(patch_visited);
    bool is_non_visited_patch_simgple =  is_simple_connected_patch(patch_non_visited);

    if(!is_visited_patch_simgple &&
       !is_non_visited_patch_simgple){
        cerr << "# [error] can not find a valid simple patch." << endl;
        return __LINE__;
      }
    if(is_visited_patch_simgple && !is_non_visited_patch_simgple){
        smaller_face_patch_with_boundary = visited_faces;
        return 0;
      }
    if(!is_visited_patch_simgple && is_non_visited_patch_simgple){
        smaller_face_patch_with_boundary = non_visited_faces;
        return 0;
      }
    if(is_visited_patch_simgple &&
       is_non_visited_patch_simgple){
        if(visited_faces.size() < non_visited_faces.size()){
            smaller_face_patch_with_boundary = visited_faces;
          }else
          smaller_face_patch_with_boundary = non_visited_faces;
        return 0;
      }
  }


  return 0;
}



/**
 * @brief
 *
 * @param linked_points
 * @param ea
 * @param outside_face_fix_axis
 * @param vector<deque<pair<size_t
 * @param edge_chain
 * @return int
 */
static int extract_chain_from_points(const vector<size_t> & linked_points,
                                     const jtf::mesh::edge2cell_adjacent & ea,
                                     const matrixst & outside_face_fix_axis,
                                     vector<deque<pair<size_t,size_t> > > & edge_chain)
{
  if(linked_points.size() < 2) {
      cerr << "# [error] there are only " << linked_points.size() << " points linked together." << endl;
      return __LINE__;
    }
  vector<pair<size_t,size_t> > possible_edges;
  for(size_t t = 0; t < linked_points.size(); ++t){
      for(size_t i = t+1; i < linked_points.size(); ++i){
          const size_t edge_idx = ea.get_edge_idx(linked_points[t],linked_points[i]);
          if( edge_idx != -1){
              const pair<size_t,size_t> & face_pair = ea.edge2cell_[edge_idx];
              if(!ea.is_boundary_edge(face_pair)){
                  if(outside_face_fix_axis[face_pair.first] ==
                     outside_face_fix_axis[face_pair.second])
                    continue;
                }
              if(linked_points[t] < linked_points[i]){
                  possible_edges.push_back(make_pair(linked_points[t],
                                                     linked_points[i]));
                }else{
                  possible_edges.push_back(make_pair(linked_points[i],
                                                     linked_points[t]));
                }
            }
        }
    }

#if 1 // check whether all points are visited
  {
    set<size_t> original_points(linked_points.begin(),linked_points.end());
    if(original_points.size() != linked_points.size())
      cerr << "# [error] input linked_points have duplicated points." << endl;
    set<size_t> linked_edge_points;
    for(size_t t = 0; t < possible_edges.size(); ++t){
        linked_edge_points.insert(possible_edges[t].first);
        linked_edge_points.insert(possible_edges[t].second);
      }
    if(linked_edge_points != original_points){
        cerr << "# [error] there are some points not visited." << endl;
        return __LINE__;
      }
  }
#endif

  //vector<deque<pair<size_t,size_t> > > edge_chain;
  edge_chain.clear();
  jtf::util::extract_chain_from_edges(possible_edges,edge_chain);
  return 0;
}

/**
 * @brief
 *
 * @param chain_type
 * @return size_t
 */
size_t get_singularity_around_axis(const deque<size_t> &chain_type)
{
  set<size_t> around_axis;
  for(size_t t = 0; t < chain_type.size(); ++t)
    around_axis.insert(axis_to_around(chain_type[t]));
  if(around_axis.size() != 1) return -1; // this singularity contain severay axis rotation
  return *around_axis.begin();
}

bool is_edge_type_accetable(const size_t & edge_type)
{
  if(is_regular_type(edge_type) || is_trivial_type(edge_type))
    return true;
  else
    return false;
}


bool is_singularity_configuration_accetable(
    const std::map<std::pair<size_t,size_t>,size_t> & face_type,
    const std::pair<size_t,size_t> & attempt_tet_pair,
    const size_t  & attempt_face_type)
{

}


/**
 * @brief         since surface type is recorded by: min |n - (f*R)colon(0)|
 *                so we should figure out where does the finial u in f*R come
 *                from f. It depends on the non-zero item in first line of R^T
 *
 * @param type
 * @return size_t
 */
size_t get_surface_fix_axis(const size_t type)
{
  const matrixd rotation_mat = trans(type_transition2(type));
  for(size_t t = 0; t < 3; ++t){
      if(fabs(rotation_mat(0,t)) > 1e-8){ // non-zero
          return t;
        }
    }
  return -1;
}

/**
 * @brief
 *
 */
struct point_fix_axis_per_edge{
  pair<size_t,size_t> edge; // point belong to this edge /**< TODO */
  vector<size_t> fix_axis; /**< TODO */
};

/**
 * @brief
 *
 */
struct point_fix_axis_per_face{
  size_t face_idx; /**< TODO */
  size_t fix_axis; /**< TODO */
};

/**
 * @brief
 *
 * @param current_point_idx
 * @param prev_point_idx
 * @param singularity_point_uvw_idx
 * @param singularity_point_uvw_ptr
 * @param temp_to_store
 * @param visited_points
 * @param has_meet_loop
 * @return int
 */
int constrate_axis_fix_ptr_deep(
    const size_t current_point_idx,
    const size_t prev_point_idx,
    const matrixst & singularity_point_uvw_idx,
    vector<set<size_t*> > &singularity_point_uvw_ptr,
    set<size_t*> & temp_to_store,
    set<size_t> &visited_points,
    bool & has_meet_loop)
{
  //  if(current_point_idx == 157)
  //    cerr << "pause";
  set<size_t*> &link_points = singularity_point_uvw_ptr[current_point_idx];
  if(link_points.empty()) return 0;
  if(link_points.size() == 2){
      if(link_points.find(const_cast<size_t*>(&singularity_point_uvw_idx[current_point_idx])) != link_points.end() &&
         link_points.find(const_cast<size_t*>(&singularity_point_uvw_idx[prev_point_idx])) != link_points.end()){ // end
          //temp_to_store.insert(link_points.begin(),link_points.end());
          link_points.insert(temp_to_store.begin(),temp_to_store.end());
          //link_points.clear();
          return 0;
        }else{
          cerr << "# [error] this point only link two points, but in wrong condition"
               << endl;
          cerr << "# [error] current point " << current_point_idx << endl;
          cerr << "# [error] prev point " << prev_point_idx << endl;

          for(set<size_t*>::const_iterator scit = link_points.begin();
              scit != link_points.end(); ++scit){
              cerr << *(*scit) << " ";
            }
          return __LINE__;
        }
    }

  const size_t visited_points_num = temp_to_store.size();
  temp_to_store.insert(link_points.begin(),link_points.end());
  if(temp_to_store.size() == visited_points_num)
    {
      link_points.clear();
      has_meet_loop = true;
      return 0; // means current point can not offer new link points
    }

  //bool has_meet_loop = false;
  for(set<size_t*>::const_iterator scit = link_points.begin();
      scit != link_points.end(); ++scit){
      if(*(*scit) == current_point_idx || *(*scit) == prev_point_idx) continue;
      const size_t prev_vis_num = visited_points.size();
      visited_points.insert(*(*scit));
      if(visited_points.size() == prev_vis_num) {
          has_meet_loop = true;
          continue; // has visited
        }
      constrate_axis_fix_ptr_deep(*(*scit),current_point_idx,singularity_point_uvw_idx,
                                  singularity_point_uvw_ptr,temp_to_store,visited_points,has_meet_loop);
      if(link_points.empty())
        return 0;
    }
  //temp_to_store.insert(link_points.begin(),link_points.end());
  link_points.clear();
  return 0;
}

/**
 * @brief
 *
 * @param point_idx
 * @param singularity_point_uvw_idx
 * @param singularity_point_uvw_ptr
 * @return int
 */
int constrate_axis_fix_ptr(const size_t point_idx,
                           const matrixst & singularity_point_uvw_idx,
                           vector<set<size_t*> > &singularity_point_uvw_ptr)
{
  assert(singularity_point_uvw_idx.size() == singularity_point_uvw_ptr.size());
  if(point_idx >= singularity_point_uvw_idx.size()) return __LINE__;
  set<size_t*> & link_points = singularity_point_uvw_ptr[point_idx];

  if(link_points.empty()) return 0;
  bool is_constrated = true;
  for(set<size_t*>::const_iterator scit = link_points.begin();
      scit != link_points.end(); ++scit){
      if(*(*scit) == point_idx) continue;
      const set<size_t*> & other_linked = singularity_point_uvw_ptr[*(*scit)];
      if(!other_linked.empty()) {
          is_constrated = false;
          break;
        }
    }

  if(is_constrated) return 0;

  bool has_meet_loop = false;
  set<size_t*> temp = link_points;
  set<size_t> visited_points;
  visited_points.insert(point_idx);
  for(set<size_t*>::const_iterator sit = link_points.begin();
      sit != link_points.end(); ++sit){
      if(*(*sit) == point_idx) continue;
      constrate_axis_fix_ptr_deep(*(*sit),point_idx,singularity_point_uvw_idx,
                                  singularity_point_uvw_ptr,temp,visited_points,has_meet_loop);
      //aux_current.insert(temp.begin(),temp.end());
    }
  //link_points.insert(aux_current.begin(),aux_current.end());
  if(!has_meet_loop)
    link_points.clear();
  else
    {
      link_points.insert(temp.begin(),temp.end());
    }

#if 0 // check
  for(set<size_t*>::const_iterator scit = link_points.begin();
      scit != link_points.end(); ++scit){
      const size_t & point_idx_ = *(*scit);
      const set<size_t*> & linked_points = singularity_point_uvw_ptr[point_idx_];
      if(point_idx_ == point_idx) continue;
      if(!linked_points.empty()) {
          cerr << "# [error] this linked_points should be empty: "
               << point_idx << " --> " << point_idx_ << endl;
        }
    }
#endif
  return 0;
}

/**
 * @brief
 *
 * @param uvw_of_point_on_surface_ptr
 * @param points_need_to_link
 * @param fix_axis
 * @return int
 */
static int link_surface_points(matrix<size_t*> &uvw_of_point_on_surface_ptr,
                               const vector<size_t> &points_need_to_link,
                               const size_t fix_axis)
{
  set<size_t> points_on_tracking;
  for(size_t t = 0; t < points_need_to_link.size(); ++t){
      points_on_tracking.insert(3*points_need_to_link[t] + fix_axis);
      size_t linked_coord = 3*points_need_to_link[t] + fix_axis;
      while(*uvw_of_point_on_surface_ptr[linked_coord] != linked_coord){
          const size_t points_on_tracking_prev_num = points_on_tracking.size();
          points_on_tracking.insert(*uvw_of_point_on_surface_ptr[linked_coord]);
          linked_coord = *uvw_of_point_on_surface_ptr[linked_coord];
        }
    }

  for(set<size_t>::const_iterator scit = points_on_tracking.begin();
      scit != points_on_tracking.end(); ++scit){
      *uvw_of_point_on_surface_ptr[*scit] = *points_on_tracking.begin();
    }
  return 0;
}

/**
 * @brief
 *
 * @param outside_face
 * @param outside_face_type
 * @param map<size_t
 * @param real_point_to_idx
 * @param idx_to_real_point
 * @param ea
 * @param uvw_of_point_on_surface_ptr
 * @return int
 */
static int constrate_axis_fix_ptr_surface(const matrixst &outside_face,
                                          const matrixst &outside_face_type,
                                          const map<size_t,size_t> &real_point_to_idx,
                                          const matrixst &idx_to_real_point,
                                          const jtf::mesh::edge2cell_adjacent & ea,
                                          matrix<size_t*> &uvw_of_point_on_surface_ptr)
{
  stack<size_t> face_stack;
  face_stack.push(0);
  vector<bool> outside_face_visited(outside_face.size(2),false);
  vector<size_t> point_uvw_on_face(3);
  size_t adjacent_face[3];

  while(!face_stack.empty()){
      const size_t face_idx = face_stack.top();
      face_stack.pop();
      if(outside_face_visited[face_idx]) continue;
      outside_face_visited[face_idx] = true;
      const size_t fix_axis = get_surface_fix_axis(outside_face_type[face_idx]);
      assert(fix_axis == 0 || fix_axis == 1 || fix_axis == 2);


      for(size_t i = 0; i < outside_face.size(1); ++i){
          map<size_t,size_t>::const_iterator mcit =
              real_point_to_idx.find(outside_face(i,face_idx));
          if(mcit == real_point_to_idx.end()){
              cerr << "# [error] can not find point " << outside_face(i,face_idx);
              return __LINE__;
            }

          point_uvw_on_face[i] = mcit->second;
          if(ea.is_boundary_edge(make_pair(outside_face(i,face_idx),
                                           outside_face((i+1)%outside_face.size(1),face_idx))))
            continue;
          const pair<size_t,size_t>  face_pair =
              ea.query(outside_face(i,face_idx),
                       outside_face((i+1)%outside_face.size(1),face_idx));
          if(face_pair.first == -1 || face_pair.second == -1){
              cerr << "# [error] strange this outside face is illegal." << endl;
              return __LINE__;
            }
          assert(face_pair.first == face_idx || face_pair.second == face_idx);
          face_stack.push(face_pair.first + face_pair.second - face_idx);
        }
      // unify the three points on one face with the same uvw axis;
      // all link to the small index point
      link_surface_points(uvw_of_point_on_surface_ptr,point_uvw_on_face,fix_axis);
    }

  if(find(outside_face_visited.begin(),outside_face_visited.end(),false)
     != outside_face_visited.end()){
      cerr << "# [error] this outside mesh is not manifold." << endl;
      return __LINE__;
    }

  // link all points and uniform them
  {
    for(size_t t = 0; t < uvw_of_point_on_surface_ptr.size(); ++t){
        const size_t &linked_point_idx = *(uvw_of_point_on_surface_ptr[t]);
        if(linked_point_idx == t) continue;
        assert(linked_point_idx < t);
        //assert(*(uvw_of_point_on_surface_ptr[linked_point_idx]) == linked_point_idx);
        const size_t & links_link = *(uvw_of_point_on_surface_ptr[linked_point_idx]);
        assert(*(uvw_of_point_on_surface_ptr[links_link]) == links_link);
        *(uvw_of_point_on_surface_ptr[t]) = links_link;
      }
  }

  return 0;
}

/**
 * @brief
 *
 * @param deque<pair<size_t
 * @param first_chain
 * @param first_chain_type
 * @param deque<pair<size_t
 * @param second_chain
 * @param second_chain_type
 * @return int
 */
static int combine_two_singualrity_with_type(
    deque<pair<size_t,size_t> > & first_chain,
    deque<size_t> & first_chain_type,
    const deque<pair<size_t,size_t> > & second_chain,
    const deque<size_t> & second_chain_type)
{
  if(first_chain.back().second == second_chain.front().first){
      first_chain.insert(first_chain.end(),second_chain.begin(),second_chain.end());
      first_chain_type.insert(first_chain_type.end(),second_chain_type.begin(),second_chain_type.end());
    }else if(first_chain.front().first == second_chain.front().first)
    {
      reverse_singularity_with_type(first_chain,first_chain_type);
      first_chain.insert(first_chain.end(),second_chain.begin(),second_chain.end());
      first_chain_type.insert(first_chain_type.end(),second_chain_type.begin(),second_chain_type.end());
    }else if(first_chain.front().first == second_chain.back().second){
      first_chain.insert(first_chain.begin(),second_chain.begin(),second_chain.end());
      first_chain_type.insert(first_chain_type.begin(),second_chain_type.begin(),second_chain_type.end());
    }else if(first_chain.back().second == second_chain.back().second){
      reverse_singularity_with_type(first_chain,first_chain_type);
      first_chain.insert(first_chain.begin(),second_chain.begin(),second_chain.end());
      first_chain_type.insert(first_chain_type.begin(),second_chain_type.begin(),second_chain_type.end());
    }

  if(!first_chain.empty()){
      for(size_t t = 0; t < first_chain.size() - 1; ++t){
          assert(first_chain[t].second == first_chain[t+1].first);
        }
    }
  return 0;
}

/**
 * @brief
 *
 * @param deque<pair<size_t
 * @param first_chain
 * @param first_chain_type
 * @param deque<pair<size_t
 * @param second_chain
 * @param second_chain_type
 * @return int
 */
static int combine_two_singualrity(
    deque<pair<size_t,size_t> > & first_chain,
    const deque<pair<size_t,size_t> > & second_chain)
{
  if(first_chain.back().second == second_chain.front().first){
      first_chain.insert(first_chain.end(),second_chain.begin(),second_chain.end());
      //first_chain_type.insert(first_chain_type.end(),second_chain_type.begin(),second_chain_type.end());
    }else if(first_chain.front().first == second_chain.front().first)
    {
      reverse_singularity(first_chain);
      first_chain.insert(first_chain.end(),second_chain.begin(),second_chain.end());
      //first_chain_type.insert(first_chain_type.end(),second_chain_type.begin(),second_chain_type.end());
    }else if(first_chain.front().first == second_chain.back().second){
      first_chain.insert(first_chain.begin(),second_chain.begin(),second_chain.end());
      //first_chain_type.insert(first_chain_type.begin(),second_chain_type.begin(),second_chain_type.end());
    }else if(first_chain.back().second == second_chain.back().second){
      reverse_singularity(first_chain);
      first_chain.insert(first_chain.begin(),second_chain.begin(),second_chain.end());
      //first_chain_type.insert(first_chain_type.begin(),second_chain_type.begin(),second_chain_type.end());
    }

  if(!first_chain.empty()){
      for(size_t t = 0; t < first_chain.size() - 1; ++t){
          assert(first_chain[t].second == first_chain[t+1].first);
        }
    }
  return 0;
}

/**
 * @brief
 *
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularities_chain
 * @param singularities_type
 * @return int
 */
int combine_singularity_chain(
    std::vector<std::deque<std::pair<size_t,size_t> > >  &singularities_chain,
    std::vector<std::deque<size_t> > &singularities_type)
{
  map<size_t,vector<size_t> > point_adjacent_chain;
  for(size_t t = 0; t < singularities_chain.size(); ++t){
      const deque<pair<size_t,size_t> > & one_chian = singularities_chain[t];
      point_adjacent_chain[one_chian.front().first].push_back(t);
      point_adjacent_chain[one_chian.back().second].push_back(t);
    }

  set<pair<size_t,size_t> > need_to_combine_chains;
  for(map<size_t,vector<size_t> >::const_iterator msvcit = point_adjacent_chain.begin();
      msvcit != point_adjacent_chain.end(); ++msvcit){
      if(msvcit->second.size() == 2) {
          if(msvcit->second.front() > msvcit->second.back()){
              need_to_combine_chains.insert(make_pair(msvcit->second.back(),
                                                      msvcit->second.front()));
            }else
            need_to_combine_chains.insert(make_pair(msvcit->second.front(),
                                                    msvcit->second.back()));
        }
    }
  vector<deque<pair<size_t,size_t> > > new_singularity_chain;
  vector<deque<size_t> > new_singularity_chain_type;
  vector<bool> remained_chain_flag(singularities_chain.size(),true);
  for(set<pair<size_t,size_t> >::const_iterator scit = need_to_combine_chains.begin();
      scit != need_to_combine_chains.end(); ++scit){
      combine_two_singualrity_with_type(singularities_chain[scit->first],
          singularities_type[scit->first],
          singularities_chain[scit->second],
          singularities_type[scit->second]);
      new_singularity_chain.push_back(singularities_chain[scit->first]);
      new_singularity_chain_type.push_back(singularities_type[scit->first]);
      remained_chain_flag[scit->first] = false;
      remained_chain_flag[scit->second] = false;
    }

  for(size_t t = 0; t < remained_chain_flag.size(); ++t){
      if(remained_chain_flag[t]){
          new_singularity_chain.push_back(singularities_chain[t]);
          new_singularity_chain_type.push_back(singularities_type[t]);
        }
    }
  if(!need_to_combine_chains.empty()){
      assert(new_singularity_chain.size() + need_to_combine_chains.size()
             == singularities_chain.size()+1);
    }
  singularities_chain = new_singularity_chain;
  singularities_type = new_singularity_chain_type;

  return 0;
}

int check_tet_rotation_type_compatible_with_jump_type(
    const std::vector<size_t> & tet_rot_type,
    const jtf::mesh::face2tet_adjacent & fa,
    const std::map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type)
{
  //const matrixd seed_matrix = eye<double>(3);
  typedef std::map<std::pair<size_t,size_t>,size_t>::const_iterator map_cit;
  matrixd rot_jump = eye<double>(3);

  for(size_t t = 0; t < fa.face2tet_.size(); ++t){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[t];
      if(fa.is_outside_face(tet_pair)) continue;
      const matrixd first = type_transition2(tet_rot_type[tet_pair.first]);
      const matrixd second = type_transition2(tet_rot_type[tet_pair.second]);
      map_cit it = inner_face_jump_type.find(tet_pair);
      if(it == inner_face_jump_type.end())
        rot_jump = eye<double>(3);
      else
        rot_jump = type_transition2(it->second);
      if(norm(first * rot_jump - second) > 1e-6){
          cerr << "# [error] jump type is not compatible, tet_pair <"
               << tet_pair.first << "," << tet_pair.second << ">: "
               << "type " << tet_rot_type[tet_pair.first] << ","
               << tet_rot_type[tet_pair.second] << endl;
          cerr << "# [error] jump type is " << type_transition1(rot_jump) << endl;
        }
    }
  return 0;
}

int check_cut_tet_inner_after_aligned(
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type)
{
  typedef boost::unordered_map<std::pair<size_t,size_t>,size_t>::const_iterator mpscit;
  for(size_t t = 0; t < fa_cut.face2tet_.size(); ++t){
      const pair<size_t,size_t> & tet_pair = fa_cut.face2tet_[t];
      if(fa_cut.is_outside_face(tet_pair)) continue;
      mpscit it = inner_face_jump_type.find(tet_pair);
      if(it == inner_face_jump_type.end()) continue;
      if(!is_trivial_type(it->second)){
          cerr << "# [error] tet pair <" << tet_pair.first << "," << tet_pair.second
               << "supposed to be identity, its type is " << it->second << endl;
        }
    }
  return 0;
}


/**
 * @brief
 *
 * @param face_group
 * @param outside_face_normal
 * @param face_normal
 * @return int
 */
static int calc_average_face_normal(const vector<size_t> & face_group,
                                    const matrixd & outside_face_normal,
                                    matrixd & face_normal)
{
  face_normal = zeros<double>(3,1);
  for(size_t t = 0; t < face_group.size(); ++t){
      face_normal += outside_face_normal(colon(),face_group[t]);
    }
  if(norm(face_normal) < 1e-8){
      cerr << " # [error] face normal degenerated" << endl;
      return __LINE__;
    }
  face_normal /= norm(face_normal);
  return 0;
}


/**
 * @brief
 *
 * @param deque<pair<size_t
 * @param ori_chain
 * @param node
 * @param vector<deque<pair<size_t
 * @param surface_chain
 * @param selected_chain_idx
 * @param ori_chain_vec
 * @param other_chain_vec
 * @return int
 */
static int find_a_loop_for_given_chain(const deque<pair<size_t,size_t> > & ori_chain,
                                       const matrixd & node,
                                       const vector<deque<pair<size_t,size_t> > > & surface_chain,
                                       const vector<size_t> & selected_chain_idx,
                                       vector<size_t> & ori_chain_vec,
                                       vector<size_t> & other_chain_vec)
{
  map<pair<size_t,size_t>,double> surface_edge_weight;
  for(size_t t = 0; t < surface_chain.size(); ++t){
      const deque<pair<size_t,size_t> >  & one_chain = surface_chain[t];
      for(size_t i = 0; i < one_chain.size(); ++i){
          surface_edge_weight[one_chain[i]] = norm(node(colon(),one_chain[i].first)
                                                   -node(colon(),one_chain[i].second));
        }
    }

  map<pair<size_t,size_t>,vector<vector<size_t> > > chain_with_loops;
  for(size_t t = 0; t < selected_chain_idx.size(); ++t){
      //      //chain_with_loops.
      map<pair<size_t,size_t>,double> surface_edge_weight_modified
          = surface_edge_weight;
      unique_ptr<vertex_connection<UNDIRECT> > vcu;
      vector<size_t> other_chain_path;
      vector<size_t> shortes_chain_vec;

      for(size_t i = 0; i < ori_chain.size(); ++i){
          shortes_chain_vec.push_back(ori_chain[i].first);
        }
      shortes_chain_vec.push_back(ori_chain.back().second);

      for(size_t i = 0; i < ori_chain.size(); ++i){
          map<pair<size_t,size_t>,double>::iterator mpdit =
              surface_edge_weight_modified.find(ori_chain[i]);
          surface_edge_weight_modified.erase(mpdit);
        }
      vcu.reset(vertex_connection<UNDIRECT>::create(surface_edge_weight_modified));
      vcu->get_shortest_path(ori_chain.front().first,
                             ori_chain.back().second,
                             other_chain_path);

      if(other_chain_path.empty()){
          cerr << "# [error] it's hard to imagine this case." << endl;
          return __LINE__;
        }

      chain_with_loops[make_pair(ori_chain.front().first,
                                 ori_chain.back().second)].push_back(shortes_chain_vec);

      chain_with_loops[make_pair(ori_chain.front().first,
                                 ori_chain.back().second)].push_back(other_chain_path);
    }

  vector<std::tuple<size_t,size_t,size_t> > priority_loop_sequence;
  typedef map<pair<size_t,size_t>,vector<vector<size_t> > >::const_iterator mpvcit;
  for(mpvcit it = chain_with_loops.begin(); it != chain_with_loops.end(); ++it){
      const pair<size_t,size_t> & chain_with_ends = it->first;
      const vector<vector<size_t> > & loops = it->second;
      size_t total_loop_length = 0;
      for(size_t t = 0 ; t < loops.size(); ++t) total_loop_length += loops[t].size() - 1;

      std::tuple<size_t,size_t,size_t> tp(total_loop_length,
                                          chain_with_ends.first,
                                          chain_with_ends.second);
      priority_loop_sequence.push_back(tp);
    }
  sort(priority_loop_sequence.begin(),priority_loop_sequence.end());

  const std::tuple<size_t,size_t,size_t> & need_to_handle_loop =
      priority_loop_sequence.front();
  mpvcit need_to_handle_loop_ptr =
      chain_with_loops.find(make_pair(get<1>(need_to_handle_loop),
                                      get<2>(need_to_handle_loop)));
  assert(need_to_handle_loop_ptr->second.size() == 2);

  ori_chain_vec = need_to_handle_loop_ptr->second.back();
  other_chain_vec = need_to_handle_loop_ptr->second.front();

  return 0;
}

/**
 * @brief
 *
 * @param node
 * @param degenerated_points
 * @param vector<deque<pair<size_t
 * @param surface_chain
 * @param selected_chain_idx
 * @param ea
 * @param outside_face_patch_idx
 * @param outside_face_normal
 * @param map<size_t
 * @param patch_idx_outside_face
 * @param orfap
 * @param outside_face
 * @param outside_face_type
 * @return int
 */
static int change_surface_type_to_remove_loop_degenerated_to_point(
    const matrixd & node,
    const vector<size_t> & degenerated_points,
    const vector<deque<pair<size_t,size_t> > > &surface_chain,
    const vector<size_t> & selected_chain_idx,
    const jtf::mesh::edge2cell_adjacent &ea,
    const matrixst &outside_face_patch_idx,
    const matrixd & outside_face_normal,
    const map<size_t,vector<size_t> > & patch_idx_outside_face,
    const jtf::mesh::one_ring_face_at_point & orfap,
    const matrixst & outside_face,
    matrixst &outside_face_type)
{
  vector<size_t> other_chain;
  vector<size_t> shortest_chain;
  if(selected_chain_idx.size() == 1){ // need to find another chain which make a loop
      const deque<pair<size_t,size_t> > & shortest_chain_edges =
          surface_chain[selected_chain_idx.front()];
      find_a_loop_for_given_chain(shortest_chain_edges,node,surface_chain,selected_chain_idx,
                                  shortest_chain,other_chain);

#if 1
      {
        vector<size_t> simple_patch;
        if(find_face_patch_with_two_closed_boundary(node,shortest_chain,
                                                    other_chain,ea,outside_face,simple_patch))
          {
            cerr  << "# [error] can not find a face patch." << endl;
            return __LINE__;
          }
        if(!is_patch_with_same_axis_fix(simple_patch,outside_face_type)){
            cerr << "# [error] this simple patch is invalid for containing several fix axes." << endl;
            return __LINE__;
          }
      }
#endif
      //    map<pair<size_t,size_t>,double> surface_edge_weight;
      //    for(size_t t = 0; t < surface_chain.size(); ++t){
      //      const deque<pair<size_t,size_t> >  & one_chain = surface_chain[t];
      //      for(size_t i = 0; i < one_chain.size(); ++i){
      //        surface_edge_weight[one_chain[i]] = norm(node(colon(),one_chain[i].first)
      //                                                 -node(colon(),one_chain[i].second));
      //      }
      //    }

      //    map<pair<size_t,size_t>,vector<vector<size_t> > > chain_with_loops;

      //    unique_ptr<vertex_connection<UNDIRECT> > vcu;
      //    vector<size_t> other_chain_path;
      //    vector<size_t> shortes_chain_vec;
      //    for(size_t t = 0; t < selected_chain_idx.size(); ++t){
      //      //chain_with_loops.
      //      map<pair<size_t,size_t>,double> surface_edge_weight_modified
      //          = surface_edge_weight;
      //      const deque<pair<size_t,size_t> > & shortest_chain =
      //          surface_chain[selected_chain_idx[t]];

      //      shortes_chain_vec.clear();

      //      for(size_t i = 0; i < shortest_chain.size(); ++i){
      //        shortes_chain_vec.push_back(shortest_chain[i].first);
      //      }
      //      shortes_chain_vec.push_back(shortest_chain.back().second);

      //      for(size_t i = 0; i < shortest_chain.size(); ++i){
      //        map<pair<size_t,size_t>,double>::iterator mpdit =
      //            surface_edge_weight_modified.find(shortest_chain[i]);
      //        surface_edge_weight_modified.erase(mpdit);
      //      }
      //      vcu.reset(vertex_connection<UNDIRECT>::create(surface_edge_weight_modified));
      //      vcu->get_shortest_path(shortest_chain.front().first,
      //                             shortest_chain.back().second,
      //                             other_chain_path);

      //      if(other_chain_path.empty()){
      //        cerr << "# [error] it's hard to imagine this case." << endl;
      //        return __LINE__;
      //      }

      //      chain_with_loops[make_pair(shortest_chain.front().first,
      //                                 shortest_chain.back().second)].push_back(shortes_chain_vec);

      //      chain_with_loops[make_pair(shortest_chain.front().first,
      //                                 shortest_chain.back().second)].push_back(other_chain_path);

      //    }

      //    vector<std::tuple<size_t,size_t,size_t> > priority_loop_sequence;
      //    typedef map<pair<size_t,size_t>,vector<vector<size_t> > >::const_iterator mpvcit;
      //    for(mpvcit it = chain_with_loops.begin(); it != chain_with_loops.end(); ++it){
      //      const pair<size_t,size_t> & chain_with_ends = it->first;
      //      const vector<vector<size_t> > & loops = it->second;
      //      size_t total_loop_length = 0;
      //      for(size_t t = 0 ; t < loops.size(); ++t) total_loop_length += loops[t].size() - 1;

      //      tuples::tuple<size_t,size_t,size_t> tp(total_loop_length,
      //                                             chain_with_ends.first,
      //                                             chain_with_ends.second);
      //      priority_loop_sequence.push_back(tp);
      //    }
      //    sort(priority_loop_sequence.begin(),priority_loop_sequence.end());

      //    const tuples::tuple<size_t,size_t,size_t> & need_to_handle_loop =
      //        priority_loop_sequence.front();
      //    mpvcit need_to_handle_loop_ptr =
      //        chain_with_loops.find(make_pair(need_to_handle_loop.get<1>(),
      //                                        need_to_handle_loop.get<2>()));
      //    assert(need_to_handle_loop_ptr->second.size() == 2);

      //    other_chain = need_to_handle_loop_ptr->second.back();
      //    shortest_chain = need_to_handle_loop_ptr->second.front();

    }else{ // if there are many chains connect two points, we just need to select two with simple connected patch

      assert(selected_chain_idx.size() > 1);
      shortest_chain.clear();
      const deque<pair<size_t,size_t> > & first_chain = surface_chain[selected_chain_idx.front()];
      for(size_t i = 0; i < first_chain.size(); ++i){
          shortest_chain.push_back(first_chain[i].first);
        }
      shortest_chain.push_back(first_chain.back().second);

      vector<size_t> simple_patch;
      for(size_t t = 1; t < selected_chain_idx.size(); ++t){
          const deque<pair<size_t,size_t> > & one_chain =
              surface_chain[selected_chain_idx[t]];
          other_chain.clear();
          for(size_t i = 0; i < one_chain.size(); ++i){
              other_chain.push_back(one_chain[i].first);
            }
          other_chain.push_back(one_chain.back().second);

          if(other_chain.front() != shortest_chain.front() &&
             other_chain.front() != shortest_chain.back()){
              continue; // these two chains can not make a loop
            }
          if(other_chain.back() != shortest_chain.front() &&
             other_chain.back() != shortest_chain.back()){
              continue; // there two chains can not make a loop
            }

          if(find_face_patch_with_two_closed_boundary(node,shortest_chain,other_chain,ea,outside_face,
                                                      simple_patch)){
              simple_patch.clear();
              continue;
            }
          else{
              if(is_patch_with_same_axis_fix(simple_patch,outside_face_type))
                break;
              else
                simple_patch.clear();
            }
        }
      if(simple_patch.empty()) {

          // there chains can not make a valid loop, needs to find a new one
          find_a_loop_for_given_chain(first_chain,node,surface_chain,selected_chain_idx,
                                      shortest_chain,other_chain);
#if 1 // check whether this patch is valid or not
          {
            vector<size_t> simple_patch;
            if(find_face_patch_with_two_closed_boundary(node,shortest_chain,
                                                        other_chain,ea,outside_face,simple_patch))
              {
                cerr  << "# [error] can not find a face patch." << endl;
                return __LINE__;
              }
            if(!is_patch_with_same_axis_fix(simple_patch,outside_face_type)){
                cerr << "# [error] this simple patch is invalid for containing several fix axes." << endl;
                return __LINE__;
              }
          }
#endif
        }
    }

  // first chain is the shortest one , which both ends must be three patches' intersection
  // the type of patch boundary with two chain can be determined by counting the adjacent
  // face type along two chains


  size_t current_loop_patch_type = -1;


  vector<size_t> simple_face_patch_with_boundary;
  { // find surface face set with is boundaried by these two chains
    if(find_face_patch_with_two_closed_boundary(node,shortest_chain,other_chain,ea,outside_face,
                                                simple_face_patch_with_boundary))
      return  __LINE__;

    if(!is_patch_with_same_axis_fix(simple_face_patch_with_boundary,outside_face_type))
      {
        cerr << "# [error] in such a simple patch, there are more than one fix axes." << endl;
        cerr << " need to find two new boundary chains." << endl;
        return __LINE__;
      }
    //    set<size_t> fix_axis_type_in_simple_patch;
    //    for(size_t t = 0; t < simple_face_patch_with_boundary.size(); ++t){
    //      fix_axis_type_in_simple_patch.insert(
    //            get_surface_fix_axis(outside_face_type(simple_face_patch_with_boundary[t])));
    //    }
    //    if(fix_axis_type_in_simple_patch.size() != 1){
    //      cerr << "# [error] in such a simple patch, there are more than one fix axes." << endl;
    //      cerr << " need to find two new boundary chains." << endl;
    //      return __LINE__;
    //    }

    current_loop_patch_type =
        outside_face_type[simple_face_patch_with_boundary.front()];
  }
#if 1
  {
    static int count_num = 0;
    matrixst simple_face_patch_with_boundary_mat(3,simple_face_patch_with_boundary.size());
    for(size_t t = 0; t < simple_face_patch_with_boundary.size(); ++t){
        simple_face_patch_with_boundary_mat = outside_face(colon(),
                                                           simple_face_patch_with_boundary[t]);
      }
    ostringstream os;
    os << count_num++;
    ofstream ofs((os.str()+"_simple_connected_patch.vtk").c_str());
    tri2vtk(ofs,&node[0],node.size(2), &simple_face_patch_with_boundary_mat[0],
        simple_face_patch_with_boundary_mat.size(2));
  }
#endif

  size_t isolated_point = 0;

  {// this part is to find the isolated points ,
    // isolated points means walk around faces of such point, you will meeet a
    // same patch several times; for example the patch sequence may be u,v,u,v
    typedef jtf::mesh::one_ring_face_at_point::p2f_type::const_iterator p2fcit;

    for(size_t t = 0; t < other_chain.size(); ++t){
        p2fcit it = orfap.p2f_.find(other_chain[t]);
        if(it == orfap.p2f_.end()) {
            cerr << "# [error] strange this point is not found" << endl;
            return __LINE__;
          }
        const vector<size_t> & around_faces = it->second;
        assert(find(around_faces.begin(),around_faces.end(),-1) == around_faces.end());
        assert(around_faces.front() == around_faces.back());
        vector<size_t> meet_patches;
        for(size_t i = 0; i < around_faces.size()-1; ++i){
            if(meet_patches.empty()) {
                meet_patches.push_back(get_surface_fix_axis(outside_face_type[around_faces[i]]));
                continue;
              }
            if(meet_patches.back() ==
               get_surface_fix_axis(outside_face_type[around_faces[i]]))
              continue;
            else
              meet_patches.push_back(get_surface_fix_axis(outside_face_type[around_faces[i]]));
          }

        if(meet_patches.front() == meet_patches.back())
          meet_patches.pop_back();

        // determin whether such a patch loop is valid
        assert(meet_patches.size() != 1); // this point should not exist on surface singularity edges
        // 2: this point is a regular surface singularity point
        // 3: this point is a regular three surface patches intesection point
        if(meet_patches.size() == 2 ||
           meet_patches.size() == 3) continue;

        isolated_point++;
        // now meet isolated points, which should be reduced into two or three patches

        vector<size_t> need_to_modify_face;
        for(size_t i = 0; i < around_faces.size(); ++i){
            if(get_surface_fix_axis(outside_face_type[around_faces[i]])
               != get_surface_fix_axis(current_loop_patch_type)){
                if(need_to_modify_face.empty()) continue;
                need_to_modify_face.push_back(around_faces[i]);
                continue;
              }else{
                if(need_to_modify_face.empty()) {
                    need_to_modify_face.push_back(around_faces[i]);
                    continue;
                  }
                if(get_surface_fix_axis(outside_face_type[need_to_modify_face.back()])
                   == get_surface_fix_axis(outside_face_type[around_faces[i]])) {
                    continue;
                  }else{
                    break; // need_to_modify_face has already sotred the faces which
                    //should be modified into current_loop_patch_type
                  }
              }
          }// end for

        for(size_t t = 0; t < need_to_modify_face.size(); ++t){
            outside_face_type[need_to_modify_face[t]] = current_loop_patch_type;
          }
      }


    if(isolated_point == 0){ // if walk along other chain can not find isolated points
        // need to check shotest_chian
        for(size_t t = 0; t < shortest_chain.size(); ++t){
            p2fcit it = orfap.p2f_.find(shortest_chain[t]);
            if(it == orfap.p2f_.end()) {
                cerr << "# [error] strange this point is not found" << endl;
                return __LINE__;
              }
            const vector<size_t> & around_faces = it->second;
            assert(find(around_faces.begin(),around_faces.end(),-1) == around_faces.end());
            assert(around_faces.front() == around_faces.back());
            vector<size_t> meet_patches;
            for(size_t i = 0; i < around_faces.size()-1; ++i){
                if(meet_patches.empty()) {
                    meet_patches.push_back(get_surface_fix_axis(outside_face_type[around_faces[i]]));
                    continue;
                  }
                if(meet_patches.back() ==
                   get_surface_fix_axis(outside_face_type[around_faces[i]]))
                  continue;
                else
                  meet_patches.push_back(get_surface_fix_axis(outside_face_type[around_faces[i]]));
              }

            if(meet_patches.front() == meet_patches.back())
              meet_patches.pop_back();

            // determin whether such a patch loop is valid
            assert(meet_patches.size() != 1); // this point should not exist on surface singularity edges
            // 2: this point is a regular surface singularity point
            // 3: this point is a regular three surface patches intesection point
            if(meet_patches.size() == 2 ||
               meet_patches.size() == 3) continue;

            isolated_point++;
            // now meet isolated points, which should be reduced into two or three patches

            vector<size_t> need_to_modify_face;
            for(size_t i = 0; i < around_faces.size(); ++i){
                if(get_surface_fix_axis(outside_face_type[around_faces[i]])
                   != get_surface_fix_axis(current_loop_patch_type)){
                    if(need_to_modify_face.empty()) continue;
                    need_to_modify_face.push_back(around_faces[i]);
                    continue;
                  }else{
                    if(need_to_modify_face.empty()) {
                        need_to_modify_face.push_back(around_faces[i]);
                        continue;
                      }
                    if(get_surface_fix_axis(outside_face_type[need_to_modify_face.back()])
                       == get_surface_fix_axis(outside_face_type[around_faces[i]])) {
                        continue;
                      }else{
                        break; // need_to_modify_face has already sotred the faces which
                        //should be modified into current_loop_patch_type
                      }
                  }
              }// end for

            for(size_t t = 0; t < need_to_modify_face.size(); ++t){
                outside_face_type[need_to_modify_face[t]] = current_loop_patch_type;
              }
          }
      }
  }

  if(isolated_point != 0) return 0; // this degeneration casued by isolated points have been removed
  // then need to address those degeneration cased by two or more patch intesection point

  { // this case is a simple patch around with severay fix axis, which means
    // this simple patch should be modified into adjacent surface type
    set<pair<size_t,size_t> > face_pair_set;

    const size_t fix_axis_of_current_patch = get_surface_fix_axis(current_loop_patch_type);
    map<size_t,matrixd > face_type_to_normal_map;
    map<size_t,size_t> face_type_with_face_idx;
    typedef map<size_t,matrixd >::iterator msmit;

    for(size_t t = 0; t < other_chain.size()-1; ++t){
        const size_t edge_idx = ea.get_edge_idx(other_chain[t],other_chain[t+1]);
        if(edge_idx == -1){
            cerr << "# [error] strange can not find such an edge." << endl;
            return __LINE__;
          }

        const pair<size_t,size_t> & face_pair = ea.edge2cell_[edge_idx];
        face_pair_set.insert(face_pair);
      }

    for(size_t t = 0; t < shortest_chain.size()-1; ++t){
        const size_t edge_idx = ea.get_edge_idx(shortest_chain[t],shortest_chain[t+1]);
        if(edge_idx == -1){
            cerr << "# [error] strange can not find such an edge." << endl;
            return __LINE__;
          }

        const pair<size_t,size_t> & face_pair = ea.edge2cell_[edge_idx];
        face_pair_set.insert(face_pair);
      }

    for(set<pair<size_t,size_t> >::const_iterator scit = face_pair_set.begin();
        scit != face_pair_set.end(); ++scit){
        msmit it_first = face_type_to_normal_map.find(
              get_surface_fix_axis(outside_face_type[scit->first]));
        if(it_first == face_type_to_normal_map.end())
          face_type_to_normal_map.insert(make_pair(get_surface_fix_axis(outside_face_type[scit->first]),
                                         outside_face_normal(colon(),scit->first)));
        else
          it_first->second += outside_face_normal(colon(),scit->first);

        face_type_with_face_idx[get_surface_fix_axis(outside_face_type[scit->first])] =
            scit->first;

        msmit it_second = face_type_to_normal_map.find(
              get_surface_fix_axis(outside_face_type[scit->second]));
        if(it_second == face_type_to_normal_map.end())
          face_type_to_normal_map.insert(make_pair(get_surface_fix_axis(outside_face_type[scit->second]),
                                         outside_face_normal(colon(),scit->second)));
        else
          it_second->second += outside_face_normal(colon(),scit->second);

        face_type_with_face_idx[get_surface_fix_axis(outside_face_type[scit->first])] =
            scit->second;
      }

    matrixd normal_of_loop_patch = zeros<double>(3,1);
    calc_average_face_normal(simple_face_patch_with_boundary,outside_face_normal,
                             normal_of_loop_patch);

    vector<pair<double,size_t> > sort_face_normal_diff;
    for(msmit it = face_type_to_normal_map.begin();
        it != face_type_to_normal_map.end(); ++it){
        it->second /= norm(it->second);
        sort_face_normal_diff.push_back(make_pair(dot(it->second,normal_of_loop_patch),
                                                  face_type_with_face_idx[it->first]));
      }
    sort(sort_face_normal_diff.begin(),sort_face_normal_diff.end());
    size_t modified_surface_type = -1;
    for(vector<pair<double,size_t> >::reverse_iterator it = sort_face_normal_diff.rbegin();
        it != sort_face_normal_diff.rend(); ++it){
        modified_surface_type = outside_face_type[it->second];
        if(modified_surface_type != current_loop_patch_type) break;
      }
    if(modified_surface_type == -1){
        cerr << "# [error] all face type around this simple patch is the same." << endl;
        return __LINE__;
      }
    assert(modified_surface_type != current_loop_patch_type);
    for(size_t t = 0; t < simple_face_patch_with_boundary.size(); ++t){
        outside_face_type[simple_face_patch_with_boundary[t]] = modified_surface_type;
      }
  }

  return 0;
}

/**
 * @brief
 *
 * @param node
 * @param map<vector<size_t>
 * @param degenerated_points
 * @param vector<deque<pair<size_t
 * @param surface_chains
 * @param outside_face_normal
 * @param outside_face_fix_axis
 * @param outside_face_with_patch_idx
 * @param map<size_t
 * @param patch_idx_outside_face
 * @param ea
 * @param orfap
 * @param outside_face
 * @param outside_face_type
 * @return int
 */
static int remove_surface_degeneration_to_point(
    const matrixd & node,
    const map<vector<size_t>,vector<size_t> > & degenerated_points,
    //const map<vector<size_t>,vector<size_t> > & all_points_linked_by_uvw,
    const vector<deque<pair<size_t,size_t> > > & surface_chains,
    const matrixd & outside_face_normal,
    const matrixst & outside_face_fix_axis,
    const matrixst &outside_face_with_patch_idx,
    const map<size_t,vector<size_t> > &patch_idx_outside_face,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::one_ring_face_at_point & orfap,
    const matrixst & outside_face,
    matrixst & outside_face_type)
{
  typedef map<vector<size_t>,vector<size_t> >::const_iterator mvvcit;
  vector<size_t> selected_edge_chain_idx;

  for(mvvcit it = degenerated_points.begin(); it != degenerated_points.end(); ++it){
      const vector<size_t> & degenerated_points_array = it->second;

      selected_edge_chain_idx.clear();
      for(size_t t = 0; t < surface_chains.size(); ++t){
          const deque<pair<size_t,size_t> > & one_chain = surface_chains[t];

          if((find(degenerated_points_array.begin(),
                   degenerated_points_array.end(),
                   one_chain.front().first) != degenerated_points_array.end()) &&
             (find(degenerated_points_array.begin(),
                   degenerated_points_array.end(),
                   one_chain.back().second) != degenerated_points_array.end()))
            {
              selected_edge_chain_idx.push_back(t);
            }
        }
      if(selected_edge_chain_idx.empty()) continue; //TODO this case needs extra address
      int rtn = change_surface_type_to_remove_loop_degenerated_to_point(
            node,
            degenerated_points_array, surface_chains,
            selected_edge_chain_idx,ea,outside_face_with_patch_idx,
            outside_face_normal,patch_idx_outside_face,orfap,outside_face,
            outside_face_type);
      if(rtn == 0) break; // successfully modify one, need to rextract surface chain

    }
  return 0;
}

/**
 * @brief
 *
 * @param deque<pair<size_t
 * @param one_chain
 * @param map<size_t
 * @param real_point_to_idx
 * @param uvw_integer_fix
 * @return bool
 */
static bool is_surface_singularity_degenerated_to_edge(
    const deque<pair<size_t,size_t> > &one_chain,
    const map<size_t,size_t> &real_point_to_idx,
    const matrixst & uvw_integer_fix)
{
  if(one_chain.front().first == one_chain.back().second){
      return true;
    }
  set<vector<size_t> > uvw_coord;
  vector<size_t> temp_coord(3);
  for(size_t t = 0; t < one_chain.size(); ++t){
      const auto it = real_point_to_idx.find(one_chain[t].first);
      if(it == real_point_to_idx.end()){
          throw std::logic_error("# [error] can not find point in real_point_to_idx");
        }
      copy(uvw_integer_fix(colon(),it->second).begin(),
           uvw_integer_fix(colon(),it->second).end(),
           temp_coord.begin());
      uvw_coord.insert(temp_coord);
    }
  const auto it = real_point_to_idx.find(one_chain.back().second);
  if(it == real_point_to_idx.end()){
      throw std::logic_error("# [error] can not find point in real_point_to_idx");
    }
  copy(uvw_integer_fix(colon(),it->second).begin(),
       uvw_integer_fix(colon(),it->second).end(),
       temp_coord.begin());
  uvw_coord.insert(temp_coord);

  if(uvw_coord.size() == 1){
      if(count((*uvw_coord.begin()).begin(),
               (*uvw_coord.begin()).end(),
               -1) == 1){
          return true;
        }
    }
  return false;
}

/**
 * @brief
 *
 * @param deque<pair<size_t
 * @param chain
 * @param vec
 * @return int
 */
static int convert_chain_to_vec(const deque<pair<size_t,size_t> > & chain,
                                vector<size_t> & vec)
{
  vec.resize(chain.size() +1);
  for(size_t t = 0; t < chain.size(); ++t){
      vec[t] = chain[t].first;
    }
  vec.back() = chain.back().second;
  return 0;
}

/**
 * @brief
 *
 * @param node
 * @param vector<deque<pair<size_t
 * @param surface_singualrity_chain
 * @param degenerated_chain_idx
 * @param outside_face_with_patch_idx
 * @param map<size_t
 * @param patch_idx_outside_face
 * @param ea
 * @param uvw_integer_fix
 * @param map<size_t
 * @param real_point_to_idx
 * @param outside_face
 * @param outside_face_type
 * @return int
 */
static int remove_surface_degeneration_to_edge(
    const matrixd & node,
    const vector<deque<pair<size_t,size_t> > > &surface_singualrity_chain,
    const vector<size_t> & degenerated_chain_idx,
    const matrixst & outside_face_with_patch_idx,
    const map<size_t,vector<size_t> > &patch_idx_outside_face,
    const jtf::mesh::edge2cell_adjacent & ea,
    const matrixst & uvw_integer_fix,
    const map<size_t,size_t> & real_point_to_idx,
    const matrixst & outside_face,
    matrixst &outside_face_type)
{

  vector<size_t> chain_loop_idx; // easy case
  // complicated cases: severay chains may introduce an degenerated loops
  map<vector<size_t>,vector<size_t> > chains_degen_loop_map;
  vector<size_t> temp_coord(3);
  for(size_t t = 0; t < degenerated_chain_idx.size(); ++t){
      const deque<pair<size_t,size_t> > & one_chain =
          surface_singualrity_chain[degenerated_chain_idx[t]];
      if(one_chain.front().first == one_chain.back().second)
        chain_loop_idx.push_back(degenerated_chain_idx[t]);
      copy(uvw_integer_fix(colon(),one_chain.front().first).begin(),
           uvw_integer_fix(colon(),one_chain.front().first).end(),
           temp_coord.begin());
      chains_degen_loop_map[temp_coord].push_back(degenerated_chain_idx[t]);
    }

  typedef map<size_t,vector<size_t> >::const_iterator msvcit;

  // remove easy cases where one chain can make a loop, just modify the simple patch to
  // remove that loop
  for(size_t t = 0; t < chain_loop_idx.size(); ++t){
      const deque<pair<size_t,size_t> > & one_chain =
          surface_singualrity_chain[chain_loop_idx[t]];
      if(one_chain.front().first  != one_chain.back().second){
          cerr << "# [error] this chain is not a pure loop degeneration case." << endl;
          return __LINE__;
        }

      const pair<size_t,size_t> face_pair = ea.query(one_chain.front().first,
                                                     one_chain.front().second);
      assert(!ea.is_boundary_edge(face_pair));
      assert(outside_face_with_patch_idx[face_pair.first] !=
          outside_face_with_patch_idx[face_pair.second]);

      msvcit it_first = patch_idx_outside_face.find(outside_face_with_patch_idx[face_pair.first]);
      msvcit it_second = patch_idx_outside_face.find(outside_face_with_patch_idx[face_pair.second]);

      msvcit it_need_to_be_modified;
      size_t changed_face_type;
      if(it_first->second.size() < it_second->second.size()){
          it_need_to_be_modified = it_first;
          changed_face_type = outside_face_type[it_second->second.front()];
        }else{
          it_need_to_be_modified = it_second;
          changed_face_type = outside_face_type[it_first->second.front()];
        }
      for(size_t i = 0; i < it_need_to_be_modified->second.size(); ++i)
        {
          outside_face_type[it_need_to_be_modified->second[i]] = changed_face_type;
        }
    }

  vector<size_t> chain0,chain1;
  vector<size_t> simple_patch;
  for(map<vector<size_t>, vector<size_t> >::const_iterator mvvcit  =
      chains_degen_loop_map.begin(); mvvcit != chains_degen_loop_map.end();
      ++mvvcit){
      const vector<size_t> & chains_idx = mvvcit->second;
      if(chains_idx.size() < 2){
          cerr << "# [error] strange edge degenerated cases, but do not have two chain with the same uvw fix" << endl;
          return __LINE__;
        }

      const deque<pair<size_t,size_t> > & chain0_ori =
          surface_singualrity_chain[chains_idx.front()];
      convert_chain_to_vec(surface_singualrity_chain[chains_idx.front()], chain0);

      for(size_t t = 1; t < chains_idx.size(); ++t){
          const deque<pair<size_t,size_t> > & chain1_ori =
              surface_singualrity_chain[chains_idx[t]];
          if(chain0_ori.front().first != chain1_ori.back().second &&
             chain0_ori.front().first != chain1_ori.front().first) continue;
          if(chain0_ori.back().second != chain1_ori.back().second &&
             chain0_ori.back().second != chain1_ori.front().first) continue;

          convert_chain_to_vec(surface_singualrity_chain[chains_idx[t]], chain1);
          if(find_face_patch_with_two_closed_boundary(node,chain0,chain1,
                                                      ea,outside_face,simple_patch))
            return __LINE__;
          itr_matrix<size_t*> simple_patch_mat(1,simple_patch.size(),&simple_patch[0]);
          matrixst simple_patch_facse = outside_face(colon(),simple_patch_mat);
          if(!is_simple_connected_patch(simple_patch_facse))
            {
              cerr << "# [error] strange the face patch is not simple connected." << endl;
              return __LINE__;
            }
          if(!is_patch_with_same_axis_fix(simple_patch,outside_face_type)) continue;
          const pair<size_t,size_t> face_pair = ea.query(chain0[0],chain0[1]);
          assert(face_pair.first != -1 && face_pair.second != -1);
          const size_t current_loop_type = outside_face_type[simple_patch.front()];
          size_t patch_type_changed_to = -1;
          if(get_surface_fix_axis(outside_face_type[face_pair.first]) !=
             get_surface_fix_axis(current_loop_type)){
              patch_type_changed_to = outside_face_type[face_pair.first];
            }else if(get_surface_fix_axis(outside_face_type[face_pair.second]) !=
                     get_surface_fix_axis(current_loop_type)){
              patch_type_changed_to = outside_face_type[face_pair.second];
            }else{
              cerr << "# [error] both face pair of chain0.front edge do not satisfy." << endl;
              return __LINE__;
            }
          assert(patch_type_changed_to != -1);
          for(size_t i = 0; i < simple_patch.size(); ++i){
              outside_face_type[simple_patch[i]] = patch_type_changed_to;
            }
          break;
        }
    }

  return 0;
}

int calc_surface_singularities_integer_constraint(
    boost::property_tree::ptree &pt,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::one_ring_face_at_point & orfap,
    const matrixst &tet,
    const matrixd &node,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularities_chain,
    const std::vector<std::deque<size_t> > &singularities_type,
    const matrixd & outside_face_normal)
{//detect surface conflits with inner singularity integer constraint
  //#define visual
  pt.put("surface_type.desc","surface type file");
  const size_t trivial_type = 9;
  matrixst outside_face_type = ones<size_t>(1,outside_face_idx.size())
      * trivial_type; //
  if(load_surface_type(pt.get<string>("surface_type.value").c_str(),
                       fa,outside_face_idx,outside_face_type))
    return __LINE__;

  //  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(outside_face));
  matrixst outside_face_fix_axis(outside_face_idx.size());
  // first step: detect all surface triangle's type
  {
    map<size_t,vector<point_fix_axis_per_face> > point_with_fix_axis_on_face;
    for(size_t t = 0 ; t < outside_face.size(2); ++t){
        const size_t & face_type = outside_face_type[t];
        size_t fix_uvw = get_surface_fix_axis(face_type);
        outside_face_fix_axis[t] = fix_uvw;
        if(fix_uvw == -1) {
            cerr << "# [error] can not find a suitable axis." << endl;
            return __LINE__;
          }
        point_fix_axis_per_face pfapf;
        pfapf.face_idx = outside_face_idx[t];
        pfapf.fix_axis = fix_uvw;
        for(size_t i = 0; i < outside_face.size(1); ++i){
            point_with_fix_axis_on_face[outside_face(i,t)].push_back(pfapf);
          }
      }

#ifdef visual
    {
      ofstream ofs("surface_with_fix_axis.vtk");
      tri2vtk(ofs,&node[0],node.size(2),&outside_face[0],outside_face.size(2));
      cell_data(ofs,&outside_face_fix_axis[0],outside_face_fix_axis.size(),"fix_axis");
    }
#endif

#ifdef  visual
    typedef  map<size_t,vector<point_fix_axis_per_face> >::const_iterator msvpcit;
    for(msvpcit mcit = point_with_fix_axis_on_face.begin();
        mcit != point_with_fix_axis_on_face.end(); ++mcit){
        cerr << "# [info] surface point " << mcit->first << " :" << endl;
        const vector<point_fix_axis_per_face> & pfapf = mcit->second;
        for(size_t t = 0; t < pfapf.size(); ++t){
            cerr << "# ------ face " << pfapf[t].face_idx << " fix ";
            if(pfapf[t].fix_axis == 0) cerr << "u " << endl;
            if(pfapf[t].fix_axis == 1) cerr << "v " << endl;
            if(pfapf[t].fix_axis == 2) cerr << "w " << endl;
          }
      }
#endif

    // add surface integer setting
    matrixst idx_to_real_point(point_with_fix_axis_on_face.size());
    matrixst uvw_of_point_on_surface(3,point_with_fix_axis_on_face.size());
    matrix<size_t*> uvw_of_point_on_surface_ptr(uvw_of_point_on_surface.size(1),
                                                uvw_of_point_on_surface.size(2));
    {
      uvw_of_point_on_surface(colon()) = colon(0,3*point_with_fix_axis_on_face.size()-1);

      map<size_t,size_t> real_point_to_idx;
      size_t idx = 0;
      for(map<size_t,vector<point_fix_axis_per_face> >::const_iterator
          msvcit = point_with_fix_axis_on_face.begin();
          msvcit != point_with_fix_axis_on_face.end(); ++msvcit, ++idx){
          idx_to_real_point[idx] = msvcit->first;
          real_point_to_idx[msvcit->first] = idx;
        }

      for(size_t t = 0; t < uvw_of_point_on_surface_ptr.size(); ++t){
          uvw_of_point_on_surface_ptr[t] = &uvw_of_point_on_surface[t];
        }


      if(constrate_axis_fix_ptr_surface(outside_face,outside_face_type,
                                        real_point_to_idx,
                                        idx_to_real_point,
                                        ea,
                                        uvw_of_point_on_surface_ptr))
        return __LINE__;

      // assign uvw integer coordinate
      matrixst uvw_integer_fix =
          -1*ones<size_t>(uvw_of_point_on_surface.size(1),
                          uvw_of_point_on_surface.size(2));
      size_t integer = 0;
      for(size_t t = 0; t < uvw_of_point_on_surface.size(); ++t){
          if(uvw_of_point_on_surface[t] == t) continue;
          if(uvw_integer_fix[uvw_of_point_on_surface[t]] == -1){
              uvw_integer_fix[uvw_of_point_on_surface[t]] = integer;
              uvw_integer_fix[t] = integer;
              ++integer;
            }else{
              uvw_integer_fix[t] = uvw_integer_fix[uvw_of_point_on_surface[t]];
            }
        }

#ifdef  visual
      {
        map<size_t,vector<size_t> > group_point;
        for(size_t t = 0; t < uvw_of_point_on_surface_ptr.size() ; ++t){
            if(*uvw_of_point_on_surface_ptr[t] == t) continue;
            group_point[*(uvw_of_point_on_surface_ptr[t])].push_back(idx_to_real_point(t/3));
          }
        size_t i = 0;
        for(map<size_t,vector<size_t> >::iterator mcit = group_point.begin();
            mcit != group_point.end(); ++mcit,++i){
            vector<size_t> & point_g = mcit->second;
            point_g.push_back(idx_to_real_point(mcit->first/3));
            ostringstream os;
            os << i;
            ofstream ofs((os.str() + "point.vtk").c_str());
            point2vtk(ofs,&node[0],node.size(2),&point_g[0],point_g.size());
          }
      }
#endif

      cerr << "# [info] start to check surface degeneration." << endl;
      { // check degeneration
        map<vector<size_t>,vector<size_t> > uvw_collect;
        vector<size_t> temp(3);
        for(size_t t = 0; t < uvw_integer_fix.size(2); ++t){
            copy(uvw_integer_fix(colon(),t).begin(),
                 uvw_integer_fix(colon(),t).end(),temp.begin());
            uvw_collect[temp].push_back(idx_to_real_point[t]);
          }

        matrixst outside_face_with_patch_idx(outside_face.size(2));
        get_surface_patch_type_mat(outside_face_with_patch_idx,outside_face_type,
                                   tet,node,outside_face,outside_face_idx,fa);

        map<size_t,vector<size_t> > patch_idx_outside_face;
        {
          for(size_t t = 0; t < outside_face_with_patch_idx.size(); ++t){
              patch_idx_outside_face[outside_face_with_patch_idx[t]].push_back(t);
            }
        }


        vector<deque<pair<size_t,size_t> > > surface_chains;
        { // extract all surface singularity chains
          set<size_t> surface_singularity_points;
          for(map<vector<size_t>,vector<size_t> >::const_iterator mvvcit =
              uvw_collect.begin();mvvcit != uvw_collect.end(); ++mvvcit){
              const vector<size_t> & uvw = mvvcit->first;
              const vector<size_t> & linked_points = mvvcit->second;
              if(count(uvw.begin(),uvw.end(),-1) < 2){ // on surface edge
                  surface_singularity_points.insert(linked_points.begin(),
                                                    linked_points.end());
                }
            }
          vector<size_t> surface_singularity_points_arrray(surface_singularity_points.size());
          copy(surface_singularity_points.begin(),
               surface_singularity_points.end(),
               surface_singularity_points_arrray.begin());
          extract_chain_from_points(surface_singularity_points_arrray,ea,outside_face_fix_axis,surface_chains);
        }

#if 1
        {
          ofstream ofs("surface_chain.vtk");
          vector<size_t> lines;
          vector<size_t> line_type;
          for(size_t t = 0; t < surface_chains.size(); ++t){
              const deque<pair<size_t,size_t> > & one_chain = surface_chains[t];
              for(size_t i = 0; i < one_chain.size(); ++i){
                  lines.push_back(one_chain[i].first);
                  lines.push_back(one_chain[i].second);
                  line_type.push_back(t);
                }
            }

          line2vtk(ofs,&node[0],node.size(2),&lines[0],lines.size()/2);
          cell_data(ofs,&line_type[0],line_type.size(),"chain_type");

          vector<deque<size_t> > surface_type(surface_chains.size());
          for(size_t t = 0; t < surface_type.size(); ++t)
            surface_type[t].resize(surface_chains[t].size(),-1);
          dump_out_singularity_chain_3("surface_chain.sc",surface_chains,surface_type);
        }

#endif
#if 1
        {
          ofstream ofs("modified_surface_type_before.vtk");
          tri2vtk(ofs,&node[0],node.size(2),&outside_face[0],outside_face.size(2));
          cell_data(ofs,&outside_face_type[0],outside_face_type.size(),"face_type");
        }
#endif
        {
          // step 0: check duplicated uvw coordinates
          map<vector<size_t>,vector<size_t> > degenerated_points;
          for(map<vector<size_t>,vector<size_t> >::const_iterator mvvcit = uvw_collect.begin();
              mvvcit != uvw_collect.end(); ++mvvcit){
              if(mvvcit->second.size() < 2) continue;
              const vector<size_t> & uvw = mvvcit->first;
              if(find(uvw.begin(),uvw.end(),-1) == uvw.end()){
                  cerr << "# [error] meet fatal degeneration: ";
                  for(size_t t = 0; t < mvvcit->second.size(); ++t){
                      cerr << mvvcit->second[t] << " ";
                    }
                  cerr << endl;
                  //remove_surface_degeneration_to_point()
                  degenerated_points.insert(*mvvcit);
                }
            }

          remove_surface_degeneration_to_point(node,degenerated_points,
                                               surface_chains,
                                               outside_face_normal,
                                               outside_face_fix_axis,
                                               outside_face_with_patch_idx,
                                               patch_idx_outside_face,ea,orfap,
                                               outside_face,
                                               outside_face_type);
#if 1
          {
            ofstream ofs("modified_surface_type1.vtk");
            tri2vtk(ofs,&node[0],node.size(2),&outside_face[0],outside_face.size(2));
            cell_data(ofs,&outside_face_type[0],outside_face_type.size(),"face_type");
          }
#endif
        } // end step 0

        { // step 1: check surface edge degeneration, all surface points which is
          // fixed with two axis should not be able to make a loop,
          // or the loop must be degenerated.
          vector<size_t> degenerate_edge_chain_idx;
          for(size_t t = 0; t < surface_chains.size(); ++t){
              const deque<pair<size_t,size_t> > & one_chain = surface_chains[t];
              if(is_surface_singularity_degenerated_to_edge(one_chain,real_point_to_idx,
                                                            uvw_integer_fix))
                {
                  degenerate_edge_chain_idx.push_back(t);
                  cerr << "# [info] meet surface edge degeneration. chain " << t << endl;
                }
              //            if(one_chain.front().first == one_chain.back().second){
              //              degenerate_edge_chain_idx.push_back(t);
              //            }
            }

          remove_surface_degeneration_to_edge(node,surface_chains,
                                              degenerate_edge_chain_idx,
                                              outside_face_with_patch_idx,
                                              patch_idx_outside_face,
                                              ea,uvw_integer_fix,real_point_to_idx,
                                              outside_face,
                                              outside_face_type);

          set<size_t> point_to_glu_chains;
          for(size_t t = 0; t < degenerate_edge_chain_idx.size(); ++t){
              const deque<pair<size_t,size_t> > & one_chain =
                  surface_chains[degenerate_edge_chain_idx[t]];
              point_to_glu_chains.insert(one_chain.front().first);
            }
          map<size_t,vector<size_t> > point_adjacent_chain_idx;
          for(size_t t  = 0; t < surface_chains.size(); ++t){
              const deque<pair<size_t,size_t> > & one_chain = surface_chains[t];
              if(find(degenerate_edge_chain_idx.begin(),
                      degenerate_edge_chain_idx.end(),
                      t) == degenerate_edge_chain_idx.end()){
                  point_adjacent_chain_idx[one_chain.front().first].push_back(t);
                  point_adjacent_chain_idx[one_chain.back().second].push_back(t);
                }
            }
          vector<deque<pair<size_t,size_t> > > combined_chain;
          vector<deque<pair<size_t,size_t> > > temp_chain_to_replace_surfac_chain;
          vector<pair<size_t,size_t> > edges_need_to_extract_chain;
          typedef map<size_t,vector<size_t> >::const_iterator msvcit;
          for(set<size_t>::const_iterator scit = point_to_glu_chains.begin();
              scit != point_to_glu_chains.end(); ++scit){
              msvcit it = point_adjacent_chain_idx.find(*scit);
              if(it == point_adjacent_chain_idx.end()){
                  continue;
                }
              if(it->second.size() == 2){
                  edges_need_to_extract_chain.clear();
                  const deque<pair<size_t,size_t> > & chain0 = surface_chains[it->second.front()];
                  const deque<pair<size_t,size_t> > & chain1 = surface_chains[it->second.back()];
                  edges_need_to_extract_chain.insert(edges_need_to_extract_chain.end(),
                                                     chain0.begin(),chain0.end());
                  edges_need_to_extract_chain.insert(edges_need_to_extract_chain.end(),
                                                     chain1.begin(),chain1.end());
                  jtf::util::extract_chain_from_edges(edges_need_to_extract_chain,combined_chain);
                  temp_chain_to_replace_surfac_chain.insert(temp_chain_to_replace_surfac_chain.end(),
                                                            combined_chain.begin(),
                                                            combined_chain.end());
                  degenerate_edge_chain_idx.push_back(it->second.front());
                  degenerate_edge_chain_idx.push_back(it->second.back());
                }
            }
          // temp_surface_chains.reserve(surface_chains.size());
          for(size_t t = 0; t < surface_chains.size(); ++t){
              if(find(degenerate_edge_chain_idx.begin(),
                      degenerate_edge_chain_idx.end(),
                      t) == degenerate_edge_chain_idx.end()){
                  temp_chain_to_replace_surfac_chain.push_back(surface_chains[t]);
                }
            }
          surface_chains = temp_chain_to_replace_surfac_chain;
        } // end step 1
#if 1
        {
          ofstream ofs("modified_surface_type0.vtk");
          tri2vtk(ofs,&node[0],node.size(2),&outside_face[0],outside_face.size(2));
          cell_data(ofs,&outside_face_type[0],outside_face_type.size(),"face_type");
        }
#endif

      }// complete surface degeneration check
    }



    pt.put("dump_surface_integer_cons.desc","file recored surface points' integer constraint under surface type.");
    dump_surface_integer_cons(pt.get<string>("dump_surface_integer_cons.value").c_str(),
                              uvw_of_point_on_surface_ptr,idx_to_real_point);
  }

  pt.put("dump_refined_surface_type.desc","dump out the refined surface type");
  dump_surface_normal_align_type(pt.get<string>("dump_refined_surface_type.value").c_str(),
                                 outside_face,outside_face_type);


  return 0;
}

/**
 * @brief
 *
 * @param current_point_idx
 * @param linked_point_ptr
 * @param singularity_point_uvw_ptr
 * @return bool
 */
static bool is_all_point_constrated(const size_t current_point_idx,
                                    const set<size_t*> &linked_point_ptr,
                                    const vector<set<size_t*> > &singularity_point_uvw_ptr)
{
  for(set<size_t*>::const_iterator scit = linked_point_ptr.begin();
      scit != linked_point_ptr.end(); ++scit){
      if(*(*scit) == current_point_idx) continue;
      const set<size_t*> & linked_deep_point_ptr = singularity_point_uvw_ptr[*(*scit)];
      if(!linked_deep_point_ptr.empty()) return false;
    }
  return true;
}

bool is_exist_inner_loop(
    const std::vector<size_t> &same_coord_points,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &singularity_edge_map_to_chain_idx,
    std::vector<deque<pair<size_t,size_t> > > & loop_chains)
{
  vector<pair<size_t,size_t> > possible_edges;
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
  for(size_t t = 0;  t < same_coord_points.size(); ++t){
      for(size_t i = t +1; i < same_coord_points.size(); ++i){
          if(same_coord_points[t] > same_coord_points[i]){
              mpscit it = singularity_edge_map_to_chain_idx.find(
                    make_pair(same_coord_points[i],same_coord_points[t]));
              if(it == singularity_edge_map_to_chain_idx.end()) continue;
              possible_edges.push_back(make_pair(same_coord_points[i],same_coord_points[t]));
            }else{
              mpscit it = singularity_edge_map_to_chain_idx.find(
                    make_pair(same_coord_points[t],same_coord_points[i]));
              if(it == singularity_edge_map_to_chain_idx.end()) continue;
              possible_edges.push_back(make_pair(same_coord_points[t],same_coord_points[i]));
            }
        }
    }// end collect possible edgs;

  vector<deque<pair<size_t,size_t> > > chains;
  jtf::util::extract_chain_from_edges(possible_edges,chains);

  if(chains.size() == 1 &&
     chains.front().front().first == chains.front().back().second){
      loop_chains = chains;
      return true;
    }


  map<size_t,set<size_t> > point_to_chain;
  for(size_t t = 0; t < chains.size(); ++t){
      const deque<pair<size_t,size_t> > & one_chain = chains[t];
      point_to_chain[one_chain.front().first].insert(t);
      point_to_chain[one_chain.back().second].insert(t);
    }

  typedef map<size_t,set<size_t> >::const_iterator msscit;
  set<size_t> possible_loop_chain;
  for(size_t t = 0; t < chains.size(); ++t){
      const deque<pair<size_t,size_t> > & one_chain = chains[t];
      msscit it_first = point_to_chain.find(one_chain.front().first);
      msscit it_second = point_to_chain.find(one_chain.back().second);
      if(it_first == point_to_chain.end() ||
         it_second == point_to_chain.end()) {
          cerr << "# [error] strange can not find the point  in point_to_chain." << endl;
          return false;
        }
      if(it_first->second.size() < 3 || it_second->second.size() < 3) continue; // regular chain not parts of a loop
      possible_loop_chain.insert(t);
    }

  loop_chains.clear();
  while(!possible_loop_chain.empty())
    {
      deque<pair<size_t,size_t> > first_chain
          = chains[*possible_loop_chain.begin()];
      possible_loop_chain.erase(possible_loop_chain.begin());

      for(set<size_t>::iterator sit = possible_loop_chain.begin();
          sit != possible_loop_chain.end(); ++sit){
          const deque<pair<size_t,size_t> > & other_chain = chains[*sit];
          if(other_chain.front().first == first_chain.front().first &&
             other_chain.back().second == first_chain.back().second){
              combine_two_singualrity(first_chain,other_chain);
              loop_chains.push_back(first_chain);
              possible_loop_chain.erase(sit);
              break;
            }
          if(other_chain.front().first == first_chain.back().second &&
             other_chain.back().second == first_chain.front().first){
              combine_two_singualrity(first_chain,other_chain);
              loop_chains.push_back(first_chain);
              possible_loop_chain.erase(sit);
              break;
            }
        } // end for
    }
  if(loop_chains.empty()) return false;
  return true;
}



template <typename container>
bool size_comp(const container& a, const container& b)
{
  return a.size() < b.size();
}

int calc_singularities_integer_constraint_per_edge(
    boost::property_tree::ptree &pt,
    boost::unordered_map<pair<size_t,size_t>,size_t> inner_face_jump_type,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const matrixst &tet,
    const matrixd &node,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularities_chain,
    const std::vector<std::deque<size_t> > &singularities_type)
{
  //#define VISUAL_INNER
  pt.put("singularity_degenerated_points.desc","store degenerated points in singularity");
  // store the fixed axis and corresponding chain idx for each singularity ends
  typedef boost::unordered_map<size_t,vector<point_fix_axis_per_edge> >::const_iterator msvcit;
  boost::unordered_map<size_t,vector<point_fix_axis_per_edge> > point_with_fix_axis;

  for(size_t t = 0; t < singularities_chain.size(); ++t){
      const deque<pair<size_t,size_t> > & chain = singularities_chain[t];
      for(size_t i = 0; i < chain.size(); ++i){
          const pair<size_t,size_t> & edge = chain[i];
          const size_t around_axis = axis_to_around(singularities_type[t][i]);
          if(around_axis == -1) {
              cerr << "# [error] around axis error." << endl;
              return __LINE__;
            }
          point_fix_axis_per_edge a;
          a.edge = edge;
          a.fix_axis.push_back((around_axis+1)%3);
          a.fix_axis.push_back((around_axis+2)%3);

          point_with_fix_axis[edge.first].push_back(a); // fix other two axis
          point_with_fix_axis[edge.second].push_back(a);
        }
    }

  {// show the integer equations
    cout << "# [info] build integer constraint for each singularity point:" << endl;
#ifdef VISUAL_INNER
    typedef map<size_t,vector<point_fix_axis_per_edge> >::const_iterator msvcit;
    for(msvcit mcit = point_with_fix_axis.begin();
        mcit != point_with_fix_axis.end(); ++mcit){
        cout << "# *** at point: " << mcit->first << endl;
        const vector<point_fix_axis_per_edge> & vpfa = mcit->second;
        for(size_t i = 0; i < vpfa.size(); ++i){
            cout << "# ++++++ edge: <" << vpfa[i].edge.first << ","
                 << vpfa[i].edge.second << "> fix: ";
            for(size_t j = 0; j < vpfa[i].fix_axis.size(); ++j){
                if(vpfa[i].fix_axis[j] == 0)
                  cout << "u ";
                if(vpfa[i].fix_axis[j] == 1)
                  cout << "v ";
                if(vpfa[i].fix_axis[j] == 2)
                  cout << "w ";
              }
            cout << endl;
          }
      }
#endif

    matrixst singularity_point_uvw_idx(3, point_with_fix_axis.size());
    singularity_point_uvw_idx(colon()) = colon(0,3 * point_with_fix_axis.size()-1);
    //matrix<size_t*> singularity_point_uvw_ptr(3,point_with_fix_axis.size());
    vector<set<size_t*> > singularity_point_uvw_ptr(3 * point_with_fix_axis.size());
    for(size_t t = 0; t < singularity_point_uvw_idx.size(); ++t)
      singularity_point_uvw_ptr[t].insert(&singularity_point_uvw_idx[t]);

    boost::unordered_map<size_t,size_t> point2idx;
    matrixst point_idx(1,point_with_fix_axis.size());
    size_t idx = 0;
    for(msvcit mcit = point_with_fix_axis.begin();
        mcit != point_with_fix_axis.end(); ++mcit,++idx){
        point2idx[mcit->first] = idx;
        point_idx[idx] = mcit->first;
      }

    // link all points with fix relation
    size_t temp = 0;
    for(msvcit mit = point_with_fix_axis.begin();
        mit != point_with_fix_axis.end(); ++mit){// for each point in singularity
        const size_t &point = mit->first;
        const vector<point_fix_axis_per_edge> & vpfa = mit->second;
        for(size_t t = 0; t < vpfa.size(); ++t){ // for adjacent point which links current point
            const point_fix_axis_per_edge & pfa = vpfa[t];
            assert(pfa.edge.first == point || pfa.edge.second == point);
            const size_t other_point = pfa.edge.first + pfa.edge.second - point;
            assert(pfa.fix_axis.size() == 2);

            for(size_t j = 0; j < pfa.fix_axis.size(); ++j){
                singularity_point_uvw_ptr[3 * point2idx[point] + pfa.fix_axis[j]].insert(
                      &(singularity_point_uvw_idx[3 * point2idx[other_point]+ pfa.fix_axis[j]]));
              }
          }
      }

    for(size_t t = 0; t < singularity_point_uvw_ptr.size(); ++t){
        constrate_axis_fix_ptr(t,singularity_point_uvw_idx,singularity_point_uvw_ptr);
      }

    for(size_t t = 0; t < singularity_point_uvw_ptr.size(); ++t){
        set<size_t*> & linked_point_ptr = singularity_point_uvw_ptr[t];
        if(linked_point_ptr.empty()) continue;
        while(!is_all_point_constrated(t,linked_point_ptr,singularity_point_uvw_ptr)){
            for(set<size_t*>::const_iterator scit = linked_point_ptr.begin();
                scit != linked_point_ptr.end(); ++scit){
                if(t == *(*scit)) continue;
                set<size_t*> & linked_deep_point_ptr = singularity_point_uvw_ptr[*(*scit)];
                if(linked_deep_point_ptr.empty()) continue;
                linked_point_ptr.insert(linked_deep_point_ptr.begin(),
                                        linked_deep_point_ptr.end());
                linked_deep_point_ptr.clear();
              }
          }
      }
#if 0 // check if this constration is enough
    for(size_t t = 0; t < singularity_point_uvw_ptr.size(); ++t){
        set<size_t*> & link_point = singularity_point_uvw_ptr[t];
        if(link_point.empty()) continue;
        while(1){
            bool is_all_point_constrated = true;
            for(set<size_t*>::const_iterator sit = link_point.begin();
                sit != link_point.end(); ++sit){
                if(*(*sit) == t) continue;
                const set<size_t*> & other_link = singularity_point_uvw_ptr[*(*sit)];
                if(other_link.empty()) continue;
                else{
                    is_all_point_constrated = false;
                    link_point.insert(other_link.begin(),other_link.end());
                    break;
                  }
              }
            if(is_all_point_constrated) break;
            cerr << "# [info] really needs to post constrate." << endl;
          }
      }
#endif

    {// set all integer
      matrixst uvw_integre = ones<size_t>(3,singularity_point_uvw_idx.size(2)) * -1;
      size_t integer = 0;
      for(size_t t = 0; t < singularity_point_uvw_ptr.size(); ++t){
          const set<size_t*> & linked_points = singularity_point_uvw_ptr[t];
          if(linked_points.empty()) continue;
          if(linked_points.size() == 1) continue; // it's varying point
          for(set<size_t*>::const_iterator scit = linked_points.begin();
              scit != linked_points.end(); ++scit){
              if(uvw_integre[*(*scit)] != -1) {
                  cerr << "# [error] point " << point2idx[*(*scit)/3] << " uvw: " <<  *(*scit)%3
                       << " has been fixed " << endl;
                }else
                uvw_integre[*(*scit)] = integer;
            }
          ++integer;
        }
#ifdef VISUAL_INNER
      for(map<size_t,size_t>::const_iterator mcit = point2idx.begin();
          mcit != point2idx.end(); ++mcit){
          cerr << "# point " << mcit->first << " uvw : ";
          cerr << uvw_integre(0,mcit->second) << " "
               << uvw_integre(1,mcit->second) << " "
               << uvw_integre(2,mcit->second) << endl;
        }
#endif
      // check whether such uvw will result in degeneration
      {
        boost::unordered_map<vector<size_t>,vector<size_t> > uvw_col;
        vector<size_t> temp(3);

        for(size_t t = 0; t < uvw_integre.size(2); ++t){
            copy(uvw_integre(colon(),t).begin(),
                 uvw_integre(colon(),t).end(),temp.begin());
            uvw_col[temp].push_back(point_idx[t]);
          }

        boost::unordered_map<pair<size_t,size_t>,size_t> singularity_edge_with_idx;
        boost::unordered_map<pair<size_t,size_t>,size_t> singularity_edge_with_type;
        for(size_t t = 0; t < singularities_chain.size(); ++t){
            const deque<pair<size_t,size_t> > & one_chain = singularities_chain[t];
            for(size_t i = 0; i < one_chain.size(); ++i){
                singularity_edge_with_type[one_chain[i]] = singularities_type[t][i];
                if(one_chain[i].first > one_chain[i].second){
                    singularity_edge_with_idx[make_pair(one_chain[i].second,
                                                        one_chain[i].first)] = t;
                  }else{
                    singularity_edge_with_idx[make_pair(one_chain[i].first,
                                                        one_chain[i].second)] = t;
                  }
              }
          }


        {
          // degeneration into points
          cerr << "# [info] check inner singularity degeneration to points." << endl;
          vector<vector<size_t> > degenerated_points;
          for(boost::unordered_map<vector<size_t>,vector<size_t> >::const_iterator mcit =
              uvw_col.begin(); mcit != uvw_col.end(); ++mcit){
              const vector<size_t> & uvw = mcit->first;
              if(find(uvw.begin(),uvw.end(),-1) != uvw.end()) continue;// -1 means varying point
              if(mcit->second.size() != 1){
                  cerr << "# [error] meet inner singularity fatal degeneration: ";
                  degenerated_points.push_back(mcit->second);
                  for(size_t t = 0; t < mcit->second.size(); ++t){
                      cerr << mcit->second[t] << " ";
                    }
                  cerr << endl;
                }
            }

          if(dump_singularity_degenerated_points(
               pt.get<string>("singularity_degenerated_points.value").c_str(),
               degenerated_points))
            return __LINE__;

          std::sort(degenerated_points.begin(),degenerated_points.end(),size_comp<vector<size_t> >);

          boost::unordered_map<size_t,set<pair<size_t,size_t> > > singularity_point_with_linked_edge;
          for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mpscit =
              singularity_edge_with_type.begin();
              mpscit != singularity_edge_with_type.end(); ++mpscit){
              singularity_point_with_linked_edge[mpscit->first.first].insert(mpscit->first);
              singularity_point_with_linked_edge[mpscit->first.second].insert(mpscit->first);
            }

          typedef boost::unordered_map<size_t,set<pair<size_t,size_t> > >::const_iterator msspit;
          bool is_remove_one_degenerated_point = false;
          for(size_t t = 0; t < degenerated_points.size(); ++t){
              const vector<size_t> & linked_points = degenerated_points[t];
              for(size_t i = 0; i < linked_points.size(); ++i){
                  msspit it = singularity_point_with_linked_edge.find(linked_points[i]);
                  if(it == singularity_point_with_linked_edge.end()){
                      cerr << "# [error] strange can not find this point in linking map." << endl;
                      return __LINE__;
                    }
                  if(it->second.size() < 3) continue; // only has 2 linked edge or one
                  //              if(!remove_singularity_degenerated_point(
                  //                   it->first, it->second,singularity_edge_with_type,
                  //                   inner_face_jump_type,tet,fa))
                  //                break;
                }// end find a linked point to be addressed
              if(is_remove_one_degenerated_point) break;
            }
        } // end degeneration into points


        { // degeneration into edges
          cerr << "# [info] check inner singularity degeneration to edges." << endl;
          vector<deque<pair<size_t,size_t> > > loops;
          for(boost::unordered_map<vector<size_t>,vector<size_t> >::const_iterator mcit =
              uvw_col.begin(); mcit != uvw_col.end(); ++mcit)
            {
              const vector<size_t> & uvw_ = mcit->first;
              const vector<size_t> & same_coord_points = mcit->second;
              vector<deque<pair<size_t,size_t> > > loop_chains;
              if(count(uvw_.begin(), uvw_.end(),-1) == 1){ // points fixed two integer
                  if(is_exist_inner_loop(same_coord_points,
                                         singularity_edge_with_idx,
                                         loop_chains))
                    {
                      assert(!loop_chains.empty());
                      loops.insert(loops.end(),loop_chains.begin(),loop_chains.end());
                      cerr << "# [error] meet inner singularity loop degeneration." << endl;
                    }
                }
            }
          //#define VISUAL_INNER
#ifdef VISUAL_INNER
          {
            ofstream ofs("singularity_edge_degenerated.vtk");
            vector<size_t> edges;
            for(size_t t = 0; t < loops.size(); ++t){
                const deque<pair<size_t,size_t> > & one_chain = loops[t];
                for(size_t i = 0; i < one_chain.size(); ++i){
                    edges.push_back(one_chain[i].first);
                    edges.push_back(one_chain[i].second);
                  }
              }
            line2vtk(ofs,&node[0],node.size(2),&edges[0],edges.size()/2);
          }
#endif

          remove_singularity_degenerated_edge(loops,inner_face_jump_type,tet,node,ortae,fa);

        }// end degenerated into edges
      }

#if 0 // check 1
      for(size_t t = 0; t < singularities_chain.size(); ++t){
          const deque<pair<size_t,size_t> > & chain = singularities_chain[t];
          for(size_t i = 0; i < chain.size(); ++i){
              const pair<size_t,size_t> & edge = chain[i];
              const size_t around_axis = axis_to_around(singularities_type[t][i]);
              if((uvw_integre((around_axis+1)%3,point2idx[edge.first]) !=
                  uvw_integre((around_axis+1)%3,point2idx[edge.second])) ||
                 (uvw_integre((around_axis+2)%3,point2idx[edge.first]) !=
                  uvw_integre((around_axis+2)%3,point2idx[edge.second]))
                 ){
                  cerr << "# [error] edge <" << edge.first << "," << edge.second
                       << "> should fix " << (around_axis+1)%3 << " "
                       << (around_axis+2)%3 << endl;
                  cerr << "# but uvw of two points: " << endl;
                  cerr << uvw_integre(0,point2idx[edge.first]) << " "
                                                               << uvw_integre(1,point2idx[edge.first]) << " "
                                                                                                       << uvw_integre(2,point2idx[edge.first]) << endl;
                  cerr << uvw_integre(0,point2idx[edge.second]) << " "
                                                                << uvw_integre(1,point2idx[edge.second]) << " "
                                                                                                         << uvw_integre(2,point2idx[edge.second]) << endl;
                }
            }
        }
#endif


#if 1 // check the uvw coordinates
      {
        for(size_t t = 0; t < singularities_chain.size(); ++t){
            const deque<pair<size_t,size_t> > & one_chain = singularities_chain[t];
            for(size_t i = 0; i < one_chain.size(); ++i){
                const size_t fix_axis = axis_to_around(singularities_type[t][i]);
                if((uvw_integre((fix_axis+1)%3,point2idx[one_chain[i].first]) !=
                    uvw_integre((fix_axis+1)%3,point2idx[one_chain[i].second])) ||
                   (uvw_integre((fix_axis+2)%3,point2idx[one_chain[i].first]) !=
                    uvw_integre((fix_axis+2)%3,point2idx[one_chain[i].second]))){
                    cerr << "# [error] uvw is not suitable on edge <"
                         << one_chain[i].first << "," << one_chain[i].second << endl;
                  }
              }
          }
      }
#endif

#ifdef VISUAL_INNER
      {
        for(size_t t = 0; t < singularity_point_uvw_ptr.size(); ++t){
            const set<size_t*> & linked_points = singularity_point_uvw_ptr[t];
            if(linked_points.empty()) continue;
            if(linked_points.size() == 1) continue; // it's varying point
            cerr << "# [info] point " << point_idx[t/3] << " ";
            if(t%3 == 0) cerr  << "u equal with: ";
            if(t%3 == 1) cerr  << "v equal with: ";
            if(t%3 == 2) cerr  << "w equal with: ";
            cerr << endl;

            for(set<size_t*>::const_iterator scit = linked_points.begin();
                scit != linked_points.end(); ++scit){
                const size_t linked_idx = *(*scit);
                cerr << "# ----- point " << point_idx[linked_idx/3] << " ";
                if(linked_idx%3 == 0) cerr  << "u ";
                if(linked_idx%3 == 1) cerr  << "v ";
                if(linked_idx%3 == 2) cerr  << "w ";
                cerr << endl;
              }
          }
      }
#endif
    }


    pt.put("dump_inner_face_jump_type.desc","dump inner face jump type file");
    dump_inner_face_jump_type(
          pt.get<string>("dump_inner_face_jump_type.value").c_str(),
          inner_face_jump_type);

    pt.put("dump_inner_singularity_integer_cons.desc","inner singularity integer constraint file");
    dump_inner_singularity_integer_cons(
          pt.get<string>("dump_inner_singularity_integer_cons.value").c_str(),
          singularity_point_uvw_ptr,point_idx);
  }

  return 0;
}

int get_tet_rot_type(
    const matrixst & cut_tet,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type,
    matrixst & tet_rot)
{
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
  stack<size_t> growing_tets;
  growing_tets.push(0);
  tet_rot.resize(cut_tet.size(2));
  itr_matrix<size_t*> tet_rot_mat(1,tet_rot.size(),&tet_rot[0]);
  tet_rot_mat = ones<size_t>(1,tet_rot.size()) * -1;
  tet_rot[0] = 9;// trivial type
  while(!growing_tets.empty()){
      const size_t tet_idx = growing_tets.top();
      growing_tets.pop();
      assert(tet_rot[tet_idx] != -1);
      for(size_t t = 0; t < 4; ++t){
          const pair<size_t,size_t> tet_pair = fa_cut.query(cut_tet(t,tet_idx),
                                                            cut_tet((t+1)%4,tet_idx),
                                                            cut_tet((t+2)%4,tet_idx));
          if(fa_cut.is_outside_face(tet_pair)) continue;
          assert(tet_pair.first == tet_idx || tet_pair.second == tet_idx);
          const size_t other_tet_idx = tet_pair.first + tet_pair.second - tet_idx;
          if(tet_rot[other_tet_idx] == -1){
              mpscit mit = inner_face_jump_type.find(make_pair(tet_idx,
                                                               other_tet_idx));
              if(mit == inner_face_jump_type.end()){
                  assert(tet_rot[tet_idx] != -1);
                  tet_rot[other_tet_idx] = tet_rot[tet_idx];
                }
              else
                tet_rot[other_tet_idx] = type_transition1(type_transition2(tet_rot[tet_idx])
                                                          * type_transition2(mit->second));
              growing_tets.push(other_tet_idx);
            }
        }
    }
  return 0;
}

int get_tet_rot_type_vec(const matrixst & cut_tet,
                         const jtf::mesh::face2tet_adjacent &fa_cut,
                         const std::map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
                         std::vector<size_t> & tet_rot)
{
  typedef map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
  stack<size_t> growing_tets;
  growing_tets.push(0);
  tet_rot.resize(cut_tet.size(2));
  itr_matrix<size_t*> tet_rot_mat(1,tet_rot.size(),&tet_rot[0]);
  tet_rot_mat = ones<size_t>(1,tet_rot.size()) * -1;
  tet_rot[0] = 9;// trivial type
  while(!growing_tets.empty()){
      const size_t tet_idx = growing_tets.top();
      growing_tets.pop();
      assert(tet_rot[tet_idx] != -1);
      for(size_t t = 0; t < 4; ++t){
          const pair<size_t,size_t> tet_pair = fa_cut.query(cut_tet(t,tet_idx),
                                                            cut_tet((t+1)%4,tet_idx),
                                                            cut_tet((t+2)%4,tet_idx));
          if(fa_cut.is_outside_face(tet_pair)) continue;
          assert(tet_pair.first == tet_idx || tet_pair.second == tet_idx);
          const size_t other_tet_idx = tet_pair.first + tet_pair.second - tet_idx;
          if(tet_rot[other_tet_idx] == -1){
              mpscit mit = inner_face_jump_type.find(make_pair(tet_idx,
                                                               other_tet_idx));
              if(mit == inner_face_jump_type.end()){
                  assert(tet_rot[tet_idx] != -1);
                  tet_rot[other_tet_idx] = tet_rot[tet_idx];
                }
              else
                tet_rot[other_tet_idx] = type_transition1(type_transition2(tet_rot[tet_idx])
                                                          * type_transition2(mit->second));
              growing_tets.push(other_tet_idx);
            }
        }
    }
  return 0;
}

int calc_tet_rotation(
    const jtf::tet_mesh &tm,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    vector<size_t> & tet_rot)
{
  tetmesh_cutter tmc(tm);
  tmc.cut(inner_face_jump_type);

  unique_ptr<jtf::mesh::face2tet_adjacent> fa_cut(jtf::mesh::face2tet_adjacent::create(tmc.cut_tm_.mesh_));
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mpscit;

  stack<size_t> growing_tets;
  growing_tets.push(0);
  tet_rot.resize(tm.tetmesh_.mesh_.size(2));
  itr_matrix<size_t*> tet_rot_mat(1,tet_rot.size(),&tet_rot[0]);
  tet_rot_mat = ones<size_t>(1,tet_rot.size()) * -1;
  tet_rot[0] = TRIVIAL_TYPE;// trivial type
  while(!growing_tets.empty()){
      const size_t tet_idx = growing_tets.top();
      growing_tets.pop();
      assert(tet_rot[tet_idx] != -1);
      for(size_t t = 0; t < 4; ++t){
          const pair<size_t,size_t> tet_pair = fa_cut->query(tmc.cut_tm_.mesh_(t,tet_idx),
                                                             tmc.cut_tm_.mesh_((t+1)%4,tet_idx),
                                                             tmc.cut_tm_.mesh_((t+2)%4,tet_idx));
          if(fa_cut->is_outside_face(tet_pair)) continue;
          assert(tet_pair.first == tet_idx || tet_pair.second == tet_idx);
          const size_t other_tet_idx = tet_pair.first + tet_pair.second - tet_idx;
          if(tet_rot[other_tet_idx] == -1){
              mpscit mit = inner_face_jump_type.find(make_pair(tet_idx,
                                                               other_tet_idx));
              if(mit == inner_face_jump_type.end()){
                  assert(tet_rot[tet_idx] != -1);
                  tet_rot[other_tet_idx] = tet_rot[tet_idx];
                }
              else
                tet_rot[other_tet_idx] = type_transition1(type_transition2(tet_rot[tet_idx])
                                                          * type_transition2(mit->second));
              growing_tets.push(other_tet_idx);
            }
        }
    }

  for(size_t t = 0; t < tet_rot.size(); ++t){
      if(tet_rot[t] > 23)
        cerr << "# [error] tet " << t << " has not been calc correctly." << endl;
    }

  return 0;
}

int load_tet_rotation(const char * filename,
                      matrixst & tet_rot)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open tet rotation file" << endl;
      return __LINE__;
    }
  string trash;
  size_t tet_num,sed_tet;
  ifs >> trash >> tet_num >> trash >> sed_tet; // "tet n root 0"
  tet_rot.resize(1,tet_num);
  for(size_t t = 0; t < tet_num; ++t){
      ifs >> tet_rot[t];
    }
  return 0;
}


int load_tet_rotation_array(const char* filename,
                            std::vector<size_t> & tet_rot)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open tet rotation file" << endl;
      return __LINE__;
    }
  string trash;
  size_t tet_num,sed_tet;
  ifs >> trash >> tet_num >> trash >> sed_tet; // "tet n root 0"
  tet_rot.resize(tet_num);
  for(size_t t = 0; t < tet_num; ++t){
      ifs >> tet_rot[t];
    }
  return 0;
}


int dump_tet_rotation_from_array(const char* filename,
                                 const vector<size_t> & tet_rot)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      cerr << "#  [error] can not open tet rot file." << endl;
      return __LINE__;
    }

  ofs << "tet " << tet_rot.size() << " root 0" << endl;
  for(size_t t = 0; t < tet_rot.size(); ++t)
    ofs << tet_rot[t] << " ";

  return 0;
}

int dump_tet_rotation_from_raw(
    const char* filename,
    const size_t *tet_rot,
    const size_t tet_num)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      cerr << "#  [error] can not open tet rot file." << endl;
      return __LINE__;
    }

  if(!tet_rot)  {
      cerr << "# [error] invalid input tet_rot vec." << endl;
      return __LINE__;
    }

  ofs << "tet " << tet_num << " root 0" << endl;
  for(size_t t = 0; t < tet_num; ++t)
    ofs << tet_rot[t] << " ";

  return 0;
}

int dump_inner_singularity_integer_cons(const char * filename,
                                        const vector<set<size_t*> > &singularity_point_uvw_ptr,
                                        const matrixst & idx_to_real_point)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      cerr << "# [info] can not open inner_singularity_integer_cons file." << endl;
      return __LINE__;
    }
  for(size_t t = 0; t < singularity_point_uvw_ptr.size(); ++t){
      const set<size_t*> & linked_points = singularity_point_uvw_ptr[t];
      if(linked_points.empty()) continue;
      if(linked_points.size() == 1) continue; // it's varying coordinate
      ofs << "point " << idx_to_real_point[t/3] << " axis " << t%3 << endl;
      ofs << linked_points.size() << endl;
      for(set<size_t*>::const_iterator scit = linked_points.begin();
          scit != linked_points.end(); ++scit){
          const size_t &other_point = *(*scit);
          assert(other_point%3 == t%3);
          ofs << idx_to_real_point[other_point/3] << " ";
        }
      ofs << endl;
    }
  return 0;
}


int dump_surface_integer_cons(const char * filename,
                              const zjucad::matrix::matrix<size_t*> & uvw_of_point_on_surface_ptr,
                              const matrixst & idx_to_real_point)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      cerr << "# [error] can not open surface integer cons file." << endl;
      return __LINE__;
    }

  map<pair<size_t,size_t>,set<size_t> > fix_points;
  typedef map<pair<size_t,size_t>,set<size_t> >::const_iterator mpscit;
  for(size_t t = 1; t < uvw_of_point_on_surface_ptr.size(); ++t){
      const size_t &linked_coordinate = *(uvw_of_point_on_surface_ptr[t]);
      if(linked_coordinate == t) continue; // varying coordinate
      const size_t fixed_axis = t%3;
      fix_points[make_pair(idx_to_real_point[linked_coordinate/3],fixed_axis)].insert(idx_to_real_point[t/3]);
    }

  for(mpscit mit = fix_points.begin(); mit != fix_points.end(); ++mit){
      ofs << "point " << mit->first.first << " axis " << mit->first.second << endl;
      ofs << mit->second.size() + 1 << endl; // link to itself
      ofs << mit->first.first << " ";
      for(set<size_t>::const_iterator scit = mit->second.begin();
          scit != mit->second.end(); ++scit){
          ofs << *scit << " ";
        }
      ofs << endl;
    }
  return 0;
}


int dump_union_integer_constraints(
    const char * filename,
    const matrixst & uvw_of_point,
    const matrixst & idx_to_real_point)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      cerr << "# [error] can not open file." << endl;
      return __LINE__;
    }

  map<pair<size_t,size_t>,set<size_t> > same_coordinate_points;
  for(size_t t = 0; t < uvw_of_point.size(); ++t){
      if(uvw_of_point[t] == -1) continue;
      const size_t point_idx = t/3;
      const size_t fix_axis = t%3;
      same_coordinate_points[make_pair(uvw_of_point[t],fix_axis)].insert(idx_to_real_point[point_idx]);
    }

  typedef map<pair<size_t,size_t>,set<size_t> >::const_iterator mpscit;
  for(mpscit it = same_coordinate_points.begin();
      it != same_coordinate_points.end(); ++it){
      const set<size_t> & group_points = it->second;
      assert(!group_points.empty());
      ofs << "point " << *group_points.begin() << " axis " << it->first.second << endl;
      ofs << group_points.size() << endl;

      for(set<size_t>::const_iterator scit = group_points.begin();
          scit != group_points.end(); ++scit){
          ofs << *scit << " ";
        }
      ofs << endl;
    }
  return 0;
}


/**
 * @brief
 *
 * @param pt
 * @param uvw_integer
 * @param idx_to_real_point
 * @param map<size_t
 * @param real_point_to_idx
 * @return int
 */
static int gather_surface_and_inner_integer_cons(
    boost::property_tree::ptree & pt,
    matrixst & uvw_integer,
    matrixst &idx_to_real_point,
    boost::unordered_map<size_t,size_t> &real_point_to_idx)
{
  pt.put("inner_singularity_integer_cons.desc","file which contains integer constraints on inner singularity edges");
  boost::unordered_map<pair<size_t,size_t>,size_t> inner_singularity_integer_with_group_idx;
  vector<set<size_t> > inner_group;

  if(load_integer_constraint_group_info(
       pt.get<string>("inner_singularity_integer_cons.value").c_str(),
       inner_singularity_integer_with_group_idx,
       inner_group)){
      cerr << "# [error] open inner_singularity_integer_cons fail." << endl;
      return __LINE__;
    }


  pt.put("surface_integer_cons.desc","file which contains integer constraints on surface");

  boost::unordered_map<pair<size_t,size_t>,size_t> surface_integer_constraint_with_group_idx;
  vector<set<size_t> > surface_group;

  if(load_integer_constraint_group_info(
       pt.get<string>("surface_integer_cons.value").c_str(),
       surface_integer_constraint_with_group_idx,
       surface_group)){
      cerr << "# [error] open surface_integer_cons fail." << endl;
      return __LINE__;
    }

  boost::unordered_map<pair<size_t,size_t>,size_t> total_inner_cons =
      surface_integer_constraint_with_group_idx;

  const size_t surface_group_num = surface_group.size();

  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::iterator mpsit;

  while(!inner_singularity_integer_with_group_idx.empty()){
      const mpscit mit = inner_singularity_integer_with_group_idx.begin();

      const size_t & inner_integer_group_idx = mit->second;
      const set<size_t> & inner_group_set = inner_group[mit->second];
      const size_t & fix_axis = mit->first.second;

      set<size_t> adjacent_surface_group;
      for(set<size_t>::const_iterator scit = inner_group_set.begin();
          scit != inner_group_set.end(); ++scit){

          mpscit mit_ = surface_integer_constraint_with_group_idx.find(
                make_pair(*scit,fix_axis));
          if(mit_ != surface_integer_constraint_with_group_idx.end())
            {
              adjacent_surface_group.insert(mit_->second);
            }
        }
      if(!adjacent_surface_group.empty() && adjacent_surface_group.size() != 1){
          cerr << "# [error] inner group has linked to different surface groups." <<endl;
          cerr << "#         means they have different integer on same axis." << endl;
          return __LINE__;
        }

      if(adjacent_surface_group.empty()) { // it's an inner group, do not reach surface
          for(set<size_t>::const_iterator scit = inner_group_set.begin();
              scit != inner_group_set.end(); ++scit){

              //TODO: why is the index way get wrong value?
              //total_inner_cons[make_pair(*scit,fix_axis)] = mit->second + surface_group_num; // assign new idx
              mpsit itt = total_inner_cons.find(make_pair(*scit,fix_axis));
              if(itt == total_inner_cons.end())
                total_inner_cons.insert(make_pair(make_pair(*scit,fix_axis),
                                                  mit->second + surface_group_num));
              else
                itt->second = mit->second + surface_group_num;

              mpsit m = inner_singularity_integer_with_group_idx.find(make_pair(*scit,fix_axis));
              assert(m != inner_singularity_integer_with_group_idx.end());
              inner_singularity_integer_with_group_idx.erase(m);
            }
        }

      if(adjacent_surface_group.size() == 1){

          for(set<size_t>::const_iterator scit = inner_group_set.begin();
              scit != inner_group_set.end(); ++scit){
              //        if(*scit == 7399){
              //          cerr << " 7399 fix axis " <<  fix_axis << " idx " << *adjacent_surface_group.begin() << endl;
              //        }
              total_inner_cons[make_pair(*scit,fix_axis)] = *adjacent_surface_group.begin(); // unifor all linked idx
              mpsit m = inner_singularity_integer_with_group_idx.find(make_pair(*scit,fix_axis));
              assert(m != inner_singularity_integer_with_group_idx.end());
              inner_singularity_integer_with_group_idx.erase(m);
            }
        }
    }

  set<size_t> point_set;
  for(mpscit mit = total_inner_cons.begin(); mit != total_inner_cons.end(); ++mit){
      point_set.insert(mit->first.first);
    }

  uvw_integer = ones<size_t>(3,point_set.size()) * -1;
  size_t point_idx = 0;
  idx_to_real_point.resize(point_set.size());
  for(set<size_t>::const_iterator scit = point_set.begin();
      scit != point_set.end(); ++scit,++point_idx){

      idx_to_real_point[point_idx] = *scit;
      real_point_to_idx[*scit] = point_idx;

      for(size_t i  = 0; i < 3; ++i) //three axes
        {
          mpscit mit_ = total_inner_cons.find(make_pair(*scit,i));
          if(mit_ == total_inner_cons.end()) continue;
          uvw_integer(i,point_idx) = mit_->second;
        }
      assert((uvw_integer(0,point_idx) != uvw_integer(1,point_idx)) ||
             (uvw_integer(0,point_idx) == -1));

      assert((uvw_integer(1,point_idx) != uvw_integer(2,point_idx)) ||
             (uvw_integer(1,point_idx) == -1));

      assert((uvw_integer(2,point_idx) != uvw_integer(0,point_idx)) ||
             (uvw_integer(2,point_idx) == -1));
    }

  return 0;
}

int check_degeneration_under_inner_singularity_and_surface(
    boost::property_tree::ptree & pt,
    jtf::mesh::face2tet_adjacent & fa,
    vector<size_t> & tet_vec,
    vector<double> & node_vec,
    jtf::mesh::one_ring_tet_at_edge & ortae,
    boost::unordered_map<pair<size_t,size_t>,size_t> inner_face_jump_type,
    const jtf::mesh::edge2cell_adjacent & ea,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const vector<deque<pair<size_t,size_t> > > &singularity_chain,
    const vector<deque<size_t> > & singularity_type)
{
  matrixst uvw_integer;
  matrixst idx_to_real_point;
  boost::unordered_map<size_t,size_t> real_point_to_idx;
  gather_surface_and_inner_integer_cons(pt,uvw_integer,
                                        idx_to_real_point,real_point_to_idx);

  pt.put("dump_union_integer_cons.desc","dump out union points integer constraint");
  dump_union_integer_constraints(
        pt.get<string>("dump_union_integer_cons.value").c_str(),
        uvw_integer, idx_to_real_point);

  pt.put("surface_type.desc","surface type file");
  const string str_outside_face_type  = pt.get<string>("surface_type.value");
  matrixst outside_face_normal_align_type;

  if(load_surface_type(str_outside_face_type.c_str(),
                       fa,outside_face_idx,outside_face_normal_align_type))
    return __LINE__;

  boost::unordered_map<vector<size_t>,vector<size_t> > uvw_integer_vec;
  vector<size_t> uvw(3);
  for(size_t t = 0; t < uvw_integer.size(2); ++t){
      copy(uvw_integer(colon(),t).begin(),
           uvw_integer(colon(),t).end(),uvw.begin());
      uvw_integer_vec[uvw].push_back(t);
    }

#if 1
  {// check degeneration points
    vector<vector<size_t> > degenerated_points;
    for(boost::unordered_map<vector<size_t>,vector<size_t> >::const_iterator mvvcit =
        uvw_integer_vec.begin(); mvvcit != uvw_integer_vec.end(); ++mvvcit){
        const vector<size_t> & uvw_ = mvvcit->first;
        if(count(uvw_.begin(),uvw_.end(),-1) == 0){
            if(mvvcit->second.size() != 1){
                cerr << "# [error] meet union degeneration points: " ;
                copy(mvvcit->second.begin(),mvvcit->second.end(),
                     ostream_iterator<size_t>(cerr," "));
                cerr << endl;
                cerr << "# ------ uvw: " << uvw_[0] << " " << uvw_[1] << " " << uvw_[2] << endl;
                degenerated_points.push_back(mvvcit->second);
              }
          }
      }

    remove_union_degeneration_points(degenerated_points, outside_face,
                                     singularity_chain, singularity_type,
                                     outside_face_normal_align_type,
                                     inner_face_jump_type);
  } // end check degeneration points

  {// check degeneatoin edges
    vector<vector<size_t> > edge_degenerated_points_vec;
    set<size_t> singularity_points;
    for(size_t t  = 0; t < singularity_chain.size(); ++t){
        const deque<pair<size_t,size_t> > & one_chain = singularity_chain[t];
        for(size_t i = 0; i < one_chain.size(); ++i){
            singularity_points.insert(one_chain[i].first);
            singularity_points.insert(one_chain[i].second);
          }
      }
    for(boost::unordered_map<vector<size_t>,vector<size_t> >::const_iterator mvvcit =
        uvw_integer_vec.begin(); mvvcit != uvw_integer_vec.end(); ++mvvcit){
        const vector<size_t> & uvw_ = mvvcit->first;
        const vector<size_t> & degenerated_points_vec = mvvcit->second;
        if(count(uvw_.begin(),uvw_.end(),-1) == 1){
            if(is_union_edge_degenerated_case(degenerated_points_vec,uvw_integer,
                                              real_point_to_idx,singularity_points,
                                              singularity_chain,singularity_type,
                                              outside_face))
              {
                edge_degenerated_points_vec.push_back(degenerated_points_vec);
              }
          }
      }
    remove_union_degeneration_edge(edge_degenerated_points_vec, outside_face,
                                   singularity_chain, singularity_type,
                                   uvw_integer,ea,fa,tet_vec,node_vec,
                                   ortae,outside_face_normal_align_type,
                                   inner_face_jump_type);
  }
#endif

  { // dump out new inner face jump type and surface type
    pt.put("dump_refined_inner_face_jump_type.desc","dump out inner face jump type file");
    const string str_inner_face_jump_type_new = pt.get<string>("dump_refined_inner_face_jump_type.value");
    pt.put("dump_refined_surface_type.desc","dump out refined surface type file");
    const string str_surface_type = pt.get<string>("dump_refined_surface_type");
    dump_surface_normal_align_type(str_surface_type.c_str(),
                                   outside_face,outside_face_normal_align_type);
    dump_inner_face_jump_type(str_inner_face_jump_type_new.c_str(),
                              inner_face_jump_type);
    cerr << "# [info] finish dump out inner face jump type and surface type" << endl;
  }
  return 0;
}


bool is_union_edge_degenerated_case(
    const std::vector<size_t> &same_coord_points_vec,
    const matrixst &uvw_integer,
    const boost::unordered_map<size_t,size_t> & real_point_to_idx,
    const std::set<size_t> &singularity_points,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_chain,
    const std::vector<std::deque<size_t> > &singularity_type,
    const matrixst & outside_face)
{
  // there must be some degenerated points located on surface, or they should be fixed         in interior step
  vector<size_t> surface_points,inner_linked_points;
  for(size_t t = 0; t < same_coord_points_vec.size(); ++t){
      if(find(singularity_points.begin(),singularity_points.end(),
              same_coord_points_vec[t]) != singularity_points.end()){
          inner_linked_points.push_back(same_coord_points_vec[t]);
        }else{
          surface_points.push_back(same_coord_points_vec[t]);
        }
    }

  if(inner_linked_points.empty()) {
      return false; // it's regular surface singularity points
    }

  for(size_t t = 0; t < inner_linked_points.size(); ++t){
      if(is_outside_vertex(inner_linked_points[t],outside_face))
        continue;
      return true; // TODO: here needs to check , whether an inner points with the
      // same integer coordinate constraint is enough
    }
  return false;
}


template <typename edge_iterator >
static double calc_projected_area_with_point_and_boundary(
    const matrixd & node,
    const size_t point_idx,
    edge_iterator begin_,
    edge_iterator end_,
    const matrixd & normal_proj)
{
  BOOST_STATIC_ASSERT((boost::is_same<typename edge_iterator::value_type,
                       pair<size_t,size_t> >::value) ||
                      (boost::is_same<typename edge_iterator::value_type,
                       const pair<size_t,size_t> >::value));

  double area = 0;
  for(edge_iterator eit = begin_; eit != end_; ++eit){
      const pair<size_t,size_t> & edge = *eit;
      area += 0.5 * fabs(dot(cross(node(colon(),edge.first) - node(colon(),point_idx),
                                   node(colon(),edge.second) - node(colon(),point_idx)),
                             normal_proj));
    }
  return area;
}

int remove_one_inner_loop_by_modify_face_type(
    const std::deque<std::pair<size_t,size_t> > & one_loop,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const matrixst &tet,
    const matrixd &node)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  matrixd xyz = eye<double>(3);

  deque<list<pair<size_t,size_t> > > loop_edge_list_queue;
  {
    list<pair<size_t,size_t> > loop_edge_list;
    for(size_t t = 0; t < one_loop.size(); ++t)
      loop_edge_list.push_back(one_loop[t]);

    loop_edge_list_queue.push_back(loop_edge_list);
  }

  while(!loop_edge_list_queue.empty()){
      list<pair<size_t,size_t> > & one_loop = loop_edge_list_queue.front();

      { // check the loop whether needs to be splited or remvoe duplicated edges
        bool is_loop_empty = false;
        while(1){
            if(one_loop.empty()) {
                loop_edge_list_queue.pop_front();
                is_loop_empty = true;
                break;
              }
            const size_t prev_loop_length = one_loop.size();
            list<pair<size_t,size_t> >::iterator lpit_prev = one_loop.begin(),
                lpit_current = one_loop.begin();
            ++lpit_current;
            while(lpit_current != one_loop.end()){
                if(lpit_prev->first == lpit_current->second &&
                   lpit_prev->second == lpit_current->first){
                    one_loop.erase(lpit_prev);
                    one_loop.erase(lpit_current++);
                    lpit_prev = lpit_current++;
                  }else{
                    ++lpit_prev;
                    ++lpit_current;
                  }
              }// end passing one loop
            if(one_loop.size() == prev_loop_length) break; // if no edge is deleted jump out
          }// end while
        if(is_loop_empty) continue;
      }

      set<size_t> possible_points;
      map<pair<size_t,size_t>,vector<size_t> > edge_with_around_points;
      {
        for(list<pair<size_t,size_t> >::const_iterator lpscit = one_loop.begin();
            lpscit  != one_loop.end(); ++lpscit){
            vector<size_t> around_points;
            find_one_ring_vertex_around_edge(tet,*lpscit,ortae,around_points);
            possible_points.insert(around_points.begin(),
                                   around_points.end());
            edge_with_around_points[*lpscit] = around_points;
          }
      }

      //vector<pair<double,size_t> > least_area_points;
      vector<std::tuple<double,double,size_t> > least_area_points;

      for(set<size_t>::const_iterator scit = possible_points.begin();
          scit != possible_points.end(); ++scit){
          double area_x = calc_projected_area_with_point_and_boundary(
                node,*scit,one_loop.begin(),one_loop.end(),xyz(colon(),0));
          double area_y = calc_projected_area_with_point_and_boundary(
                node,*scit,one_loop.begin(),one_loop.end(),xyz(colon(),1));
          double area_z = calc_projected_area_with_point_and_boundary(
                node,*scit,one_loop.begin(),one_loop.end(),xyz(colon(),2));
          //least_area_points.push_back(make_pair(area_x * area_y * area_z,*scit));
          least_area_points.push_back(std::make_tuple(area_x * area_y * area_z,
                                                      area_x + area_y + area_z,
                                                      *scit));
        }

      sort(least_area_points.begin(),least_area_points.end());
      size_t selected_point = get<2>(least_area_points.front());

      for(map<pair<size_t,size_t>,vector<size_t> >::const_iterator mpvcit =
          edge_with_around_points.begin(); mpvcit !=  edge_with_around_points.end();
          ++mpvcit){
          const vector<size_t> & around_points = mpvcit->second;
          const pair<size_t,size_t> & edge = mpvcit->first;
          if(find(around_points.begin(), around_points.end(),
                  selected_point) != around_points.end())
            {
              const size_t face_idx = fa.get_face_idx(edge.first,edge.second,selected_point);
              if(face_idx == -1){
                  cerr << "# [error] strange can not find face. "
                       << edge.first << edge.second << selected_point << endl;
                  return __LINE__;
                }
              const pair<size_t,size_t> & tet_pair = fa.face2tet_[face_idx];
              oecit it = ortae.e2t_.find(edge);
              if(it ==  ortae.e2t_.end())
                it = ortae.e2t_.find(make_pair(edge.second, edge.first));

              modify_face_jump_type_at_given_tets_edge(tet_pair.first,tet_pair.second,
                                                       inner_face_jump_type,it->second);


              {// update the one_loop
                list<pair<size_t,size_t> >::iterator lit =
                    find(one_loop.begin(),one_loop.end(), edge); //one_loop.find(edge);
                if(lit == one_loop.end()){
                    cerr << "# [error] strange can not find this edge in loop." << endl;
                    return __LINE__;
                  }
                *lit = make_pair(edge.first,selected_point);
                ++lit;
                one_loop.insert(lit,make_pair(selected_point,edge.second));
              }
              break;
            }// end if
        }// end for
    }

  return 0;
}


int remove_singularity_degenerated_edge(
    const std::vector<std::deque<std::pair<size_t,size_t> > > &loops,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    const matrixst &tet,
    const matrixd &node,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::face2tet_adjacent &fa)
{
  for(size_t t = 0; t < loops.size(); ++t){
      const deque<pair<size_t,size_t> > & one_loop = loops[t];
      remove_one_inner_loop_by_modify_face_type(one_loop,inner_face_jump_type,fa,ortae,tet,node);
    }
  return 0;
}

int remove_union_degeneration_points(
    const std::vector<std::vector<size_t> > & degenerated_points,
    const matrixst & outside_face,
    const vector<deque<pair<size_t,size_t> > > & singularity_chain,
    const vector<deque<size_t> > & singularity_type,
    matrixst & outside_face_normal_align_type,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type)
{
  cerr << "# [error] function is not ready." << endl;


  return 0;
}


int remove_union_degeneration_edge(
    const std::vector<std::vector<size_t> > & degenerated_points_vec,
    const matrixst & outside_face,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_chain,
    const std::vector<std::deque<size_t> > & singularity_type,
    const matrixst & uvw_integer,
    const jtf::mesh::edge2cell_adjacent & ea,
    jtf::mesh::face2tet_adjacent & fa,
    std::vector<size_t> & new_tet,
    std::vector<double> & new_node,
    jtf::mesh::one_ring_tet_at_edge & ortae,
    matrixst &outside_face_normal_align_type,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type)
{
  cerr << "# [error] remove union degeneration edge function is not finished yet." << endl;

  {
    set<size_t> near_surface_singularity_chain;
    for(size_t t = 0; t < degenerated_points_vec.size(); ++t){
        const vector<size_t> & degenerated_point = degenerated_points_vec[t];
        vector<size_t> linked_near_surface_chains;
        if(is_near_surface_degenerated_case(degenerated_point,
                                            linked_near_surface_chains)){
            assert(!linked_near_surface_chains.empty());
            near_surface_singularity_chain.insert(linked_near_surface_chains.begin(),
                                                  linked_near_surface_chains.end());
          }
      }

    itr_matrix<size_t*> tet_mat(4,new_tet.size()/4,&new_tet[0]);
    itr_matrix<double*> node_mat(3,new_node.size()/3,&new_node[0]);

    map<pair<size_t,size_t>,double> edge_weight;
    for(size_t t = 0; t < ea.edges_.size(); ++t){
        const pair<size_t,size_t> &edge = ea.edges_[t];
        edge_weight[edge] = norm(node_mat(colon(),edge.first) -
                                 node_mat(colon(),edge.second));
      }
    unique_ptr<vertex_connection<UNDIRECT> > vc(
          vertex_connection<UNDIRECT>::create(edge_weight));

    for(set<size_t>::const_iterator scit = near_surface_singularity_chain.begin();
        scit != near_surface_singularity_chain.end(); ++scit){
        const deque<pair<size_t,size_t> > & one_chain =
            singularity_chain[*scit];
        modify_face_type_to_remove_near_surface(tet_mat,node_mat,outside_face,ea,
                                                fa,*vc,ortae,one_chain,inner_face_jump_type);
      }

  }

  return 0;
}

int extract_inner_face_type_from_zyz(
    const jtf::mesh::face2tet_adjacent & fa,
    const matrixd & zyz,
    boost::unordered_map<std::pair<size_t,size_t>, size_t> & inner_face_jump_type)
{
  matrix<matrixd > frame(zyz.size(2));
  for(size_t t = 0; t < zyz.size(2); ++t){
      frame[t].resize(3,3);
      zyz_angle_2_rotation_matrix1(&zyz(0,t), &frame[t][0]);
    }

  extract_inner_face_type_from_frame(fa, frame, inner_face_jump_type);
  return 0;
}


int extract_inner_face_type_from_frame(
    const jtf::mesh::face2tet_adjacent & fa,
    const zjucad::matrix::matrix<matrixd > & frame,
    boost::unordered_map<std::pair<size_t,size_t>, size_t> & inner_face_jump_type)
{
  matrixd rot = eye<double>(3);
  for(size_t t = 0; t < fa.face2tet_.size(); ++t){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[t];
      if(fa.is_outside_face(tet_pair)) continue;
      get_best_alignment(&frame[tet_pair.first][0], &frame[tet_pair.second][0],
          &rot[0]);
      if(is_trivial_type(type_transition1(rot))) continue;
      inner_face_jump_type[tet_pair] = type_transition1(rot);
      inner_face_jump_type[make_pair(tet_pair.second, tet_pair.first)]
          = type_transition1(trans(rot));
    }
  return 0;
}


int convert_surface_rotation_to_axis_type(const size_t &type)
{
  assert(type < 24);
  matrixd basis = eye<double>(3);
  matrixd direction = type_transition2(type)(colon(),0);
  vector<pair<double,size_t> > err;
  for(size_t ti = 0; ti < 6; ++ti){
      double sign = 1;
      if(ti % 2 == 1) sign = -1;
      err.push_back(make_pair(norm(direction - sign * basis(colon(), ti/2)), ti));
    }
  sort(err.begin(), err.end());
  return err.front().second;
}

int detect_topology_graph_valid(
    const vector<size_t> &rot_type,
    const matrixst & cut_tet,
    const matrixst & cut_outside_face,
    const matrixst & cut_outside_face_idx,
    const matrixst & cut_tet2tet,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixd & node,
    const matrixd & cut_node,
    const vector<pair<size_t,size_t> > & jump_face_vec,
    const matrixst &  face_pair_cut,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    singularity_graph & g,
    const vector<bool> & is_surface,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &tet_pair2rot_idx,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & face_pair_to_rot_type_idx,
    const boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
    bool no_surface)
{
  cerr << "# [info] start to test_graph." << endl;

  cerr << "# [info] total face number " << is_surface.size() << endl;

  size_t fi = 0;
  for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator bucit =
      face_pair_to_rot_type_idx.begin(); bucit != face_pair_to_rot_type_idx.end();
      ++bucit, ++fi){
      if(fi % 1000 == 0){
          cerr << "# [info] add faces number " << fi << endl;
        }
      g.insert_transition(bucit->first,cut_tet,cut_tet2tet,
                          fa_cut, ortae,rot_type,face_pair_to_rot_type_idx,
                          tet_pair2rot_idx,2); // not check unready_edge utill the last
    }
  if(!no_surface){
      cerr << "# [info] Surface face adding." << endl;
      for(boost::unordered_map<size_t,size_t>::const_iterator bucit =
          surface_idx_to_rot_idx.begin(); bucit != surface_idx_to_rot_idx.end();
          ++bucit, ++fi){
          if(fi % 1000 == 0){
              cerr << "# [info] add faces number " << fi << endl;
            }
          g.insert_surface_transition(
                rot_type[bucit->second], fa_cut, cut_tet2tet, bucit->first,
              rot_type,tet_pair2rot_idx );
        }
    }

  g.check_unready_edges(ortae, fa_cut, cut_tet, cut_tet2tet, rot_type,
                        tet_pair2rot_idx);

  if(g.is_valid_with_info(cut_tet2tet,rot_type, ortae,surface_idx_to_rot_idx)){
      cerr << "# [info] valid graph." << endl;
    }else
    cerr << "# [error] invalid graph." << endl;

  matrixd param;
  if(load_wyz_param_file("param", cut_tet.size(2),param)){
      cerr << "# [param type] original node" << endl;
      param = cut_node;
    }else{
      cerr << "# [param type] param node" << endl;
    }

  //#endif

#ifdef debug
  {
    vector<size_t> outside_face_vec;
    vector<size_t> outside_face_type_vec;

    vector<size_t> rot_fi;
    for(size_t fi = 0; fi < face_pair_cut.size(); ++fi){
        if(face_pair_cut[fi] == -1){

            outside_face_vec.insert(
                  outside_face_vec.end(),
                  cut_outside_face(colon(),fi).begin(),
                  cut_outside_face(colon(),fi).end());

            const auto it = surface_idx_to_rot_idx.find(cut_outside_face_idx[fi]);
            if(it == surface_idx_to_rot_idx.end()){
                throw std::logic_error("can not find surface_idx in surface_idx_to_rot_idx");
              }
            outside_face_type_vec.push_back(rot_type[it->second]);
            //cerr << "rot_idx " << surface_idx_to_rot_idx[cut_outside_face_idx[fi]] << endl;
          }
      }

    ofstream ofs("surface_patch.vtk");

    tri2vtk(ofs, &param[0], param.size(2),&outside_face_vec[0],
        outside_face_vec.size()/3);
    cell_data(ofs, &outside_face_type_vec[0], outside_face_type_vec.size(),
        "patch");
  }
#endif

  //#ifdef debug
  //  {
  //    ofstream ofs("jump_faces.vtk");
  //    tri2vtk(ofs, &param[0], param.size(2),&jump_faces[0], jump_faces.size()/3);
  //  }
  //#endif

#ifdef debug
  {
    //    ofstream ofs("face_pair_jump_transition");
    //    for(size_t fi = 0; fi < jump_face_vec.size(); ++fi){
    //      const pair<size_t,size_t> & face_pair = jump_face_vec[fi];
    //      const vector<size_t> & f0 = fa_cut.faces_[face_pair.first];
    //      const vector<size_t> & f1 = fa_cut.faces_[face_pair.second];
    //      ofs << "# face pair " << endl;
    //      ofs << "# ------" << f0[0] << " " << f0[1] << " " << f0[2] << endl;
    //      ofs << "# ------" << f1[0] << " " << f1[1] << " " << f1[2] << endl;
    //      ofs << "# ------ rot type " << rot_type[fi] << endl;
    //      ofs << "# ------ gnode u,v,w " ;
    //      for(size_t t = 0; t < 3; ++t){
    //        if(g.gnode_flag_[fi * 3 + t])
    //          ofs << g.gnode_[fi * 3 +t] << " ";
    //        else
    //          ofs << "unknown ";
    //      }
    //      ofs << endl;
    //    }

    //    //ofstream ofs_compund("compound_edge.vtk");
    //    dump_singularity_to_vtk("compound_edge.vtk", param, compound_edges);
    //    //dump_singularity_to_vtk("loop_edge.vtk", node, loop_edges);

    vector<size_t> loop_points;
    vector<size_t> loop_points_idx;
    ofstream ofs_loop_file("loop_points");
    for(size_t ci = 0; ci < g.loop_edges_.size(); ++ci){
        const deque<pair<size_t,size_t> > & chain = g.loop_edges_[ci];
        boost::unordered_set<size_t> point_set;
        for(size_t ei = 0; ei < chain.size(); ++ei){
            point_set.insert(chain[ei].first);
            point_set.insert(chain[ei].second);
          }
        loop_points.insert(loop_points.end(), point_set.begin(), point_set.end());

        {
          ofs_loop_file << "loops " << ci << " " << point_set.size() << endl;
          for(boost::unordered_set<size_t>::const_iterator bucit = point_set.begin();
              bucit != point_set.end(); ++bucit){
              ofs_loop_file << *bucit << " ";
            }
          ofs_loop_file << endl;
        }

        for(size_t i = 0; i < point_set.size(); ++i){
            loop_points_idx.push_back(ci);
          }
      }

    vector<size_t> loop_edge;
    vector<size_t> edge_type;
    for(size_t ci = 0; ci < g.loop_edges_.size(); ++ci){
        const deque<pair<size_t,size_t> > & chain = g.loop_edges_[ci];
        for(size_t i = 0; i < chain.size(); ++i){
            loop_edge.push_back(chain[i].first);
            loop_edge.push_back(chain[i].second);
            edge_type.push_back(i);
          }
      }

    //    {// dump out degenerated edgs
    //      ofstream ofs("degenerated_edges");
    //      ofs << loop_edge.size()/2 << endl;
    //      for(size_t ei = 0; ei < loop_edge.size()/2; ++ei){
    //        ofs << loop_edge[ei * 2 + 0] << " " << loop_edge[ei * 2 + 1] << endl;
    //      }

    //    }
    ofstream ofs_loop_edge("loop_edge.vtk");
    line2vtk(ofs_loop_edge, &param[0], param.size(2), &loop_edge[0], loop_edge.size()/2);
    cell_data(ofs_loop_edge, &edge_type[0], edge_type.size(), "idx");

    ofstream ofs_loop("loop_points.vtk");
    point2vtk(ofs_loop, &param[0], param.size(2), &loop_points[0], loop_points.size());
    cell_data(ofs_loop, &loop_points_idx[0], loop_points_idx.size(), "idx");
  }
#endif

  vector<pair<size_t,size_t> > g_unknonw_jump_face_pair;
  {
    // draw the gnode which is not 0 or unknown
    set<pair<size_t,size_t> > g_unknonw_jump_face_pair_set;
    boost::unordered_map<size_t,size_t> g_not_zero_faces;
    for(size_t gi = 0; gi < g.gnode_.size(); ++gi){
        if(!g.gnode_flag_[gi]){
            const pair<size_t,size_t> face_pair = g.get_jump_face_from_gnode_idx(gi);
            g_not_zero_faces.insert(make_pair(face_pair.first,0));
            g_not_zero_faces.insert(make_pair(face_pair.second,0));
            g_unknonw_jump_face_pair_set.insert(face_pair);
          }else{
            if(fabs(g.gnode_[gi]) > 1e-6){
                const pair<size_t,size_t> face_pair = g.get_jump_face_from_gnode_idx(gi);
                g_not_zero_faces.insert(make_pair(face_pair.first,1));
                g_not_zero_faces.insert(make_pair(face_pair.second,1));
              }
          }
      }
    g_unknonw_jump_face_pair.resize(g_unknonw_jump_face_pair_set.size());
    copy(g_unknonw_jump_face_pair_set.begin(), g_unknonw_jump_face_pair_set.end(),
         g_unknonw_jump_face_pair.begin());
    vector<size_t> g_not_zeros_face_vec;
    vector<size_t> g_not_zeros_face_vec_type;
    for(boost::unordered_map<size_t,size_t>::const_iterator bucit = g_not_zero_faces.begin();
        bucit != g_not_zero_faces.end(); ++bucit){
        const vector<size_t> & face_vec = fa_cut.faces_[bucit->first];
        g_not_zeros_face_vec.insert(g_not_zeros_face_vec.end(),
                                    face_vec.begin(), face_vec.end());
        g_not_zeros_face_vec_type.push_back(bucit->second);
      }
    ofstream ofs("g_not_zeros_faces.vtk");
    tri2vtk(ofs, &cut_node[0], cut_node.size(2),
        &g_not_zeros_face_vec[0], g_not_zeros_face_vec.size()/3);
    cell_data(ofs, &g_not_zeros_face_vec_type[0],
        g_not_zeros_face_vec_type.size(), "0_unknown_1_not_zero");

    ofstream ofs_g_unknown("g_unknown_jump_face");
    ofs_g_unknown << g_unknonw_jump_face_pair.size() << endl;
    for(size_t fi = 0; fi < g_unknonw_jump_face_pair.size(); ++fi){
        ofs_g_unknown << g_unknonw_jump_face_pair[fi].first << " "
                      << g_unknonw_jump_face_pair[fi].second << endl;
      }
  }

  //#ifdef debug
  boost::unordered_set<size_t> restricted_nodes;
  {
    vector<deque<pair<size_t,size_t> > > restrict_edges(1);
    vector<deque<size_t> > restrict_edges_type(1);
    vector<pair<size_t,size_t> > restrict_edges_cut;
    vector<size_t> restrict_edges_cut_type;

    unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
          jtf::mesh::edge2cell_adjacent::create(cut_outside_face));
    for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
        const pair<size_t,size_t> & one_edge_cut = ea->edges_[ei];
        vector<size_t> grouped_axes;
        for(size_t ai = 0; ai < 3; ++ai){
            pair<size_t,size_t> node_edge(
                  g.get_node_idx_from_point_idx(one_edge_cut.first, ai),
                  g.get_node_idx_from_point_idx(one_edge_cut.second, ai));
            if(g.fnode_[node_edge.first] == g.fnode_[node_edge.second])
              grouped_axes.push_back(ai);
          }
        if(grouped_axes.size() < 2) continue;
        if(grouped_axes.size() == 2){
            const size_t axis_need_to_build_path =
                (0+1+2) - (std::accumulate(grouped_axes.begin(), grouped_axes.end(),0));
            restrict_edges[0].push_back(one_edge_cut);
            restrict_edges_type[0].push_back(axis_need_to_build_path*3);

            restrict_edges_cut.push_back(one_edge_cut);
            restrict_edges_cut_type.push_back(axis_need_to_build_path);

            for(size_t ai = 0; ai < 3; ++ai){
                if(ai != axis_need_to_build_path){
                    // attention should add all fnode belong to this group
                    const size_t fnode_idx_0 =
                        g.get_node_idx_from_point_idx(one_edge_cut.first, ai);
                    const size_t fnode_idx_1 =
                        g.get_node_idx_from_point_idx(one_edge_cut.second, ai);

                    restricted_nodes.insert(g.groups_.at(g.fnode_.at(fnode_idx_0)).begin(),
                                            g.groups_.at(g.fnode_.at(fnode_idx_0)).end());
                    restricted_nodes.insert(g.groups_.at(g.fnode_.at(fnode_idx_1)).begin(),
                                            g.groups_.at(g.fnode_.at(fnode_idx_1)).end());
                  }
              }
          }else{
            restrict_edges_cut.push_back(one_edge_cut);
            restrict_edges_cut_type.push_back(-1);

            restrict_edges[0].push_back(one_edge_cut);
            restrict_edges_type[0].push_back(9);
          }
      }

    dump_restricted_edge("restricted_edges", restrict_edges_cut,
                         restrict_edges_cut_type);

    dump_singularity_chain_to_vtk_2("restricted_edge.vtk", param,
                                    restrict_edges,restrict_edges_type);
  }

  // to gather restrict nodes according to the surface type
  {
    if(!no_surface){
        for(size_t fi = 0; fi < cut_outside_face.size(2); ++fi){
            const size_t & surface_idx = cut_outside_face_idx[fi];

            boost::unordered_map<size_t,size_t>::const_iterator bcit =
                surface_idx_to_rot_idx.find(surface_idx);
            if(bcit == surface_idx_to_rot_idx.end())  continue;
            assert(rot_type.size() > bcit->second);
            assert(rot_type[bcit->second] < 3);
            for(size_t pi = 0; pi < 3; ++pi){
                // attention should add all fnode belong to this group
                const size_t fnode_idx =
                    g.get_node_idx_from_point_idx(cut_outside_face(pi,fi),
                                                  rot_type[bcit->second]);

                restricted_nodes.insert(g.groups_.at(g.fnode_.at(fnode_idx)).begin(),
                                        g.groups_.at(g.fnode_.at(fnode_idx)).end());
                //        restricted_nodes.insert(
                //              g.get_node_idx_from_point_idx(cut_outside_face(pi,fi), rot_type[bcit->second]));
              }
          }
      }
  }

  dump_out_param_config(
        "group_file", "equation_file", "chain_file","cut_tet.tet",
        cut_tet,cut_node,cut_tet2tet,rot_type, fa_cut,g,
        restricted_nodes,face_pair_to_rot_type_idx,g_unknonw_jump_face_pair);

#ifdef debug
  {
    ofstream ofs("fnode_param");
    for(size_t i = 0; i < g.fnode_.size()/3; ++i){
        const pair<size_t,size_t> fnode_idx = g.get_point_idx_from_node_idx(i * 3);
        ofs << fnode_idx.first << " " << g.fnode_[i * 3 + 0] << " "
            << g.fnode_[i * 3 + 1] << " " << g.fnode_[i * 3 + 2] << endl;
      }
  }
#endif
#ifdef debug
  {
    vector<deque<pair<size_t,size_t> > > chain_vec = g.chains_;
    for(size_t ci = 0; ci < chain_vec.size(); ++ci){
        for(size_t ei = 0; ei < chain_vec[ci].size(); ++ei){
            chain_vec[ci][ei].first
                = g.get_point_idx_from_node_idx(chain_vec[ci][ei].first%g.fnode_.size()).first;
            chain_vec[ci][ei].second
                = g.get_point_idx_from_node_idx(chain_vec[ci][ei].second%g.fnode_.size()).first;
          }
      }
    dump_singularity_to_vtk("chains.vtk", param, chain_vec);
  }
#endif

  return 0;
}

static bool is_fnode_same(const std::vector<std::stack<size_t> > & fnode0,
                          const std::vector<std::stack<size_t> > & fnode1)
{
  if(fnode0.size() != fnode1.size()){
      cerr << "# [info] different fnode size " << fnode0.size() << ","
           << fnode1.size() << endl;
      return false;
    }
  for(size_t fi = 0; fi < fnode0.size(); ++fi){
      if(fnode0[fi].size() != fnode1[fi].size()){
          cerr << "# [info] different group size of " << fi << endl;
        }
    }
  return true;
}

int aggressive_assign_transition_with_type(
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_type,
    boost::unordered_map<size_t,size_t> & surface_type,
    const bool no_surface,
    const zjucad::matrix::matrix<matrixd> * frame_ptr)
{
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(orig_tet));
  unique_ptr<jtf::mesh::face2tet_adjacent> fa_cut(jtf::mesh::face2tet_adjacent::create(cut_tet));
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;

  if(!fa.get() || !fa_cut.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent" << endl;
      return __LINE__;
    }

  matrixst face_pair_cut;
  analysis_transition_raw(orig_tet, node, *fa, *fa_cut, cut_tet,face_pair_cut);


  searching_strategy::orig_check_graph(
        orig_tet, cut_tet, node, *fa, *fa_cut,  face_pair_cut,
        surface_type, face_type, no_surface, frame_ptr);
  //   searching_strategy::searching_solutions_by_using_disk(
  //   orig_tet, cut_tet, node, *fa, *fa_cut,  face_pair_cut, surface_type, face_type, frame_ptr);

  //  searching_strategy::searching_solutions_by_using_memory_binary_search_group(
  //        orig_tet, cut_tet, node, *fa, *fa_cut,  face_pair_cut,
  //        surface_type, face_type, frame_ptr);

  //  searching_strategy::searching_solutions_by_using_memory_binary_search_group(
  //        orig_tet, cut_tet, node, *fa, *fa_cut,  face_pair_cut,
  //        surface_type, face_type, frame_ptr);

  //  searching_strategy::searching_solutions_by_using_memory(
  //        orig_tet, cut_tet, node, *fa, *fa_cut,  face_pair_cut, surface_type, face_type, frame_ptr);

  //  searching_strategy::searching_solutions_by_minimal_error_search_group(
  //        *fa, *fa_cut, orig_tet, cut_tet, node, face_pair_cut, *frame_ptr, face_type);

  return 0;
}

int gradually_fix_reliable_constraints(
    const matrixst & cut_tet,
    matrixd &node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    boost::unordered_map<size_t,size_t> & surface_type)
{
  cerr << "# [error] gradually fix reliable constraints is not finished." << endl;
  return __LINE__;
}

int find_patch_for_each_near_miss(
    const boost::unordered_set<pair<size_t,size_t> > & all_singularity_edges,
    const deque<pair<size_t,size_t> > & one_chain,
    const boost::unordered_set<size_t> & g_unknow_face_set,
    const jtf::mesh::one_ring_face_at_edge & orfae,
    const jtf::mesh::edge2cell_adjacent & ea,
    const jtf::mesh::face2tet_adjacent & fa,
    boost::unordered_set<size_t> & patch)
{
  list<pair<size_t,size_t> > one_chain_list;
  for(size_t ei = 0; ei < one_chain.size(); ++ei){
      if(one_chain[ei].first > one_chain[ei].second)
        one_chain_list.push_back(
              make_pair(one_chain[ei].second,one_chain[ei].first));
      else
        one_chain_list.push_back(one_chain[ei]);
    }

  patch.clear();
  pair<size_t,size_t> first_multi_face_edge(-1,-1);
  while(!one_chain_list.empty()){
      list<pair<size_t,size_t> >::iterator lit = one_chain_list.begin();
      const pair<size_t,size_t> & edge = *lit;

      { // check whether this edge is surface edge
        const size_t edge_idx = ea.get_edge_idx(edge.first, edge.second);
        if(edge_idx != -1){
            one_chain_list.pop_front();
            continue;
          }
      }

      jtf::mesh::one_ring_face_at_edge::edge2_face_type::const_iterator
          ecit = orfae.e2f_.find(edge);
      if(ecit == orfae.e2f_.end())
        ecit = orfae.e2f_.find(make_pair(edge.second, edge.first));
      if(ecit == orfae.e2f_.end()){
          cerr << "# [error] strange can not find faces arounding edge: "
               << edge.first << " " << edge.second << endl;
          return __LINE__;
        }
      const vector<size_t> & face_around = ecit->second;
      vector<size_t> possible_faces;
      for(size_t i = 0; i < face_around.size(); ++i){
          boost::unordered_set<size_t>::const_iterator buscit =
              g_unknow_face_set.find(face_around[i]);
          if(buscit != g_unknow_face_set.end() &&
             (patch.find(face_around[i]) == patch.end()))
            possible_faces.push_back(face_around[i]);
        }

      size_t face_idx = -1;
      if(possible_faces.size() > 1){
          if(first_multi_face_edge == edge)
            face_idx = possible_faces.front();
          else{
              if(first_multi_face_edge.first == -1 &&
                 first_multi_face_edge.second == -1)
                first_multi_face_edge = edge;
              one_chain_list.push_back(edge);
              one_chain_list.pop_front();
              continue;
            }
        }else{
          first_multi_face_edge = make_pair(-1,-1);
          if(possible_faces.size() == 1){
              face_idx = possible_faces.front();
            }else if(possible_faces.empty()){
              cerr << "# [error] strange can not find possible faces." << endl;
              return __LINE__;
            }
        }

      patch.insert(face_idx);
      {
        const vector<size_t> & face_vec = fa.faces_[face_idx];
        assert(find(face_vec.begin(), face_vec.end(), edge.first) !=
            face_vec.end());
        assert(find(face_vec.begin(), face_vec.end(), edge.second) !=
            face_vec.end());
        const size_t other_point =
            std::accumulate(face_vec.begin(), face_vec.end(), static_cast<size_t>(0))
            - edge.first - edge.second;
        vector<pair<size_t,size_t> > candidate_edges;
        candidate_edges.push_back(make_pair(other_point, edge.first));
        candidate_edges.push_back(make_pair(other_point, edge.second));
        if(one_chain_list.front().first == first_multi_face_edge.first &&
           one_chain_list.front().second == first_multi_face_edge.second){
            first_multi_face_edge = make_pair(-1,-1);
          }
        one_chain_list.pop_front();
        for(size_t ei = 0; ei < candidate_edges.size(); ++ei){
            if(candidate_edges[ei].first > candidate_edges[ei].second){
                swap(candidate_edges[ei].first, candidate_edges[ei].second);
              }
            // check if this edge is already in this chain list, delete it
            // or, should check whether this edge is a singularity edge,
            // if yes, can not include this one to the chain
            // if no, add it to this chain list
            list<pair<size_t,size_t> >::iterator lit =
                find(one_chain_list.begin(), one_chain_list.end(), candidate_edges[ei]);
            if(lit != one_chain_list.end()){
                one_chain_list.erase(lit);
                continue;
              }else{
                boost::unordered_set<pair<size_t,size_t> >::const_iterator
                    buscit = all_singularity_edges.find(candidate_edges[ei]);
                if(buscit != all_singularity_edges.end()){
                    continue;
                  }else{
                    // if reach surface, then stop
                    const size_t edge_idx = ea.get_edge_idx(candidate_edges[ei].first,
                                                            candidate_edges[ei].second);
                    if(edge_idx != -1) continue;
                    one_chain_list.push_front(candidate_edges[ei]);
                  }
              }
          }
      }
    }

  return 0;
}

int remove_near_miss_by_minimal_cut_new(
    const matrixst &cut_tet,
    const matrixst &tet,
    const matrixd & cut_node,
    const matrixd & node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    boost::unordered_map<size_t,size_t> &surface_type,
    const vector<pair<size_t,size_t> > & g_unknown_face_pair,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const vector<size_t> & fnode)
{
  matrixst cut_tet2tet(max(cut_tet)+1);
  cut_tet2tet(cut_tet) = tet(colon());

  jtf::mesh::one_ring_tet_at_edge ortae;
  ortae.add_tets(tet, fa);
  ortae.sort_into_loop(tet, node);
  matrixst outside_face, outside_face_cut, outside_face_cut_idx;

  get_outside_face(fa, outside_face);
  get_outside_face(fa_cut, outside_face_cut);
  get_outside_face_idx(fa_cut, outside_face_cut_idx);

  boost::unordered_set<pair<size_t,size_t> > singularity_edges_set;
  vector<deque<pair<size_t,size_t> > > singularity_chains;
  vector<deque<size_t> > singularity_type;

  jtf::mesh::meshes tm_mesh;
  tm_mesh.mesh_ = tet;
  tm_mesh.node_ = node;
  jtf::tet_mesh tm(tm_mesh);
  singularity_extractor se(tm);

  {
    se.extract(inner_face_jump_type, singularity_chains,singularity_type);

    for(size_t ci = 0; ci < singularity_chains.size(); ++ci){
        const deque<pair<size_t,size_t> > & one_chain = singularity_chains[ci];
        for(size_t ei = 0; ei < one_chain.size(); ++ei){
            const pair<size_t,size_t> & one_edge = one_chain[ei];
            if(one_edge.first > one_edge.second)
              singularity_edges_set.insert(make_pair(one_edge.second, one_edge.first));
            else
              singularity_edges_set.insert(one_edge);
          }
      }

    { // visual
      dump_singularity_chain_to_vtk_2("before_remove_near_miss.vtk", node,
                                      singularity_chains,  singularity_type);
    }
  }

  // tuple<size_t,size_t,size_t>: p0,p1,type(0:u, 1:v, 2:w, 3:compound, -1:free)
  boost::unordered_map<pair<size_t,size_t>,
      vector<std::tuple<size_t,size_t,size_t> > > orig_edges2cut_edges;
  {
    jtf::mesh::one_ring_tet_at_edge ortae_cut;
    ortae_cut.add_tets(cut_tet, fa_cut);
    vector<size_t> fixed_axis;
    size_t type = -1;
    for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator ecit
        = ortae_cut.e2t_.begin(); ecit != ortae_cut.e2t_.end(); ++ecit){
        const pair<size_t,size_t> & cut_edge = ecit->first;

        fixed_axis.clear();
        for(size_t i = 0; i < 3; ++i){
            assert(cut_edge.first * 3 + i < fnode.size());
            assert(cut_edge.second * 3 + i < fnode.size());

            if(fnode[cut_edge.first * 3 + i] ==
               fnode[cut_edge.second * 3 + i])
              fixed_axis.push_back(i);
          }
        if(fixed_axis.empty()) type = -1; // free
        else if(fixed_axis.size() == 3) type= 3; // compound
        else if(fixed_axis.size() == 2) type = 3 - fixed_axis[0] - fixed_axis[1];
        else if(fixed_axis.size() == 1) type = -1;

        pair<size_t,size_t> orig_edge(cut_tet2tet[cut_edge.first],
            cut_tet2tet[cut_edge.second]);

        //      if(orig_edge.first == 1857 || orig_edge.second == 1857){
        //        if(orig_edge.first == 10081 || orig_edge.second == 10081){
        //          cerr << fnode[6261] << " " << fnode[3228] << endl;
        //          cerr << endl;
        //        }
        //      }

        if(orig_edge.first > orig_edge.second){
            orig_edges2cut_edges[
                make_pair(orig_edge.second,orig_edge.first)].push_back(
                  std::make_tuple(cut_edge.second, cut_edge.first, type));
          }else
          orig_edges2cut_edges[orig_edge].push_back(
                std::make_tuple(cut_edge.first,
                                cut_edge.second, type));
      }
  }

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea_cut(
        jtf::mesh::edge2cell_adjacent::create(outside_face_cut));
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));

  if(!ea_cut.get() || !ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }

  jtf::mesh::one_ring_face_at_point orfap;
  orfap.add_all_faces(outside_face_cut, *ea_cut);
  orfap.sort_int_loop(outside_face_cut, cut_node);

  jtf::mesh::one_ring_face_at_edge orfae;
  orfae.build(ortae, fa, tet);

  boost::unordered_set<size_t> g_unknown_face_in_original_set;
  for(size_t fi = 0; fi < g_unknown_face_pair.size(); ++fi){
      const vector<size_t> & face_vec = fa_cut.faces_[g_unknown_face_pair[fi].first];
      const size_t face_original =
          fa.get_face_idx(cut_tet2tet[face_vec[0]],cut_tet2tet[face_vec[1]],
          cut_tet2tet[face_vec[2]]);
      if(face_original == -1){
          cerr << "# [error] can not find face " << g_unknown_face_pair[fi].first
               << " in original tetmesh." << endl;
          return __LINE__;
        }
      g_unknown_face_in_original_set.insert(face_original);
    }

  jtf::tetmesh::check_manifold(g_unknown_face_in_original_set.begin(),
                               g_unknown_face_in_original_set.end(), fa);


  // adjust the order by singularity chain size, the longer the fronter
  vector<pair<size_t,size_t> > singularity_chain_with_size_order;
  deque<size_t> singularity_chain_queue;
  {
    for(size_t ci = 0; ci < singularity_chains.size(); ++ci)
      singularity_chain_with_size_order.push_back(
            make_pair(singularity_chains[ci].size(), ci));

    sort(singularity_chain_with_size_order.begin(),
         singularity_chain_with_size_order.end());

    for(size_t ci = singularity_chain_with_size_order.size() -1; ci != -1; --ci)
      singularity_chain_queue.push_back(singularity_chain_with_size_order[ci].second);
  }

  // this
  size_t queue_size = singularity_chain_queue.size();
  size_t queue_size_keeping_time = 0;
  while(!singularity_chain_queue.empty()){
      const deque<pair<size_t,size_t> > & one_chain =
          singularity_chains[singularity_chain_queue.front()];
      if(is_near_miss_chain(one_chain, orig_edges2cut_edges, surface_type,
                            orfap, outside_face_cut_idx, outside_face)){

          boost::unordered_set<size_t> patch, patch_temp;
          find_patch_for_each_near_miss(
                singularity_edges_set, one_chain, g_unknown_face_in_original_set,
                orfae, *ea, fa, patch);

          for(boost::unordered_set<size_t>::const_iterator cit = patch_temp.begin();
              cit != patch_temp.end(); ++cit){
              const size_t &face_idx = *cit;
              boost::unordered_set<size_t>::iterator sit =
                  g_unknown_face_in_original_set.find(face_idx);
              if(sit != g_unknown_face_in_original_set.end())
                g_unknown_face_in_original_set.erase(sit);
            }

          patch_temp = patch;
          if(remove_near_miss_chain_by_modifying_g_unknown_faces_new(
               one_chain, fa, ortae, orfae, *ea, singularity_edges_set,
               g_unknown_face_in_original_set,
               patch, inner_face_jump_type, outside_face, node)){
              singularity_chain_queue.push_back(singularity_chain_queue.front());
              singularity_chain_queue.pop_front();
              g_unknown_face_in_original_set.insert(patch_temp.begin(), patch_temp.end());
            }else{
              g_unknown_face_in_original_set.insert(patch.begin(), patch.end());
              singularity_chain_queue.pop_front();

              // if this chain has one_point inside the tet
              // we should reextract all singularity chains because the near_miss
              // removing will change the topology of singularity graph
              if(find(outside_face.begin(), outside_face.end(), one_chain.front().first)
                 == outside_face.end() ||
                 find(outside_face.begin(), outside_face.end(), one_chain.back().second)
                 == outside_face.end()){
                  se.extract(inner_face_jump_type, singularity_chains, singularity_type);

                  singularity_chain_queue.clear();
                  singularity_chain_with_size_order.clear();
                  for(size_t ci = 0; ci < singularity_chains.size(); ++ci)
                    singularity_chain_with_size_order.push_back(
                          make_pair(singularity_chains[ci].size(), ci));

                  sort(singularity_chain_with_size_order.begin(),
                       singularity_chain_with_size_order.end());

                  for(size_t ci = singularity_chain_with_size_order.size() -1;
                      ci != -1; --ci)
                    singularity_chain_queue.push_back(
                          singularity_chain_with_size_order[ci].second);
                  continue;
                }
            }
        }else
        singularity_chain_queue.pop_front();

      if(singularity_chain_queue.size() == queue_size)
        ++queue_size_keeping_time ;
      else{
          queue_size_keeping_time = 0;
          queue_size = singularity_chain_queue.size();
        }
      if(queue_size_keeping_time == queue_size)
        break;
    }

  { // visual
    ostringstream os;
    os << num++;
    vector<deque<pair<size_t,size_t> > > singularity_chain;
    vector<deque<size_t> > singularity_type;
    se.extract(inner_face_jump_type,singularity_chain, singularity_type);

    string file_name = os.str() + "_singularity.vtk";
    dump_singularity_chain_to_vtk_2(file_name.c_str(), node,
                                    singularity_chain,  singularity_type);

    vector<size_t> left_g_unknown_faces;
    for(boost::unordered_set<size_t>::const_iterator buscit =
        g_unknown_face_in_original_set.begin(); buscit != g_unknown_face_in_original_set.end();
        ++buscit){
        const vector<size_t> & face_vec = fa.faces_[*buscit];
        left_g_unknown_faces.insert(left_g_unknown_faces.end(), face_vec.begin(), face_vec.end());
      }
    ofstream ofs((os.str()+"_left_g_unknown_face.vtk").c_str());
    tri2vtk(ofs, &node[0], node.size(2), &left_g_unknown_faces[0], left_g_unknown_faces.size()/3);
  }
  return 0;
}

int remove_near_miss_chain_by_modifying_g_unknown_faces_new(
    const deque<pair<size_t,size_t> > & one_chain,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::one_ring_face_at_edge & orfae,
    const jtf::mesh::edge2cell_adjacent & ea, // surface edge
    const boost::unordered_set<pair<size_t,size_t> > & all_singularity_edges,
    const boost::unordered_set<size_t> & g_unknown_face_in_original_without_patch,
    boost::unordered_set<size_t> & g_unknown_face_in_patch,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type,
    const matrixst & outside_face,
    const matrixd & node)
{
  boost::unordered_map<pair<size_t,size_t>,size_t> inner_face_jump_type_temp =
      inner_face_jump_type;

  list<pair<size_t,size_t> > one_chain_list;
  for(size_t ei = 0; ei < one_chain.size(); ++ei){
      const pair<size_t,size_t> & edge = one_chain[ei];
      if(edge.first > edge.second)
        one_chain_list.push_back(make_pair(edge.second, edge.first));
      else
        one_chain_list.push_back(edge);
    }
  pair<size_t,size_t> first_edge(-1,-1); // to avoid dead cycle

  matrixd rot = eye<double>(3);
  matrixd rot_begin = eye<double>(3);
  while(!one_chain_list.empty()){
      list<pair<size_t,size_t> >::iterator lit = one_chain_list.begin();

      // check if surface edge
      {
        const size_t edge_idx = ea.get_edge_idx(lit->first, lit->second);
        if(edge_idx != -1){
            one_chain_list.erase(lit);
            continue;
          }
      }

#if 0
      {
        ostringstream os;
        os << num++;
        vector<deque<pair<size_t,size_t> > > singularity_chain;
        vector<deque<size_t> > singularity_type;
        find_singularities_use_face_type(ortae,inner_face_jump_type, outside_face,
                                         singularity_chain, singularity_type);
        string file_name = os.str() + "_singularity.vtk";
        dump_singularity_chain_to_vtk_2(file_name.c_str(), node,
                                        singularity_chain,  singularity_type);

        vector<size_t> left_g_unknown_faces;
        for(boost::unordered_set<size_t>::const_iterator buscit =
            g_unknown_face_in_patch.begin(); buscit != g_unknown_face_in_patch.end();
            ++buscit){
            const vector<size_t> & face_vec = fa.faces_[*buscit];
            left_g_unknown_faces.insert(left_g_unknown_faces.end(), face_vec.begin(), face_vec.end());
          }

        ofstream ofs((os.str()+"_left_g_unknown_face.vtk").c_str());
        tri2vtk(ofs, &node[0], node.size(2), &left_g_unknown_faces[0], left_g_unknown_faces.size()/3);
      }
#endif

      jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator ecit
          = ortae.e2t_.find(*lit);
      if(ecit == ortae.e2t_.end())
        ecit = ortae.e2t_.find(make_pair(lit->second, lit->first));
      if(ecit == ortae.e2t_.end()){
          cerr << "# [error] can not find edge " << lit->first << " "
               << lit->second << " in one_ring_tet_at_edge " << endl;
          return __LINE__;
        }

      jtf::mesh::one_ring_face_at_edge::edge2_face_type::const_iterator
          cit = orfae.e2f_.find(*lit);
      if(cit == orfae.e2f_.end())
        cit = orfae.e2f_.find(make_pair(lit->second, lit->first));
      if(cit == orfae.e2f_.end()){
          cerr << "# [error] can not find edge " << lit->first << " "
               << lit->second << " in one_ring_face_at_edge " << endl;
          return __LINE__;
        }
      const vector<size_t> &around_faces = cit->second;
      size_t face_idx = -1;
      vector<size_t> face_idx_candidate;
      for(size_t fi = 0; fi < around_faces.size(); ++fi){
          boost::unordered_set<size_t>::iterator buit =
              g_unknown_face_in_patch.find(around_faces[fi]);
          if(buit != g_unknown_face_in_patch.end()){
              face_idx_candidate.push_back(around_faces[fi]);
            }
        }
      // here choose crease first, if all singularity edges are not crease
      // pick the first one
      if(face_idx_candidate.size() > 1) {
          if(one_chain_list.front() == first_edge){
              face_idx = face_idx_candidate.front();
              boost::unordered_set<size_t>::iterator buit =
                  g_unknown_face_in_patch.find(face_idx);
              assert(buit != g_unknown_face_in_patch.end());
              g_unknown_face_in_patch.erase(buit);
              first_edge = make_pair(-1,-1);
            }else{
              if(first_edge.first == -1 && first_edge.second == -1)
                first_edge = one_chain_list.front();
              one_chain_list.push_back(one_chain_list.front());
              one_chain_list.pop_front();
              continue;
            }
        }else{
          if(one_chain_list.front() == first_edge)
            first_edge = make_pair(-1,-1);
          if(face_idx_candidate.empty())
            face_idx = -1;
          else{
              face_idx = face_idx_candidate.front();
              boost::unordered_set<size_t>::iterator buit =
                  g_unknown_face_in_patch.find(face_idx);
              assert(buit != g_unknown_face_in_patch.end());
              g_unknown_face_in_patch.erase(buit);
            }
        }
      if(face_idx == -1){
          // can not find other face to go, then such edge must be identity or has
          // other candidate face in other patches. // check
          const vector<size_t> & tet_loop = ecit->second;
          rot = eye<double>(3);
          for(size_t ti = 0; ti < tet_loop.size() -1; ++ti) {
              boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit
                  = inner_face_jump_type.find(make_pair(tet_loop[ti], tet_loop[ti + 1]));
              if(cit == inner_face_jump_type.end()) continue;
              rot = temp(rot * type_transition2(cit->second));
            }
          if(type_transition1(rot) != TRIVIAL_TYPE){
              bool have_other_candidate_face = false;
              for(size_t fi = 0; fi < around_faces.size(); ++fi){
                  boost::unordered_set<size_t>::iterator buit =
                      g_unknown_face_in_original_without_patch.find(around_faces[fi]);
                  if(buit != g_unknown_face_in_original_without_patch.end()){
                      one_chain_list.erase(lit);
                      have_other_candidate_face = true;
                      break;
                    }
                }
              if(have_other_candidate_face) continue;

              cerr << "# [error] edge " << lit->first << " " << lit->second
                   << " is not trivial type while having no face to go." << endl;
              inner_face_jump_type = inner_face_jump_type_temp;

              return __LINE__;
            }else{
              one_chain_list.erase(lit);
              continue;
            }
        }

      // find the face which is in g_unkown_face_set, and need to find a suitable
      // rotation to ensure that this edge type is identity.
      {
        const pair<size_t,size_t> & tet_pair = fa.face2tet_[face_idx];
        assert(!fa.is_outside_face(tet_pair));

        vector<size_t> tet_loop = ecit->second;

        { // rotate loop at tet_pair
          tet_loop.pop_back();
          assert(find(tet_loop.begin(), tet_loop.end(), tet_pair.first) !=
              tet_loop.end());
          assert(find(tet_loop.begin(), tet_loop.end(), tet_pair.second) !=
              tet_loop.end());

          size_t ti = 0;
          for(; ti < tet_loop.size(); ++ti){
              if(tet_loop[ti] == tet_pair.first || tet_loop[ti] == tet_pair.second)
                break;
            }
          if(ti == 0 && (tet_loop.back() == tet_pair.first
                         ||tet_loop.back() == tet_pair.second)){
              tet_loop.insert(tet_loop.begin(), tet_loop.back());
              //tet_loop.pop_back();
            }else{
              std::rotate(tet_loop.begin(), tet_loop.begin() + ti, tet_loop.end());
              tet_loop.push_back(tet_loop.front());
            }
        }

        {
          // T = R01 * R12 * R23 * R30
          // T^-1 * T = (T^-1 * R01) * R12 * R23 * R30
          rot = eye<double>(3);
          rot_begin = eye<double>(3);
          for(size_t ti = 0; ti < tet_loop.size() -1; ++ti) {
              boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit
                  = inner_face_jump_type.find(make_pair(tet_loop[ti], tet_loop[ti + 1]));
              if(cit == inner_face_jump_type.end()) continue;
              if(ti == 0)
                rot_begin = type_transition2(cit->second);
              rot = temp(rot * type_transition2(cit->second));
            }

          // this edge is already regular, just remove this edge
          // and recovery the g_unknown_face and continue
          if(type_transition1(rot) == TRIVIAL_TYPE){
              one_chain_list.erase(lit);
              g_unknown_face_in_patch.insert(face_idx);
              continue;
            }

          rot_begin = temp(trans(rot) * rot_begin);
          inner_face_jump_type[make_pair(tet_loop[0], tet_loop[1])]
              = type_transition1(rot_begin);
          inner_face_jump_type[make_pair(tet_loop[1], tet_loop[0])]
              = type_transition1(trans(rot_begin));
        }

        const vector<size_t> & face_vec = fa.faces_[face_idx];
        const size_t other_point =
            std::accumulate(face_vec.begin(), face_vec.end(), static_cast<size_t>(0))
            - lit->first - lit->second;
        pair<size_t,size_t> edge = *lit;

        vector<pair<size_t,size_t> > candidate_edges(2);
        candidate_edges[0] = make_pair(other_point,edge.second);
        candidate_edges[1] = make_pair(edge.first, other_point);

        one_chain_list.erase(lit);

        for(size_t ei = 0; ei < candidate_edges.size(); ++ei){
            if(candidate_edges[ei].first > candidate_edges[ei].second){
                swap(candidate_edges[ei].first, candidate_edges[ei].second);
              }
            // check if this edge is already in this chain list, delete it
            // or, should check whether this edge is a singularity edge,
            // if yes, can not include this one to the chain
            // if no, add it to this chain list
            list<pair<size_t,size_t> >::iterator lit =
                find(one_chain_list.begin(), one_chain_list.end(), candidate_edges[ei]);
            if(lit != one_chain_list.end()){
                one_chain_list.erase(lit); // is this needed??
                continue;
              }else{
                boost::unordered_set<pair<size_t,size_t> >::const_iterator
                    buscit = all_singularity_edges.find(candidate_edges[ei]);
                if(buscit != all_singularity_edges.end()){
                    continue;
                  }else{
                    // if reach surface, then stop
                    const size_t edge_idx = ea.get_edge_idx(candidate_edges[ei].first,
                                                            candidate_edges[ei].second);
                    if(edge_idx != -1) continue;
                    one_chain_list.push_front(candidate_edges[ei]);
                  }
              }
          }
      }
    }
  return 0;
}

bool is_near_miss_chain(
    const deque<pair<size_t,size_t> > & chain_orig,
    const boost::unordered_map<pair<size_t,size_t>, vector<
    std::tuple<size_t,size_t,size_t> > > &orig_edges2cut_edges,
    const boost::unordered_map<size_t,size_t> & surface_cut_type,
    const jtf::mesh::one_ring_face_at_point &orfap,
    const matrixst & outside_face_cut_idx,
    const matrixst & outside_face)
{
  pair<size_t,size_t> front_edge = chain_orig.front();
  pair<size_t,size_t> back_edge = chain_orig.back();

  vector<bool> is_reverse_edge(2, false);
  if(front_edge.first > front_edge.second){
      is_reverse_edge[0] = true;
      swap(front_edge.first, front_edge.second);
    }

  if(back_edge.first > back_edge.second) {
      is_reverse_edge[1] = true;
      swap(back_edge.first, back_edge.second);
    }

  typedef boost::unordered_map<pair<size_t,size_t>, vector<
      std::tuple<size_t,size_t,size_t> > > orig_edge2_cut_edge_type;

  for(size_t i = 0; i < 2; ++i){
      if(i == 0){
          if(!is_reverse_edge[i] &&
             find(outside_face.begin(), outside_face.end(), front_edge.first)
             == outside_face.end())
            continue;
          if(is_reverse_edge[i] &&
             find(outside_face.begin(), outside_face.end(), front_edge.second)
             == outside_face.end())
            continue;
        }else{
          if(!is_reverse_edge[i] &&
             find(outside_face.begin(), outside_face.end(), back_edge.second)
             == outside_face.end())
            continue;
          if(is_reverse_edge[i] &&
             find(outside_face.begin(), outside_face.end(), back_edge.first)
             == outside_face.end())
            continue;
        }

      orig_edge2_cut_edge_type::const_iterator ocit =
          orig_edges2cut_edges.find((i == 0?front_edge:back_edge));
      if(ocit == orig_edges2cut_edges.end()){
          cerr << "# [error] can not find edge " << __LINE__ << endl;
          return __LINE__;
        }

      for(size_t ei = 0; ei < ocit->second.size(); ++ei){
          const size_t type = get<2>(ocit->second[ei]);
          if(type == -1) continue; // free edge
          pair<size_t,size_t> cut_edge(get<0>(ocit->second[ei]),
                                       get<1>(ocit->second[ei]));
          jtf::mesh::one_ring_face_at_point::p2f_type::const_iterator orfapcit;
          if(i == 0){// chain begin should check the faces arounding first point
              if(is_reverse_edge[i]){
                  orfapcit =  orfap.p2f_.find(cut_edge.second);
                }else
                orfapcit = orfap.p2f_.find(cut_edge.first);
            }else{// chain end should check the faces arounding back point
              if(is_reverse_edge[i]){
                  orfapcit =  orfap.p2f_.find(cut_edge.first);
                }else
                orfapcit = orfap.p2f_.find(cut_edge.second);
            }
          if(orfapcit == orfap.p2f_.end()){
              cerr << "# [error] can not find one_ring_face for vertex." << endl;
              return __LINE__;
            }

          const vector<size_t> & around_faces = orfapcit->second;
          for(size_t fi = 0; fi < around_faces.size(); ++fi){
              const size_t &face_idx = outside_face_cut_idx[around_faces[fi]];
              boost::unordered_map<size_t,size_t>::const_iterator bumcit =
                  surface_cut_type.find(face_idx);
              if(bumcit != surface_cut_type.end()){
                  if(bumcit->second != type)
                    return true;
                }
            }
        }
    }
  return false;
}


int remove_near_miss_by_minimal_cut_old(
    const matrixst &cut_tet,
    const matrixst &tet,
    const matrixd & cut_node,
    const matrixd & node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    boost::unordered_map<size_t,size_t> &surface_type,
    const vector<pair<size_t,size_t> > & g_unknown_face_pair,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const vector<size_t> & fnode)
{
  matrixst cut_tet2tet(max(cut_tet)+1);
  cut_tet2tet(cut_tet) = tet(colon());

  jtf::mesh::one_ring_tet_at_edge ortae;
  ortae.add_tets(tet, fa);
  ortae.sort_into_loop(tet, node);
  matrixst outside_face, outside_face_cut, outside_face_cut_idx;
  //vector<deque<pair<size_t,size_t> > > possible_near_miss;
  get_outside_face(fa, outside_face);
  get_outside_face(fa_cut, outside_face_cut);
  get_outside_face_idx(fa_cut, outside_face_cut_idx);

  jtf::mesh::meshes tm_mesh;
  tm_mesh.mesh_ = tet;
  tm_mesh.node_ = node;
  jtf::tet_mesh tm(tm_mesh);
  singularity_extractor se(tm);

  boost::unordered_set<pair<size_t,size_t> > singularity_edges_set;
  vector<deque<pair<size_t,size_t> > > singularity_chains;
  {

    vector<deque<size_t> > singularity_type;
    se.extract(inner_face_jump_type, singularity_chains, singularity_type);

    for(size_t ci = 0; ci < singularity_chains.size(); ++ci){
        const deque<pair<size_t,size_t> > & one_chain = singularity_chains[ci];
        for(size_t ei = 0; ei < one_chain.size(); ++ei){
            const pair<size_t,size_t> & one_edge = one_chain[ei];
            if(one_edge.first > one_edge.second)
              singularity_edges_set.insert(make_pair(one_edge.second, one_edge.first));
            else
              singularity_edges_set.insert(one_edge);
          }
      }

    { // visual
      dump_singularity_chain_to_vtk_2("before_remove_near_miss.vtk", node,
                                      singularity_chains,  singularity_type);
    }

    //    for(size_t ci = 0; ci < singularity_chains.size(); ++ci){
    //      // a near miss should reach surface at both ends
    //      if(find(outside_face.begin(), outside_face.end(),
    //              singularity_chains[ci].front().first) != outside_face.end()
    //         && find(outside_face.begin(), outside_face.end(),
    //                 singularity_chains[ci].back().second) != outside_face.end())
    //      {
    //        possible_near_miss.push_back(singularity_chains[ci]);
    //      }
    //    }
  }

  // tuple<size_t,size_t,size_t>: p0,p1,type(0:u, 1:v, 2:w, 3:compound, -1:free)
  boost::unordered_map<pair<size_t,size_t>,
      vector<std::tuple<size_t,size_t,size_t> > > orig_edges2cut_edges;
  {
    jtf::mesh::one_ring_tet_at_edge ortae_cut;
    ortae_cut.add_tets(cut_tet, fa_cut);
    vector<size_t> fixed_axis;
    size_t type = -1;
    for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator ecit
        = ortae_cut.e2t_.begin(); ecit != ortae_cut.e2t_.end(); ++ecit){
        const pair<size_t,size_t> & cut_edge = ecit->first;
        fixed_axis.clear();
        for(size_t i = 0; i < 3; ++i){
            assert(cut_edge.first * 3 + i < fnode.size());
            assert(cut_edge.second * 3 + i < fnode.size());

            if(fnode[cut_edge.first * 3 + i] ==
               fnode[cut_edge.second * 3 + i])
              fixed_axis.push_back(i);
          }
        if(fixed_axis.empty()) type = -1; // free
        else if(fixed_axis.size() == 3) type= 3; // compound
        else if(fixed_axis.size() == 2) type = 3 - fixed_axis[0] - fixed_axis[1];
        else if(fixed_axis.size() == 1) type = -1;

        pair<size_t,size_t> orig_edge(cut_tet2tet[cut_edge.first],
            cut_tet2tet[cut_edge.second]);
        if(orig_edge.first > orig_edge.second){
            orig_edges2cut_edges[
                make_pair(orig_edge.second,orig_edge.first)].push_back(
                  std::make_tuple(cut_edge.second, cut_edge.first, type));
          }else
          orig_edges2cut_edges[orig_edge].push_back(
                std::make_tuple(cut_edge.first,
                                cut_edge.second, type));
      }
  }

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea_cut(
        jtf::mesh::edge2cell_adjacent::create(outside_face_cut));
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));

  if(!ea_cut.get() || !ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }

  jtf::mesh::one_ring_face_at_point orfap;
  orfap.add_all_faces(outside_face_cut, *ea_cut);
  orfap.sort_int_loop(outside_face_cut, cut_node);

  jtf::mesh::one_ring_face_at_edge orfae;
  orfae.build(ortae, fa, tet);

  boost::unordered_set<size_t> g_unknown_face_in_original_set;
  for(size_t fi = 0; fi < g_unknown_face_pair.size(); ++fi){
      const vector<size_t> & face_vec = fa_cut.faces_[g_unknown_face_pair[fi].first];
      const size_t face_original =
          fa.get_face_idx(cut_tet2tet[face_vec[0]],cut_tet2tet[face_vec[1]],
          cut_tet2tet[face_vec[2]]);
      if(face_original == -1){
          cerr << "# [error] can not find face " << g_unknown_face_pair[fi].first
               << " in original tetmesh." << endl;
          return __LINE__;
        }
      g_unknown_face_in_original_set.insert(face_original);
    }


  deque<size_t> singularity_chain_queue;
  for(size_t ci = 0; ci < singularity_chains.size(); ++ci)
    singularity_chain_queue.push_back(ci);
  //for(size_t ci = 0; ci < singularity_chains.size(); ++ci){
  size_t queue_size = singularity_chain_queue.size();
  size_t queue_size_keeping_time = 0;
  while(!singularity_chain_queue.empty()){
      const deque<pair<size_t,size_t> > & one_chain =
          singularity_chains[singularity_chain_queue.front()];
      if(is_near_miss_chain(one_chain, orig_edges2cut_edges, surface_type,
                            orfap, outside_face_cut_idx, outside_face)){
          if(remove_near_miss_chain_by_modifying_g_unknown_faces_old(
               one_chain, fa, ortae, orfae,
               *ea, g_unknown_face_in_original_set,
               inner_face_jump_type, outside_face, node)){
              singularity_chain_queue.push_back(singularity_chain_queue.front());
              singularity_chain_queue.pop_front();
            }else{
              singularity_chain_queue.pop_front();
            }
        }else
        singularity_chain_queue.pop_front();

      if(singularity_chain_queue.size() == queue_size)
        ++queue_size_keeping_time ;
      else{
          queue_size_keeping_time = 0;
          queue_size = singularity_chain_queue.size();
        }
      if(queue_size_keeping_time == queue_size)
        break;
    }

  { // visual
    ostringstream os;
    os << num++;
    vector<deque<pair<size_t,size_t> > > singularity_chain;
    vector<deque<size_t> > singularity_type;
    se.extract(inner_face_jump_type, singularity_chain, singularity_type);

    string file_name = os.str() + "_singularity.vtk";
    dump_singularity_chain_to_vtk_2(file_name.c_str(), node,
                                    singularity_chain,  singularity_type);

    vector<size_t> left_g_unknown_faces;
    for(boost::unordered_set<size_t>::const_iterator buscit =
        g_unknown_face_in_original_set.begin(); buscit != g_unknown_face_in_original_set.end();
        ++buscit){
        const vector<size_t> & face_vec = fa.faces_[*buscit];
        left_g_unknown_faces.insert(left_g_unknown_faces.end(), face_vec.begin(), face_vec.end());
      }
    ofstream ofs((os.str()+"_left_g_unknown_face.vtk").c_str());
    tri2vtk(ofs, &node[0], node.size(2), &left_g_unknown_faces[0], left_g_unknown_faces.size()/3);
  }
  return 0;
}

int remove_near_miss_chain_by_modifying_g_unknown_faces_old(
    const deque<pair<size_t,size_t> > & one_chain,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::one_ring_face_at_edge & orfae,
    const jtf::mesh::edge2cell_adjacent & ea, // surface edge
    boost::unordered_set<size_t> & g_unknown_face_in_original,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type,
    const matrixst & outside_face,
    const matrixd & node)
{
  boost::unordered_set<size_t>  g_unknown_face_in_original_temp =
      g_unknown_face_in_original;

  boost::unordered_map<pair<size_t,size_t>,size_t> inner_face_jump_type_temp =
      inner_face_jump_type;

  list<pair<size_t,size_t> > one_chain_list;
  one_chain_list.insert(one_chain_list.end(), one_chain.begin(), one_chain.end());

  pair<size_t,size_t> first_edge(-1,-1); // to avoid dead cycle

  matrixd rot = eye<double>(3);
  matrixd rot_begin = eye<double>(3);
  while(!one_chain_list.empty()){
      list<pair<size_t,size_t> >::iterator lit = one_chain_list.begin();

      // check if surface edge
      {
        const size_t edge_idx = ea.get_edge_idx(lit->first, lit->second);
        if(edge_idx != -1){
            one_chain_list.erase(lit);
            continue;
          }
      }

#if 0
      {
        ostringstream os;
        os << num++;
        vector<deque<pair<size_t,size_t> > > singularity_chain;
        vector<deque<size_t> > singularity_type;
        find_singularities_use_face_type(ortae,inner_face_jump_type, outside_face,
                                         singularity_chain, singularity_type);
        string file_name = os.str() + "_singularity.vtk";
        dump_singularity_chain_to_vtk_2(file_name.c_str(), node,
                                        singularity_chain,  singularity_type);

        vector<size_t> left_g_unknown_faces;
        for(boost::unordered_set<size_t>::const_iterator buscit =
            g_unknown_face_in_original.begin(); buscit != g_unknown_face_in_original.end();
            ++buscit){
            const vector<size_t> & face_vec = fa.faces_[*buscit];
            left_g_unknown_faces.insert(left_g_unknown_faces.end(), face_vec.begin(), face_vec.end());
          }
        ofstream ofs((os.str()+"_left_g_unknown_face.vtk").c_str());
        tri2vtk(ofs, &node[0], node.size(2), &left_g_unknown_faces[0], left_g_unknown_faces.size()/3);
      }
#endif

      jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator ecit
          = ortae.e2t_.find(*lit);
      if(ecit == ortae.e2t_.end())
        ecit = ortae.e2t_.find(make_pair(lit->second, lit->first));
      if(ecit == ortae.e2t_.end()){
          cerr << "# [error] can not find edge " << lit->first << " "
               << lit->second << " in one_ring_tet_at_edge " << endl;
          return __LINE__;
        }

      jtf::mesh::one_ring_face_at_edge::edge2_face_type::const_iterator
          cit = orfae.e2f_.find(*lit);
      if(cit == orfae.e2f_.end())
        cit = orfae.e2f_.find(make_pair(lit->second, lit->first));
      if(cit == orfae.e2f_.end()){
          cerr << "# [error] can not find edge " << lit->first << " "
               << lit->second << " in one_ring_face_at_edge " << endl;
          return __LINE__;
        }
      const vector<size_t> &around_faces = cit->second;
      size_t face_idx = -1;
      vector<size_t> face_idx_candidate;
      for(size_t fi = 0; fi < around_faces.size(); ++fi){
          boost::unordered_set<size_t>::iterator buit =
              g_unknown_face_in_original.find(around_faces[fi]);
          if(buit != g_unknown_face_in_original.end()){
              //g_unknown_face_in_original.erase(buit);
              face_idx_candidate.push_back(around_faces[fi]);
              //break;
            }
        }
      // here choose crease first, if all singularity edges are not crease
      // pick the first one
      if(face_idx_candidate.size() > 1) {
          if(one_chain_list.front() == first_edge){
              face_idx = face_idx_candidate.front();
              boost::unordered_set<size_t>::iterator buit =
                  g_unknown_face_in_original.find(face_idx);
              assert(buit != g_unknown_face_in_original.end());
              g_unknown_face_in_original.erase(buit);
              first_edge = make_pair(-1,-1);
            }else{
              if(first_edge.first == -1 && first_edge.second == -1)
                first_edge = one_chain_list.front();
              one_chain_list.push_back(one_chain_list.front());
              one_chain_list.pop_front();
              continue;
            }
        }else{
          if(one_chain_list.front() == first_edge)
            first_edge = make_pair(-1,-1);
          if(face_idx_candidate.empty())
            face_idx = -1;
          else{
              face_idx = face_idx_candidate.front();
              boost::unordered_set<size_t>::iterator buit =
                  g_unknown_face_in_original.find(face_idx);
              assert(buit != g_unknown_face_in_original.end());
              g_unknown_face_in_original.erase(buit);
            }
        }
      if(face_idx == -1){
          // can not find other face to go, then such edge must be identity.
          // check
          const vector<size_t> & tet_loop = ecit->second;
          rot = eye<double>(3);
          for(size_t ti = 0; ti < tet_loop.size() -1; ++ti) {
              boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit
                  = inner_face_jump_type.find(make_pair(tet_loop[ti], tet_loop[ti + 1]));
              if(cit == inner_face_jump_type.end()) continue;
              rot = temp(rot * type_transition2(cit->second));
            }
          if(type_transition1(rot) != TRIVIAL_TYPE){
              cerr << "# [error] edge " << lit->first << " " << lit->second
                   << " is not trivial type while having no face to go." << endl;
              inner_face_jump_type = inner_face_jump_type_temp;
              g_unknown_face_in_original = g_unknown_face_in_original_temp;
              return __LINE__;
            }else{
              one_chain_list.erase(lit);
              continue;
            }
        }

      // find the face which is in g_unkown_face_set, and need to find a suitable
      // rotation to ensure that this edge type is identity.
      {
        const pair<size_t,size_t> & tet_pair = fa.face2tet_[face_idx];
        assert(!fa.is_outside_face(tet_pair));

        vector<size_t> tet_loop = ecit->second;

        { // rotate loop at tet_pair
          tet_loop.pop_back();
          assert(find(tet_loop.begin(), tet_loop.end(), tet_pair.first) !=
              tet_loop.end());
          assert(find(tet_loop.begin(), tet_loop.end(), tet_pair.second) !=
              tet_loop.end());

          size_t ti = 0;
          for(; ti < tet_loop.size(); ++ti){
              if(tet_loop[ti] == tet_pair.first || tet_loop[ti] == tet_pair.second)
                break;
            }
          if(ti == 0 && (tet_loop.back() == tet_pair.first
                         ||tet_loop.back() == tet_pair.second)){
              tet_loop.insert(tet_loop.begin(), tet_loop.back());
              //tet_loop.pop_back();
            }else{
              std::rotate(tet_loop.begin(), tet_loop.begin() + ti, tet_loop.end());
              tet_loop.push_back(tet_loop.front());
            }
        }

        {
          // T = R01 * R12 * R23 * R30
          // T^-1 * T = (T^-1 * R01) * R12 * R23 * R30
          rot = eye<double>(3);
          rot_begin = eye<double>(3);
          for(size_t ti = 0; ti < tet_loop.size() -1; ++ti) {
              boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit
                  = inner_face_jump_type.find(make_pair(tet_loop[ti], tet_loop[ti + 1]));
              if(cit == inner_face_jump_type.end()) continue;
              if(ti == 0)
                rot_begin = type_transition2(cit->second);
              rot = temp(rot * type_transition2(cit->second));
            }

          // this edge is already regular, just remove this edge
          // and recovery the g_unknown_face and continue
          if(type_transition1(rot) == TRIVIAL_TYPE){
              one_chain_list.erase(lit);
              g_unknown_face_in_original.insert(face_idx);
              continue;
            }

          rot_begin = temp(trans(rot) * rot_begin);
          inner_face_jump_type[make_pair(tet_loop[0], tet_loop[1])]
              = type_transition1(rot_begin);
          inner_face_jump_type[make_pair(tet_loop[1], tet_loop[0])]
              = type_transition1(trans(rot_begin));
        }

        const vector<size_t> & face_vec = fa.faces_[face_idx];
        const size_t other_point =
            std::accumulate(face_vec.begin(), face_vec.end(), static_cast<size_t>(0))
            - lit->first - lit->second;
        pair<size_t,size_t> edge = *lit;

        //      { // check whether the added two new edges are singularity edges
        //        // if they are, means we meet with other singularity edge, this chain
        //        // should recovery the inital state
        //        // if no, move on
        //        pair<size_t,size_t> edge_0(other_point, edge.second);
        //        if(edge_0.first > edge_0.second)
        //          swap(edge_0.first, edge_0.second);
        //        pair<size_t,size_t> edge_1(other_point, edge.first);
        //        if(edge_1.first > edge_1.second)
        //          swap(edge_1.first, edge_1.second);

        //        boost::unordered_set<pair<size_t,size_t> >::iterator it_0 =
        //            singularity_edges_set.find(edge_0);
        //        boost::unordered_set<pair<size_t,size_t> >::iterator it_1 =
        //            singularity_edges_set.find(edge_1);
        //        if(it_0 == singularity_edges_set.end() &&
        //           it_1 == singularity_edges_set.end()){
        //          singularity_edges_set.insert(edge_0);
        //          singularity_edges_set.insert(edge_1);
        //        }else{
        //          // meet other singularity chain, which means this chain should
        //          // not be processed, left it what it was
        //          inner_face_jump_type = inner_face_jump_type_temp;
        //          singularity_edges_set = singularity_edges_set_temp;
        //          return __LINE__;
        //        }
        //        if(edge.first > edge.second)
        //          swap(edge.first, edge.second);

        //        singularity_edges_set.insert(edge);
        //      }

        one_chain_list.erase(lit);
        one_chain_list.insert(one_chain_list.begin(),
                              make_pair(other_point,edge.second));
        one_chain_list.insert(one_chain_list.begin(),
                              make_pair(edge.first, other_point));
      }
    }
  return 0;
}


int iterately_remove_surface_degenerated_patch(
    matrixst &tet, matrixd &node,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type,
    boost::unordered_map<size_t,size_t> & surface_type)
{
  const size_t iteration_num = 10;
  std::vector<std::pair<size_t,size_t> > g_unknown_face_pair;
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }

  matrixst cut_tet = tet;
  matrixst cut_tet2tet(max(cut_tet)+1);
  cut_tet2tet(cut_tet) = tet(colon());
  matrixd cut_node = node;
  matrixst outside_face_idx, outside_face;
  boost::unordered_set<size_t> loop_points;
  boost::unordered_map<pair<size_t,size_t>,size_t> restricted_edges, restricted_edges_orig;

  for(size_t i = 0; i < iteration_num; ++i){
      if(!inner_face_jump_type.empty()){
          cut_tet.resize(0);
          cut_tetmesh(tet, node, *fa, inner_face_jump_type,
                      surface_type, true, cut_tet);
          if(cut_tet.size(2) != tet.size(2)){
              cerr << "# [error] cut_tet is broken." << endl;
              return __LINE__;
            }
          cut_tet2tet.resize(max(cut_tet)+1);
          cut_tet2tet(cut_tet) = tet(colon());
        }else{
          cut_tet = tet;
        }
      aggressive_assign_transition_with_type(
            tet, cut_tet, node, inner_face_jump_type, surface_type,false);


      if(load_loop_points("loop_points", loop_points))
        return __LINE__;

      if(loop_points.empty()){
          break;
        }

      if(load_g_unknown_face("g_unknown_jump_face", g_unknown_face_pair))
        return __LINE__;

      // there must be degeneration patches
      if(load_restricted_edge("restricted_edges",restricted_edges))
        return __LINE__;

      if(!restricted_edges.empty()){
          for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator
              cit = restricted_edges.begin(); cit != restricted_edges.end(); ++cit){
              pair<size_t,size_t> edge = cit->first;
              edge.first = cut_tet2tet[edge.first];
              edge.second = cut_tet2tet[edge.second];
              if(edge.first > edge.second)
                swap(edge.first, edge.second);
              restricted_edges_orig[edge] = cit->second;
            }
        }

      //    int rtn = collapse_degenerated_edges(tet, node, *fa, inner_face_jump_type,
      //                                         surface_type, restricted_edges);
      //    if(rtn == 1) continue; // collapse several degenerated edges

      get_outside_face_idx(*fa, outside_face_idx);
      get_outside_face(*fa, outside_face);

      unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
            jtf::mesh::edge2cell_adjacent::create(outside_face));
      if(!ea.get()){
          cerr << "# [error] can not build edge2cell_adjacent." << endl;
          return __LINE__;
        }

      //    remove_surface_critcial_points(
      //          tet, node, outside_face, *fa, *ea, surface_type);

      straighten_patch_boundary(outside_face, outside_face_idx, *ea, surface_type);

      remove_dgenerated_patch_by_modify_type(
            tet, cut_tet,node, *fa, outside_face, outside_face_idx, *ea,
            g_unknown_face_pair, surface_type);

      //    // remove surface zigzag by splitting.
      //    remove_surface_zigzag_by_restricted_edges_from_graph(
      //          tet, node, outside_face_idx, *fa, *ea, restricted_edges_orig, surface_type,
      //          tet, node, inner_face_jump_type);

      cerr << "# [info] finish remove surface_zigzag." << endl;

      cerr << "# [info] finish iteration " << i << endl;
    }

  // if there is one edge does not exist on any patch boundary, but both of
  // it's ending points are on patch boundary, this edge should be splitted
  // or it will result in degeneration
  cerr << "# [info] remove face_degeneration_by_splitting." << endl;
  remove_face_degeneration_by_splitting(*fa, surface_type, node,
                                        tet, inner_face_jump_type);

  //remove_degenerated_patch_and_break_through(tet, node, surface_type);

  if(!inner_face_jump_type.empty()){
      cut_tet.resize(0);
      cut_tetmesh(tet, node, *fa, inner_face_jump_type,
                  surface_type, true, cut_tet);
      cut_tet2tet.resize(max(cut_tet)+1);
      cut_tet2tet(cut_tet) = tet(colon());
    }else
    cut_tet = tet;
  aggressive_assign_transition_with_type(
        tet, cut_tet, node, inner_face_jump_type, surface_type,false);

  cerr << "# [info] ++++++++++++++++type issue+++++++++++++++++++++" << endl;
  cerr << "# [info] boundary issue number " << boundary_issue_num << endl;
  cerr << "# [info] degree issue number " << degree_issue_num << endl;
  cerr << "# [info] finial dump out type setting." << endl;
  return 1;
}
