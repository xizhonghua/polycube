#include <sstream>
#include <numeric>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/lambda/lambda.hpp>

#include <jtflib/mesh/util.h>

#include "../common/transition_type.h"
#include "../common/timer.h"
#include "../tetmesh/util.h"
#include "../common/vtk.h"
#include "../common/transition.h"
#include "../common/visualize_tool.h"
#include "../numeric/util.h"
#include "../equation_graph/equation_graph.h"
#include "../equation_graph/util.h"
#include "solution_searching.h"
#include "topology_analysis.h"
#include "hex_param.h"


using namespace std;
using namespace zjucad::matrix;

static bool is_plane_degenerated(
    const vector<size_t> & rot_type,
    const boost::unordered_map<size_t,size_t> &surface_idx_to_rot_idx)
{
  size_t type = -1;
  for(boost::unordered_map<size_t,size_t>::const_iterator cit
      = surface_idx_to_rot_idx.begin(); cit != surface_idx_to_rot_idx.end();
      ++cit){
      const size_t &idx = cit->second;
      if(rot_type[idx] == -1)
        return false;
      if(type == -1){
          type = rot_type[idx];
          continue;
        }else{
          if(type != rot_type[idx])
            return false;
        }
    }
  return true;
}


int searching_strategy::searching_solutions_by_using_disk(
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
    const boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const zjucad::matrix::matrix<matrixd> * frame_ptr)
{
  boost::unordered_map<size_t,size_t> surface_type_cut;
  matrixd cut_node;
  matrixst cut_tet2tet,outside_face, outside_face_cut,outside_face_idx_in_cut, outside_face_idx;
  vector<pair<size_t,size_t> > jump_face_vec;
  jtf::mesh::one_ring_tet_at_edge  ortae;
  jtf::mesh::one_ring_tet_at_edge ortae_cut;
  prev_process(orig_tet,node,cut_tet, fa, fa_cut,face_pair_cut, surface_type, cut_node,
               cut_tet2tet, outside_face, outside_face_cut, outside_face_idx_in_cut, outside_face_idx,
               jump_face_vec, surface_type_cut, ortae, ortae_cut);

  vector<size_t> rot_type(jump_face_vec.size() + surface_type_cut.size(),-1);
  matrixst type_candidates(24,rot_type.size());
  vector<bool> is_surface(rot_type.size(), false);

  boost::unordered_map<size_t,size_t> surface_idx_to_rot_idx;
  boost::unordered_map<size_t,size_t> surface_idx_with_rot_idx;
  boost::unordered_map<std::pair<size_t,size_t>,size_t> face_pair_to_rot_idx;
  boost::unordered_map<size_t,std::pair<size_t,size_t> > face_pair_with_rot_idx;
  boost::unordered_map<std::pair<size_t,size_t>,size_t> tet_pair_to_rot_idx;

  adjust_face_and_type_order_to_speed_up(
        jump_face_vec, inner_face_jump_type, surface_type_cut, fa_cut, cut_tet2tet,
        type_candidates, surface_idx_to_rot_idx,surface_idx_with_rot_idx,
        face_pair_to_rot_idx,face_pair_with_rot_idx,tet_pair_to_rot_idx,
        is_surface,frame_ptr);

  unique_ptr<singularity_graph> g(
        singularity_graph::create(ortae_cut,ortae,cut_tet2tet,outside_face_cut,
                                  face_pair_cut,outside_face_idx_in_cut,
                                  surface_idx_to_rot_idx,face_pair_to_rot_idx));
  if(!g.get()){
      cerr << "# [error] can not build singularity graph." << endl;
      return __LINE__;
    }

  vector<size_t> tried_type_idx(rot_type.size(), -1);
  // start to check the question tree, and try to figure out the answer
  //step_state ss;
  size_t fi = 0;
  size_t searching_steps = 0;
  const size_t fn = rot_type.size();
  const string iter_data = "iter_data_";
  double write_time = 0.0, read_time = 0.0, check_time = 0.0, trans_time = 0.0;
  timer tr;
  bool is_valid;
  deque<size_t> check_numbers;
  while(fi < fn) {
      // cerr << "# [info] checking face " << fi << ", tried "  << tried_type_idx[fi] << endl;
      size_t ri;
      for(ri = tried_type_idx[fi]+1; ri < 24; ++ri) {
          ++searching_steps;
          is_valid = true;
          const size_t &next_type = type_candidates(ri, fi);

          if(next_type == -1) break;  // index = -1 meas the surface type is tested over

          long_number_add(check_numbers, 1);

          ostringstream oos_fi, oos_ri;
          oos_fi << fi;
          oos_ri << ri;
          const string iter_data_file = iter_data + oos_fi.str() + "_" + oos_ri.str();

          rot_type[fi] = next_type;
          tried_type_idx[fi] = ri;

          if(is_surface[fi] &&
             is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
            continue;

          tr.start();
          g->save_state(iter_data_file);
          tr.finish();
          write_time += tr.result_c();

          assert(!g->is_group_info_broken());

          tr.start();
          if(is_surface[fi]) { // the third type = -1 means this face is surface
              g->insert_surface_transition(next_type, fa_cut, cut_tet2tet,
                                           surface_idx_with_rot_idx[fi], rot_type,
                                           tet_pair_to_rot_idx);
            }else {
              if(g->insert_transition(face_pair_with_rot_idx[fi],cut_tet,cut_tet2tet,
                                      fa_cut, ortae,rot_type,face_pair_to_rot_idx,
                                      tet_pair_to_rot_idx,1))
                is_valid = false;
            }
          tr.finish();
          trans_time += tr.result_c();

          assert(!g->is_group_info_broken());

          tr.start();
          if(is_valid)
            is_valid = g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                             surface_idx_to_rot_idx,1);
          tr.finish();
          check_time += tr.result_c();

          if(is_valid)      break;

          {
            assert(!is_valid);
            rot_type[fi] = -1;// try_rot_type.second;
            tr.start();
            g->load_state(iter_data_file);
            tr.finish();
            read_time += tr.result_c();

            assert(!g->is_group_info_broken());
            assert(g->is_valid_with_info(cut_tet2tet,rot_type, ortae,
                                         surface_idx_to_rot_idx));

            // to save the disk space, here we remvoe this file
            boost::filesystem::path p(iter_data_file.c_str());
            assert(boost::filesystem::exists(p));
            boost::filesystem::remove(p);
            if(is_valid)
              cerr << "# shit happens" << endl;
          }
        }
      // revert to the previous cut_face setting.
      if(( !is_surface[fi] && ri == 24) || (is_surface[fi] && ri == 3)) {
          if(fi == 0) {
              break; // can not find one solution
            }

          ostringstream oos_fi,oos_ti;
          oos_fi << fi - 1;
          oos_ti << tried_type_idx[fi-1];
          const string iter_data_file = iter_data + oos_fi.str() + "_" + oos_ti.str();

          tried_type_idx[fi] = -1;
          rot_type[fi] = -1;
          --fi;

          tr.start();
          g->load_state(iter_data_file);
          tr.finish();
          read_time += tr.result_c();

          assert(!g->is_group_info_broken());
          assert(g->is_valid_with_info(cut_tet2tet,rot_type, ortae,
                                       surface_idx_to_rot_idx));
          // to save the disk space, here we remvoe this file
          boost::filesystem::path p(iter_data_file.c_str());
          assert(boost::filesystem::exists(p));
          boost::filesystem::remove(p);
          rot_type[fi] = -1;
        }else{
          ++fi;
        }

    }

  cerr << "# write time " << write_time/1e6 << endl;
  cerr << "# read time " << read_time/1e6 << endl;
  cerr << "# trans time " << trans_time/1e6 << endl;
  cerr << "# check time " << check_time/1e6 << endl;
  cerr << "# searching steps " << searching_steps << endl;
  cerr << "# iterate numbers ";
  print_long_number(cerr, check_numbers);
  if(fi == 0) {
      cerr << "cannot find a solution" << endl;
      return __LINE__;
    } else {
      cerr << "find a solution. " << endl;

      post_process(rot_type,cut_node, node, outside_face, fa_cut,ortae,
                   surface_idx_to_rot_idx, tet_pair_to_rot_idx, face_pair_to_rot_idx,
                   inner_face_jump_type);

      simple_check(rot_type, is_surface, jump_face_vec, cut_tet, cut_node, node,
                   fa_cut, cut_tet2tet, ortae, ortae_cut, outside_face_cut,
                   outside_face_idx_in_cut, face_pair_cut, surface_idx_to_rot_idx,
                   face_pair_to_rot_idx,tet_pair_to_rot_idx);
      return 0;
    }
}

int searching_strategy::adjust_face_and_type_order_to_speed_up(
    const std::vector<std::pair<size_t,size_t> > & jump_face_vec,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const boost::unordered_map<size_t,size_t> & surface_type_cut,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & cut_tet2tet,
    matrixst & type_candidates,
    boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
    boost::unordered_map<size_t,size_t> & surface_idx_with_rot_idx,
    boost::unordered_map<std::pair<size_t,size_t>,size_t > & face_pair_to_rot_idx,
    boost::unordered_map<size_t, pair<size_t,size_t> > & face_pair_with_rot_idx,
    boost::unordered_map<pair<size_t,size_t>,size_t> & tet_pair_to_rot_idx,
    vector<bool> & is_surface,
    const matrix<matrixd> * frame_ptr)
{
  face_pair_to_rot_idx.clear();
  const size_t type_num = jump_face_vec.size() + surface_type_cut.size();
  type_candidates = ones<size_t>(24, type_num) * (-1);
  is_surface.resize(type_num);
  // this data record the index in faces_mat before reorder to surface index
  boost::unordered_map<size_t,size_t> surface_cut_with_special_idx;
  vector<size_t> reorder_face_idx;
  { // reorder the inner faces and surface faces
    matrixst faces(3, jump_face_vec.size() + surface_type_cut.size());
    size_t fi = 0;
    for(; fi < jump_face_vec.size(); ++fi){
        const vector<size_t> & face_of_cut = fa_cut.faces_[jump_face_vec[fi].first];
        for(size_t pi = 0; pi < 3; ++pi)
          faces(pi,fi) = cut_tet2tet[face_of_cut[pi]];
      }
    for(boost::unordered_map<size_t,size_t>::const_iterator bumcit
        = surface_type_cut.begin(); bumcit != surface_type_cut.end();
        ++bumcit, ++fi){
        const vector<size_t> & face_of_cut = fa_cut.faces_[bumcit->first];
        for(size_t pi = 0; pi < 3; ++pi)
          faces(pi,fi) = cut_tet2tet[face_of_cut[pi]];
        surface_cut_with_special_idx[fi] = bumcit->first;
      }

    boost::unordered_map<pair<size_t,size_t>, vector<size_t> > edge2face;
    for(size_t fi = 0; fi < faces.size(2); ++fi){
        for(size_t pi = 0; pi < faces.size(1); ++pi){
            pair<size_t,size_t> one_edge(faces(pi,fi), faces((pi+1)%3,fi));
            if(one_edge.first > one_edge.second)
              swap(one_edge.first, one_edge.second);
            edge2face[one_edge].push_back(fi);
          }
      }

#define separate_mst
#ifdef  whole_mst
    { // whole mst // there is a bug inside, should not be used anymore
      vector<bool> visited_face(faces.size(2), false);
      cerr << visited_face.size() << endl;
      stack<size_t> temp;
      temp.push(0);
      visited_face[0] = true;
      while(!temp.empty()){
          const size_t face_idx_in_mat = temp.top();
          temp.pop();
          //assert(visited_face[face_idx_in_mat] == false);
          //visited_face[face_idx_in_mat] = true;
          reorder_face_idx.push_back(face_idx_in_mat);
          for(size_t pi = 0; pi < 3; ++pi){
              pair<size_t,size_t> one_edge(faces(pi, face_idx_in_mat),
                                           faces((pi+1)%3, face_idx_in_mat));
              if(one_edge.first > one_edge.second)
                swap(one_edge.first, one_edge.second);
              boost::unordered_map<pair<size_t,size_t>,vector<size_t> >::const_iterator
                  bcit = edge2face.find(one_edge);
              assert(bcit != edge2face.end());
              const vector<size_t> & linked_faces = bcit->second;
              for(size_t li = 0; li < linked_faces.size(); ++li){
                  if(visited_face[linked_faces[li]]) continue;
                  temp.push(linked_faces[li]);
                  visited_face[linked_faces[li]] = true;
                }
            }
        }
      assert(find(visited_face.begin(), visited_face.end(), false)
             == visited_face.end());
    }
#endif

#ifdef separate_mst
    { // separate mst
      vector<bool> visited_face(faces.size(), false);
      stack<size_t> inner_mst;
      size_t non_visited_face_idx = 0;
      while(non_visited_face_idx < jump_face_vec.size()){
          assert(inner_mst.empty());
          inner_mst.push(non_visited_face_idx);
          visited_face[non_visited_face_idx] = true;
          while(!inner_mst.empty()){
              const size_t current_inner_face = inner_mst.top();
              reorder_face_idx.push_back(current_inner_face);
              inner_mst.pop();
              for(size_t pi = 0; pi < 3; ++pi){
                  pair<size_t,size_t> one_edge(faces(pi, current_inner_face),
                                               faces((pi+1)%3,current_inner_face));
                  if(one_edge.first > one_edge.second)
                    swap(one_edge.first, one_edge.second);
                  boost::unordered_map<pair<size_t,size_t>,vector<size_t> >::const_iterator
                      bcit = edge2face.find(one_edge);
                  assert(bcit != edge2face.end());
                  const vector<size_t> & linked_faces = bcit->second;
                  for(size_t li = 0; li < linked_faces.size(); ++li){
                      if(visited_face[linked_faces[li]]) continue;
                      if(linked_faces[li] < jump_face_vec.size()){ // inner face
                          inner_mst.push(linked_faces[li]);
                          visited_face[linked_faces[li]] = true;
                        }
                    }
                }
            }
          non_visited_face_idx = find(visited_face.begin(), visited_face.end(), false)
              - visited_face.begin();
        }

      assert(reorder_face_idx.size() == jump_face_vec.size());
      stack<size_t> out_mst;
      while(non_visited_face_idx < type_num){
          out_mst.push(non_visited_face_idx);
          visited_face[non_visited_face_idx] = true;

          while(!out_mst.empty()){
              const size_t current_inner_face = out_mst.top();
              reorder_face_idx.push_back(current_inner_face);
              out_mst.pop();
              for(size_t pi = 0; pi < 3; ++pi){
                  pair<size_t,size_t> one_edge(faces(pi, current_inner_face),
                                               faces((pi+1)%3,current_inner_face));
                  if(one_edge.first > one_edge.second)
                    swap(one_edge.first, one_edge.second);
                  boost::unordered_map<pair<size_t,size_t>,vector<size_t> >::const_iterator
                      bcit = edge2face.find(one_edge);
                  assert(bcit != edge2face.end());
                  const vector<size_t> & linked_faces = bcit->second;
                  for(size_t li = 0; li < linked_faces.size(); ++li){
                      if(visited_face[linked_faces[li]]) continue;
                      if(linked_faces[li] > jump_face_vec.size()){ // inner face
                          out_mst.push(linked_faces[li]);
                          visited_face[linked_faces[li]] = true;
                        }
                    }
                }
            }
          vector<bool>::const_iterator cit =
              find(visited_face.begin(), visited_face.end(), false);
          if(cit == visited_face.end())
            break;
          else
            non_visited_face_idx = cit  - visited_face.begin();
        }

      assert(reorder_face_idx.size() == type_num);
    }
#endif

#ifdef depth_priority
    {
      vector<bool> visited_face(faces.size(), false);
      deque<size_t> inner_deque;
      size_t non_visited_face_idx = 0;
      while(non_visited_face_idx < jump_face_vec.size()){
          assert(inner_deque.empty());
          inner_deque.push_back(non_visited_face_idx);
          visited_face[non_visited_face_idx] = true;
          while(!inner_deque.empty()){
              const size_t current_inner_face = inner_deque.front();
              reorder_face_idx.push_back(current_inner_face);
              inner_deque.pop_front();
              for(size_t pi = 0; pi < 3; ++pi){
                  pair<size_t,size_t> one_edge(faces(pi, current_inner_face),
                                               faces((pi+1)%3,current_inner_face));
                  if(one_edge.first > one_edge.second)
                    swap(one_edge.first, one_edge.second);
                  boost::unordered_map<pair<size_t,size_t>,vector<size_t> >::const_iterator
                      bcit = edge2face.find(one_edge);
                  assert(bcit != edge2face.end());
                  const vector<size_t> & linked_faces = bcit->second;
                  for(size_t li = 0; li < linked_faces.size(); ++li){
                      if(visited_face[linked_faces[li]]) continue;
                      if(linked_faces[li] < jump_face_vec.size()){ // inner face
                          inner_deque.push_back(linked_faces[li]);
                          visited_face[linked_faces[li]] = true;
                        }
                    }
                }
            }
          non_visited_face_idx = find(visited_face.begin(), visited_face.end(), false)
              - visited_face.begin();
        }

      assert(reorder_face_idx.size() == jump_face_vec.size());
      deque<size_t> out_deq;
      while(non_visited_face_idx < type_num){
          out_deq.push_back(non_visited_face_idx);
          visited_face[non_visited_face_idx] = true;

          while(!out_deq.empty()){
              const size_t current_inner_face = out_deq.front();
              reorder_face_idx.push_back(current_inner_face);
              out_deq.pop_front();
              for(size_t pi = 0; pi < 3; ++pi){
                  pair<size_t,size_t> one_edge(faces(pi, current_inner_face),
                                               faces((pi+1)%3,current_inner_face));
                  if(one_edge.first > one_edge.second)
                    swap(one_edge.first, one_edge.second);
                  boost::unordered_map<pair<size_t,size_t>,vector<size_t> >::const_iterator
                      bcit = edge2face.find(one_edge);
                  assert(bcit != edge2face.end());
                  const vector<size_t> & linked_faces = bcit->second;
                  for(size_t li = 0; li < linked_faces.size(); ++li){
                      if(visited_face[linked_faces[li]]) continue;
                      if(linked_faces[li] > jump_face_vec.size()){ // inner face
                          out_deq.push_back(linked_faces[li]);
                          visited_face[linked_faces[li]] = true;
                        }
                    }
                }
            }
          vector<bool>::const_iterator cit =
              find(visited_face.begin(), visited_face.end(), false);
          if(cit == visited_face.end())
            break;
          else
            non_visited_face_idx = cit  - visited_face.begin();
        }

      assert(reorder_face_idx.size() == type_num);
    }
#endif

#ifdef reliability
    {
      assert(frame_ptr);
      vector<pair<double,size_t> > inner_face_ordered, rot24(24);
      for(size_t fi = 0; fi < jump_face_vec.size(); ++fi){
          const pair<size_t,size_t> & face_pair = jump_face_vec[fi];
          assert(face_pair.first != -1 && face_pair.second != -1);
          const pair<size_t,size_t> & tet_pair_s = fa_cut.face2tet_[face_pair.first];
          const pair<size_t,size_t> & tet_pair_t = fa_cut.face2tet_[face_pair.second];
          assert(fa_cut.is_outside_face(tet_pair_s));
          assert(fa_cut.is_outside_face(tet_pair_t));
          pair<size_t,size_t> tet_pair(
                (tet_pair_s.first == -1?tet_pair_s.second:tet_pair_s.first),
                (tet_pair_t.first == -1?tet_pair_t.second:tet_pair_t.first));
          for(size_t ri = 0; ri < 24; ++ri){
              rot24[ri] = make_pair(norm((*frame_ptr)[tet_pair.first] * type_transition2(ri)
                                    - (*frame_ptr)[tet_pair.second]),ri);
            }
          sort(rot24.begin(), rot24.end());
          inner_face_ordered.push_back(make_pair( rot24[0].first / rot24[1].first,fi));
        }
      sort(inner_face_ordered.begin(), inner_face_ordered.end());

      for(size_t fi = 0; fi < inner_face_ordered.size(); ++fi)
        reorder_face_idx.push_back(inner_face_ordered[fi].second);

      assert(reorder_face_idx.size() == jump_face_vec.size());

      vector<bool> visited_face(faces.size(), false);
      std::fill(visited_face.begin(), visited_face.begin()+jump_face_vec.size(), true);
      deque<size_t> out_deq;
      size_t non_visited_face_idx = jump_face_vec.size();

      while(non_visited_face_idx < type_num){
          out_deq.push_back(non_visited_face_idx);
          visited_face[non_visited_face_idx] = true;

          while(!out_deq.empty()){
              const size_t current_inner_face = out_deq.front();
              reorder_face_idx.push_back(current_inner_face);
              out_deq.pop_front();
              for(size_t pi = 0; pi < 3; ++pi){
                  pair<size_t,size_t> one_edge(faces(pi, current_inner_face),
                                               faces((pi+1)%3,current_inner_face));
                  if(one_edge.first > one_edge.second)
                    swap(one_edge.first, one_edge.second);
                  boost::unordered_map<pair<size_t,size_t>,vector<size_t> >::const_iterator
                      bcit = edge2face.find(one_edge);
                  assert(bcit != edge2face.end());
                  const vector<size_t> & linked_faces = bcit->second;
                  for(size_t li = 0; li < linked_faces.size(); ++li){
                      if(visited_face[linked_faces[li]]) continue;
                      if(linked_faces[li] > jump_face_vec.size()){ // outside face
                          out_deq.push_back(linked_faces[li]);
                          visited_face[linked_faces[li]] = true;
                        }
                    }
                }
            }
          vector<bool>::const_iterator cit =
              find(visited_face.begin(), visited_face.end(), false);
          if(cit == visited_face.end())
            break;
          else
            non_visited_face_idx = cit  - visited_face.begin();
        }
      size_t total_num =
          std::accumulate(reorder_face_idx.begin(),
                          reorder_face_idx.end(), static_cast<size_t>(0));
      if(2*total_num != type_num * (type_num -1)){
          cerr << "# shit happens: total_num " << total_num
               << " face num " << type_num << endl;
          sort(reorder_face_idx.begin(), reorder_face_idx.end());
          for(size_t i = 0; i < reorder_face_idx.size(); ++i){
              reorder_face_idx[i] -= i;
              if(reorder_face_idx[i] != 0)
                cerr << "# [error] meet shit " << endl;
            }

          return __LINE__;
        }
    }
#endif
  }

  vector<pair<double, size_t> > left_order;
  matrixd rot = eye<double>(3);
  for(size_t fi = 0; fi < reorder_face_idx.size(); ++fi){
      if(reorder_face_idx[fi] < jump_face_vec.size()){ // inner face
          is_surface[fi] = false;
          const pair<size_t,size_t> & face_pair = jump_face_vec[reorder_face_idx[fi]];
          face_pair_to_rot_idx[face_pair] = fi;
          face_pair_with_rot_idx[fi] = face_pair;
          assert(face_pair.first != -1 && face_pair.second != -1);
          const pair<size_t,size_t> & tet_pair_s = fa_cut.face2tet_[face_pair.first];
          const pair<size_t,size_t> & tet_pair_t = fa_cut.face2tet_[face_pair.second];
          assert(fa_cut.is_outside_face(tet_pair_s));
          assert(fa_cut.is_outside_face(tet_pair_t));
          pair<size_t,size_t> tet_pair(
                (tet_pair_s.first == -1?tet_pair_s.second:tet_pair_s.first),
                (tet_pair_t.first == -1?tet_pair_t.second:tet_pair_t.first));
          tet_pair_to_rot_idx[tet_pair] = fi;
          boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator bumpcit
              = inner_face_jump_type.find(tet_pair);
          if(bumpcit == inner_face_jump_type.end()){
              type_candidates(0, fi) = TRIVIAL_TYPE;
            }else
            type_candidates(0, fi) = bumpcit->second;
          if(frame_ptr){
              left_order.clear();
              for(size_t i = 0; i < 24; ++i){
                  if(i == type_candidates(0,fi)) continue;
                  get_best_alignment(&(*frame_ptr)[tet_pair.first][0],
                      &(*frame_ptr)[tet_pair.second][0],
                      &rot[0]);
                  left_order.push_back(make_pair(norm(rot), i));
                }
              sort(left_order.begin(), left_order.end());
              assert(left_order.size() == 23);
              for(size_t i = 0; i < 23; ++i){
                  type_candidates(i+1, fi) = left_order[i].second;
                }
            }else{
              for(size_t i = 0; i < 23; ++i)
                type_candidates(i+1, fi) = (type_candidates(i,fi) + 1)%24;
            }
          assert(find(type_candidates(colon(),fi).begin(),
                      type_candidates(colon(),fi).end(),-1) ==
                 type_candidates(colon(),fi).end());
        }else{ // out surface
          is_surface[fi] = true;
          boost::unordered_map<size_t,size_t>::const_iterator cit =
              surface_cut_with_special_idx.find(reorder_face_idx[fi]);
          assert(cit != surface_cut_with_special_idx.end());
          surface_idx_to_rot_idx[cit->second] = fi;
          surface_idx_with_rot_idx[fi] = cit->second;

          //      static int a = 0;
          //      if(a++ == 0){
          //        type_candidates(0,fi) = (surface_type_cut[cit->second]/2+rand())%3;
          //        cerr << "# [info] attention, the first surface type is random." << endl;
          //      }
          //      else
          const auto it = surface_type_cut.find(cit->second);
          if(it == surface_type_cut.end())
            cerr << "# [error] can not find cit->second in surface_type_cut." << endl;
          else
          type_candidates(0,fi) = it->second%3;
          type_candidates(1,fi) = (type_candidates(0,fi)+1)%3;
          type_candidates(2,fi) = (type_candidates(1,fi)+1)%3;
          assert(type_candidates(0,fi) < 3);
        }
    }

  return 0;
}

int searching_strategy::searching_solutions_by_using_memory(
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
    const boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const zjucad::matrix::matrix<matrixd> * frame_ptr)
{
  boost::unordered_map<size_t,size_t> surface_type_cut;
  matrixd cut_node;
  matrixst cut_tet2tet,outside_face, outside_face_cut,outside_face_idx_in_cut, outside_face_idx;
  vector<pair<size_t,size_t> > jump_face_vec;
  jtf::mesh::one_ring_tet_at_edge  ortae;
  jtf::mesh::one_ring_tet_at_edge ortae_cut;
  prev_process(orig_tet,node,cut_tet, fa, fa_cut,face_pair_cut, surface_type, cut_node,
               cut_tet2tet, outside_face, outside_face_cut, outside_face_idx_in_cut, outside_face_idx,
               jump_face_vec, surface_type_cut, ortae, ortae_cut);

  vector<size_t> rot_type(jump_face_vec.size() + surface_type_cut.size(),-1);
  matrixst type_candidates(24,rot_type.size());
  vector<bool> is_surface(rot_type.size(), false);

  boost::unordered_map<size_t,size_t> surface_idx_to_rot_idx;
  boost::unordered_map<size_t,size_t> surface_idx_with_rot_idx;
  boost::unordered_map<pair<size_t,size_t>,size_t> face_pair_to_rot_idx;
  boost::unordered_map<pair<size_t,size_t>,size_t> tet_pair_to_rot_idx;
  boost::unordered_map<size_t,pair<size_t,size_t> > face_pair_with_rot_idx;

  adjust_face_and_type_order_to_speed_up(
        jump_face_vec, inner_face_jump_type, surface_type_cut, fa_cut, cut_tet2tet,
        type_candidates, surface_idx_to_rot_idx,surface_idx_with_rot_idx,
        face_pair_to_rot_idx,face_pair_with_rot_idx,tet_pair_to_rot_idx,
        is_surface,frame_ptr);

  unique_ptr<singularity_graph> g(
        singularity_graph::create(ortae_cut,ortae,cut_tet2tet,outside_face_cut,
                                  face_pair_cut,outside_face_idx_in_cut,
                                  surface_idx_to_rot_idx,face_pair_to_rot_idx));
  if(!g.get()){
      cerr << "# [error] can not build singularity graph." << endl;
      return __LINE__;
    }

  vector<size_t> tried_type_idx(rot_type.size(), -1);
  // start to check the question tree, and try to figure out the answer
  size_t fi = 0;
  size_t searching_steps = 0;
  const size_t fn = rot_type.size();
  stack<state_each_step> ss;
  state_each_step ses;
  double write_time = 0.0, read_time = 0.0, check_time = 0.0, trans_time = 0.0;
  timer tr;
  bool is_valid;
  deque<size_t> check_numbers;
  while(fi < fn) {
      // cerr << "# [info] checking face " << fi << ", tried "  << tried_type_idx[fi] << endl;
      size_t ri;
      for(ri = tried_type_idx[fi]+1; ri < 24; ++ri) {
          ++searching_steps;
          is_valid = true;
          const size_t &next_type = type_candidates(ri, fi);

          if(next_type == -1) break;  // index = -1 meas the surface type is tested over

          long_number_add(check_numbers, 1);
          rot_type[fi] = next_type;
          tried_type_idx[fi] = ri;

          if(is_surface[fi] &&
             is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
            continue;

          tr.start();
          g->save_state_mem(ses);
          ss.push(ses);
          tr.finish();
          write_time += tr.result_c();
          //g->get_step_state(ss);

          assert(!g->is_group_info_broken());

          tr.start();
          if(is_surface[fi]) { // the third type = -1 means this face is surface
              g->insert_surface_transition(next_type, fa_cut, cut_tet2tet,
                                           surface_idx_with_rot_idx[fi], rot_type,
                                           tet_pair_to_rot_idx);
            }else {
              if(g->insert_transition(face_pair_with_rot_idx[fi],cut_tet,cut_tet2tet,
                                      fa_cut, ortae,rot_type,face_pair_to_rot_idx,
                                      tet_pair_to_rot_idx,1))
                is_valid = false;
            }
          tr.finish();
          trans_time += tr.result_c();

          assert(!g->is_group_info_broken());

          tr.start();
          if(is_valid)
            is_valid = g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                             surface_idx_to_rot_idx,1);
          tr.finish();
          check_time += tr.result_c();

          if(is_valid)      break;

          {
            assert(!is_valid);
            rot_type[fi] = -1;// try_rot_type.second;
            tr.start();
            g->load_state_mem(ss.top());
            ss.pop();
            tr.finish();
            read_time += tr.result_c();

            assert(!g->is_group_info_broken());
            assert(g->is_valid_with_info(cut_tet2tet,rot_type, ortae,
                                         surface_idx_to_rot_idx));
          }
        }
      // revert to the previous cut_face setting.
      if(( !is_surface[fi] && ri == 24) || (is_surface[fi] && ri == 3)) {
          if(fi == 0) {
              break; // can not find one solution
            }

          tried_type_idx[fi] = -1;
          rot_type[fi] = -1;
          --fi;

          tr.start();
          g->load_state_mem(ss.top());
          ss.pop();
          tr.finish();
          read_time += tr.result_c();

          assert(!g->is_group_info_broken());
          assert(g->is_valid_with_info(cut_tet2tet,rot_type, ortae,
                                       surface_idx_to_rot_idx));
          rot_type[fi] = -1;
        }else{
          ++fi;
        }

    }

  cerr << "# write time " << write_time/1e6 << endl;
  cerr << "# read time " << read_time/1e6 << endl;
  cerr << "# trans time " << trans_time/1e6 << endl;
  cerr << "# check time " << check_time/1e6 << endl;
  cerr << "# searching steps " << searching_steps << endl;
  cerr << "# iterate numbers ";
  print_long_number(cerr, check_numbers);

  if(fi == 0) {
      cerr << "cannot find a solution" << endl;
      return __LINE__;
    } else {
      cerr << "find a solution. " << endl;
      post_process(rot_type,cut_node, node, outside_face, fa_cut, ortae,
                   surface_idx_to_rot_idx, tet_pair_to_rot_idx,face_pair_to_rot_idx,
                   inner_face_jump_type);

      simple_check(rot_type, is_surface, jump_face_vec, cut_tet, cut_node, node,
                   fa_cut, cut_tet2tet, ortae, ortae_cut, outside_face_cut,
                   outside_face_idx_in_cut, face_pair_cut, surface_idx_to_rot_idx,
                   face_pair_to_rot_idx,tet_pair_to_rot_idx);
      return 0;
    }
  return 0;
}

int searching_strategy::searching_solutions_by_using_memory_binary_search_group_polycube(
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
    const boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const zjucad::matrix::matrix<matrixd> * frame_ptr)
{
  boost::unordered_map<size_t,size_t> surface_type_cut;
  matrixd cut_node;
  matrixst cut_tet2tet,outside_face, outside_face_cut,
      outside_face_idx_in_cut, outside_face_idx;
  vector<pair<size_t,size_t> > jump_face_vec;
  jtf::mesh::one_ring_tet_at_edge  ortae;
  jtf::mesh::one_ring_tet_at_edge ortae_cut;

  prev_process(
        orig_tet,node,cut_tet, fa, fa_cut,face_pair_cut, surface_type, cut_node,
        cut_tet2tet, outside_face, outside_face_cut, outside_face_idx_in_cut,
        outside_face_idx, jump_face_vec, surface_type_cut, ortae, ortae_cut);

  vector<vector<pair<size_t,size_t> > > jump_face_groups;
  vector<vector<size_t> > surface_groups;
  vector<pair<size_t,double> > group_type_with_frame_reliablity,
      group_type_with_surface_reliablity;
  const double frame_reliablity_threshold = 0.5;

  // this function is used to group faces according to the type and reliablity

  boost::unordered_map<pair<size_t,size_t>,size_t>  face_pair_to_rot_idx;
  boost::unordered_map<pair<size_t,size_t>,size_t>  tet_pair_to_rot_idx;
  boost::unordered_map<size_t,size_t>  surface_idx_to_rot_idx;
  matrixst type_candidates;
  vector<bool> is_surface;

  group_faces_according_face_reliablity(
        cut_tet,cut_node,fa_cut, outside_face_cut,outside_face_idx_in_cut,
        face_pair_cut, *frame_ptr, jump_face_groups, surface_groups,
        face_pair_to_rot_idx, tet_pair_to_rot_idx, surface_idx_to_rot_idx,
        type_candidates, is_surface, frame_reliablity_threshold);

  unique_ptr<singularity_graph> g(
        singularity_graph::create(
          ortae_cut, ortae, cut_tet2tet, outside_face_cut,  face_pair_cut,
          outside_face_idx_in_cut, surface_idx_to_rot_idx, face_pair_to_rot_idx));

  if(!g.get()){
      cerr << "# [info] build singularity graph failed." << endl;
      return __LINE__;
    }

  vector<size_t> rot_type(jump_face_groups.size() + surface_groups.size() ,-1);
  vector<size_t> tried_type_idx(rot_type.size(), -1);

  std::fill(rot_type.begin(), rot_type.begin() + jump_face_groups.size(), TRIVIAL_TYPE);

  // start to check the question tree, and try to figure out the answer
  size_t searching_steps = 0;
  const size_t fn = rot_type.size();
  deque<state_each_step> ss;
  double write_time = 0.0, read_time = 0.0, check_time = 0.0, trans_time = 0.0;
  timer tr;
  bool is_valid;

  state_each_step ses;
  deque<size_t> check_numbers;


  size_t fi = jump_face_groups.size();
  size_t windows_size = fn - jump_face_groups.size();
  {
    is_valid = true;
    for(size_t gi = 0; gi < jump_face_groups.size(); ++gi){
        if(g->insert_transition_group(
             jump_face_groups[gi],cut_tet,cut_tet2tet,
             fa_cut, ortae,rot_type,face_pair_to_rot_idx,
             tet_pair_to_rot_idx,1))
          is_valid = false;
        g->save_state_mem(ses);
        ss.push_back(ses);
      }
    if(is_valid)
      is_valid = g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                       surface_idx_to_rot_idx,1);
    if(!is_valid){
        cerr << "# [error] strange invalid inner face setting." << endl;
        return __LINE__;
      }
    cerr << "# [info] finish polycube_zyz setting." << endl;

  }

  while(fi < fn) {
      is_valid = true;

      assert(g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                   surface_idx_to_rot_idx,1));

      size_t fii = fi;
      long_number_add(check_numbers,1);
      for(; fii < fi + windows_size; ++fii){
          //cerr << fii << " " << fn << endl;
          tried_type_idx[fii] += 1;

          assert(tried_type_idx[fii] < 24);

          const size_t &next_type = type_candidates(tried_type_idx[fii], fii);
          rot_type[fii] = next_type;
          tr.start();
          if(is_surface[fii]) { // the third type = -1 means this face is surface
              g->insert_surface_transition_group(
                    next_type, fa_cut, cut_tet2tet,
                    surface_groups[fii - jump_face_groups.size()], rot_type,
                  tet_pair_to_rot_idx);
            }else {
              if(g->insert_transition_group(
                   jump_face_groups[fii],cut_tet,cut_tet2tet,
                   fa_cut, ortae,rot_type,face_pair_to_rot_idx,
                   tet_pair_to_rot_idx,1))
                is_valid = false;
            }
          tr.finish();
          trans_time += tr.result_c();
          tr.start();
          g->save_state_mem(ses);
          ss.push_back(ses);
          tr.finish();
          write_time += tr.result_c();
          if(!is_valid) break;
        }

      if(is_valid){
          if(is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
            is_valid = false;
        }

      tr.start();
      if(is_valid)
        is_valid = g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                         surface_idx_to_rot_idx,1);
      tr.finish();
      check_time += tr.result_c();

      if(is_valid){
          fi += windows_size;
          windows_size = ((fn - fi)/2 > 0?(fn-fi)/2:1);
          continue;
        }

      {
        assert(!is_valid);
        // recover the states
        size_t begin = fi;
        // to ensure the range is [begin, end),
        if(fii != fi + windows_size)
          ++fii;
        size_t end = fii;

        vector<size_t> rot_type_bkp;
        while(1){
            size_t mid_face_idx = (begin + end)/2;
            tr.start();
            g->load_state_mem(ss[mid_face_idx]);
            tr.finish();
            read_time += tr.result_c();
            rot_type_bkp = rot_type;
            if(fii > mid_face_idx + 1)
              std::fill(rot_type_bkp.begin() + mid_face_idx+1,
                        rot_type_bkp.begin() + fii,  -1);
            is_valid = true;
            tr.start();
            is_valid = g->is_valid_with_info(cut_tet2tet, rot_type_bkp, ortae,
                                             surface_idx_to_rot_idx,1);
            if(is_valid){
                if(is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
                  is_valid = false;
              }
            tr.finish();
            check_time += tr.result_c();
            if(is_valid){
                if(mid_face_idx == begin){
                    rot_type = rot_type_bkp;
                    break;
                  }
                begin = mid_face_idx;
                continue;
              }else{
                if(fii > mid_face_idx + 1)
                  std::fill(rot_type.begin() + mid_face_idx+1,
                            rot_type.begin() + fii,  -1);
                if(end == mid_face_idx + 1)
                  end = mid_face_idx;
                else
                  end = mid_face_idx+1;
                if(end == begin){
                    rot_type[end] = -1;
                    tr.start();
                    if(end == 0){
                        g->clear_state_mem();
                      }else
                      g->load_state_mem(ss[end-1]);
                    tr.finish();
                    read_time += tr.result_c();
                    //cerr << "# [strange] it's supposed not to be here." << endl;
                    break;
                  }
              }
          }

        ss.resize(end);

        for(size_t ti = end + 1; ti < fii; ++ti)
          tried_type_idx[ti] = -1;

        fi = end;

        while(1){
            size_t ri;
            for(ri = tried_type_idx[fi]+1; ri < 24; ++ri){
                is_valid = true;
                if(is_surface[fi] && ri == 3) break;
                rot_type[fi] = type_candidates(ri, fi);
                tried_type_idx[fi] = ri;

                tr.start();
                if(is_surface[fi] &&
                   is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
                  continue;
                tr.finish();
                check_time += tr.result_c();
                tr.start();
                if(is_surface[fi]) { // the third type = -1 means this face is surface
                    g->insert_surface_transition_group(
                          rot_type[fi], fa_cut, cut_tet2tet,
                          surface_groups[fi - jump_face_groups.size()], rot_type,
                        tet_pair_to_rot_idx);
                  }else {
                    if(g->insert_transition_group(
                         jump_face_groups[fi],cut_tet,cut_tet2tet,
                         fa_cut, ortae,rot_type,face_pair_to_rot_idx,
                         tet_pair_to_rot_idx,1))
                      is_valid = false;
                  }
                tr.finish();
                trans_time += tr.result_c();
                tr.start();
                if(is_valid)
                  is_valid = g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                                   surface_idx_to_rot_idx,1);
                tr.finish();
                check_time += tr.result_c();
                if(is_valid){
                    tr.start();
                    g->save_state_mem(ses);
                    ss.push_back(ses);
                    tr.finish();
                    write_time += tr.result_c();
                    break;
                  }
                {
                  assert(!is_valid);
                  rot_type[fi] = -1;
                  tr.start();
                  if(ss.empty())
                    g->clear_state_mem();
                  else
                    g->load_state_mem(ss.back());

                  tr.finish();
                  read_time += tr.result_c();

                  assert(g->is_valid_with_info(cut_tet2tet,rot_type, ortae,
                                               surface_idx_to_rot_idx));
                }
              }
            if((!is_surface[fi] && ri == 24) || (is_surface[fi] && ri == 3)){
                if(fi == jump_face_groups.size()){
                    cerr << "# [error] cannot find a solution" << endl;
                    return __LINE__;
                  }

                tried_type_idx[fi] = -1;
                rot_type[fi] = -1;
                --fi;

                ss.pop_back();
                tr.start();
                if(ss.empty())
                  g->clear_state_mem();
                else
                  g->load_state_mem(ss.back());
                tr.finish();
                read_time += tr.result_c();
                rot_type[fi] = -1;
              }else{
                ++fi;
                windows_size = ((fn - fi)/2 > 0?(fn-fi)/2:1);
                break;
              }
          }
      }
    }

  cerr << "# write time " << write_time/1e6 << endl;
  cerr << "# read time " << read_time/1e6 << endl;
  cerr << "# trans time " << trans_time/1e6 << endl;
  cerr << "# check time " << check_time/1e6 << endl;
  cerr << "# searching steps " << searching_steps << endl;
  cerr << "# iterate numbers ";
  print_long_number(cerr, check_numbers);

  {
    cerr << "find a solution. " << endl;
    post_process(rot_type,cut_node, node, outside_face, fa_cut, ortae,
                 surface_idx_to_rot_idx, tet_pair_to_rot_idx, face_pair_to_rot_idx,
                 inner_face_jump_type);

    if(find(rot_type.begin(), rot_type.end(), -1) != rot_type.end()){
        cerr << "# [error] strange, there is rot type that has not been set" << endl;
        return __LINE__;
      }
    simple_check(rot_type, is_surface, jump_face_vec, cut_tet, cut_node, node,
                 fa_cut, cut_tet2tet, ortae, ortae_cut, outside_face_cut,
                 outside_face_idx_in_cut, face_pair_cut, surface_idx_to_rot_idx,
                 face_pair_to_rot_idx, tet_pair_to_rot_idx);
    return 0;
  }

  return 0;
}

int searching_strategy::searching_solutions_by_using_memory_binary_search_group(
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
    const boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const zjucad::matrix::matrix<matrixd> * frame_ptr)
{
  boost::unordered_map<size_t,size_t> surface_type_cut;
  matrixd cut_node;
  matrixst cut_tet2tet,outside_face, outside_face_cut,
      outside_face_idx_in_cut, outside_face_idx;
  vector<pair<size_t,size_t> > jump_face_vec;
  jtf::mesh::one_ring_tet_at_edge  ortae;
  jtf::mesh::one_ring_tet_at_edge ortae_cut;

  prev_process(
        orig_tet,node,cut_tet, fa, fa_cut,face_pair_cut, surface_type, cut_node,
        cut_tet2tet, outside_face, outside_face_cut, outside_face_idx_in_cut,
        outside_face_idx, jump_face_vec, surface_type_cut, ortae, ortae_cut);

  vector<vector<pair<size_t,size_t> > > jump_face_groups;
  vector<vector<size_t> > surface_groups;
  vector<pair<size_t,double> > group_type_with_frame_reliablity,
      group_type_with_surface_reliablity;
  const double frame_reliablity_threshold = 0.1;

  // this function is used to group faces according to the type and reliablity

  boost::unordered_map<pair<size_t,size_t>,size_t>  face_pair_to_rot_idx;
  boost::unordered_map<pair<size_t,size_t>,size_t>  tet_pair_to_rot_idx;
  boost::unordered_map<size_t,size_t>  surface_idx_to_rot_idx;
  matrixst type_candidates;
  vector<bool> is_surface;

  group_faces_according_face_reliablity(
        cut_tet,cut_node,fa_cut, outside_face_cut,outside_face_idx_in_cut,
        face_pair_cut, *frame_ptr, jump_face_groups, surface_groups,
        face_pair_to_rot_idx, tet_pair_to_rot_idx, surface_idx_to_rot_idx,
        type_candidates, is_surface, frame_reliablity_threshold);

  unique_ptr<singularity_graph> g(
        singularity_graph::create(
          ortae_cut, ortae, cut_tet2tet, outside_face_cut,  face_pair_cut,
          outside_face_idx_in_cut, surface_idx_to_rot_idx, face_pair_to_rot_idx));

  if(!g.get()){
      cerr << "# [info] build singularity graph failed." << endl;
      return __LINE__;
    }

  vector<size_t> rot_type(jump_face_groups.size() + surface_groups.size() ,-1);
  vector<size_t> tried_type_idx(rot_type.size(), -1);

  // start to check the question tree, and try to figure out the answer
  size_t searching_steps = 0;
  const size_t fn = rot_type.size();
  deque<state_each_step> ss;
  double write_time = 0.0, read_time = 0.0, check_time = 0.0, trans_time = 0.0;
  timer tr;
  bool is_valid;
  size_t windows_size = fn;
  state_each_step ses;
  deque<size_t> check_numbers;

#define dump_path
  size_t fi = 0;
  //vector<std::tuple<size_t,size_t,size_t> > time_line;
  //size_t step = -1;
  while(fi < fn) {
      //++step;
      is_valid = true;

      assert(g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                   surface_idx_to_rot_idx,1));

      size_t fii = fi;
      long_number_add(check_numbers,1);


      for(; fii < fi + windows_size; ++fii){
          //cerr << fii << " " << fn << endl;
          tried_type_idx[fii] += 1;

          //time_line.push_back(std::make_tuple(fii, tried_type_idx[fii], step) );

          assert(tried_type_idx[fii] < 24);

          const size_t &next_type = type_candidates(tried_type_idx[fii], fii);
          rot_type[fii] = next_type;
          tr.start();
          if(is_surface[fii]) { // the third type = -1 means this face is surface
              g->insert_surface_transition_group(
                    next_type, fa_cut, cut_tet2tet,
                    surface_groups[fii - jump_face_groups.size()], rot_type,
                  tet_pair_to_rot_idx);
            }else {
              if(g->insert_transition_group(
                   jump_face_groups[fii],cut_tet,cut_tet2tet,
                   fa_cut, ortae,rot_type,face_pair_to_rot_idx,
                   tet_pair_to_rot_idx,1))
                is_valid = false;
            }
          tr.finish();
          trans_time += tr.result_c();
          tr.start();
          g->save_state_mem(ses);
          ss.push_back(ses);
          tr.finish();
          write_time += tr.result_c();
          if(!is_valid) break;
        }

      if(is_valid){
          if(is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
            is_valid = false;
        }

      tr.start();
      if(is_valid)
        is_valid = g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                         surface_idx_to_rot_idx,1);
      tr.finish();
      check_time += tr.result_c();

      if(is_valid){
#ifdef dump_path
          {
            for(size_t i = fi ; i < fi + windows_size; ++i){
                if(is_surface[i]){
                    const vector<size_t> & vec = surface_groups[i - jump_face_groups.size()];
                    for(size_t fjj = 0; fjj < vec.size(); ++fjj){
                        cout << "+ " << vec[fjj] << " "
                             << type_candidates(tried_type_idx[i], i) << endl;
                      }
                  }else{
                    const vector<pair<size_t,size_t> > & vec = jump_face_groups[i];
                    for(size_t fjj = 0; fjj < vec.size(); ++fjj){
                        cout << "+ " << vec[fjj].first << " "
                             << type_candidates(tried_type_idx[i], i) << endl;
                      }
                  }
              }
          }
#endif
          fi += windows_size;
          windows_size = ((fn - fi)/2 > 0?(fn-fi)/2:1);

          continue;
        }

      {
        assert(!is_valid);
        // recover the states
        size_t begin = fi;
        // to ensure the range is [begin, end),
        if(fii != fi + windows_size)
          ++fii;
        size_t end = fii;

        vector<size_t> rot_type_bkp;
        while(1){
            size_t mid_face_idx = (begin + end)/2;
            tr.start();
            g->load_state_mem(ss[mid_face_idx]);
            tr.finish();
            read_time += tr.result_c();
            rot_type_bkp = rot_type;
            if(fii > mid_face_idx + 1)
              std::fill(rot_type_bkp.begin() + mid_face_idx+1,
                        rot_type_bkp.begin() + fii,  -1);
            is_valid = true;
            tr.start();
            is_valid = g->is_valid_with_info(cut_tet2tet, rot_type_bkp, ortae,
                                             surface_idx_to_rot_idx,1);
            if(is_valid){
                if(is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
                  is_valid = false;
              }
            tr.finish();
            check_time += tr.result_c();
            if(is_valid){
                if(mid_face_idx == begin){
                    rot_type = rot_type_bkp;
                    break;
                  }
                begin = mid_face_idx;
                continue;
              }else{
                if(fii > mid_face_idx + 1)
                  std::fill(rot_type.begin() + mid_face_idx+1,
                            rot_type.begin() + fii,  -1);
                if(end == mid_face_idx + 1)
                  end = mid_face_idx;
                else
                  end = mid_face_idx+1;
                if(end == begin){
                    rot_type[end] = -1;
                    tr.start();
                    if(end == 0){
                        g->clear_state_mem();
                      }else
                      g->load_state_mem(ss[end-1]);
                    tr.finish();
                    read_time += tr.result_c();
                    //cerr << "# [strange] it's supposed not to be here." << endl;
                    break;
                  }
              }
          }

        ss.resize(end);
#ifdef dump_path
        {
          for(size_t i = fi; i < end + 1; ++i){
              if(is_surface[i]){
                  const vector<size_t> & vec = surface_groups[i - jump_face_groups.size()];
                  for(size_t fjj = 0; fjj < vec.size(); ++fjj){
                      cout << "+ " << vec[fjj] << " "
                           << type_candidates(tried_type_idx[i], i) << endl;
                    }
                }else{
                  const vector<pair<size_t,size_t> > & vec = jump_face_groups[i];
                  for(size_t fjj = 0; fjj < vec.size(); ++fjj){
                      cout << "+ " << vec[fjj].first << " "
                           << type_candidates(tried_type_idx[i], i) << endl;
                    }
                }
            }
          if(is_surface[end])
            cout << "- " << surface_groups[end - jump_face_groups.size()].size() << endl;
          else
            cout << "- " << jump_face_groups[end].size() << endl;
        }
#endif
        for(size_t ti = end + 1; ti < fii; ++ti)
          tried_type_idx[ti] = -1;

        fi = end;

        while(1){
            size_t ri;
            for(ri = tried_type_idx[fi]+1; ri < 24; ++ri){
                // ++step;
                //time_line.push_back(std::make_tuple(fi, ri, step) );
                is_valid = true;
                if(is_surface[fi] && ri == 3) break;
                rot_type[fi] = type_candidates(ri, fi);
                tried_type_idx[fi] = ri;

                tr.start();
                if(is_surface[fi] &&
                   is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
                  continue;
                tr.finish();
                check_time += tr.result_c();
                tr.start();
                if(is_surface[fi]) { // the third type = -1 means this face is surface
                    g->insert_surface_transition_group(
                          rot_type[fi], fa_cut, cut_tet2tet,
                          surface_groups[fi - jump_face_groups.size()], rot_type,
                        tet_pair_to_rot_idx);
#ifdef dump_path
                    {
                      const vector<size_t> & vec =  surface_groups[fi - jump_face_groups.size()];
                      for(size_t fjj = 0; fjj < vec.size(); ++fjj){
                          cout << "+ " << vec[fjj] << " " << rot_type[fi] << endl;
                        }
                    }
#endif
                  }else {
                    if(g->insert_transition_group(
                         jump_face_groups[fi],cut_tet,cut_tet2tet,
                         fa_cut, ortae,rot_type,face_pair_to_rot_idx,
                         tet_pair_to_rot_idx,1))
                      is_valid = false;
#ifdef dump_path
                    {
                      const vector<pair<size_t,size_t> > & vec =  jump_face_groups[fi];
                      for(size_t fjj = 0; fjj < vec.size(); ++fjj){
                          cout << "+ " << vec[fjj].first << " " << rot_type[fi] << endl;
                        }
                    }
#endif
                  }
                tr.finish();
                trans_time += tr.result_c();
                tr.start();
                if(is_valid)
                  is_valid = g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                                   surface_idx_to_rot_idx,1);
                tr.finish();
                check_time += tr.result_c();
                if(is_valid){
                    tr.start();
                    g->save_state_mem(ses);
                    ss.push_back(ses);
                    tr.finish();
                    write_time += tr.result_c();
                    break;
                  }
                {
#ifdef dump_path
                  {
                    if(is_surface[fi]){
                        cout << "- " << surface_groups[fi - jump_face_groups.size()].size() << endl;
                      }else{
                        cout << "- " << jump_face_groups[fi].size() << endl;
                      }
                  }
#endif
                  assert(!is_valid);
                  rot_type[fi] = -1;
                  tr.start();
                  if(ss.empty())
                    g->clear_state_mem();
                  else
                    g->load_state_mem(ss.back());

                  tr.finish();
                  read_time += tr.result_c();

                  assert(g->is_valid_with_info(cut_tet2tet,rot_type, ortae,
                                               surface_idx_to_rot_idx));
                }
              }
            if((!is_surface[fi] && ri == 24) || (is_surface[fi] && ri == 3)){
                if(fi == 0){
                    cerr << "# [error] cannot find a solution" << endl;
                    return __LINE__;
                  }

                tried_type_idx[fi] = -1;
                rot_type[fi] = -1;
                --fi;
#ifdef dump_path
                {
                  if(is_surface[fi]){
                      cout << "- " << surface_groups[fi - jump_face_groups.size()].size() << endl;
                    }else{
                      cout << "- " << jump_face_groups[fi].size() << endl;
                    }
                }
#endif
                ss.pop_back();
                tr.start();
                if(ss.empty())
                  g->clear_state_mem();
                else
                  g->load_state_mem(ss.back());
                tr.finish();
                read_time += tr.result_c();
                rot_type[fi] = -1;
              }else{
                ++fi;
                windows_size = ((fn - fi)/2 > 0?(fn-fi)/2:1);
                break;
              }
          }
      }
    }

  // cerr << "[info] start to draw time lines" << endl;
  //  save_2d_point("time_line", "time_line",time_line);
  // save_3d_point("time_line", "time_line",time_line);
  //  {
  //    ofstream ofs("time_line");
  //    for(size_t pi = 0; pi < time_line.size(); ++pi){
  //      ofs << time_line[pi].first << " " << time_line[pi].second << endl;
  //    }
  //   }
  cerr << "# write time " << write_time/1e6 << endl;
  cerr << "# read time " << read_time/1e6 << endl;
  cerr << "# trans time " << trans_time/1e6 << endl;
  cerr << "# check time " << check_time/1e6 << endl;
  cerr << "# searching steps " << searching_steps << endl;
  cerr << "# iterate numbers ";
  print_long_number(cerr, check_numbers);

  {
    cerr << "find a solution. " << endl;
    post_process(rot_type,cut_node, node, outside_face, fa_cut, ortae,
                 surface_idx_to_rot_idx, tet_pair_to_rot_idx, face_pair_to_rot_idx,
                 inner_face_jump_type);

    if(find(rot_type.begin(), rot_type.end(), -1) != rot_type.end()){
        cerr << "# [error] strange, there is rot type that has not been set" << endl;
        return __LINE__;
      }
    simple_check(rot_type, is_surface, jump_face_vec, cut_tet, cut_node, node,
                 fa_cut, cut_tet2tet, ortae, ortae_cut, outside_face_cut,
                 outside_face_idx_in_cut, face_pair_cut, surface_idx_to_rot_idx,
                 face_pair_to_rot_idx, tet_pair_to_rot_idx);
    return 0;
  }

  return 0;
}

int searching_strategy::searching_solutions_by_using_memory_binary_search(
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
    const boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const zjucad::matrix::matrix<matrixd> * frame_ptr)
{
  boost::unordered_map<size_t,size_t> surface_type_cut;
  matrixd cut_node;
  matrixst cut_tet2tet,outside_face, outside_face_cut,outside_face_idx_in_cut, outside_face_idx;
  vector<pair<size_t,size_t> > jump_face_vec;
  jtf::mesh::one_ring_tet_at_edge  ortae;
  jtf::mesh::one_ring_tet_at_edge ortae_cut;

  prev_process(orig_tet,node,cut_tet, fa, fa_cut,face_pair_cut, surface_type, cut_node,
               cut_tet2tet, outside_face, outside_face_cut, outside_face_idx_in_cut, outside_face_idx,
               jump_face_vec, surface_type_cut, ortae, ortae_cut);

  vector<size_t> rot_type(jump_face_vec.size() + surface_type_cut.size(),-1);
  matrixst type_candidates(24,rot_type.size());
  vector<bool> is_surface(rot_type.size(), false);

  boost::unordered_map<size_t,size_t> surface_idx_to_rot_idx;
  boost::unordered_map<size_t,size_t> surface_idx_with_rot_idx;
  boost::unordered_map<pair<size_t,size_t>,size_t> face_pair_to_rot_idx;
  boost::unordered_map<pair<size_t,size_t>,size_t> tet_pair_to_rot_idx;
  boost::unordered_map<size_t,pair<size_t,size_t> > face_pair_with_rot_idx;

  adjust_face_and_type_order_to_speed_up(
        jump_face_vec, inner_face_jump_type, surface_type_cut, fa_cut, cut_tet2tet,
        type_candidates, surface_idx_to_rot_idx,surface_idx_with_rot_idx,
        face_pair_to_rot_idx,face_pair_with_rot_idx,tet_pair_to_rot_idx,
        is_surface,frame_ptr);

  unique_ptr<singularity_graph> g(
        singularity_graph::create(ortae_cut,ortae,cut_tet2tet,outside_face_cut,
                                  face_pair_cut,outside_face_idx_in_cut,
                                  surface_idx_to_rot_idx,face_pair_to_rot_idx));
  if(!g.get()){
      cerr << "# [error] can not build singularity graph." << endl;
      return __LINE__;
    }

  vector<size_t> tried_type_idx(rot_type.size(), -1);
  // start to check the question tree, and try to figure out the answer
  size_t fi = 0;
  size_t searching_steps = 0;
  const size_t fn = rot_type.size();
  deque<state_each_step> ss;
  double write_time = 0.0, read_time = 0.0, check_time = 0.0, trans_time = 0.0;
  timer tr;
  bool is_valid;
  size_t windows_size = fn;
  state_each_step ses;
  deque<size_t> check_numbers;
  while(fi < fn) {
      is_valid = true;

      assert(g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                   surface_idx_to_rot_idx,1));

      size_t fii = fi;
      long_number_add(check_numbers,1);
      for(; fii < fi + windows_size; ++fii){

          tried_type_idx[fii] += 1;

          assert(tried_type_idx[fii] < 24);

          const size_t &next_type = type_candidates(tried_type_idx[fii], fii);
          rot_type[fii] = next_type;
          tr.start();
          if(is_surface[fii]) { // the third type = -1 means this face is surface
              g->insert_surface_transition(next_type, fa_cut, cut_tet2tet,
                                           surface_idx_with_rot_idx[fii], rot_type,
                                           tet_pair_to_rot_idx);
            }else {
              if(g->insert_transition(face_pair_with_rot_idx[fii],cut_tet,cut_tet2tet,
                                      fa_cut, ortae,rot_type,face_pair_to_rot_idx,
                                      tet_pair_to_rot_idx,1)){
                  is_valid = false;
                }
            }
          tr.finish();
          trans_time += tr.result_c();
          tr.start();
          g->save_state_mem(ses);
          ss.push_back(ses);
          tr.finish();
          write_time += tr.result_c();
          if(!is_valid) break;
        }

      if(is_valid){
          if(is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
            is_valid = false;
        }

      tr.start();
      if(is_valid)
        is_valid = g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                         surface_idx_to_rot_idx,1);
      tr.finish();
      check_time += tr.result_c();

      if(is_valid){
          fi += windows_size;
          windows_size = ((fn - fi)/2 > 0?(fn-fi)/2:1);
          continue;
        }

      {
        assert(!is_valid);
        // recover the states
        size_t begin = fi;
        // to ensure the range is [begin, end),
        if(fii != fi + windows_size)
          ++fii;
        size_t end = fii;

        vector<size_t> rot_type_bkp;
        while(1){
            size_t mid_face_idx = (begin + end)/2;
            tr.start();
            g->load_state_mem(ss[mid_face_idx]);
            tr.finish();
            read_time += tr.result_c();
            rot_type_bkp = rot_type;
            if(fii > mid_face_idx + 1)
              std::fill(rot_type_bkp.begin() + mid_face_idx+1,
                        rot_type_bkp.begin() + fii,  -1);
            is_valid = true;
            tr.start();
            is_valid = g->is_valid_with_info(cut_tet2tet, rot_type_bkp, ortae,
                                             surface_idx_to_rot_idx,1);
            if(is_valid){
                if(is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
                  is_valid = false;
              }
            tr.finish();
            check_time += tr.result_c();
            if(is_valid){
                if(mid_face_idx == begin){
                    rot_type = rot_type_bkp;
                    break;
                  }
                begin = mid_face_idx;
                continue;
              }else{
                if(fii > mid_face_idx + 1)
                  std::fill(rot_type.begin() + mid_face_idx+1,
                            rot_type.begin() + fii,  -1);
                if(end == mid_face_idx + 1)
                  end = mid_face_idx;
                else
                  end = mid_face_idx+1;
                if(end == begin){
                    rot_type[end] = -1;
                    tr.start();
                    if(end == 0){
                        g->clear_state_mem();
                      }else
                      g->load_state_mem(ss[end-1]);
                    tr.finish();
                    read_time += tr.result_c();
                    //cerr << "# [strange] it's supposed not to be here." << endl;
                    break;
                  }
              }
          }

        ss.resize(end);

        for(size_t ti = end + 1; ti < fii; ++ti)
          tried_type_idx[ti] = -1;

        fi = end;

        while(1){
            size_t ri;
            for(ri = tried_type_idx[fi]+1; ri < 24; ++ri){
                is_valid = true;
                if(is_surface[fi] && ri == 3) break;
                rot_type[fi] = type_candidates(ri, fi);
                tried_type_idx[fi] = ri;

                tr.start();
                if(is_surface[fi] &&
                   is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
                  continue;
                tr.finish();
                check_time += tr.result_c();
                tr.start();
                if(is_surface[fi]) { // the third type = -1 means this face is surface
                    g->insert_surface_transition(rot_type[fi], fa_cut, cut_tet2tet,
                                                 surface_idx_with_rot_idx[fi], rot_type,
                                                 tet_pair_to_rot_idx);
                  }else {
                    if(g->insert_transition(face_pair_with_rot_idx[fi],cut_tet,cut_tet2tet,
                                            fa_cut, ortae,rot_type,face_pair_to_rot_idx,
                                            tet_pair_to_rot_idx,1))
                      is_valid = false;
                  }
                tr.finish();
                trans_time += tr.result_c();
                tr.start();
                if(is_valid)
                  is_valid = g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                                   surface_idx_to_rot_idx,1);
                tr.finish();
                check_time += tr.result_c();
                if(is_valid){
                    tr.start();
                    g->save_state_mem(ses);
                    ss.push_back(ses);
                    tr.finish();
                    write_time += tr.result_c();
                    break;
                  }
                {
                  assert(!is_valid);
                  rot_type[fi] = -1;
                  tr.start();
                  if(ss.empty())
                    g->clear_state_mem();
                  else
                    g->load_state_mem(ss.back());

                  tr.finish();
                  read_time += tr.result_c();

                  assert(g->is_valid_with_info(cut_tet2tet,rot_type, ortae,
                                               surface_idx_to_rot_idx));
                }
              }
            if((!is_surface[fi] && ri == 24) || (is_surface[fi] && ri == 3)){
                if(fi == 0){
                    cerr << "# [error] cannot find a solution" << endl;
                    return __LINE__;
                  }

                tried_type_idx[fi] = -1;
                rot_type[fi] = -1;
                --fi;

                ss.pop_back();
                tr.start();
                if(ss.empty())
                  g->clear_state_mem();
                else
                  g->load_state_mem(ss.back());
                tr.finish();
                read_time += tr.result_c();
                rot_type[fi] = -1;
              }else{
                ++fi;
                windows_size = ((fn - fi)/2 > 0?(fn-fi)/2:1);
                break;
              }
          }
      }
    }

  cerr << "# write time " << write_time/1e6 << endl;
  cerr << "# read time " << read_time/1e6 << endl;
  cerr << "# trans time " << trans_time/1e6 << endl;
  cerr << "# check time " << check_time/1e6 << endl;
  cerr << "# searching steps " << searching_steps << endl;
  cerr << "# iterate numbers ";
  print_long_number(cerr, check_numbers);

  {
    cerr << "find a solution. " << endl;
    post_process(rot_type,cut_node, node, outside_face, fa_cut, ortae,
                 surface_idx_to_rot_idx, tet_pair_to_rot_idx, face_pair_to_rot_idx,
                 inner_face_jump_type);

    assert(find(rot_type.begin(), rot_type.end(), -1) == rot_type.end());
    simple_check(rot_type, is_surface, jump_face_vec, cut_tet, cut_node, node,
                 fa_cut, cut_tet2tet, ortae, ortae_cut, outside_face_cut,
                 outside_face_idx_in_cut, face_pair_cut, surface_idx_to_rot_idx,
                 face_pair_to_rot_idx, tet_pair_to_rot_idx);
    return 0;
  }
  return 0;
}


int searching_strategy::searching_solutions_by_minimal_error(
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const matrixst & face_pair_cut,
    const matrix<matrixd> & frame,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type)
{
  boost::unordered_map<size_t,size_t> surface_type_cut;
  matrixd cut_node;
  matrixst cut_tet2tet,outside_face, outside_face_cut,outside_face_cut_idx, outside_face_idx;
  boost::unordered_map<size_t,size_t> surface_type;
  vector<pair<size_t,size_t> > jump_face_vec;
  jtf::mesh::one_ring_tet_at_edge  ortae;
  jtf::mesh::one_ring_tet_at_edge ortae_cut;
  prev_process(orig_tet,node,cut_tet, fa, fa_cut,face_pair_cut, surface_type, cut_node,
               cut_tet2tet, outside_face, outside_face_cut, outside_face_cut_idx, outside_face_idx,
               jump_face_vec, surface_type_cut, ortae, ortae_cut);


  assert(frame.size() == cut_tet.size(2));
  size_t surface_num = 0;
  for(size_t fi = 0; fi < face_pair_cut.size(); ++fi)
    if(face_pair_cut[fi] == -1)   ++surface_num;

  const size_t total_face_num = jump_face_vec.size() + surface_num;
  vector<deque<pair<double,size_t> > > diff(total_face_num);
  vector<bool> is_surface(total_face_num, false);

  boost::unordered_map<size_t,size_t> surface_idx_to_rot_idx;
  boost::unordered_map<size_t,size_t> surface_idx_with_rot_idx;
  boost::unordered_map<std::pair<size_t,size_t>,size_t> face_pair_to_rot_idx;
  boost::unordered_map<size_t,std::pair<size_t,size_t> > face_pair_with_rot_idx;
  boost::unordered_map<std::pair<size_t,size_t>,size_t> tet_pair_to_rot_idx;

  // set the inner face jump type
  for(size_t fi = 0; fi < jump_face_vec.size(); ++fi){
      const pair<size_t,size_t> & face_pair = jump_face_vec[fi];
      const pair<size_t,size_t> & tet_pair_0 = fa_cut.face2tet_[face_pair.first];
      const pair<size_t,size_t> & tet_pair_1 = fa_cut.face2tet_[face_pair.second];

      pair<size_t,size_t> tet_pair(
            tet_pair_0.first == -1?tet_pair_0.second:tet_pair_0.first,
            tet_pair_1.first == -1?tet_pair_1.second:tet_pair_1.first);

      assert(fa_cut.is_outside_face(tet_pair_0));
      assert(fa_cut.is_outside_face(tet_pair_1));

      face_pair_to_rot_idx[face_pair] = fi;
      face_pair_with_rot_idx[fi] = face_pair;
      tet_pair_to_rot_idx[tet_pair] = fi;

      for(size_t ri = 0; ri < 24; ++ri){
          diff[fi].push_back(
                make_pair(norm(frame[tet_pair.first] * type_transition2(ri)
                          - frame[tet_pair.second]), ri));
        }
      diff[fi].push_back(make_pair(std::numeric_limits<double>::infinity(),-1));
      sort(diff[fi].begin(), diff[fi].end());
    }

  // set the surface type
  size_t fi = jump_face_vec.size();
  matrixd face_normal = zeros<double>(3,1);
  vector<pair<double,size_t> > axis_difference(3);
  for(size_t fii = 0; fii < face_pair_cut.size(); ++fii){
      if(face_pair_cut[fii] == -1){
          is_surface[fi] = true;
          surface_idx_to_rot_idx[outside_face_cut_idx[fii]] = fi;
          surface_idx_with_rot_idx[fi] = outside_face_cut_idx[fii];
          jtf::mesh::cal_face_normal(outside_face_cut(colon(),fii),
                                     cut_node,face_normal);
          jtf::tetmesh::orient_one_face_normal_outside_tetmesh(
                cut_tet, cut_node,outside_face_cut(colon(),fii),
                outside_face_cut_idx[fii], fa_cut,face_normal);

          const pair<size_t,size_t> &tet_pair =
              fa_cut.face2tet_[outside_face_cut_idx[fii]];
          assert(fa_cut.is_outside_face(tet_pair));
          const size_t & tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);

          for(size_t ai = 0; ai < 3; ++ai){
              axis_difference[ai] = make_pair(dot(frame[tet_idx](colon(),ai),face_normal),ai);
              const double dd = dot(-1 * frame[tet_idx](colon(),ai),face_normal);
              if(dd > axis_difference[ai].first) {
                  // this direction is more suitable, but in diff we should record the
                  // difference, the smaller the better, so we should use cross to
                  // choose it
                  axis_difference[ai] =
                      make_pair(1.0 - dd,ai);
                }else{
                  axis_difference[ai] =
                      make_pair(1.0 - axis_difference[ai].first,ai);
                }
            }
          sort(axis_difference.begin(), axis_difference.end());

          for(size_t adi = 0; adi < 3; ++adi)
            diff[fi].push_back(axis_difference[adi]);
          diff[fi].push_back(make_pair(std::numeric_limits<double>::infinity(),-1));
          ++fi;
        }
    }

  unique_ptr<singularity_graph> g(
        singularity_graph::create(
          ortae_cut, ortae, cut_tet2tet, outside_face_cut,  face_pair_cut,
          outside_face_cut_idx, surface_idx_to_rot_idx, face_pair_to_rot_idx));

  //////////////////////////////////////////////////////////////////////////////
  {// debug
    vector<size_t> outside_face_vec;
    vector<size_t> outside_face_type;
    for(boost::unordered_map<size_t,size_t>::const_iterator bumcit =
        surface_idx_to_rot_idx.begin(); bumcit != surface_idx_to_rot_idx.end();
        ++bumcit){
        const size_t &face_idx = bumcit->first;
        const size_t &face_type = diff[bumcit->second].front().second;

        const vector<size_t> &face_vec= fa_cut.faces_[face_idx];
        outside_face_vec.insert(outside_face_vec.end(), face_vec.begin(), face_vec.end());
        outside_face_type.push_back(face_type);
      }
    ofstream ofs("surface_type_in_min_error.vtk");
    tri2vtk(ofs, &cut_node[0],cut_node.size(2), &outside_face_vec[0],
        outside_face_vec.size()/3);
    cell_data(ofs, &outside_face_type[0], outside_face_type.size(), "idx");
  }

  //////////////////////////////////////////////////////////////////////////////

  vector<bool> visited_face_flag(total_face_num, false);
  vector<size_t> rot_type(total_face_num, -1);
  deque<size_t> check_numbers;
  vector<state_each_step> ss;
  ss.reserve(total_face_num);
  state_each_step ses;
  vector<pair<size_t,size_t> > face_type_stack;
  face_type_stack.reserve(total_face_num);
  size_t searching_steps = 0;
  timer tr;
  double write_time = 0.0, read_time = 0.0, trans_time = 0.0, check_time = 0.0;
  while(face_type_stack.size() < total_face_num ){
      ++searching_steps;
      // difference, face_idx, type_idx
      vector<std::tuple<double,size_t,size_t> > candidate_faces;
      for(size_t fi = 0; fi < total_face_num; ++fi) {
          if(!visited_face_flag[fi]){
              candidate_faces.push_back(
                    std::make_tuple(diff[fi].front().first,fi,
                                    diff[fi].front().second));
            }
        }
      assert(candidate_faces.size() + face_type_stack.size() == total_face_num);

      const std::tuple<double,size_t,size_t> & best_choose
          = *(min_element(candidate_faces.begin(), candidate_faces.end()));

      face_type_stack.push_back(make_pair(get<1>(best_choose), get<2>(best_choose)));
      tr.start();
      g->save_state_mem(ses);
      ss.push_back(ses);
      tr.finish();
      write_time += tr.result_c();

      rot_type[get<1>(best_choose)] = get<2>(best_choose);
      long_number_add(check_numbers, 1);
      bool is_valid = true;
      tr.start();
      if(is_surface[get<1>(best_choose)]) { // the third type = -1 means this face is surface
          g->insert_surface_transition(
                get<2>(best_choose), fa_cut, cut_tet2tet,
                surface_idx_with_rot_idx[get<1>(best_choose)], rot_type,
              tet_pair_to_rot_idx);
        }else {
          if(g->insert_transition(
               jump_face_vec[get<1>(best_choose)],cut_tet,cut_tet2tet,
               fa_cut, ortae,rot_type,face_pair_to_rot_idx,
               tet_pair_to_rot_idx,1))
            is_valid = false;
        }
      tr.finish();
      trans_time += tr.result_c();

      tr.start();
      if(is_valid){
          is_valid = g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                           surface_idx_to_rot_idx,1);
        }
      tr.finish();
      check_time += tr.result_c();

      if(is_valid){
          visited_face_flag[get<1>(best_choose)] = true;
          continue;
        }


      {
        assert(is_valid == false);
        // have tested all candidates of this face, need to go back
        while(!face_type_stack.empty()){
            const pair<size_t,size_t> & latest_choose = face_type_stack.back();
            if(diff[latest_choose.first][1].second == -1){
                // if the next type choose of face latest_choose.first
                // ---diff[latest_choose.first][1].second--- is -1, which means all type
                // of this face has been tried, need to go back

                sort(diff[latest_choose.first].begin(), diff[latest_choose.first].end());
                visited_face_flag[latest_choose.first] = false;
                rot_type[latest_choose.first] = -1;
                face_type_stack.pop_back();
                ss.pop_back();
              }else{
                diff[latest_choose.first].push_back(diff[latest_choose.first].front());
                diff[latest_choose.first].pop_front();
                visited_face_flag[latest_choose.first] = false;
                rot_type[latest_choose.first] = -1;
                face_type_stack.pop_back();
                tr.start();
                g->load_state_mem(ss.back());
                tr.finish();
                read_time += tr.result_c();

                ss.pop_back();
                break;
              }
          }
        if(face_type_stack.empty())
          break;
      }
    }

  {
    ofstream ofs("debug_face_type");
    for(size_t ti = 0; ti < face_type_stack.size(); ++ti)
      ofs << face_type_stack[ti].first << " "
          << face_type_stack[ti].second << endl;
  }
  cerr << "# write time " << write_time/1e6 << endl;
  cerr << "# read time " << read_time/1e6 << endl;
  cerr << "# trans time " << trans_time/1e6 << endl;
  cerr << "# check time " << check_time/1e6 << endl;
  cerr << "# searching steps " << searching_steps << endl;
  cerr << "# iterate numbers ";
  print_long_number(cerr, check_numbers);

  if(face_type_stack.empty())
    cerr << "# [error] can not find a solution." << endl;
  else{
      cerr << "# [info] find a solution." << endl;

      post_process(rot_type, cut_node, node, outside_face, fa_cut,ortae,
                   surface_idx_to_rot_idx,tet_pair_to_rot_idx,
                   face_pair_to_rot_idx, inner_face_jump_type);

      simple_check(
            rot_type, is_surface, jump_face_vec, cut_tet, cut_node, node,
            fa_cut, cut_tet2tet, ortae, ortae_cut, outside_face_cut,
            outside_face_cut_idx, face_pair_cut, surface_idx_to_rot_idx,
            face_pair_to_rot_idx, tet_pair_to_rot_idx);
    }
  return 0;
}


int searching_strategy::searching_solutions_by_minimal_error_binary_search(
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const matrixst & face_pair_cut,
    const matrix<matrixd> & frame,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type)
{
  boost::unordered_map<size_t,size_t> surface_type_cut;
  matrixd cut_node;
  matrixst cut_tet2tet,outside_face, outside_face_cut,outside_face_cut_idx, outside_face_idx;
  boost::unordered_map<size_t,size_t> surface_type;
  vector<pair<size_t,size_t> > jump_face_vec;
  jtf::mesh::one_ring_tet_at_edge  ortae;
  jtf::mesh::one_ring_tet_at_edge ortae_cut;
  prev_process(orig_tet,node,cut_tet, fa, fa_cut,face_pair_cut, surface_type, cut_node,
               cut_tet2tet, outside_face, outside_face_cut, outside_face_cut_idx, outside_face_idx,
               jump_face_vec, surface_type_cut, ortae, ortae_cut);

  assert(frame.size() == cut_tet.size(2));
  size_t surface_num = 0;
  for(size_t fi = 0; fi < face_pair_cut.size(); ++fi)
    if(face_pair_cut[fi] == -1)   ++surface_num;

  const size_t total_face_num = jump_face_vec.size() + surface_num;
  vector<deque<pair<double,size_t> > > diff(total_face_num);
  vector<bool> is_surface(total_face_num, false);

  boost::unordered_map<size_t,size_t> surface_idx_to_rot_idx;
  boost::unordered_map<size_t,size_t> surface_idx_with_rot_idx;
  boost::unordered_map<std::pair<size_t,size_t>,size_t> face_pair_to_rot_idx;
  boost::unordered_map<size_t,std::pair<size_t,size_t> > face_pair_with_rot_idx;
  boost::unordered_map<std::pair<size_t,size_t>,size_t> tet_pair_to_rot_idx;

  // set the inner face jump type
  for(size_t fi = 0; fi < jump_face_vec.size(); ++fi){
      const pair<size_t,size_t> & face_pair = jump_face_vec[fi];
      const pair<size_t,size_t> & tet_pair_0 = fa_cut.face2tet_[face_pair.first];
      const pair<size_t,size_t> & tet_pair_1 = fa_cut.face2tet_[face_pair.second];

      pair<size_t,size_t> tet_pair(
            tet_pair_0.first == -1?tet_pair_0.second:tet_pair_0.first,
            tet_pair_1.first == -1?tet_pair_1.second:tet_pair_1.first);

      assert(fa_cut.is_outside_face(tet_pair_0));
      assert(fa_cut.is_outside_face(tet_pair_1));

      face_pair_to_rot_idx[face_pair] = fi;
      face_pair_with_rot_idx[fi] = face_pair;
      tet_pair_to_rot_idx[tet_pair] = fi;

      for(size_t ri = 0; ri < 24; ++ri){
          diff[fi].push_back(
                make_pair(norm(frame[tet_pair.first] * type_transition2(ri)
                          - frame[tet_pair.second]), ri));
        }
      diff[fi].push_back(make_pair(std::numeric_limits<double>::infinity(),-1));
      sort(diff[fi].begin(), diff[fi].end());
    }

  // set the surface type
  size_t fi = jump_face_vec.size();
  matrixd face_normal = zeros<double>(3,1);
  vector<pair<double,size_t> > axis_difference(3);
  for(size_t fii = 0; fii < face_pair_cut.size(); ++fii){
      if(face_pair_cut[fii] == -1){
          is_surface[fi] = true;
          surface_idx_to_rot_idx[outside_face_cut_idx[fii]] = fi;
          surface_idx_with_rot_idx[fi] = outside_face_cut_idx[fii];
          jtf::mesh::cal_face_normal(outside_face_cut(colon(),fii),
                                     cut_node,face_normal);
          jtf::tetmesh::orient_one_face_normal_outside_tetmesh(
                cut_tet, cut_node,outside_face_cut(colon(),fii),
                outside_face_cut_idx[fii], fa_cut,face_normal);

          const pair<size_t,size_t> &tet_pair =
              fa_cut.face2tet_[outside_face_cut_idx[fii]];
          assert(fa_cut.is_outside_face(tet_pair));
          const size_t & tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);

          for(size_t ai = 0; ai < 3; ++ai){
              axis_difference[ai] = make_pair(dot(frame[tet_idx](colon(),ai),face_normal),ai);
              const double dd = dot(-1 * frame[tet_idx](colon(),ai),face_normal);
              if(dd > axis_difference[ai].first) {
                  // this direction is more suitable, but in diff we should record the
                  // difference, the smaller the better, so we should use cross to
                  // choose it
                  axis_difference[ai] =
                      make_pair(1.0 - dd,ai);
                }else{
                  axis_difference[ai] =
                      make_pair(1.0 - axis_difference[ai].first,ai);
                }
            }
          sort(axis_difference.begin(), axis_difference.end());

          for(size_t adi = 0; adi < 3; ++adi)
            diff[fi].push_back(axis_difference[adi]);
          diff[fi].push_back(make_pair(std::numeric_limits<double>::infinity(),-1));
          ++fi;
        }
    }

  unique_ptr<singularity_graph> g(
        singularity_graph::create(
          ortae_cut, ortae, cut_tet2tet, outside_face_cut,  face_pair_cut,
          outside_face_cut_idx, surface_idx_to_rot_idx, face_pair_to_rot_idx));
  if(!g.get()){
      cerr << "# [error] can not build singularity graph" << endl;
      return __LINE__;
    }

  //////////////////////////////////////////////////////////////////////////////
  {// debug
    vector<size_t> outside_face_vec;
    vector<size_t> outside_face_type;
    for(boost::unordered_map<size_t,size_t>::const_iterator bumcit =
        surface_idx_to_rot_idx.begin(); bumcit != surface_idx_to_rot_idx.end();
        ++bumcit){
        const size_t &face_idx = bumcit->first;
        const size_t &face_type = diff[bumcit->second].front().second;

        const vector<size_t> &face_vec= fa_cut.faces_[face_idx];
        outside_face_vec.insert(outside_face_vec.end(), face_vec.begin(), face_vec.end());
        outside_face_type.push_back(face_type);
      }
    ofstream ofs("surface_type_in_min_error.vtk");
    tri2vtk(ofs, &cut_node[0],cut_node.size(2), &outside_face_vec[0],
        outside_face_vec.size()/3);
    cell_data(ofs, &outside_face_type[0], outside_face_type.size(), "idx");
  }

  //////////////////////////////////////////////////////////////////////////////

  vector<bool> visited_face_flag(total_face_num, false);
  vector<size_t> rot_type(total_face_num, -1);
  deque<size_t> check_number;
  vector<state_each_step> ss;
  ss.reserve(total_face_num);
  state_each_step ses;
  vector<pair<size_t,size_t> > face_type_stack;
  face_type_stack.reserve(total_face_num);
  size_t searching_steps = 0;
  timer tr;
  size_t windows_size = total_face_num;
  //set<size_t> invalid_face;
  double write_time = 0.0, read_time = 0.0, trans_time = 0.0, check_time = 0.0;
  while(face_type_stack.size() < total_face_num ){
      ++searching_steps;

      // difference, face_idx, type_idx
      size_t begin_ = face_type_stack.size();

      vector<std::tuple<double,size_t,size_t> > candidate_faces;
      for(size_t fi = 0; fi < total_face_num; ++fi) {
          if(!visited_face_flag[fi]){
              candidate_faces.push_back(
                    std::make_tuple(diff[fi].front().first,fi,
                                    diff[fi].front().second));
            }
        }

      sort(candidate_faces.begin(), candidate_faces.end());
      size_t ci = 0;
      bool is_valid = true;
      long_number_add(check_number,1);
      for(; ci < candidate_faces.size(); ++ci){
          face_type_stack.push_back(make_pair(get<1>(candidate_faces[ci]),
                                              get<2>(candidate_faces[ci])));

          rot_type[get<1>(candidate_faces[ci])] = get<2>(candidate_faces[ci]);
          visited_face_flag[get<1>(candidate_faces[ci])] = true;
          tr.start();
          if(is_surface[get<1>(candidate_faces[ci])]) {
              g->insert_surface_transition(
                    get<2>(candidate_faces[ci]), fa_cut, cut_tet2tet,
                    surface_idx_with_rot_idx[get<1>(candidate_faces[ci])], rot_type,
                  tet_pair_to_rot_idx);
            }else {
              if(g->insert_transition(
                   jump_face_vec[get<1>(candidate_faces[ci])],cut_tet,cut_tet2tet,
                   fa_cut, ortae,rot_type,face_pair_to_rot_idx,
                   tet_pair_to_rot_idx,1))
                is_valid = false;
            }
          tr.finish();
          trans_time += tr.result_c();

          tr.start();
          g->save_state_mem(ses);
          ss.push_back(ses);
          tr.finish();
          write_time += tr.result_c();

          if(!is_valid) break;
        }

      if(is_valid)
        is_valid = g->is_valid_with_info(cut_tet2tet, rot_type, ortae,
                                         surface_idx_to_rot_idx,1);

      if(is_valid){
          //      for(size_t ti = 0; ti < candidate_faces.size(); ++ti){
          //        assert(!visited_face_flag[candidate_faces[ti].get<1>()]);
          //        visited_face_flag[candidate_faces[ti].get<1>()] = true;
          //      }
          //invalid_face.clear();
          windows_size = (total_face_num - face_type_stack.size())/2>0?
                (total_face_num - face_type_stack.size())/2:1;
          continue;
        }

      { // check the graph
        if(ci != candidate_faces.size())
          ++ci; // to ensure the [) range

        vector<size_t> rot_type_bkp;
        vector<bool> visited_flag_bkp;
        size_t check_begin = begin_;
        size_t check_end = ci+begin_;
        bool is_valid;
        while(1){
            size_t mid_face_idx = (check_begin + check_end)/2;
            g->load_state_mem(ss[mid_face_idx]);
            rot_type_bkp = rot_type;
            visited_flag_bkp = visited_face_flag;
            // recovery the rot_type
            for(size_t ti = mid_face_idx + 1; ti < ci+begin_; ++ti){
                rot_type_bkp[get<1>(candidate_faces[ti-begin_])] = -1;
                visited_flag_bkp[get<1>(candidate_faces[ti-begin_])] = false;
              }
            is_valid = true;
            is_valid = g->is_valid_with_info(cut_tet2tet, rot_type_bkp, ortae,
                                             surface_idx_to_rot_idx,1);
            if(is_valid){
                if(is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
                  is_valid = false;
              }
            if(is_valid){
                if(mid_face_idx == check_begin){
                    rot_type = rot_type_bkp;
                    visited_face_flag = visited_flag_bkp;
                    break;
                  }
                check_begin = mid_face_idx;
                continue;
              }else{
                rot_type = rot_type_bkp;
                visited_face_flag = visited_flag_bkp;
                if(check_end == mid_face_idx+1)
                  check_end = mid_face_idx;
                else
                  check_end = mid_face_idx + 1;
                if(check_end == check_begin){
                    rot_type[get<1>(candidate_faces[check_end - begin_])] = -1;
                    visited_face_flag[get<1>(candidate_faces[check_end - begin_])] = false;

                    g->load_state_mem(ss[check_end-1]);
                    break;
                  }
              }
          }

        face_type_stack.resize(check_end+1);
        ss.resize(check_end);

        //      for(size_t ti = check_end + 1; ti < ci+begin_; ++ti){
        //        deque<pair<double,size_t> > & one_deque =
        //            diff[candidate_faces[ti-begin_].get<1>()];
        //        if(invalid_face.find(candidate_faces[ti-begin_].get<1>()) ==
        //           invalid_face.end()){
        //          sort(one_deque.begin(), one_deque.end());
        //        }else{
        //          one_deque.push_back(one_deque.front());
        //          one_deque.pop_front();
        //          if(one_deque.front().second == -1)
        //            sort(one_deque.begin(), one_deque.end());
        //        }
        //      }

      }

      {
        // have tested all candidates of this face, need to go back
        while(!face_type_stack.empty()){
            const pair<size_t,size_t> & latest_choose = face_type_stack.back();
            if(diff[latest_choose.first][1].second == -1){
                // if the next type choose of face latest_choose.first
                // ---diff[latest_choose.first][1].second--- is -1, which means all type
                // of this face has been tried, need to go back

                sort(diff[latest_choose.first].begin(), diff[latest_choose.first].end());
                visited_face_flag[latest_choose.first] = false;
                rot_type[latest_choose.first] = -1;
                face_type_stack.pop_back();
                ss.pop_back();
              }else{
                //invalid_face.insert(latest_choose.first);
                diff[latest_choose.first].push_back(diff[latest_choose.first].front());
                diff[latest_choose.first].pop_front();
                visited_face_flag[latest_choose.first] = false;
                rot_type[latest_choose.first] = -1;
                face_type_stack.pop_back();
                tr.start();
                g->load_state_mem(ss.back());
                tr.finish();
                read_time += tr.result_c();

                //ss.pop_back();
                break;
              }
          }
        if(face_type_stack.empty())
          break;
      }
    }


  cerr << "# write time " << write_time/1e6 << endl;
  cerr << "# read time " << read_time/1e6 << endl;
  cerr << "# trans time " << trans_time/1e6 << endl;
  cerr << "# check time " << check_time/1e6 << endl;
  cerr << "# searching steps " << searching_steps << endl;
  cerr << "# iterate numbers " ;
  print_long_number(cerr, check_number);

  if(face_type_stack.empty())
    cerr << "# [error] can not find a solution." << endl;
  else{
      cerr << "# [info] find a solution." << endl;

      post_process(rot_type, cut_node, node, outside_face, fa_cut,ortae,
                   surface_idx_to_rot_idx,tet_pair_to_rot_idx,
                   face_pair_to_rot_idx, inner_face_jump_type);

      simple_check(
            rot_type, is_surface, jump_face_vec, cut_tet, cut_node, node,
            fa_cut, cut_tet2tet, ortae, ortae_cut, outside_face_cut,
            outside_face_cut_idx, face_pair_cut, surface_idx_to_rot_idx,
            face_pair_to_rot_idx, tet_pair_to_rot_idx);
    }
  return 0;

}

int searching_strategy::post_process(
    const vector<size_t> & rot_type,
    const matrixd & cut_node,
    const matrixd & node,
    const matrixst &outside_face,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair_to_rot_idx,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_pair_to_rot_idx,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_jump_type)
{
  inner_jump_type.clear();
  assert(find(rot_type.begin(), rot_type.end(), -1) == rot_type.end());
  for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator
      bumpscit = tet_pair_to_rot_idx.begin(); bumpscit != tet_pair_to_rot_idx.end();
      ++bumpscit){
      const pair<size_t,size_t> & tet_pair = bumpscit->first;
      const size_t &type  = rot_type[bumpscit->second];
      inner_jump_type[tet_pair] = type;
      inner_jump_type[make_pair(tet_pair.second, tet_pair.first)] =
          get_trans_type(type);
    }
  dump_inner_face_jump_type("inner_face_jump_type_after_searcing", inner_jump_type);

  { // visual
//    std::vector<std::deque<std::pair<size_t,size_t> > > chain_list;
//    std::vector<deque<size_t> > singularities_type;
//    find_singularities_use_face_type(ortae, inner_jump_type, outside_face, chain_list,
//                                     singularities_type);
//    dump_singularity_chain_to_vtk_2("singularity_after_lazy_determin2.vtk",
//                                    node, chain_list, singularities_type);
    vector<size_t> outside_face_vec;
    vector<size_t> outside_face_type;
    //for(size_t fi = jump_face_vec.size(); fi < rot_type.size(); ++fi){
    for(boost::unordered_map<size_t,size_t>::const_iterator bumcit =
        surface_idx_to_rot_idx.begin(); bumcit != surface_idx_to_rot_idx.end();
        ++bumcit){
        const size_t face_idx = bumcit->first;
        const vector<size_t> & face_cut = fa_cut.faces_[face_idx];
        for(size_t i = 0; i < face_cut.size(); ++i){
            outside_face_vec.push_back(face_cut[i]);
          }
        outside_face_type.push_back(rot_type[bumcit->second]);
      }

    ofstream ofs_surface("surface_axis_type_after_searcing");
    dump_surface_axis_type(ofs_surface ,&outside_face_vec[0],
        outside_face_vec.size()/3,&outside_face_type[0]);

    ofstream ofs("surface_type_after_lazy_determin2.vtk");
    tri2vtk(ofs, &cut_node[0], cut_node.size(2), &outside_face_vec[0],
        outside_face_vec.size()/3);
    cell_data(ofs, &outside_face_type[0], outside_face_type.size(), "surface_patch");
  }

  {
    vector<size_t> faces;
    vector<size_t> face_type;

    for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator bumcit
        = face_pair_to_rot_idx.begin(); bumcit != face_pair_to_rot_idx.end(); ++bumcit){
        const pair<size_t,size_t> & face_pair = bumcit->first;
        if(rot_type[bumcit->second] == TRIVIAL_TYPE) continue;

        faces.insert(faces.end(), fa_cut.faces_[face_pair.first].begin(),
            fa_cut.faces_[face_pair.first].end());
        face_type.push_back(bumcit->second);
      }

    ofstream ofs("cut_face_in_solution.vtk");
    tri2vtk(ofs, &cut_node[0], cut_node.size(2), &faces[0], faces.size()/3);
    cell_data(ofs, &face_type[0], face_type.size(), "face_type");
  }
  return 0;
}

int searching_strategy::simple_check(
    const std::vector<size_t> & rot_type,
    const std::vector<bool> & is_surface,
    const std::vector<std::pair<size_t,size_t> > & jump_face_vec,
    const matrixst & cut_tet,
    const matrixd & cut_node,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst &cut_tet2tet,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::one_ring_tet_at_edge & ortae_cut,
    const matrixst & outside_face_cut,
    const matrixst & outside_face_idx_in_cut,
    const matrixst & face_pair_cut,
    const boost::unordered_map<size_t,size_t> &surface_idx_to_rot_idx,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &face_pair_to_rot_idx,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &tet_pair_to_rot_idx,
    bool no_surface)
{
  cerr << "# [info] simple check." << endl;
  unique_ptr<singularity_graph> g_test(
        singularity_graph::create(
          ortae_cut, ortae, cut_tet2tet,outside_face_cut, face_pair_cut,
          outside_face_idx_in_cut, surface_idx_to_rot_idx, face_pair_to_rot_idx,
          no_surface));
  if(!g_test.get()){
      cerr << "# [error] can not build singularity graph." << endl;
      return __LINE__;
    }
  detect_topology_graph_valid(
        rot_type, cut_tet, outside_face_cut, outside_face_idx_in_cut,
        cut_tet2tet, fa_cut, node, cut_node, jump_face_vec,face_pair_cut,
        ortae, *g_test, is_surface, tet_pair_to_rot_idx, face_pair_to_rot_idx,
        surface_idx_to_rot_idx, no_surface);

  return 0;
}

int searching_strategy::prev_process(
    const matrixst & orig_tet,
    const matrixd & orig_node,
    const matrixst &cut_tet,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
    const boost::unordered_map<size_t,size_t> & surface_type,
    matrixd &cut_node,
    matrixst & cut_tet2tet,
    matrixst & outside_face,
    matrixst & outside_face_cut,
    matrixst & outside_face_idx_in_cut,
    matrixst & outside_face_idx,
    std::vector<std::pair<size_t,size_t> > &jump_face_vec,
    boost::unordered_map<size_t,size_t> &surface_type_cut,
    jtf::mesh::one_ring_tet_at_edge  &ortae,
    jtf::mesh::one_ring_tet_at_edge  &ortae_cut)
{
  cerr << "# [info] start prev_process " << endl;
  cut_tet2tet.resize(max(cut_tet)+1,1);
  cut_tet2tet(cut_tet) = orig_tet(colon());

  get_outside_face_idx(fa_cut, outside_face_idx_in_cut);
  get_outside_face(fa_cut, outside_face_cut);

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea_cut(
        jtf::mesh::edge2cell_adjacent::create(outside_face_cut));

  get_outside_face(fa, outside_face);
  get_outside_face_idx(fa, outside_face_idx);


  // notice that: the surface_type stores the original face_idx with its type
  // we need to convert it to cut face_idx with its type

  boost::unordered_map<size_t,size_t> surface_idx_cut2orig;
  for(size_t fi = 0; fi < face_pair_cut.size(); ++fi){
      if(face_pair_cut[fi] == -1){
          const vector<size_t> & face_cut = fa_cut.faces_[outside_face_idx_in_cut[fi]];
          const size_t face_idx_orig =
              fa.get_face_idx(cut_tet2tet[face_cut[0]], cut_tet2tet[face_cut[1]],
              cut_tet2tet[face_cut[2]]);
          surface_idx_cut2orig[outside_face_idx_in_cut[fi]] = face_idx_orig;
        }
    }

  // assert(surface_idx_cut2orig.size() == surface_type.size());

  for(boost::unordered_map<size_t,size_t>::const_iterator bumcit
      = surface_idx_cut2orig.begin(); bumcit != surface_idx_cut2orig.end(); ++bumcit){
      boost::unordered_map<size_t,size_t>::const_iterator it =
          surface_type.find(bumcit->second);
      // assert(it != surface_type.end()); // in some searching strategy, the surface type is not needed,
      if(it != surface_type.end()){
          surface_type_cut[bumcit->first] = it->second;
        }
      else{
          if(surface_type.empty()) continue;
          else{
              cerr << "# [error] strange can not find face " << bumcit->first << endl;
              return __LINE__;
            }
        }
    }

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));

  if(!ea_cut.get() || !ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent" << endl;
      return __LINE__;
    }

  //  { // debug
  //    for(size_t ei = 0; ei < ea->edge2cell_.size(); ++ei){
  //      if(ea->is_boundary_edge(ea->edge2cell_[ei]))
  //        cerr << ea->edges_[ei].first << " " << ea->edges_[ei].second << endl;
  //    }
  //  }


  {
    boost::unordered_set<pair<size_t,size_t> > jump_face_set;
    for(size_t fi = 0; fi < face_pair_cut.size(); ++fi){
        if(face_pair_cut[fi] != -1){
            pair<size_t,size_t> face_pair(outside_face_idx_in_cut[fi],
                                          outside_face_idx_in_cut[face_pair_cut[fi]]);
            if(face_pair.first > face_pair.second)
              swap(face_pair.first, face_pair.second);
            jump_face_set.insert(face_pair);
          }
      }
    jump_face_vec.resize(jump_face_set.size());
    copy(jump_face_set.begin(), jump_face_set.end(), jump_face_vec.begin());
  }

  {
    ortae.add_tets(orig_tet,fa);
    ortae.sort_into_loop(orig_tet, orig_node);

    ortae_cut.add_tets(cut_tet, fa_cut);
  }

  cut_node = zeros<double>(3,cut_tet2tet.size());
  for(size_t pi = 0; pi < cut_node.size(2); ++pi){
      cut_node(colon(), pi) = orig_node(colon(), cut_tet2tet[pi]);
    }
  cerr << "# [info] end prev_process " << endl;

  return 0;
}

int searching_strategy::group_faces(
    const matrixst & cut_tet,
    const matrixd & cut_node,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & outside_face_cut,
    const matrixst & outside_face_cut_idx,
    const matrixst & face_pair_cut,
    const matrix<matrixd> & frame,
    vector<vector<pair<size_t,size_t> > > &jump_face_groups,
    vector<vector<size_t> > &surface_groups,
    vector<pair<size_t, double> > &group_type_with_frame_reliablity,
    vector<pair<size_t, double> > &group_type_with_surface_reliablity,
    const double frame_reliablity_threshold)
{
  jump_face_groups.clear();
  surface_groups.clear();

  group_type_with_frame_reliablity.clear();
  group_type_with_surface_reliablity.clear();

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face_cut));

  if(!ea.get()){
      cerr << "# [error] can not create edge2cell_adjacent" << endl;
      return __LINE__;
    }

  boost::unordered_map<size_t,size_t> outside_face_idx2idx_in_vec;
  {
    for(size_t fi = 0; fi < outside_face_cut_idx.size(); ++fi)
      outside_face_idx2idx_in_vec[outside_face_cut_idx[fi]] = fi;
  }

  vector<bool> visited_face_flag(outside_face_cut.size(), false);

  boost::unordered_map<size_t,size_t> jump_face_pair;
  boost::unordered_map<size_t, size_t> surface_to_idx_in_vec;
  {
    size_t jump_face_num = 0, surface_num = 0;
    for(size_t fi = 0; fi < face_pair_cut.size(); ++fi){
        if(face_pair_cut[fi] != -1){
            ++jump_face_num;
            jump_face_pair[outside_face_cut_idx[fi]] =
                outside_face_cut_idx[face_pair_cut[fi]];
          }else{
            ++surface_num;
            surface_to_idx_in_vec[outside_face_cut_idx[fi]] = fi;
          }
      }
    cerr << "# [info] jump_face_pair number " << jump_face_num/2 << endl;
    cerr << "# [info] surface number " << surface_num << endl;
  }

  size_t type = -1;
  double frame_reliablity = 0;
  vector<pair<size_t,size_t> > one_inner_group;
  vector<size_t> one_surface_group;
  for(size_t fi = 0; fi < face_pair_cut.size(); ++fi){
      if(face_pair_cut[fi] != -1 && !visited_face_flag[fi]){
          grow_mst_for_jump_faces(
                fa_cut,make_pair(outside_face_cut_idx[fi],
                                 outside_face_cut_idx[face_pair_cut[fi]]),
              outside_face_idx2idx_in_vec,
              frame_reliablity_threshold, *ea, frame, outside_face_cut_idx,
              jump_face_pair, one_inner_group, visited_face_flag, type, frame_reliablity);

          assert(!one_inner_group.empty());

          group_type_with_frame_reliablity.push_back(
                make_pair(type, frame_reliablity));
          jump_face_groups.push_back(one_inner_group);
          continue;
        }
      if(face_pair_cut[fi] == -1 && !visited_face_flag[fi]){
          grow_mst_for_surface_faces(
                outside_face_cut_idx[fi], frame_reliablity_threshold, *ea, frame,
                cut_tet, cut_node, fa_cut, outside_face_cut_idx,
                surface_to_idx_in_vec, visited_face_flag, one_surface_group, type,
                frame_reliablity);
          assert(!one_surface_group.empty());
          group_type_with_surface_reliablity.push_back(
                make_pair(type, frame_reliablity));
          surface_groups.push_back(one_surface_group);
          continue;
        }
    }

  {// recheck
    for(size_t fi = 0; fi < face_pair_cut.size(); ++fi){
        if(visited_face_flag[fi] == false)
          cerr << "# [error] face " << fi << " has not been visited." << endl;
      }
  }

  cerr << "# [info] inner faces are clustered into " << jump_face_groups.size()
       << " groups " << endl;
  cerr << "# [info] surface are clustered into " << surface_groups.size() << endl;

  {// visual
    vector<size_t> face_grouped;
    vector<size_t> face_patch;
    for(size_t gi = 0; gi < jump_face_groups.size(); ++gi){
        const vector<pair<size_t,size_t> > & group_ = jump_face_groups[gi];
        for(size_t fi = 0; fi < group_.size(); ++fi){
            const vector<size_t> & face_vec = fa_cut.faces_[group_[fi].first];
            face_grouped.insert(face_grouped.end(), face_vec.begin(), face_vec.end());
            face_patch.push_back(gi);
          }
      }
    ofstream ofs("jump_face_grouped.vtk");
    tri2vtk(ofs, &cut_node[0], cut_node.size(2), &face_grouped[0], face_grouped.size()/3);
    cell_data(ofs, &face_patch[0], face_patch.size(), "patch");
  }

  {
    vector<size_t> surface_face_grouped;
    vector<size_t> surface_face_patch;
    for(size_t gi = 0; gi < surface_groups.size(); ++gi){
        const vector<size_t> & group_ = surface_groups[gi];
        for(size_t fi = 0; fi < group_.size(); ++fi){
            const vector<size_t> & face_vec = fa_cut.faces_[group_[fi]];
            surface_face_grouped.insert(surface_face_grouped.end(),
                                        face_vec.begin(), face_vec.end());
            surface_face_patch.push_back(gi);
          }
      }
    ofstream ofs("surface_face_grouped.vtk");
    tri2vtk(ofs, &cut_node[0], cut_node.size(2), &surface_face_grouped[0], surface_face_grouped.size()/3);
    cell_data(ofs, &surface_face_patch[0], surface_face_patch.size(), "patch");
  }
  return 0;
}

int searching_strategy::group_faces_according_face_reliablity_for_minimal_error(
    const matrixst & cut_tet,
    const matrixd & cut_node,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & outside_face_cut,
    const matrixst & outside_face_cut_idx,
    const matrixst & face_pair_cut,
    const zjucad::matrix::matrix<matrixd> & frame,
    std::vector<std::vector<std::pair<size_t,size_t> > > &jump_face_groups,
    std::vector<std::vector<size_t> > &surface_groups,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_pair_to_rot_idx,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair_to_rot_idx,
    boost::unordered_map<size_t,size_t> & surface_to_rot_idx,
    std::vector<std::deque<std::pair<double,size_t> > > & diff,
    std::vector<bool> & is_surface,
    const double frame_reliablity_threshold)
{
  jump_face_groups.clear();
  surface_groups.clear();

  face_pair_to_rot_idx.clear();
  tet_pair_to_rot_idx.clear();
  surface_to_rot_idx.clear();

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face_cut));

  if(!ea.get()){
      cerr << "# [error] can not create edge2cell_adjacent" << endl;
      return __LINE__;
    }

  boost::unordered_map<size_t,size_t> outside_face_idx2idx_in_vec;
  {
    for(size_t fi = 0; fi < outside_face_cut_idx.size(); ++fi)
      outside_face_idx2idx_in_vec[outside_face_cut_idx[fi]] = fi;
  }

  vector<bool> visited_face_flag(outside_face_cut.size(), false);

  boost::unordered_map<size_t,size_t> jump_face_pair;
  boost::unordered_map<size_t, size_t> surface_to_idx_in_vec;
  {
    size_t jump_face_num = 0, surface_num = 0;
    for(size_t fi = 0; fi < face_pair_cut.size(); ++fi){
        if(face_pair_cut[fi] != -1){
            ++jump_face_num;
            jump_face_pair[outside_face_cut_idx[fi]] =
                outside_face_cut_idx[face_pair_cut[fi]];
          }else{
            ++surface_num;
            surface_to_idx_in_vec[outside_face_cut_idx[fi]] = fi;
          }
      }
    cerr << "# [info] jump_face_pair number " << jump_face_num/2 << endl;
    cerr << "# [info] surface number " << surface_num << endl;
  }

  vector<std::tuple<double,size_t,size_t> > group_type_with_frame_reliablity,
      group_type_with_surface_reliablity;

  size_t type = -1;
  double frame_reliablity = 0;
  vector<pair<size_t,size_t> > one_inner_group;
  vector<size_t> one_surface_group;
  for(size_t fi = 0; fi < face_pair_cut.size(); ++fi){
      if(face_pair_cut[fi] != -1 && !visited_face_flag[fi]){
          grow_mst_for_jump_faces(
                fa_cut,make_pair(outside_face_cut_idx[fi],
                                 outside_face_cut_idx[face_pair_cut[fi]]),
              outside_face_idx2idx_in_vec,
              frame_reliablity_threshold, *ea, frame, outside_face_cut_idx,
              jump_face_pair, one_inner_group, visited_face_flag, type, frame_reliablity);

          assert(!one_inner_group.empty());

          group_type_with_frame_reliablity.push_back(
                std::make_tuple(frame_reliablity,type,jump_face_groups.size()));
          jump_face_groups.push_back(one_inner_group);
          continue;
        }
      if(face_pair_cut[fi] == -1 && !visited_face_flag[fi]){
          grow_mst_for_surface_faces(
                outside_face_cut_idx[fi], frame_reliablity_threshold, *ea, frame,
                cut_tet, cut_node, fa_cut, outside_face_cut_idx,
                surface_to_idx_in_vec, visited_face_flag, one_surface_group, type,
                frame_reliablity);
          assert(!one_surface_group.empty());
          group_type_with_surface_reliablity.push_back(
                std::make_tuple(frame_reliablity,type,surface_groups.size()));
          surface_groups.push_back(one_surface_group);
          continue;
        }
    }

  //////////////////////////////////////////////////////////////////////////////
  {
    ///// gather the type_candidates
    diff.resize(jump_face_groups.size() + surface_groups.size());
    is_surface.resize(diff.size(), true);
    std::fill(is_surface.begin(), is_surface.begin() + jump_face_groups.size(), false);
    vector<pair<double,size_t> > jump_frame_candidates(24);
    vector<pair<double,size_t> > surface_type_candidates(3);

    for(size_t gi = 0; gi < jump_face_groups.size(); ++gi){
        cal_frame_difference_from_face_idx(
              jump_face_groups[gi].front().first,
              jump_face_groups[gi].front().second,
              fa_cut, frame, jump_frame_candidates);
        diff[gi].resize(24);
        copy(jump_frame_candidates.begin(), jump_frame_candidates.end(),
             diff[gi].begin());
        diff[gi].push_back(make_pair(std::numeric_limits<double>::infinity(),-1));
      }

    for(size_t gi = 0; gi < surface_groups.size(); ++gi){
        cal_frame_to_normal_difference(
              surface_groups[gi].front(), cut_tet, cut_node, fa_cut,
              frame, surface_type_candidates);
        diff[gi+jump_face_groups.size()].resize(3);
        copy(surface_type_candidates.begin(), surface_type_candidates.end(),
             diff[gi+jump_face_groups.size()].begin());
        diff[gi+jump_face_groups.size()].push_back(make_pair(std::numeric_limits<double>::infinity(),-1));
      }
  }

  // record the surface idx to rot_idx and inner jump face to rot_idx
  {
    for(size_t gi= 0; gi < jump_face_groups.size(); ++gi){
        const vector<pair<size_t,size_t> > & one_group = jump_face_groups[gi];
        for(size_t fi = 0; fi < one_group.size(); ++fi){
            face_pair_to_rot_idx[one_group[fi]] = gi;
            const pair<size_t,size_t> & face_pair = one_group[fi];
            const pair<size_t,size_t>& tet_pair_first =
                fa_cut.face2tet_[face_pair.first];
            const pair<size_t,size_t>& tet_pair_second =
                fa_cut.face2tet_[face_pair.second];
            pair<size_t,size_t> tet_pair(
                  tet_pair_first.first == -1?tet_pair_first.second:tet_pair_first.first,
                  tet_pair_second.first == -1?tet_pair_second.second:tet_pair_second.first);
            tet_pair_to_rot_idx[tet_pair] = gi;
          }
      }

    for(size_t gi = 0; gi < surface_groups.size(); ++gi){
        const vector<size_t> & one_group = surface_groups[gi];
        for(size_t fi = 0; fi < one_group.size(); ++fi){
            surface_to_rot_idx[one_group[fi]] = gi + jump_face_groups.size();
          }
      }
  }
  return 0;
}

int searching_strategy::group_faces_according_face_reliablity(
    const matrixst & cut_tet,
    const matrixd & cut_node,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & outside_face_cut,
    const matrixst & outside_face_cut_idx,
    const matrixst & face_pair_cut,
    const zjucad::matrix::matrix<matrixd> & frame,
    std::vector<std::vector<std::pair<size_t,size_t> > > &jump_face_groups,
    std::vector<std::vector<size_t> > &surface_groups,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_pair_to_rot_idx,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair_to_rot_idx,
    boost::unordered_map<size_t,size_t> & surface_to_rot_idx,
    matrixst & type_candidates,
    std::vector<bool> & is_surface,
    const double frame_reliablity_threshold)
{
  jump_face_groups.clear();
  surface_groups.clear();

  face_pair_to_rot_idx.clear();
  tet_pair_to_rot_idx.clear();
  surface_to_rot_idx.clear();

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face_cut));

  if(!ea.get()){
      cerr << "# [error] can not create edge2cell_adjacent" << endl;
      return __LINE__;
    }

  boost::unordered_map<size_t,size_t> outside_face_idx2idx_in_vec;
  {
    for(size_t fi = 0; fi < outside_face_cut_idx.size(); ++fi)
      outside_face_idx2idx_in_vec[outside_face_cut_idx[fi]] = fi;
  }

  vector<bool> visited_face_flag(outside_face_cut.size(), false);

  boost::unordered_map<size_t,size_t> jump_face_pair;
  boost::unordered_map<size_t, size_t> surface_to_idx_in_vec;
  {
    size_t jump_face_num = 0, surface_num = 0;
    for(size_t fi = 0; fi < face_pair_cut.size(); ++fi){
        if(face_pair_cut[fi] != -1){
            ++jump_face_num;
            jump_face_pair[outside_face_cut_idx[fi]] =
                outside_face_cut_idx[face_pair_cut[fi]];
          }else{
            ++surface_num;
            surface_to_idx_in_vec[outside_face_cut_idx[fi]] = fi;
          }
      }
    cerr << "# [info] jump_face_pair number " << jump_face_num/2 << endl;
    cerr << "# [info] surface number " << surface_num << endl;
  }

  // error, type, idx
  vector<std::tuple<double,size_t,size_t> > group_type_with_frame_reliablity,
      group_type_with_surface_reliablity;

  size_t type = -1;
  double frame_reliablity = 0;
  vector<pair<size_t,size_t> > one_inner_group;
  vector<size_t> one_surface_group;
  for(size_t fi = 0; fi < face_pair_cut.size(); ++fi){
      if(face_pair_cut[fi] != -1 && !visited_face_flag[fi]){
          grow_mst_for_jump_faces(
                fa_cut,make_pair(outside_face_cut_idx[fi],
                                 outside_face_cut_idx[face_pair_cut[fi]]),
              outside_face_idx2idx_in_vec,
              frame_reliablity_threshold, *ea, frame, outside_face_cut_idx,
              jump_face_pair, one_inner_group, visited_face_flag, type, frame_reliablity);

          assert(!one_inner_group.empty());

          group_type_with_frame_reliablity.push_back(
                std::make_tuple(frame_reliablity,type,jump_face_groups.size()));
          jump_face_groups.push_back(one_inner_group);
          continue;
        }
      if(face_pair_cut[fi] == -1 && !visited_face_flag[fi]){
          grow_mst_for_surface_faces(
                outside_face_cut_idx[fi], frame_reliablity_threshold, *ea, frame,
                cut_tet, cut_node, fa_cut, outside_face_cut_idx,
                surface_to_idx_in_vec, visited_face_flag, one_surface_group, type,
                frame_reliablity);
          assert(!one_surface_group.empty());
          group_type_with_surface_reliablity.push_back(
                std::make_tuple(frame_reliablity,type,surface_groups.size()));
          surface_groups.push_back(one_surface_group);
          continue;
        }
    }

  {// recheck
    for(size_t fi = 0; fi < face_pair_cut.size(); ++fi){
        if(visited_face_flag[fi] == false)
          cerr << "# [error] face " << fi << " has not been visited." << endl;
      }
  }

  cerr << "# [info] inner faces are clustered into " << jump_face_groups.size()
       << " groups " << endl;
  cerr << "# [info] surface are clustered into " << surface_groups.size() << endl;

  ////////////////////////////////////////////////////////////////////////////
  // rearange the order
  sort(group_type_with_frame_reliablity.begin(),
       group_type_with_frame_reliablity.end());
  sort(group_type_with_surface_reliablity.begin(),
       group_type_with_surface_reliablity.end());

  {
    vector<vector<pair<size_t,size_t> > > jump_face_groups_temp;
    vector<vector<size_t> > surface_groups_temp;
    jump_face_groups_temp.reserve(jump_face_groups.size());
    surface_groups_temp.reserve(surface_groups.size());
    for(size_t gi = 0; gi < group_type_with_frame_reliablity.size(); ++gi){
        jump_face_groups_temp.push_back(
              jump_face_groups[get<2>(group_type_with_frame_reliablity[gi])]);
      }

    jump_face_groups = jump_face_groups_temp;

    for(size_t gi = 0; gi < group_type_with_surface_reliablity.size(); ++gi){
        surface_groups_temp.push_back(
              surface_groups[get<2>(group_type_with_surface_reliablity[gi])]);
      }

    surface_groups = surface_groups_temp;

    ///// gather the type_candidates
    type_candidates = ones<size_t>(24, jump_face_groups.size() +
                                   surface_groups.size()) * -1;

    vector<pair<double,size_t> > jump_frame_candidates(24);
    vector<pair<double,size_t> > surface_type_candidates(3);

    for(size_t gi = 0; gi < jump_face_groups.size(); ++gi){
        cal_frame_difference_from_face_idx(
              jump_face_groups[gi].front().first,jump_face_groups[gi].front().second,
              fa_cut, frame, jump_frame_candidates);
        for(size_t ri = 0; ri < 24; ++ri){
            type_candidates(ri, gi) = jump_frame_candidates[ri].second;
          }
      }

    for(size_t gi = 0; gi < surface_groups.size(); ++gi){
        cal_frame_to_normal_difference(
              surface_groups[gi].front(), cut_tet, cut_node, fa_cut, frame,
              surface_type_candidates);
        for(size_t ri =0; ri < surface_type_candidates.size(); ++ri){
            type_candidates(ri, gi + jump_face_groups.size()) =
                surface_type_candidates[ri].second;
          }
      }
    //    cerr << type_candidates << endl;

    ///// construct the mapping
    is_surface.resize(jump_face_groups.size() + surface_groups.size(), true);
    std::fill(is_surface.begin(), is_surface.begin() + jump_face_groups.size(),false);

    for(size_t gi = 0; gi < jump_face_groups.size(); ++gi){
        const vector<pair<size_t,size_t> > & one_group = jump_face_groups[gi];
        for(size_t fi = 0; fi < one_group.size(); ++fi){
            const pair<size_t,size_t> & face_pair = one_group[fi];
            const pair<size_t,size_t> & tet_pair_first =
                fa_cut.face2tet_[face_pair.first] ;
            const pair<size_t,size_t> & tet_pair_second =
                fa_cut.face2tet_[face_pair.second];

            assert(fa_cut.is_outside_face(tet_pair_first) &&
                   fa_cut.is_outside_face(tet_pair_second));

            const pair<size_t,size_t> tet_pair(
                  (tet_pair_first.first == -1?
                     tet_pair_first.second: tet_pair_first.first),
                  (tet_pair_second.first == -1?
                     tet_pair_second.second: tet_pair_second.first));

            face_pair_to_rot_idx[one_group[fi]] = gi;
            tet_pair_to_rot_idx[tet_pair] = gi;
          }
      }
    for(size_t gi = 0; gi < surface_groups.size(); ++gi){
        const vector<size_t> & one_group = surface_groups[gi];
        for(size_t fi = 0; fi < one_group.size(); ++fi){
            surface_to_rot_idx[one_group[fi]] = gi + jump_face_groups.size();
          }
      }
  }

  {// visual
    vector<size_t> face_grouped;
    vector<size_t> face_patch;
    for(size_t gi = 0; gi < jump_face_groups.size(); ++gi){
        const vector<pair<size_t,size_t> > & group_ = jump_face_groups[gi];
        for(size_t fi = 0; fi < group_.size(); ++fi){
            const vector<size_t> & face_vec = fa_cut.faces_[group_[fi].first];
            face_grouped.insert(face_grouped.end(), face_vec.begin(), face_vec.end());
            face_patch.push_back(gi);
          }
      }
    ofstream ofs("jump_face_grouped.vtk");
    tri2vtk(ofs, &cut_node[0], cut_node.size(2), &face_grouped[0], face_grouped.size()/3);
    cell_data(ofs, &face_patch[0], face_patch.size(), "patch");
  }

  {
    vector<size_t> surface_face_grouped;
    vector<size_t> surface_face_patch;
    for(size_t gi = 0; gi < surface_groups.size(); ++gi){
        const vector<size_t> & group_ = surface_groups[gi];
        for(size_t fi = 0; fi < group_.size(); ++fi){
            const vector<size_t> & face_vec = fa_cut.faces_[group_[fi]];
            surface_face_grouped.insert(surface_face_grouped.end(),
                                        face_vec.begin(), face_vec.end());
            surface_face_patch.push_back(gi);
          }
      }
    ofstream ofs("surface_face_grouped.vtk");
    tri2vtk(ofs, &cut_node[0], cut_node.size(2), &surface_face_grouped[0], surface_face_grouped.size()/3);
    cell_data(ofs, &surface_face_patch[0], surface_face_patch.size(), "patch");
  }
  return 0;
}

int searching_strategy::cal_frame_difference_from_face_idx(
    const size_t & face_from,
    const size_t & face_to,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const zjucad::matrix::matrix<matrixd> & frame,
    std::vector<std::pair<double,size_t> > & frame_difference)
{
  assert(face_from < fa_cut.face2tet_.size() && face_to < fa_cut.face2tet_.size());

  const pair<size_t,size_t> & tet_pair_first = fa_cut.face2tet_[face_from];

  const pair<size_t,size_t> & tet_pair_second =  fa_cut.face2tet_[face_to];

  if(!fa_cut.is_outside_face(tet_pair_first)||
     !fa_cut.is_outside_face(tet_pair_second)) {
      cerr << "# [error] input face is not jump faces" << endl;
      return __LINE__;
    }

  cal_frame_difference(
        (tet_pair_first.first == -1? tet_pair_first.second:tet_pair_first.first),
        (tet_pair_second.first == -1? tet_pair_second.second:tet_pair_second.first),
        frame, frame_difference);

  return 0;
}

int searching_strategy::cal_frame_difference(
    const size_t & tet_from,
    const size_t & tet_to,
    const zjucad::matrix::matrix<matrixd> & frame,
    std::vector<std::pair<double,size_t> > & frame_difference)
{
  assert(tet_from < frame.size() && tet_to < frame.size());
  if(frame_difference.size() != 24)
    frame_difference.resize(24);

  for(size_t ri = 0; ri < 24; ++ri){
      frame_difference[ri] = make_pair(norm(frame[tet_from]* type_transition2(ri)
                                            - frame[tet_to]), ri);
    }
  sort(frame_difference.begin(), frame_difference.end());
  return 0;
}

int searching_strategy::cal_frame_to_normal_difference(
    const size_t & face_idx,
    const matrixst & cut_tet,
    const matrixd & cut_node,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const zjucad::matrix::matrix<matrixd> & frame,
    std::vector<std::pair<double, size_t> > & frame_to_normal_difference)
{
  frame_to_normal_difference.clear();
  assert(face_idx < fa_cut.faces_.size());
  const vector<size_t> & face_vec = fa_cut.faces_.at(face_idx);
  matrixst face_mat(3,1);
  copy(face_vec.begin(), face_vec.end(), face_mat.begin());
  matrixd face_normal = zeros<double>(3,1);

  jtf::mesh::cal_face_normal(face_mat, cut_node, face_normal);
  jtf::tetmesh::orient_one_face_normal_outside_tetmesh(
        cut_tet, cut_node, face_mat, face_idx, fa_cut,face_normal);

  const pair<size_t,size_t> & tet_pair = fa_cut.face2tet_[face_idx];
  assert(fa_cut.is_outside_face(tet_pair));
  const size_t & tet_idx = (tet_pair.first == -1? tet_pair.second: tet_pair.first);
  for(size_t ai = 0; ai < 3; ++ai){
      double error_0 = 1- dot(face_normal, frame[tet_idx](colon(),ai));
      double error_1 = 1- dot(face_normal, -1*frame[tet_idx](colon(),ai));
      frame_to_normal_difference.push_back(make_pair(min(error_0,error_1),ai));
    }

  sort(frame_to_normal_difference.begin(), frame_to_normal_difference.end());
  return 0;
}

int searching_strategy::grow_mst_for_jump_faces(
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const pair<size_t,size_t> & inner_face_pair,
    const boost::unordered_map<size_t,size_t> & outside_face_idx2idx_in_vec,
    const double &frame_reliablity_threshold,
    const jtf::mesh::edge2cell_adjacent & ea,
    const matrix<matrixd> & frame,
    const matrixst & outside_face_idx_cut,
    const boost::unordered_map<size_t,size_t> & jump_face_pair,
    std::vector<std::pair<size_t,size_t> > &one_jump_face_group,
    std::vector<bool> &visited_face_flag,
    size_t & type,
    double & frame_reliablity)
{
  one_jump_face_group.clear();
  boost::unordered_map<size_t,size_t>::const_iterator bumcit =
      outside_face_idx2idx_in_vec.find(inner_face_pair.first);
  if(bumcit == outside_face_idx2idx_in_vec.end()){
      cerr << "# [error] can not find face " << inner_face_pair.first
           << " in outside_face_cut" << endl;
      return __LINE__;
    }
  assert(bumcit->second < visited_face_flag.size());
  visited_face_flag.at(bumcit->second) = true;

  bumcit = outside_face_idx2idx_in_vec.find(inner_face_pair.second);
  if(bumcit == outside_face_idx2idx_in_vec.end()){
      cerr << "# [error] can not find face " << inner_face_pair.second
           << " in outside_face_cut" << endl;
      return __LINE__;
    }
  assert(bumcit->second < visited_face_flag.size());
  visited_face_flag.at(bumcit->second) = true;

  vector<pair<double, size_t> > frame_difference(24);

  cal_frame_difference_from_face_idx(
        inner_face_pair.first, inner_face_pair.second,fa_cut,
        frame, frame_difference);

  type = frame_difference.front().second;
  frame_reliablity = frame_difference.front().first / frame_difference.at(1).first;

  if(frame_difference.at(0).first / frame_difference.at(1).first >
     frame_reliablity_threshold){
      one_jump_face_group.push_back(inner_face_pair);
      return 0;
    }

  const size_t rotation_type = frame_difference.front().second;
  one_jump_face_group.push_back(inner_face_pair);

  stack<size_t> face_stack;
  face_stack.push(inner_face_pair.first);

  while(!face_stack.empty()){
      const size_t face_idx = face_stack.top();
      face_stack.pop();

      const vector<size_t> & face_vec = fa_cut.faces_[face_idx];

      for(size_t pi = 0; pi < face_vec.size(); ++pi){
          const size_t edge_idx =
              ea.get_edge_idx(face_vec[pi], face_vec[(pi+1)%face_vec.size()]);
          const pair<size_t,size_t> & face_pair = ea.edge2cell_.at(edge_idx);
          const size_t other_face_idx =
              (outside_face_idx_cut[face_pair.first] == face_idx?
                outside_face_idx_cut[face_pair.second]
              :outside_face_idx_cut[face_pair.first]);

          boost::unordered_map<size_t,size_t>::const_iterator bucit =
              outside_face_idx2idx_in_vec.find(other_face_idx);
          assert(bucit->second < visited_face_flag.size());
          if(visited_face_flag[bucit->second]) continue;

          // find a linked face
          boost::unordered_map<size_t,size_t>::const_iterator linked_face_pair_it
              = jump_face_pair.find(other_face_idx);
          if(linked_face_pair_it == jump_face_pair.end()){
              continue; // it's a original surface
            }
          const size_t & face_pair_for_linked_face = linked_face_pair_it->second;
          boost::unordered_map<size_t,size_t>::const_iterator lf2idx =
              outside_face_idx2idx_in_vec.find(face_pair_for_linked_face);
          if(lf2idx == outside_face_idx2idx_in_vec.end()){
              cerr << "# [error] strange can not find face " << face_pair_for_linked_face
                   << " in outside_face_idx2idx_in_vec." << endl;
              return __LINE__;
            }
          assert(lf2idx->second < visited_face_flag.size());
          if(visited_face_flag.at(lf2idx->second) == true){
              cerr << "# [error] strange: face " << face_pair_for_linked_face
                   << " at " << lf2idx->second << " in vec should not be visited."
                   << endl;
              return __LINE__;
            }

          // at this step, we need to check whether the face pair
          // (other_face_idx, face_pair_for_linked_face) has the same type and
          // whether the frame reliablity is correct
          cal_frame_difference_from_face_idx(
                other_face_idx, face_pair_for_linked_face,
                fa_cut, frame, frame_difference);

          if(frame_difference.front().second  != type) continue;
          if(frame_difference.front().first / frame_difference.at(1).first >
             frame_reliablity_threshold) continue;

          // find a new face pair which can be grouped
          one_jump_face_group.push_back(make_pair(other_face_idx, face_pair_for_linked_face));
          visited_face_flag.at(bucit->second) = true;
          visited_face_flag.at(lf2idx->second) = true;
          face_stack.push(other_face_idx);
        }
    }
  return 0;
}

int searching_strategy::grow_mst_for_surface_faces(
    const size_t & surface_idx,
    const double & frame_reliablity_threshold,
    const jtf::mesh::edge2cell_adjacent &ea,
    const matrix<matrixd> & frame,
    const matrixst & cut_tet,
    const matrixd & cut_node,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & outside_face_cut_idx,
    const boost::unordered_map<size_t,size_t> &surface_to_idx_in_vec,
    std::vector<bool> & visited_face_flag,
    std::vector<size_t> &one_surface_group,
    size_t &type,
    double &frame_reliablity)
{
  one_surface_group.clear();

  vector<pair<double,size_t> > frame_to_normal_difference;
  cal_frame_to_normal_difference(surface_idx,cut_tet, cut_node, fa_cut,frame,
                                 frame_to_normal_difference);
  type = frame_to_normal_difference.front().second;
  frame_reliablity = frame_to_normal_difference.front().first;

  boost::unordered_map<size_t,size_t>::const_iterator bucit =
      surface_to_idx_in_vec.find(surface_idx);

  if(bucit == surface_to_idx_in_vec.end()){
      cerr << "# [error] strange can not find surface idx " << surface_idx
           << endl;
      return __LINE__;
    }

  visited_face_flag[bucit->second] = true;
  one_surface_group.push_back(surface_idx);

  if(frame_reliablity > frame_reliablity_threshold){
      return 0;
    }

  stack<size_t> face_stack;
  face_stack.push(surface_idx);

  while(!face_stack.empty()){
      const size_t face_idx = face_stack.top();
      face_stack.pop();

      const vector<size_t> & face_vec = fa_cut.faces_[face_idx];
      for(size_t pi = 0; pi < face_vec.size(); ++pi){
          const size_t edge_idx =
              ea.get_edge_idx(face_vec[pi], face_vec[(pi+1)%face_vec.size()]);
          assert(edge_idx < ea.edges_.size());
          const pair<size_t,size_t> & tri_pair = ea.edge2cell_.at(edge_idx);
          const size_t other_face_idx =
              (outside_face_cut_idx[tri_pair.first] == face_idx?
                outside_face_cut_idx[tri_pair.second]
              : outside_face_cut_idx[tri_pair.first]);

          boost::unordered_map<size_t,size_t>::const_iterator other_face_idx_2_vec_idx
              = surface_to_idx_in_vec.find(other_face_idx);
          if(other_face_idx_2_vec_idx == surface_to_idx_in_vec.end()){
              continue; // it's an inner face
            }

          if(visited_face_flag.at(other_face_idx_2_vec_idx->second)) continue;

          cal_frame_to_normal_difference(other_face_idx,cut_tet, cut_node, fa_cut,frame,
                                         frame_to_normal_difference);
          if(frame_to_normal_difference.front().second != type) continue;
          if(frame_to_normal_difference.front().first /
             frame_to_normal_difference.at(1).first > frame_reliablity_threshold)
            continue;

          one_surface_group.push_back(other_face_idx);
          visited_face_flag.at(other_face_idx_2_vec_idx->second) = true;
        }
    }
  return 0;
}

int searching_strategy::orig_check_graph(
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
    const boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const bool no_surface,
    const zjucad::matrix::matrix<matrixd> * frame_ptr)
{
  boost::unordered_map<size_t,size_t> surface_type_cut;
  matrixd cut_node;
  matrixst cut_tet2tet,outside_face, outside_face_cut,
      outside_face_idx_in_cut, outside_face_idx;
  vector<pair<size_t,size_t> > jump_face_vec;
  jtf::mesh::one_ring_tet_at_edge  ortae;
  jtf::mesh::one_ring_tet_at_edge ortae_cut;

  prev_process(
        orig_tet,node,cut_tet, fa, fa_cut,face_pair_cut, surface_type, cut_node,
        cut_tet2tet, outside_face, outside_face_cut, outside_face_idx_in_cut,
        outside_face_idx, jump_face_vec, surface_type_cut, ortae, ortae_cut);

  vector<size_t> rot_type(jump_face_vec.size() + surface_type_cut.size(),-1);
  matrixst type_candidates(24,rot_type.size());
  vector<bool> is_surface(rot_type.size(), false);

  boost::unordered_map<size_t,size_t> surface_idx_to_rot_idx;
  boost::unordered_map<size_t,size_t> surface_idx_with_rot_idx;
  boost::unordered_map<pair<size_t,size_t>,size_t> face_pair_to_rot_idx;
  boost::unordered_map<pair<size_t,size_t>,size_t> tet_pair_to_rot_idx;
  boost::unordered_map<size_t,pair<size_t,size_t> > face_pair_with_rot_idx;

  adjust_face_and_type_order_to_speed_up(
        jump_face_vec, inner_face_jump_type, surface_type_cut, fa_cut, cut_tet2tet,
        type_candidates, surface_idx_to_rot_idx,surface_idx_with_rot_idx,
        face_pair_to_rot_idx,face_pair_with_rot_idx,tet_pair_to_rot_idx,
        is_surface,frame_ptr);

  copy(type_candidates(0,colon()).begin(), type_candidates(0,colon()).end(), rot_type.begin());


  {
    ofstream ofs("inner_face_jump_type.debug");
    for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit
        = tet_pair_to_rot_idx.begin(); cit != tet_pair_to_rot_idx.end(); ++cit){
        if(is_trivial_type(rot_type[cit->second])) continue;
        ofs << cit->first.first << " " << cit->first.second
            << " " << rot_type[cit->second] << endl;
        ofs << cit->first.second << " " << cit->first.first << " "
            << get_trans_type(rot_type[cit->second]) << endl;
      }
    ofstream ofs_surf("surface_type.debug");
    for(boost::unordered_map<size_t,size_t>::const_iterator cit
        = surface_idx_to_rot_idx.begin(); cit != surface_idx_to_rot_idx.end(); ++cit){
        ofs_surf << cit->first << " " << rot_type[cit->second] << endl;
      }
  }

  {// visual rotation face
    vector<size_t> rot_face;
    vector<size_t> rot_face_type;
    for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit
        = face_pair_to_rot_idx.begin(); cit != face_pair_to_rot_idx.end(); ++cit)
      {
        if(rot_type[cit->second] != TRIVIAL_TYPE){
            const vector<size_t> &face_vec = fa_cut.faces_[cit->first.first];
            rot_face.insert(rot_face.end(), face_vec.begin(), face_vec.end());
            rot_face_type.push_back(rot_type[cit->second]);
          }
      }
    ofstream ofs("rot_face.vtk");
    tri2vtk(ofs, &cut_node[0], cut_node.size(2), &rot_face[0], rot_face.size()/3);
    cell_data(ofs, &rot_face_type[0], rot_face_type.size(), "rot_type");
  }

  simple_check(rot_type, is_surface, jump_face_vec, cut_tet, cut_node, node,
               fa_cut, cut_tet2tet, ortae, ortae_cut, outside_face_cut,
               outside_face_idx_in_cut, face_pair_cut, surface_idx_to_rot_idx,
               face_pair_to_rot_idx,tet_pair_to_rot_idx,no_surface);
  return 0;
}

int searching_strategy::searching_solutions_by_using_memory_group(
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
    const boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const zjucad::matrix::matrix<matrixd> * frame_ptr)
{
  boost::unordered_map<size_t,size_t> surface_type_cut;
  matrixd cut_node;
  matrixst cut_tet2tet,outside_face, outside_face_cut,
      outside_face_idx_in_cut, outside_face_idx;
  vector<pair<size_t,size_t> > jump_face_vec;
  jtf::mesh::one_ring_tet_at_edge  ortae;
  jtf::mesh::one_ring_tet_at_edge ortae_cut;

  prev_process(
        orig_tet,node,cut_tet, fa, fa_cut,face_pair_cut, surface_type, cut_node,
        cut_tet2tet, outside_face, outside_face_cut, outside_face_idx_in_cut,
        outside_face_idx, jump_face_vec, surface_type_cut, ortae, ortae_cut);

  vector<vector<pair<size_t,size_t> > > jump_face_groups;
  vector<vector<size_t> > surface_groups;
  vector<pair<size_t,double> > group_type_with_frame_reliablity,
      group_type_with_surface_reliablity;
  const double frame_reliablity_threshold = 0.1;

  // this function is used to group faces according to the type and reliablity

  boost::unordered_map<pair<size_t,size_t>,size_t>  face_pair_to_rot_idx;
  boost::unordered_map<pair<size_t,size_t>,size_t>  tet_pair_to_rot_idx;
  boost::unordered_map<size_t,size_t>  surface_idx_to_rot_idx;
  matrixst type_candidates;
  vector<bool> is_surface;

  group_faces_according_face_reliablity(
        cut_tet,cut_node,fa_cut, outside_face_cut,outside_face_idx_in_cut,
        face_pair_cut, *frame_ptr, jump_face_groups, surface_groups,
        face_pair_to_rot_idx, tet_pair_to_rot_idx, surface_idx_to_rot_idx,
        type_candidates, is_surface, frame_reliablity_threshold);

  unique_ptr<singularity_graph> g(
        singularity_graph::create(
          ortae_cut, ortae, cut_tet2tet, outside_face_cut,  face_pair_cut,
          outside_face_idx_in_cut, surface_idx_to_rot_idx, face_pair_to_rot_idx));

  if(!g.get()){
      cerr << "# [info] build singularity graph failed." << endl;
      return __LINE__;
    }

  vector<size_t> rot_type(jump_face_groups.size() + surface_groups.size() ,-1);
  vector<size_t> tried_type_idx(rot_type.size(), -1);

  // start to check the question tree, and try to figure out the answer
  size_t fi = 0;
  size_t searching_steps = 0;
  const size_t fn = rot_type.size();
  stack<state_each_step> ss;
  state_each_step ses;
  double write_time = 0.0, read_time = 0.0, check_time = 0.0, trans_time = 0.0;
  timer tr;
  bool is_valid;
  deque<size_t> check_numbers;
  while(fi < fn) {
      // cerr << "# [info] checking face " << fi << ", tried "  << tried_type_idx[fi] << endl;
      size_t ri;
      for(ri = tried_type_idx[fi]+1; ri < 24; ++ri) {
          ++searching_steps;
          is_valid = true;
          const size_t &next_type = type_candidates(ri, fi);

          if(next_type == -1) break;  // index = -1 meas the surface type is tested over

          long_number_add(check_numbers, 1);
          rot_type[fi] = next_type;
          tried_type_idx[fi] = ri;

          if(is_surface[fi] &&
             is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
            continue;

          tr.start();
          g->save_state_mem(ses);
          ss.push(ses);
          tr.finish();
          write_time += tr.result_c();
          //g->get_step_state(ss);

          assert(!g->is_group_info_broken());

          tr.start();
          if(is_surface[fi]) { // the third type = -1 means this face is surface
              g->insert_surface_transition_group(
                    next_type, fa_cut, cut_tet2tet,
                    surface_groups[fi - jump_face_groups.size()], rot_type,
                  tet_pair_to_rot_idx);
            }else {
              if(g->insert_transition_group(
                   jump_face_groups[fi],cut_tet,cut_tet2tet,
                   fa_cut, ortae,rot_type,face_pair_to_rot_idx,
                   tet_pair_to_rot_idx,1))
                is_valid = false;
            }
          tr.finish();
          trans_time += tr.result_c();

          assert(!g->is_group_info_broken());

          tr.start();
          if(is_valid)
            is_valid = g->is_valid_with_info(cut_tet2tet,rot_type,ortae,
                                             surface_idx_to_rot_idx,1);
          tr.finish();
          check_time += tr.result_c();

          if(is_valid)      break;

          {
            assert(!is_valid);
            rot_type[fi] = -1;// try_rot_type.second;
            tr.start();
            g->load_state_mem(ss.top());
            ss.pop();
            tr.finish();
            read_time += tr.result_c();

            assert(!g->is_group_info_broken());
            assert(g->is_valid_with_info(cut_tet2tet,rot_type, ortae,
                                         surface_idx_to_rot_idx));
          }
        }
      // revert to the previous cut_face setting.
      if(( !is_surface[fi] && ri == 24) || (is_surface[fi] && ri == 3)) {
          if(fi == 0) {
              break; // can not find one solution
            }

          tried_type_idx[fi] = -1;
          rot_type[fi] = -1;
          --fi;

          tr.start();
          g->load_state_mem(ss.top());
          ss.pop();
          tr.finish();
          read_time += tr.result_c();

          assert(!g->is_group_info_broken());
          assert(g->is_valid_with_info(cut_tet2tet,rot_type, ortae,
                                       surface_idx_to_rot_idx));
          rot_type[fi] = -1;
        }else{
          ++fi;
        }
    }


  cerr << "# write time " << write_time/1e6 << endl;
  cerr << "# read time " << read_time/1e6 << endl;
  cerr << "# trans time " << trans_time/1e6 << endl;
  cerr << "# check time " << check_time/1e6 << endl;
  cerr << "# searching steps " << searching_steps << endl;
  cerr << "# iterate numbers ";
  print_long_number(cerr, check_numbers);

  if(fi == 0) {
      cerr << "cannot find a solution" << endl;
      return __LINE__;
    } else {
      cerr << "find a solution. " << endl;
      post_process(rot_type,cut_node, node, outside_face, fa_cut, ortae,
                   surface_idx_to_rot_idx, tet_pair_to_rot_idx,face_pair_to_rot_idx,
                   inner_face_jump_type);

      simple_check(rot_type, is_surface, jump_face_vec, cut_tet, cut_node, node,
                   fa_cut, cut_tet2tet, ortae, ortae_cut, outside_face_cut,
                   outside_face_idx_in_cut, face_pair_cut, surface_idx_to_rot_idx,
                   face_pair_to_rot_idx,tet_pair_to_rot_idx);
      return 0;
    }

  return 0;
}

int searching_strategy::searching_solutions_by_minimal_error_search_group(
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const matrixst & face_pair_cut,
    const zjucad::matrix::matrix<matrixd> & frame,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type)
{
  matrixd cut_node;
  matrixst outside_face_cut, outside_face_cut_idx;
  matrixst cut_tet2tet, outside_face;
  vector<pair<size_t,size_t> > jump_face_vec;
  jtf::mesh::one_ring_tet_at_edge  ortae;
  jtf::mesh::one_ring_tet_at_edge ortae_cut;
  {
    boost::unordered_map<size_t,size_t> surface_type_cut;
    matrixst outside_face_idx;
    boost::unordered_map<size_t,size_t> surface_type;


    prev_process(orig_tet,node,cut_tet, fa, fa_cut,face_pair_cut, surface_type, cut_node,
                 cut_tet2tet, outside_face, outside_face_cut, outside_face_cut_idx, outside_face_idx,
                 jump_face_vec, surface_type_cut, ortae, ortae_cut);
  }


  //group_faces_according_face_reliablity()
  assert(frame.size() == cut_tet.size(2));
  size_t surface_num = 0;
  for(size_t fi = 0; fi < face_pair_cut.size(); ++fi)
    if(face_pair_cut[fi] == -1)   ++surface_num;


  boost::unordered_map<size_t,size_t> surface_idx_to_rot_idx;
  boost::unordered_map<std::pair<size_t,size_t>,size_t> face_pair_to_rot_idx;
  boost::unordered_map<std::pair<size_t,size_t>,size_t> tet_pair_to_rot_idx;

  const double frame_reliablity_threshold = 0.1;

  vector<vector<pair<size_t,size_t> > > jump_face_groups;
  vector<vector<size_t> > surface_groups;

  vector<deque<pair<double,size_t> > > diff;
  vector<bool> is_surface;

  group_faces_according_face_reliablity_for_minimal_error(
        cut_tet, cut_node, fa_cut, outside_face_cut, outside_face_cut_idx,
        face_pair_cut, frame, jump_face_groups, surface_groups,
        face_pair_to_rot_idx, tet_pair_to_rot_idx, surface_idx_to_rot_idx, diff,
        is_surface, frame_reliablity_threshold);

  unique_ptr<singularity_graph> g(
        singularity_graph::create(
          ortae_cut, ortae, cut_tet2tet, outside_face_cut,  face_pair_cut,
          outside_face_cut_idx, surface_idx_to_rot_idx, face_pair_to_rot_idx));
  if(!g.get()){
      cerr << "# [error] can not build singularity graph" << endl;
      return __LINE__;
    }

  const size_t total_face_num = is_surface.size();
  vector<bool> visited_face_flag(total_face_num, false);
  vector<size_t> rot_type(total_face_num, -1);
  vector<state_each_step> ss;
  ss.reserve(rot_type.size());
  state_each_step ses;
  double write_time = 0.0, read_time = 0.0, check_time = 0.0, trans_time = 0.0;
  vector<pair<size_t,size_t> > face_type_stack;
  size_t windows_size = total_face_num;
  timer tr;
  while(face_type_stack.size() < total_face_num){
      size_t begin_ = face_type_stack.size();
      vector<std::tuple<double, size_t, size_t> > candidate_faces;
      for(size_t fi = 0; fi < total_face_num; ++fi){
          if(!visited_face_flag[fi]){
              candidate_faces.push_back(
                    std::make_tuple(diff[fi].front().first,fi,
                                    diff[fi].front().second));
            }
        }

      sort(candidate_faces.begin(), candidate_faces.end());
      size_t ci = 0;
      bool is_valid = true;
      for(; ci < candidate_faces.size(); ++ci){
          face_type_stack.push_back(make_pair(get<1>(candidate_faces[ci]),
                                              get<2>(candidate_faces[ci])));

          rot_type[get<1>(candidate_faces[ci])] = get<2>(candidate_faces[ci]);
          visited_face_flag[get<1>(candidate_faces[ci])] = true;
          tr.start();
          if(is_surface[get<1>(candidate_faces[ci])]) {
              g->insert_surface_transition_group(
                    get<2>(candidate_faces[ci]), fa_cut, cut_tet2tet,
                    surface_groups[get<1>(candidate_faces[ci]) - jump_face_groups.size()], rot_type,
                  tet_pair_to_rot_idx);
            }else {
              if(g->insert_transition_group(
                   jump_face_groups[get<1>(candidate_faces[ci])],cut_tet,cut_tet2tet,
                   fa_cut, ortae,rot_type,face_pair_to_rot_idx,
                   tet_pair_to_rot_idx,1))
                is_valid = false;
            }
          tr.finish();
          trans_time += tr.result_c();

          tr.start();
          g->save_state_mem(ses);
          ss.push_back(ses);
          tr.finish();
          write_time += tr.result_c();

          if(!is_valid) break;
        }

      if(is_valid)
        is_valid = g->is_valid_with_info(cut_tet2tet, rot_type, ortae,
                                         surface_idx_to_rot_idx,1);

      if(is_valid){
#ifdef dump_path
          {
            for(size_t i = 0; i < candidate_faces.size(); ++i){
                if(is_surface[get<1>(candidate_faces[i])]) {
                    const vector<size_t> & vec =
                        surface_groups[get<1>(candidate_faces[i]) - jump_face_groups.size()];
                    for(size_t j = 0; j < vec.size(); ++j)
                      cout << "+ " << vec[j] << " "
                           << get<2>(candidate_faces[i]) << endl;
                  }else{
                    const vector<pair<size_t,size_t> > & vec =
                        jump_face_groups[get<1>(candidate_faces[i])];
                    for(size_t j = 0; j < vec.size(); ++j)
                      cout << "+ " << vec[j].first << " "
                           << get<2>(candidate_faces[i]) << endl;
                  }
              }
          }
#endif
          windows_size = (total_face_num - face_type_stack.size())/2>0?
                (total_face_num - face_type_stack.size())/2:1;
          continue;
        }

      { // check the graph
        if(ci != candidate_faces.size())
          ++ci; // to ensure the [) range

        vector<size_t> rot_type_bkp;
        vector<bool> visited_flag_bkp;
        size_t check_begin = begin_;
        size_t check_end = ci+begin_;
        bool is_valid;
        while(1){
            size_t mid_face_idx = (check_begin + check_end)/2;
            g->load_state_mem(ss[mid_face_idx]);
            rot_type_bkp = rot_type;
            visited_flag_bkp = visited_face_flag;
            // recovery the rot_type
            for(size_t ti = mid_face_idx + 1; ti < ci+begin_; ++ti){
                rot_type_bkp[get<1>(candidate_faces[ti-begin_])] = -1;
                visited_flag_bkp[get<1>(candidate_faces[ti-begin_])] = false;
              }
            is_valid = true;
            is_valid = g->is_valid_with_info(cut_tet2tet, rot_type_bkp, ortae,
                                             surface_idx_to_rot_idx,1);
            if(is_valid){
                if(is_plane_degenerated(rot_type, surface_idx_to_rot_idx))
                  is_valid = false;
              }
            if(is_valid){
                if(mid_face_idx == check_begin){
                    rot_type = rot_type_bkp;
                    visited_face_flag = visited_flag_bkp;
                    break;
                  }
                check_begin = mid_face_idx;
                continue;
              }else{
                rot_type = rot_type_bkp;
                visited_face_flag = visited_flag_bkp;
                if(check_end == mid_face_idx+1)
                  check_end = mid_face_idx;
                else
                  check_end = mid_face_idx + 1;
                if(check_end == check_begin){
                    rot_type[get<1>(candidate_faces[check_end - begin_])] = -1;
                    visited_face_flag[get<1>(candidate_faces[check_end - begin_])] = false;

                    g->load_state_mem(ss[check_end-1]);
                    break;
                  }
              }
          }

        face_type_stack.resize(check_end+1);
        ss.resize(check_end);
#ifdef dump_path
        {
          for(size_t i = 0; i < check_end - begin_+1; ++i){
              if(is_surface[get<1>(candidate_faces[i])]) {
                  const vector<size_t> & vec =
                      surface_groups[get<1>(candidate_faces[i]) - jump_face_groups.size()];
                  for(size_t j = 0; j < vec.size(); ++j)
                    cout << "+ " << vec[j] << " "
                         << get<2>(candidate_faces[i]) << endl;
                }else{
                  const vector<pair<size_t,size_t> > & vec =
                      jump_face_groups[get<1>(candidate_faces[i])];
                  for(size_t j = 0; j < vec.size(); ++j)
                    cout << "+ " << vec[j].first << " "
                         << get<2>(candidate_faces[i]) << endl;
                }
            }
          if(is_surface[get<1>(candidate_faces[check_end - begin_])]){
              cout << "- " << surface_groups[get<1>(candidate_faces[check_end - begin_]) - jump_face_groups.size()].size() << endl;
            }else{
              cout << "- " << jump_face_groups[get<1>(candidate_faces[check_end - begin_])].size() << endl;
            }
        }
#endif
      }

      {
        // have tested all candidates of this face, need to go back
        while(!face_type_stack.empty()){
            const pair<size_t,size_t> & latest_choose = face_type_stack.back();
#ifdef dump_path
            {
              if(is_surface[latest_choose.first]){
                  cerr << "- " << surface_groups[latest_choose.first
                          - jump_face_groups.size()].size() << endl;
                }else{
                  cerr << "- " << jump_face_groups[latest_choose.first].size() << endl;
                }
            }
#endif
            if(diff[latest_choose.first][1].second == -1){
                // if the next type choose of face latest_choose.first
                // ---diff[latest_choose.first][1].second--- is -1, which means all type
                // of this face has been tried, need to go back

                sort(diff[latest_choose.first].begin(), diff[latest_choose.first].end());
                visited_face_flag[latest_choose.first] = false;
                rot_type[latest_choose.first] = -1;
                face_type_stack.pop_back();
                ss.pop_back();
              }else{
                //invalid_face.insert(latest_choose.first);
                diff[latest_choose.first].push_back(diff[latest_choose.first].front());
                diff[latest_choose.first].pop_front();
                visited_face_flag[latest_choose.first] = false;
                rot_type[latest_choose.first] = -1;
                face_type_stack.pop_back();
                tr.start();
                g->load_state_mem(ss.back());
                tr.finish();
                read_time += tr.result_c();

                //ss.pop_back();
                break;
              }
          }
        if(face_type_stack.empty())
          break;
      }
    }

  cerr << "# write time " << write_time/1e6 << endl;
  cerr << "# read time " << read_time/1e6 << endl;
  cerr << "# trans time " << trans_time/1e6 << endl;
  cerr << "# check time " << check_time/1e6 << endl;
  //cerr << "# searching steps " << searching_steps << endl;
  //cerr << "# iterate numbers " ;
  //print_long_number(cerr, check_number);

  if(face_type_stack.empty())
    cerr << "# [error] can not find a solution." << endl;
  else{
      cerr << "# [info] find a solution." << endl;

      post_process(rot_type, cut_node, node, outside_face, fa_cut,ortae,
                   surface_idx_to_rot_idx,tet_pair_to_rot_idx,
                   face_pair_to_rot_idx, inner_face_jump_type);

      simple_check(
            rot_type, is_surface, jump_face_vec, cut_tet, cut_node, node,
            fa_cut, cut_tet2tet, ortae, ortae_cut, outside_face_cut,
            outside_face_cut_idx, face_pair_cut, surface_idx_to_rot_idx,
            face_pair_to_rot_idx, tet_pair_to_rot_idx);
    }
  return 0;
}

//static int searching_strategy::searching_solutions_by_minimal_error_binary_search_group(
//    const jtf::mesh::face2tet_adjacent & fa,
//    const jtf::mesh::face2tet_adjacent & fa_cut,
//    const matrixst & orig_tet,
//    const matrixst & cut_tet,
//    const matrixd & node,
//    const matrixst & face_pair_cut,
//    const zjucad::matrix::matrix<matrixd> & frame,
//    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type)
//{

//  return 0;
//}
