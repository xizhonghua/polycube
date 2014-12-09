#include <map>
#include <vector>
#include <memory>
#include <set>
#include <deque>
#include <iostream>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <fstream>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_set.hpp>

#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <jtflib/util/util.h>

#include "hex_param.h"
#include "topology_analysis.h"
#include "find_singularities.h"
#include "../tetmesh/tetmesh.h"
#include "../common/transition.h"
#include "../common/zyz.h"
#include "../common/vtk.h"
#include "../common/util.h"
#include "../common/visualize_tool.h"
#include "../common/transition_type.h"
#include "common.h"
#include "../spherical_harmonics/rot_cubic_f_SH.h"
#include "global_alignment.h"
#include <jtflib/util/container_operation.h>
using namespace std;
using namespace zjucad::matrix;



void singularity_extractor::extract(const zjucad::matrix::matrix<double> &zyz,
                                    std::vector<std::deque<std::pair<size_t, size_t> > > &chain,
                                    std::vector<std::deque<size_t> > &singularities_type) const
{
  matrix<matrix<double> > frame;
  zyz2frame(zyz, frame);
  extract(frame, chain, singularities_type);
}

void singularity_extractor::extract(const zjucad::matrix::matrix<zjucad::matrix::matrix<double> > &frame,
                                    std::vector<std::deque<std::pair<size_t, size_t> > > &chain,
                                    std::vector<std::deque<size_t> > &singularities_type) const
{
  chain.clear();
  singularities_type.clear();

  boost::unordered_set<pair<size_t,size_t> > singularities_segment;
  boost::unordered_map<pair<size_t,size_t>,size_t> singularity_type;
  matrixd transition = eye<double>(3);
  matrixd shuffle_rot(3,3);

  // find out all the singulatiy segments
  for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator eti = tm_.ortae_.e2t_.begin();
      eti != tm_.ortae_.e2t_.end(); ++eti ) {
      const vector<size_t> &loop = eti->second;
      const pair<size_t,size_t> & edge = eti->first;
      if(loop.front() != loop.back()) // error format or surface edge or surface edge
        continue;
      if(loop.front() == -1) // open loop
        continue;

      transition = eye<double>(3);

      for(size_t i = 0; i < loop.size()-1; ++i) {
          get_best_alignment(&frame[loop[i]][0], &frame[loop[i+1]][0], &shuffle_rot[0]);
          transition = temp(transition * shuffle_rot);
        }

      // label the singularity edge
      if(norm(transition - eye<double>(3)) > 1e-8) {
          singularities_segment.insert(edge);
          const size_t type = type_transition1(transition);
          singularity_type[edge] = type;
          singularity_type[make_pair(edge.second, edge.first)] = get_trans_type(type);
        }
    }

  extract_chain_from_edges_with_outside_points(singularities_segment,tm_.outside_face_, chain);

  singularities_type.resize(chain.size());
  for(size_t ci = 0; ci < chain.size(); ++ci){
      const std::deque<std::pair<size_t, size_t> > & one_chain = chain[ci];
      for(const auto & one_edge : one_chain){
          auto it = singularity_type.find(one_edge);
          if(it == singularity_type.end())
            throw std::logic_error("strange , can not find singularity type.");
          singularities_type[ci].push_back(it->second);
        }
    }
}


void singularity_extractor::extract(const boost::unordered_map<std::pair<size_t, size_t>, size_t> &inner_type,
                                    std::vector<std::deque<std::pair<size_t, size_t> > > &chain,
                                    std::vector<std::deque<size_t> > &singularities_type,
                                    bool with_unknown_face_type) const
{
  chain.clear();
  singularities_type.clear();
  matrixd rot = eye<double>(3);
  boost::unordered_set<pair<size_t,size_t> > singularities_segment;
  boost::unordered_map<pair<size_t,size_t>,size_t> singularity_type;

  size_t t;
  for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator it = tm_.ortae_.e2t_.begin();
      it != tm_.ortae_.e2t_.end(); ++it){
      const vector<size_t> &loop = it->second;
      if(loop.front() == -1) continue; // open loop
      const pair<size_t,size_t> &edge = it->first;

      assert(loop.front() == loop.back());
      rot = eye<double>(3);
      for(t = 0; t < loop.size() -1; ++t){
          auto it_ = inner_type.find(make_pair(loop[t],loop[t+1]));
          if((it_ == inner_type.end()) && with_unknown_face_type) break;
          if((it_ == inner_type.end()) && !with_unknown_face_type) continue;
          rot = temp(rot * type_transition2(it_->second));
        }
      if(t != loop.size() - 1) {// breaked
          continue;}
      if(fabs(norm(rot - eye<double>(3))) > 1e-8) {
          singularities_segment.insert(it->first);
          const size_t type = type_transition1(rot);
          singularity_type[edge] = type;
          singularity_type[make_pair(edge.second, edge.first)] = get_trans_type(type);
        }
    }

  extract_chain_from_edges_with_outside_points(singularities_segment,tm_.outside_face_, chain);

  singularities_type.resize(chain.size());
  for(size_t ci = 0; ci < chain.size(); ++ci){
      const std::deque<std::pair<size_t, size_t> > & one_chain = chain[ci];
      for(const auto & one_edge : one_chain){
          auto it = singularity_type.find(one_edge);
          if(it == singularity_type.end())
            throw std::logic_error("strange , can not find singularity type.");
          singularities_type[ci].push_back(it->second);
        }
    }
}


// It's an inefficient data structure, for fast implementation currently.
typedef map<pair<size_t, size_t>, vector<size_t> > edge2face; // edge2face[(e0, e1)] = {f0,f1}


void make_e2f(const jtf::mesh::face2tet_adjacent &fa, const matrixst &face_idx, edge2face &e2f)
{
  for(size_t fi = 0; fi < face_idx.size(); ++fi) {
      const size_t surf_face_id = face_idx[fi];
      const vector<size_t> &face = fa.faces_[surf_face_id];
      for(int vi = 0; vi < 3; ++vi)
        e2f[make_edge(face[vi], face[(vi+1)%3])].push_back(surf_face_id);
    }
}

size_t get_edge_type_with_part_face_type_map(
    const std::vector<size_t> & tet_seq,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_type )
{
  typedef  boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
  matrixd rot = eye<double>(3);
  for(size_t t = 0; t < tet_seq.size()-1; ++t){
      mpscit it = face_type.find(make_pair(tet_seq[t], tet_seq[t+1]));
      if(it == face_type.end()) continue;
      rot = temp(rot * type_transition2(it->second));
    }
  return type_transition1(rot);
}

size_t get_edge_type_with_part_face_type_map(
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const std::pair<size_t,size_t> & edge,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_type )
{
  typedef  jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator ecit;
  ecit it = ortae.e2t_.find(edge);
  bool need_reverse = false;
  if(it == ortae.e2t_.end()){
      it = ortae.e2t_.find(make_pair(edge.second, edge.first));
      if(it == ortae.e2t_.end())
        return -1;
      need_reverse = true;
    }
  const vector<size_t> & loop = it->second;
  if(!ortae.is_inner_edge(loop)) return -1;
  assert(loop.front() == loop.back());

  size_t rtn = get_edge_type_with_part_face_type_map(loop, face_type);
  if(need_reverse) return get_trans_type(rtn);
  else return rtn;
}

size_t get_edge_type_with_full_face_type_map(
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const std::pair<size_t,size_t> & edge,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_type)
{
  typedef  jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator ecit;
  ecit it = ortae.e2t_.find(edge);
  bool need_reverse = false;
  if(it == ortae.e2t_.end()){
      it = ortae.e2t_.find(make_pair(edge.second, edge.first));
      if(it == ortae.e2t_.end())
        return -1;
      need_reverse = true;
    }
  const vector<size_t> & loop = it->second;
  if(!ortae.is_inner_edge(loop)) return -1;
  assert(loop.front() == loop.back());

  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
  matrixd rot = eye<double>(3);
  for(size_t t = 0; t < loop.size()-1; ++t){
      mpscit it = face_type.find(make_pair(loop[t], loop[t+1]));
      if(it == face_type.end()) return -1;
      rot = temp(rot * type_transition2(it->second));
    }
  if(need_reverse) trans(rot);
  return type_transition1(rot);
}

int extract_chain_from_edges_with_outside_points(
    const boost::unordered_set<pair<size_t,size_t> > & segments,
    const matrixst &outside_face,
    vector<deque<pair<size_t,size_t> > > &chain_list)
{
  chain_list.clear();
  boost::unordered_set<size_t> outside_points(outside_face.begin(), outside_face.end());
  boost::unordered_map<size_t,set<size_t> > point_count;
  boost::unordered_set<pair<size_t,size_t> > segments_temp = segments;

  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit = segments_temp.begin();
      cit != segments_temp.end(); ++cit){
      const pair<size_t,size_t> & edge = *cit;
      point_count[edge.first].insert(edge.second);
      point_count[edge.second].insert(edge.first);
    }

  while(!segments_temp.empty()){
      deque<pair<size_t,size_t> > one_chain;
      one_chain.push_back(*segments_temp.begin());
      segments_temp.erase(segments_temp.begin());
      while(1){
          const size_t chain_length = one_chain.size();
          for(boost::unordered_set<pair<size_t,size_t> >::iterator it = segments_temp.begin();
              it != segments_temp.end(); ){
              const pair<size_t,size_t> & edge = *it;
              if(outside_points.find(one_chain.front().first) == outside_points.end()){ // not surface point
                  boost::unordered_map<size_t,set<size_t> >::const_iterator it_first =
                      point_count.find(one_chain.front().first);
                  assert(it_first != point_count.end());
                  if(it_first->second.size() == 2) {
                      if(edge.first == one_chain.front().first){
                          one_chain.push_front(make_pair(edge.second, edge.first));
                          segments_temp.erase(it++);
                          continue;
                        } else if(edge.second == one_chain.front().first){
                          one_chain.push_front(edge);
                          segments_temp.erase(it++);
                          continue;
                        }
                    }
                }
              if(outside_points.find(one_chain.back().second) == outside_points.end()){
                  boost::unordered_map<size_t,set<size_t> >::const_iterator it_second =
                      point_count.find(one_chain.back().second);
                  assert(it_second != point_count.end());
                  if(it_second->second.size() == 2){
                      if(edge.first == one_chain.back().second){
                          one_chain.push_back(edge);
                          segments_temp.erase(it++);
                          continue;
                        }
                      else if(edge.second == one_chain.back().second){
                          one_chain.push_back(make_pair(edge.second, edge.first));
                          segments_temp.erase(it++);
                          continue;
                        }
                    }
                }
              ++it;
            }
          if(one_chain.size() == chain_length){
              break;
            }
        }
      chain_list.push_back(one_chain);
    }
  cerr << "# [info] chain num: " << chain_list.size() << endl;
  return 0;
}


int find_singularities_global_align_with_type(
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const matrixst & tet_rot_type,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t > &inner_face_jump,
    const matrixst &outside_face,
    std::vector<std::deque<std::pair<size_t,size_t> > > &chain_list,
    std::vector<deque<size_t> > &singularities_type)
{
  chain_list.clear();
  singularities_type.clear();
  matrixd rot = eye<double>(3);
  boost::unordered_set<pair<size_t,size_t> > singularities_segment;

  typedef boost::unordered_map<pair<size_t,size_t>,size_t >::const_iterator mcit;
  for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator it = ortae.e2t_.begin();
      it != ortae.e2t_.end(); ++it){
      const vector<size_t> &loop = it->second;
      if(loop.size() == 2) continue; // sharp edge
      if(loop.front() == -1) continue; // open loop
      const pair<size_t,size_t> &edge = it->first;


      assert(loop.front() == loop.back());
      rot = eye<double>(3);
      for(size_t t = 0; t < loop.size() -1; ++t){
          mcit it_ = inner_face_jump.find(make_pair(loop[t],loop[t+1]));
          if(it_ == inner_face_jump.end()) continue; // means identity
          rot = temp(rot * type_transition2(it_->second));
        }

      rot = temp((type_transition2(tet_rot_type[loop.front()])) * rot);
      rot = temp(rot * trans(type_transition2(tet_rot_type[loop.front()])));

      if(fabs(norm(rot - eye<double>(3))) > 1e-8) {
          singularities_segment.insert(it->first);
        }
    }

  extract_chain_from_edges_with_outside_points(singularities_segment,outside_face,chain_list);
  singularities_type.resize(chain_list.size());

  for(size_t t = 0; t < chain_list.size(); ++t){
      const deque<pair<size_t,size_t> > &single_chain = chain_list[t];
      for(size_t i = 0; i < single_chain.size(); ++i){
          jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator it =
              ortae.e2t_.find(single_chain[i]);
          if(it == ortae.e2t_.end())
            it =ortae.e2t_.find(make_pair(single_chain[i].second,single_chain[i].first));
          if(it == ortae.e2t_.end()) {
              cerr << "# [error] can not find this edge: " << single_chain[i].first
                   << "-->" << single_chain[i].second << endl;
              return __LINE__;
            }
          const vector<size_t> &loop = it->second;
          rot = eye<double>(3);
          for(size_t j = 0; j < loop.size() -1; ++j){
              mcit it_ = inner_face_jump.find(make_pair(loop[j],loop[j+1]));
              if(it_ == inner_face_jump.end()) continue; // means identity
              rot = temp(rot * type_transition2(it_->second));
            }
          rot = temp((type_transition2(tet_rot_type[loop.front()])) * rot);
          rot = temp(rot * trans(type_transition2(tet_rot_type[loop.front()])));

          singularities_type[t].push_back(type_transition1(rot));
        }
    }
  return 0;
}


void cal_singularity_adj_frame(const matrixst & tet,
                               const matrixd & node,
                               const size_t tet_idx,
                               const matrixd &frame_rot,
                               matrixd & cross_frame_edge)
{
  assert(frame_rot.size(2) == 3 && frame_rot.size(1) == 3);
  assert(cross_frame_edge.size(1) == 3 && cross_frame_edge.size(2) == 7)    ;

  matrixd bary_center = zeros<double>(3,1);

  for(size_t t = 0 ; t < 4; ++t)
    bary_center += node(colon(),tet(t,tet_idx));

  bary_center /= 4.0;
  //static const matrixd I = eye<double>(3);
  const double len = norm(bary_center - node(colon(),tet(0,tet_idx))) / 2.0;
  cross_frame_edge(colon(),0) = bary_center + len * frame_rot(colon(),0);
  cross_frame_edge(colon(),1) = bary_center - len * frame_rot(colon(),0);
  cross_frame_edge(colon(),2) = bary_center + len * frame_rot(colon(),1);
  cross_frame_edge(colon(),3) = bary_center - len * frame_rot(colon(),1);
  cross_frame_edge(colon(),4) = bary_center + len * frame_rot(colon(),2);
  cross_frame_edge(colon(),5) = bary_center - len * frame_rot(colon(),2);
  cross_frame_edge(colon(),6) = bary_center ;
}


void draw_singularity_adj_frame_full(ofstream &ofs,
                                     const matrixst & tet,
                                     const matrixd & node,
                                     const matrix<matrixd > &frame_rot,
                                     const vector<vector<size_t> >singularities_tet_loop
                                     )
{
  size_t adj_tet_num = 0;
  for(size_t t = 0; t < singularities_tet_loop.size(); ++t)
    for(size_t p = 0; p < singularities_tet_loop[t].size(); ++p)
      adj_tet_num += 1;
  matrix<matrixd > singularity_adj_frame(adj_tet_num);


  for(size_t t = 0,i = 0; t < singularities_tet_loop.size(); ++t)
    for(size_t p = 0; p < singularities_tet_loop[t].size(); ++p, ++i)
      {
        const size_t &tet_idx = singularities_tet_loop[t][p];
        singularity_adj_frame[i].resize(3,7);
        cal_singularity_adj_frame(tet,node,tet_idx,frame_rot[tet_idx],singularity_adj_frame[i]);
      }
  matrixst edges(2,adj_tet_num * 6);
  for(size_t t = 0; t < adj_tet_num ; ++t)
    for(size_t p = 0; p < 6; ++p)
      {
        edges(0, t * 6 + p) = t * 7 + p;
        edges(1, t * 6 + p) = t * 7 + 6;
      }

  const double color_table[6][4] = {{255,0,0,255},
                                    {255,160,160,255},

                                    {0,255,0,255},
                                    {160,255,160,255},

                                    {0,0,255,255},
                                    {160,160,255,255}};
  itr_matrix<const double *> color(4,6,&color_table[0][0]);

  matrixd color_edges(4,adj_tet_num * 6);
  for(size_t t = 0; t < adj_tet_num; ++t)
    {
      for(size_t p = 0; p < 6; ++p)
        color_edges(colon(), t * 6 + p) = color(colon(),p);
    }
  color_edges /= 255.0;
  matrixd saf(3, adj_tet_num * 7);
  for(size_t t = 0; t < adj_tet_num; ++t)
    {
      for(size_t p = 0; p < 7; ++p)
        saf(colon(),t * 7 + p) = singularity_adj_frame[t](colon(),p);
    }
  line2vtk(ofs,&saf[0], adj_tet_num * 7, &edges[0], adj_tet_num * 6 );
  cell_data_rgba(ofs,&color_edges[0],adj_tet_num * 6, "frame_color");
}

int find_surface_singularity_edge_using_normal_align_type(
    const matrixst & outside_face,
    const vector<size_t> & outside_face_align_type,
    const jtf::mesh::edge2cell_adjacent & e2t,
    std::vector<std::deque<std::pair<size_t,size_t> > > & singularity_chain,
    std::vector<std::deque<size_t> > & singularity_type)
{
  vector<pair<size_t,size_t> > singularity_edges;

  boost::unordered_map<pair<size_t,size_t>,size_t> edge2type;

  for(size_t ei = 0; ei < e2t.edge2cell_.size(); ++ei){
      const pair<size_t,size_t> & edge = e2t.edges_[ei];
      const pair<size_t,size_t> & two_faces = e2t.edge2cell_[ei];
      if(e2t.is_boundary_edge(two_faces)) continue;
      if(outside_face_align_type[two_faces.first]
         != outside_face_align_type[two_faces.second]){
          singularity_edges.push_back(edge);
          // since the type is u,-u,v,-v,w,-w
          const size_t free_type = (0+1+2) -
              //                               (outside_face_align_type[two_faces.second] +
              //                                outside_face_align_type[two_faces.first]);

              (axis_to_around(outside_face_align_type[two_faces.second]) +
              axis_to_around(outside_face_align_type[two_faces.first]));
          assert(free_type < 3);
          if(edge.first > edge.second)
            edge2type[make_pair(edge.second, edge.first)] = free_type;
          else
            edge2type[edge] = free_type;
        }
    }

  jtf::util::extract_chain_from_edges(singularity_edges, singularity_chain);
  singularity_type.resize(singularity_chain.size());

  for(size_t chi = 0; chi < singularity_chain.size(); ++chi){
      singularity_type[chi].resize(singularity_chain[chi].size());
      for(size_t ei = 0; ei < singularity_chain[chi].size(); ++ei){
          pair<size_t,size_t> one_edge = singularity_chain[chi][ei];

          if(one_edge.first > one_edge.second)
            swap(one_edge.first, one_edge.second);

          boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator
              bumcit = edge2type.find(one_edge);
          if(bumcit == edge2type.end()){
              cerr << "# [error] strange can not find edge pair <" << one_edge.first
                   << "," << one_edge.second << endl;
              return __LINE__;
            }
          singularity_type[chi][ei] = bumcit->second;
        }
    }
  return 0;

}

void draw_singularity_adj_frame_sparse(ofstream &ofs,
                                       const matrixst & tet,
                                       const matrixd & node,
                                       const matrix<matrixd > &frame_rot,
                                       const vector<vector<size_t> >singularities_tet_loop
                                       )
{
  const size_t adj_tet_num = singularities_tet_loop.size();
  matrix<matrixd > singularity_adj_frame(adj_tet_num);


  for(size_t t = 0; t < singularities_tet_loop.size(); ++t)
    {
      const size_t &tet_idx = singularities_tet_loop[t][0];
      singularity_adj_frame[t].resize(3,7);
      cal_singularity_adj_frame(tet,node,tet_idx,frame_rot[tet_idx],singularity_adj_frame[t]);
    }
  matrixst edges(2,adj_tet_num * 6);
  for(size_t t = 0; t < adj_tet_num ; ++t)
    for(size_t p = 0; p < 6; ++p)
      {
        edges(0, t * 6 + p) = t * 7 + p;
        edges(1, t * 6 + p) = t * 7 + 6;
      }

  const double color_table[6][4] = {{255,0,0,255},
                                    {255,160,160,255},

                                    {0,255,0,255},
                                    {160,255,160,255},

                                    {0,0,255,255},
                                    {160,160,255,255}};
  itr_matrix<const double *> color(4,6,&color_table[0][0]);

  matrixd color_edges(4,adj_tet_num * 6);
  for(size_t t = 0; t < adj_tet_num; ++t)
    {
      for(size_t p = 0; p < 6; ++p)
        color_edges(colon(), t * 6 + p) = color(colon(),p);
    }
  color_edges /= 255.0;
  matrixd saf(3, adj_tet_num * 7);
  for(size_t t = 0; t < adj_tet_num; ++t)
    {
      for(size_t p = 0; p < 7; ++p)
        saf(colon(),t * 7 + p) = singularity_adj_frame[t](colon(),p);
    }
  line2vtk(ofs,&saf[0], adj_tet_num * 7, &edges[0], adj_tet_num * 6 );
  cell_data_rgba(ofs,&color_edges[0],adj_tet_num * 6, "frame_color");
}

static void replace_vertex(const vector<size_t> &vertex_map,
                           matrixst &new_tet)
{
  for(size_t t = 0; t < vertex_map.size(); ++t){
      if(vertex_map[t] != t){
          replace(new_tet.begin(),new_tet.end(),t,vertex_map[t]);
        }
    }
}

static void check_tet_face_valid(const matrixst &new_tet,
                                 const matrixd *node = 0)
{
  map<vector<size_t>,vector<size_t> > map_face_v_idx;
  for(size_t t = 0; t < new_tet.size(2); ++t)
    {
      size_t face[3];
      vector<size_t> face_v(3);
      for(size_t i = 0; i < 4; ++i){
          for(size_t j = 0; j < 3; ++j)
            face[j] = new_tet((j + i)%4,t);

          copy(face,face+3,face_v.begin());
          sort(face_v.begin(),face_v.end());
          map_face_v_idx[face_v].push_back(t);
        }
    }

  typedef map<vector<size_t>,vector<size_t> >::iterator mit;
  size_t single_face_num = 0;
  vector<size_t> error_tet;
  error_tet.reserve(24);
  //static int num = 0;
  vector<size_t> face(3);
  for(mit it = map_face_v_idx.begin();
      it !=map_face_v_idx.end(); ++it ) {
      if(it->second.size() > 2){
          cerr << "# strang face: ";
          copy(it->first.begin(), it->first.end(), ostream_iterator<size_t>(cerr, " "));
          copy(it->first.begin(), it->first.end(),face.begin());
          cerr << endl;
          for(size_t i = 0; i < it->second.size(); ++i) {
              cerr << "# in tet: " << it->second[i];
              const matrixst & tet_= new_tet(colon(),it->second[i]);
              //copy(tet_.begin(), tet_.end(), error_tet.begin() + num);
              for(size_t t = 0; t < 4; ++t)
                error_tet.push_back(tet_[t]);

              cerr << tet_ ;
            }
          break;
        }
      if(it->second.size() == 1)
        ++single_face_num ;
    }
  cerr << "# error tet size " << error_tet.size() << endl;
  ofstream ofs("error_tet.vtk");
  tet2vtk(ofs, &(*node)[0],node->size(2),&error_tet[0], error_tet.size()/4);
  ofstream ofs_f("error_tet_face.vtk");
  tri2vtk(ofs_f, &(*node)[0],node->size(2),&face[0],1);
}

static void check_tet_valid(const matrixst &merged_tet,
                            vector<bool> &is_valid_tet)
{
  map<vector<size_t>,size_t> map_;
  typedef map<vector<size_t>,size_t>::iterator mit;
  vector<size_t> tmp(4);
  for(size_t t = 0; t < merged_tet.size(2); ++t){
      if(is_valid_tet[t]){
          copy(merged_tet(colon(),t).begin(),merged_tet(colon(),t).end(),tmp.begin());
          sort(tmp.begin(),tmp.end());
          mit i = map_.find(tmp);
          if(i != map_.end())
            is_valid_tet[t] = false;
          else
            map_[tmp] = t;
        }
    }
}

#include <sstream>

static void merge_black_edge(const size_t from,
                             const size_t to,
                             vector<bool> &is_valid_tet,
                             matrixst &merged_tet,
                             const matrixd *node = 0)
{
  replace(merged_tet.begin(),merged_tet.end(),from,to);

  for(size_t t = 0 ;t < merged_tet.size(2); ++t){
      if(is_valid_tet[t]){
          set<size_t> tmp;
          for(size_t i = 0; i < 4; ++i)
            tmp.insert(merged_tet(i,t));
          if(tmp.size() != 4) // degenerated tet
            is_valid_tet[t] = false;
        }
    }

#if 1
  size_t num = 0;
  for(size_t t = 0; t < is_valid_tet.size(); ++t)
    if(is_valid_tet[t]) ++num;
  matrixst new_tet(4,num);
  for(size_t t = 0,i = 0; t < is_valid_tet.size(); ++t){
      if(is_valid_tet[t])
        new_tet(colon(),i++) = merged_tet(colon(),t);
    }
  check_tet_face_valid(new_tet,node);
  //    std::stringstream Num;
  //    Num << from;
  //    string tet_name = "new_tet_";
  //    tet_name += Num.str();
  //    tet_name += ".vtk";

  //    ofstream ofs(tet_name.c_str());
  //    tet2vtk(ofs,&(*node)[0],node->size(2),&new_tet[0],new_tet.size(2));
#endif
  //check_tet_valid(merged_tet,is_valid_tet);
}


static void merge_black_chain(const deque<pair<size_t,size_t> > &dqp,
                              const size_t end,
                              vector<bool> &is_valid_tet,
                              matrixst &merged_tet,
                              const matrixd &node)
{
  if(end == dqp.front().first){
      for(size_t t = dqp.size() - 1; t  != -1; --t)
        merge_black_edge(dqp[t].second,dqp[t].first,is_valid_tet,merged_tet);
    }else if(end == dqp.back().second) {
      for(size_t t = 0; t < dqp.size(); ++t)
        merge_black_edge(dqp[t].first,dqp[t].second,is_valid_tet,merged_tet,&node);
    }
  else{
      size_t t;
      for(t = 0; dqp[t].first != end; ++t) {
          merge_black_edge(dqp[t].first,dqp[t].second,is_valid_tet,merged_tet);
        }
      for(size_t i = dqp.size() - 1; i >= t ; --i){
          merge_black_edge(dqp[i].second,dqp[i].first,is_valid_tet,merged_tet);
        }
    }
}


int remove_black_lines_and_rebuild(const matrixst &tet,
                                   const matrixd & node,
                                   const jtf::mesh::face2tet_adjacent &fa,
                                   const vector<size_t> & singularity_type,
                                   const vector<pair<size_t, size_t> > &singularities,
                                   const size_t black_type,
                                   matrixst &new_tet,
                                   matrixst &map_from_new_to_old)
{
  assert(singularity_type.size() == singularities.size());

  matrixst outside_face;
  get_outside_face(fa,outside_face);

  set<pair<size_t,size_t> > need_to_merge;

  for(size_t t = 0 ; t < singularity_type.size(); ++t){
      if(singularity_type[t] == black_type){
          need_to_merge.insert(singularities[t]);
        }
    }
  if(need_to_merge.empty())
    return 1;
  vector<deque<pair<size_t,size_t> > > merge_edge;

  // after reconstruct, each black chain will not have surface vertex interior
  // if there are outside vertex, these vertex will only exist on the end of chain

  while(!need_to_merge.empty())
    {
      deque<pair<size_t,size_t> > black_chain;
      set<pair<size_t,size_t> >::iterator sit = need_to_merge.begin();
      black_chain.push_back(*sit);
      need_to_merge.erase(sit);
      while(1)
        {
          size_t black_chain_size = black_chain.size();
          for(set<pair<size_t,size_t> >::iterator i = need_to_merge.begin();
              i != need_to_merge.end(); ++i)
            {
              const pair<size_t,size_t> & bedge = *i;
              if(!is_outside_vertex(black_chain.back().second,outside_face)
                 && (bedge.first == black_chain.back().second))
                {
                  black_chain.push_back(bedge);
                  need_to_merge.erase(i);
                }
              else if (!is_outside_vertex(black_chain.front().first,outside_face)
                       && (bedge.second == black_chain.front().first)) {
                  black_chain.push_front(bedge);
                  need_to_merge.erase(i);
                }
              else if(is_outside_vertex(black_chain.back().second,outside_face)
                      && is_outside_vertex(black_chain.front().first,outside_face))
                break;
            }
          if(black_chain.size() == black_chain_size)
            break;
        }
      merge_edge.push_back(black_chain);
    }
  cerr << "# find " << merge_edge.size() << " black line chain." << endl;


  matrixst merged_tet = tet;
  vector<size_t> vertex_map(node.size(2));
  vector<bool> is_valid_tet(tet.size(2),true);
  for(size_t t = 0; t < node.size(2); ++t) vertex_map[t] = t;
  for(size_t t = 0; t < merge_edge.size(); ++t) {
      cerr << "# black_chain "<< t << " has " << merge_edge[t].size() << " edges." << endl;
      //copy(merge_edge[t].begin(),merge_edge[t].end(),ostream_iterator<size_t>(cerr," "));
#if 1 // debug dump out edge
      for(size_t i = 0; i < merge_edge[t].size(); ++i)
        cerr << merge_edge[t][i].first << " ";
      cerr << merge_edge[t].back().second << endl;
      cerr << endl;
#endif

#if 0 // check
      if(merge_edge[t].size() > 1)
        {
          size_t tmp = 0;
          for(size_t i = 0; i < merge_edge[t].size(); ++i){
              if(i == 0)
                tmp += merge_edge[t][i].second;
              else if(i == merge_edge[t].size() - 1) {
                  tmp -= merge_edge[t][i].first;
                } else {
                  tmp -= merge_edge[t][i].first;
                  tmp += merge_edge[t][i].second;
                }
            }
          if(tmp != 0) {
              cerr << "# strange black line chain." << endl;
              return __LINE__;
            }
        }
#endif
      // size_t begin = -1;
      size_t end = -1;

      const deque<pair<size_t,size_t> > &dqp = merge_edge[t];
      if(is_outside_vertex(dqp.front().first,outside_face)
         && is_outside_vertex(dqp.back().second,outside_face) ) {
          end = dqp.front().second;
          // begin = merge_edge[t].back.second;
        }else if(is_outside_vertex(dqp.front().first,outside_face)){
          end = dqp.front().first;
          //begin = merge_edge[t].front.first;
        }else if(is_outside_vertex(dqp.back().second,outside_face))
        end = dqp.back().second;
      else
        end = dqp.back().second;
      //#define new_recursive
#ifndef new_recursive
      for(size_t t = 0; t < dqp.size(); ++t){
          vertex_map[dqp[t].first] = end;
          vertex_map[dqp[t].second] = end;
        }
#endif
#ifdef new_recursive
      merge_black_chain(dqp,end,is_valid_tet,merged_tet,node);
#endif
    }


#ifndef new_recursive
  // recursively adjust the vertex_map
  for(size_t i = 0; i < vertex_map.size(); ++i){
      if(vertex_map[i] != i){
          while(vertex_map[vertex_map[i]] != vertex_map[i])
            vertex_map[i] = vertex_map[vertex_map[i]];
        }
    }


  replace_vertex(vertex_map,merged_tet);
#endif
#if 0
  for(size_t t = 0; t < merge_edge.size(); ++t) {
      for(size_t i = 0; i < merge_edge[t].size(); ++i){
          const pair<size_t,size_t> &edge = merge_edge[t][i];
          if(find(merged_tet.begin(),merged_tet.end(),edge.first) != merged_tet.end())
            cerr << "# error" << endl;
          if(find(merged_tet.begin(),merged_tet.end(),edge.second) != merged_tet.end())
            cerr << "# error" << endl;
        }
    }
#endif
#ifndef new_recursive
  map<vector<size_t>,size_t> map_tet_v_idx;

  vector<size_t> map_;
  map_.reserve(tet.size(2));
  size_t degenerate_tet_num = 0;
  for(size_t t = 0; t < merged_tet.size(2); ++t)
    {
      set<size_t> tmp_;
      for(size_t k = 0; k < 4; ++k) {
          tmp_.insert(merged_tet(k,t));

        }
      if(tmp_.size() == 4) {
          vector<size_t> tmp(4);
          copy(merged_tet(colon(),t).begin(),merged_tet(colon(),t).end(),tmp.begin());
          sort(tmp.begin(),tmp.end());

          if(map_tet_v_idx.find(tmp) == map_tet_v_idx.end())
            {
              map_tet_v_idx[tmp] = t;
              map_.push_back(t);
            }
        }else
        ++degenerate_tet_num;
    }
  cerr << "# degenerate_tet_num " << degenerate_tet_num << endl;
  //swap(map,map_from_new_to_old);
  //if(!map_from_new_to_old.empty()) map_from_new_to_old.clear();
  map_from_new_to_old.resize(map_.size());
  copy(map_.begin(),map_.end(),map_from_new_to_old.begin());
  //matrixst merged_tet(4,map.size());

  new_tet.resize(4,map_.size());
  for(size_t t = 0; t < map_.size(); ++t) {
      new_tet(colon(),t) = merged_tet(colon(),map_[t]);
    }
#endif

#ifdef new_recursive
  vector<size_t> map_;
  map_.reserve(tet.size(2));
  for(size_t i  = 0; i < is_valid_tet.size(); ++i)
    if(is_valid_tet[i]) map_.push_back(i);

  map_from_new_to_old.resize(map_.size());
  copy(map_.begin(),map_.end(),map_from_new_to_old.begin());
  //matrixst merged_tet(4,map.size());

  new_tet.resize(4,map_.size());
  for(size_t t = 0; t < map_.size(); ++t) {
      new_tet(colon(),t) = merged_tet(colon(),map_[t]);
    }
#endif

#if 0 // debug
  ofstream ofs("new_tet.vtk");
  tet2vtk(ofs, &node[0],node.size(2),&new_tet[0], new_tet.size(2));
  // ofstream ofs_("input_tet.vtk");
  // tet2vtk(ofs_, &node[0],node.size(2),&tet[0], tet.size(2));
#endif

#if 1 // check the face to tet
  map<vector<size_t>,vector<size_t> > map_face_v_idx;
  for(size_t t = 0; t < new_tet.size(2); ++t)
    {
      size_t face[3];
      vector<size_t> face_v(3);
      for(size_t i = 0; i < 4; ++i){
          for(size_t j = 0; j < 3; ++j)
            face[j] = new_tet((j + i)%4,t);

          copy(face,face+3,face_v.begin());
          sort(face_v.begin(),face_v.end());
          map_face_v_idx[face_v].push_back(t);
        }
    }

  typedef map<vector<size_t>,vector<size_t> >::iterator mit;
  size_t single_face_num = 0;
  vector<size_t> error_tet;
  error_tet.reserve(24);
  for(mit it = map_face_v_idx.begin();
      it !=map_face_v_idx.end(); ++it ) {
      if(it->second.size() > 2){

          cerr << "# strang face: ";
          copy(it->first.begin(), it->first.end(), ostream_iterator<size_t>(cerr, " "));
          cerr << endl;
          for(size_t i = 0; i < it->second.size(); ++i) {
              cerr << "# in tet: " << it->second[i];
              const matrixst & tet_= new_tet(colon(),it->second[i]);
              //copy(tet_.begin(), tet_.end(), error_tet.begin() + num);
              for(size_t t = 0; t < 4; ++t)
                error_tet.push_back(tet_[t]);
              cerr << tet_ ;
            }
        }
      if(it->second.size() == 1)
        ++single_face_num ;
    }
  cerr << "# error tet size " << error_tet.size() << endl;
  ofstream ofs("error_tet.vtk");
  tet2vtk(ofs, &node[0],node.size(2),&error_tet[0], error_tet.size()/4);
#endif
#if 0 //check the merged_tet
  for(size_t t = 0; t < new_tet.size(2); ++t)
    {
      set<size_t> tmp_;
      for(size_t k = 0; k < 4; ++k) {
          tmp_.insert(new_tet(k,t));
        }
      if(tmp_.size() != 4) {
          cerr << "# error tet " << t << endl;
          cerr << new_tet(colon(),t);
        }
    }
#endif
#if 0  // check whether the new_tet is correct
  unique_ptr<jtf::mesh::face2tet_adjacent> fa_new(jtf::mesh::face2tet_adjacent::create(new_tet));
  matrixst out_side_face_new;
  get_outside_face_idx(*fa_new,out_side_face_new);
  if(out_side_face_new.size() != single_face_num)
    {
      cerr << "# incorrect tet " << endl;
      cerr << "# outsideface " << out_side_face_new << endl;
      cerr << "# single face " << single_face_num << endl;
    }
#endif

  return 0;
}


int check_singularity_edge_satisify_cubecover(const vector<deque<pair<size_t,size_t> > > &chain_list,
                                              const matrixst &outside_face)
{
  map<size_t,size_t> ends;
  for(size_t t = 0; t < chain_list.size(); ++t){
      if(chain_list[t].front().first == chain_list[t].back().second)
        continue;
      if(ends.find(chain_list[t].front().first) == ends.end())
        ends[chain_list[t].front().first] = 1;
      else
        ends[chain_list[t].front().first] += 1;

      if(ends.find(chain_list[t].back().second) == ends.end())
        ends[chain_list[t].back().second] = 1;
      else
        ends[chain_list[t].back().second] += 1;
    }

  bool is_satisfy_cubecover = true;
  for(size_t t = 0; t < chain_list.size(); ++t){
      if(chain_list[t].front().first == chain_list[t].back().second)
        continue;
      if(!is_outside_vertex(chain_list[t].front().first,outside_face )
         && ends[chain_list[t].front().first] == 1)
        {
          is_satisfy_cubecover = false;
          cerr << "# find inner singularity end " << chain_list[t].front().first << endl;
        }
      if(!is_outside_vertex(chain_list[t].back().second,outside_face )
         && ends[chain_list[t].back().second] == 1)
        {
          is_satisfy_cubecover = false;
          cerr << "# find inner singularity end " << chain_list[t].back().second << endl;
        }
    }
  if(is_satisfy_cubecover)
    cerr << "# satisfy the cubecover command." << endl;
  return 0;
}
