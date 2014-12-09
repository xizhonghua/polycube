#include "singularity_adjustment.h"
#include "../common/vtk.h"
#include "../common/IO.h"
#include "../common/util.h"
#include "../common/timer.h"
#include "common.h"
#include "../common/visualize_tool.h"
#include "../common/transition_type.h"
#include "topology_operation.h"
#include "topology_analysis.h"
#include "../numeric/util.h"
#include "../tetmesh/util.h"
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <numeric>

using namespace std;
using namespace zjucad::matrix;
using namespace boost::tuples;

int relabel_to_converge_singularities(
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::face2tet_adjacent& fa,
    std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    std::vector<std::deque<size_t> > &singularity_type,
    std::map<std::pair<size_t,size_t>, size_t> &inner_face_jump_type)
{
  cerr << "# [error] not finished." << endl;
  return __LINE__;
}

int relabel_singularity_chain_by_modify_face_jump_type(
    matrixst &tet,
    matrixd &node,
    matrixst &outside_face,
    matrixst &outside_face_idx,
    matrix<matrixd > & frame,
    jtf::mesh::face2tet_adjacent &fa,
    jtf::mesh::one_ring_tet_at_edge &ortae,
    vector<deque<pair<size_t,size_t> > > &singularity_edges,
    vector<deque<size_t> > &singularity_type,
    boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type)
{
  // this step may remove all zigzag on each chain
  relabel_zigzag_by_modify_face_jump_type(tet,node,outside_face,fa,ortae,
                                          singularity_edges,singularity_type,
                                          inner_face_jump_type);


  cerr << "# [info] finish removing zigzag." << endl;
#if 1 // visual
  {
    jtf::mesh::meshes tm_mesh;
    tm_mesh.mesh_ = tet;
    tm_mesh.node_ = node;
    jtf::tet_mesh tm(tm_mesh);
    singularity_extractor se(tm);
    vector<deque<pair<size_t,size_t> > > singularity_edges;
    std::vector<deque<size_t> > singularities_type;
    se.extract(inner_face_jump_type, singularity_edges, singularities_type);

    dump_singularity_chain_to_vtk_2("after_remove_zigzag.vtk",node, singularity_edges,
                                    singularity_type);
  }
#endif
  {
    vector<size_t> tet_array(tet.size());
    copy(tet.begin(), tet.end(), tet_array.begin());
    vector<double> node_array(node.size());
    copy(node.begin(), node.end(), node_array.begin());
    vector<matrixd > frame_array(frame.size());
    copy(frame.begin(), frame.end(), frame_array.begin());

    timer tr;
    tr.start();
    relabel_remove_black_lines_by_modify_face_jump_type_and_splitting_fans(
          tet_array,node_array,frame_array,fa,outside_face,ortae,singularity_edges,
          singularity_type,inner_face_jump_type);

    tr.finish();
    cerr << "# [compound time] cost " << tr.result()
         << "ms to remove inner compound edges" << endl;

    get_outside_face(fa, outside_face);
    get_outside_face_idx(fa, outside_face_idx);
    itr_matrix<size_t*> new_tet(4, tet_array.size() / 4, &tet_array[0]);
    tet = new_tet;
    itr_matrix<double*> new_node(3, node_array.size() / 3, &node_array[0]);
    node = new_node;
    frame.resize(frame_array.size());
    copy(frame_array.begin(), frame_array.end(), frame.begin());
  }
  cerr << "# finish remove black_lines." << endl;

  {
    relabel_zigzag_by_modify_face_jump_type(tet,node,outside_face,fa,ortae,
                                            singularity_edges,singularity_type,
                                            inner_face_jump_type);

    cerr << "# [info] removing zigzag again." << endl;
  }

#if 1 //test

  jtf::mesh::meshes tm_mesh;
  tm_mesh.mesh_ = tet;
  tm_mesh.node_ = node;
  jtf::tet_mesh tm(tm_mesh);
  singularity_extractor se(tm);
  //  vector<deque<pair<size_t,size_t> > > singularity_edges;
  std::vector<deque<size_t> > singularities_type;
  se.extract(inner_face_jump_type, singularity_edges, singularities_type);
  \
  dump_singularity_chain_to_vtk_2("recalculated_after_modification.vtk",
                                  node,singularity_edges, singularity_type);

  dump_singularity_to_cylinder("recalculated_after_modification.obj",
                               node,singularity_edges,0.002);

  //  vector<size_t > singularity_edges_;
  //  extract_singulairty_edges(ortae,inner_face_jump_type,singularity_edges_);
  //  ofstream ofs("recalculated_after_modification.vtk");
  //  line2vtk(ofs,&node[0],node.size(2),&singularity_edges_[0],singularity_edges_.size()/2);
#endif
  return 0;
}

int relabel_remove_black_lines_by_modify_face_jump_type_and_splitting_fans(
    std::vector<size_t> &tet_array,
    std::vector<double> &node_array,
    std::vector<matrixd > & frame_array,
    jtf::mesh::face2tet_adjacent &fa,
    const matrixst &outside_face,
    jtf::mesh::one_ring_tet_at_edge &ortae,
    std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    std::vector<std::deque<size_t> > &singularity_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type)
{
  boost::unordered_map<size_t,boost::unordered_set<size_t> > vertex_count; // first: vertex_idx, second: the singulairty_chain_idx list
  vector<size_t> black_lines;
  for(size_t t = 0; t < singularity_edges.size(); ++t) {
      const deque<pair<size_t,size_t> > & singularity = singularity_edges[t];
      vertex_count[singularity.front().first].insert(t);
      vertex_count[singularity.back().second].insert(t);

      if(is_black_line_new(singularity_type[t][0])) black_lines.push_back(t);
    }

  cerr << "# [info] black chain num: " << black_lines.size() << endl;
#define visual 1
#if visual
  vector<size_t> black_line_other_path;
  vector<size_t> remove_faces;
#endif

#if visual
  {
    vector<size_t> black_edges;
    for(size_t bi = 0; bi < black_lines.size(); ++bi){
        const deque<pair<size_t,size_t> > & one_chain =
            singularity_edges[black_lines[bi]];
        for(size_t ei = 0; ei < one_chain.size(); ++ei){
            black_edges.push_back(one_chain[ei].first);
            black_edges.push_back(one_chain[ei].second);
          }
      }
    cerr << "# [compound num] compound edge num " << black_edges.size() / 2 << endl;
    ofstream ofs("compound_singularity.vtk");
    line2vtk(ofs, &node_array[0], node_array.size()/3,
        &black_edges[0], black_edges.size()/2);
  }
#endif
  // for each black_chain
  vector<bool> need_to_be_removed_chain(singularity_edges.size(),false);
  vector<deque<pair<size_t,size_t> > > need_to_be_insert_chain;

  for(size_t t = 0; t < black_lines.size(); ++t){
      cerr << "# [info] processing black line " << t << endl;

      itr_matrix<size_t*> tet(4, tet_array.size() / 4, &tet_array[0]);
      itr_matrix<double *> node(3, node_array.size()/3, &node_array[0]);

      deque<pair<size_t,size_t> > & black_chain = singularity_edges[black_lines[t]];
      size_t begin_vertex = black_chain.front().first;
      size_t end_vertex = black_chain.back().second;
      if(is_outside_vertex(begin_vertex,outside_face) &&
         vertex_count[begin_vertex].size() == 1){
          swap(begin_vertex,end_vertex);
          reverse_singularity(black_chain);
        }

      size_t other_chain_idx = -1;

      // to find how many chain linked at this begin_vertex
      // if linked_chain_num > 2, it's normal case, we do it as the paper said
      // if linked_chain_num = 2, it's a near surface black line, both ends reach the surface, but since there is another chain, we can still process this
      // if linked_chain_num = 1, it's a near surface black line, no other chain touch it.
      size_t linked_chain_num = vertex_count[begin_vertex].size();
      if(linked_chain_num == 1) {
          cerr << "# [error] both of ends reach surface, "
               << "and it can be regarded as a near surface lines." << endl;
          continue;
        }
      // to find a chain that linked black line chain
      for(boost::unordered_set<size_t>::const_iterator sci = vertex_count[begin_vertex].begin();
          sci != vertex_count[begin_vertex].end(); ++sci)
        if(*sci != black_lines[t]){ // not the black line itself
            other_chain_idx = *sci ;
            break;
          }

      // to get the one edge away vertex which locates on other chain
      const deque<pair<size_t,size_t> > & other_chain = singularity_edges[other_chain_idx];

      size_t begin_ = -1;
      if(other_chain.front().first == begin_vertex) begin_ = other_chain.front().second;
      else
        if(other_chain.back().second == begin_vertex) begin_ = other_chain.back().first;

      deque<size_t> path;
      // conversion

      deque<size_t> black_chain_;
      black_chain_.push_back(black_chain[0].first);
      for(size_t k = 0; k < black_chain.size(); ++k){
          black_chain_.push_back(black_chain[k].second);
        }
      //find_shortest_path_new(begin_,end_vertex,tet,node,black_chain_,path);
      find_triangle_fan_path(begin_, tet, node, black_chain_, ortae, fa, path);

      assert(path.front() == begin_ &&
             path.back() == black_chain_.back());

      if(split_fans_between_two_path(black_chain_, path, tet_array, node_array,
                                     frame_array, fa, ortae, inner_face_jump_type))
        return __LINE__;
      iterative_modify_face_type(fa,black_chain_,path,ortae,inner_face_jump_type,remove_faces);

      {// adjust the singularity chain
        deque<pair<size_t,size_t> > & other_chain = singularity_edges[other_chain_idx];
        if(other_chain.front().first == begin_vertex){
            other_chain.pop_front();
            for(size_t k = 0; k < path.size() - 1; ++k){
                other_chain.push_front(make_pair(path[k+1],path[k]));
              }
          }else if (other_chain.back().second == begin_vertex) {
            other_chain.pop_back();
            for(size_t k = 0; k < path.size() - 1; ++k){
                other_chain.push_back(make_pair(path[k],path[k+1]));
              }
          }

        // need to detect the other_chain, if parts of this chain is on surface, we should split it.
        // and remove those edges which are on surface
        vector<deque<pair<size_t,size_t> > > chain_after_split;
        bool is_other_chain_changed = false;
        for(size_t k = 0; k < other_chain.size(); ++k){
            const pair<size_t,size_t> & each_edge = other_chain[k];
            if(!is_outside_edge(each_edge.first,each_edge.second,outside_face)) {
                if(chain_after_split.empty()) {
                    deque<pair<size_t,size_t> > empty_deque;
                    chain_after_split.push_back(empty_deque);
                  }
                chain_after_split.back().push_back(each_edge);
              }else{
                is_other_chain_changed = true;
                if(!chain_after_split.empty() && !chain_after_split.back().empty()){
                    deque<pair<size_t,size_t> > empty_deque;
                    chain_after_split.push_back(empty_deque);
                  }
              }
          }
        if(is_other_chain_changed) {
            need_to_be_removed_chain[other_chain_idx] = true;
            if(!chain_after_split.empty()) {
                if(chain_after_split.size() == 1 && !chain_after_split[0].empty()) //to ensure there is singularity edge
                  need_to_be_insert_chain.insert(
                        need_to_be_insert_chain.end(),chain_after_split.begin(),
                        chain_after_split.end());
              }
          }

        ////////////////////////////////////////////////////////////////////////////////////////////
        if(linked_chain_num > 2){ // need to merge the black line to 3rd chain
            size_t third_chain_idx = -1;
            for(boost::unordered_set<size_t>::const_iterator sci = vertex_count[begin_vertex].begin();
                sci != vertex_count[begin_vertex].end(); ++sci){
                if(*sci != other_chain_idx && * sci != black_lines[t]) {
                    third_chain_idx = *sci;
                    break;
                  }
              }
            if(third_chain_idx == -1){
                cerr << "# strange, can not find 3rd chain, but this vertex has at least three chain linked." << endl;
                return __LINE__;
              }

            deque<pair<size_t,size_t> > & third_chain = singularity_edges[third_chain_idx];
            if(third_chain.front().first == begin_vertex){
                for(size_t k = 0; k < black_chain.size(); ++k){
                    third_chain.push_front(make_pair(black_chain[k].second,black_chain[k].first));
                  }
              }else if (third_chain.back().second == begin_vertex) {
                for(size_t k = 0; k < black_chain.size(); ++k){
                    third_chain.push_back(black_chain[k]);
                  }
              }

            // need to detect the third_chain, if parts of this chain is on surface, we should split it.
            chain_after_split.clear();
            bool is_third_chain_changed = false;
            for(size_t k = 0; k < third_chain.size(); ++k){
                const pair<size_t,size_t> & each_edge = third_chain[k];
                if(!is_outside_edge(each_edge.first,each_edge.second,outside_face)) {
                    if(chain_after_split.empty()) {
                        deque<pair<size_t,size_t> > empty_deque;
                        chain_after_split.push_back(empty_deque);
                      }
                    chain_after_split.back().push_back(each_edge);
                  }else{
                    is_third_chain_changed = true;
                    if(!chain_after_split.empty() && !chain_after_split.back().empty())
                      {
                        deque<pair<size_t,size_t> > empty_deque_;
                        chain_after_split.push_back(empty_deque_);
                      }
                  }
              }
            if(is_third_chain_changed){
                need_to_be_removed_chain[third_chain_idx] = true;
                if(!chain_after_split.empty())
                  {
                    if(chain_after_split.size() == 1 && !chain_after_split[0].empty())
                      need_to_be_insert_chain.insert(need_to_be_insert_chain.end(),chain_after_split.begin(),chain_after_split.end());
                  }
              }
            /////////////////////////////////////////

            need_to_be_removed_chain[black_lines[t]] = true;
          }
      }
    }// end for

#if 0
  {
    itr_matrix<size_t*> tet(4, tet_array.size() / 4, &tet_array[0]);
    itr_matrix<double *> node(3, node_array.size()/3, &node_array[0]);
    ofstream ofs("path.vtk");
    ofstream ofs_face("faces.vtk");
    line2vtk(ofs,&node[0],node.size(2),&black_line_other_path[0],black_line_other_path.size()/2);
    tri2vtk(ofs_face,&node[0],node.size(2),&remove_faces[0],remove_faces.size()/3);
  }
#endif

  // assemble the singularities edges
  vector<deque<pair<size_t,size_t> > > new_singularity_edges;
  new_singularity_edges.reserve(singularity_edges.size());
  for(size_t t = 0; t < need_to_be_removed_chain.size(); ++t){
      if(!need_to_be_removed_chain[t]){
          new_singularity_edges.push_back(singularity_edges[t]);
        }
    }
  new_singularity_edges.insert(new_singularity_edges.end(),
                               need_to_be_insert_chain.begin(),
                               need_to_be_insert_chain.end());
  swap(new_singularity_edges,singularity_edges);

#if 1
  for(size_t i = 0; i < singularity_edges.size(); ++i){
      if(singularity_edges[i].size() == 0)
        cerr << "# [error] this singularity is empty." << endl;
    }
#endif
  return 0;
}

int relabel_remove_black_lines_by_modify_face_jump_type(
    const matrixst &tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent &fa,
    const matrixst &outside_face,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    vector<deque<pair<size_t,size_t> > > &singularity_edges,
    vector<deque<size_t> > &singularity_type,
    boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type)
{
  map<size_t,set<size_t> > vertex_count; // first: vertex_idx, second: the singulairty_chain_idx list
  vector<size_t> black_lines;
  for(size_t t = 0; t < singularity_edges.size(); ++t)
    {
      const deque<pair<size_t,size_t> > & singularity = singularity_edges[t];
      vertex_count[singularity.front().first].insert(t);
      vertex_count[singularity.back().second].insert(t);

      if(is_black_line_new(singularity_type[t][0])) black_lines.push_back(t);
    }
  cerr << "# [info] black chain num: " << black_lines.size() << endl;
#define visual 1
#if visual
  vector<size_t> black_line_other_path;
  vector<size_t> remove_faces;
#endif
  // for each black_chain
  vector<bool> need_to_be_removed_chain(singularity_edges.size(),false);
  vector<deque<pair<size_t,size_t> > > need_to_be_insert_chain;

  for(size_t t = 0; t < black_lines.size(); ++t){
      cerr << "# [info] processing black line " << t << endl;
      deque<pair<size_t,size_t> > & black_chain = singularity_edges[black_lines[t]];
      size_t begin_vertex = black_chain.front().first;
      size_t end_vertex = black_chain.back().second;
      if(is_outside_vertex(begin_vertex,outside_face) &&
         vertex_count[begin_vertex].size() == 1){
          swap(begin_vertex,end_vertex);
          deque<pair<size_t,size_t> > temp_chain;
          for(size_t k = black_chain.size() - 1; k != -1; --k){
              temp_chain.push_back(make_pair(black_chain[k].second,black_chain[k].first));
            }
          swap(temp_chain,black_chain);
        }

      size_t other_chain_idx = -1;

      // to find how many chain linked at this begin_vertex
      // if linked_chain_num > 2, it's normal case, we do it as the paper said
      // if linked_chain_num = 2, it's a near surface black line, both ends reach the surface, but since there is another chain, we can still process this
      // if linked_chain_num = 1, it's a near surface black line, no other chain touch it.
      size_t linked_chain_num = vertex_count[begin_vertex].size();
      if(linked_chain_num == 1) {
          cerr << "# both of ends reach surface, and it can be regarded as a near surface lines." << endl;
          continue;
        }
      // to find a chain that linked black line chain
      for(set<size_t>::const_iterator sci = vertex_count[begin_vertex].begin();
          sci != vertex_count[begin_vertex].end(); ++sci)
        if(*sci != black_lines[t]){ // not the black line itself
            other_chain_idx = *sci ;
            break;
          }

      // to get the one edge away vertex which locates on other chain
      const deque<pair<size_t,size_t> > & other_chain = singularity_edges[other_chain_idx];

#if 0 // test to show the result
      static size_t count = 0;
      if(count == 0 )
        {
          vector<deque<pair<size_t,size_t> > > test_;
          test_.push_back(other_chain);
          dump_singularity_to_vtk("other_chain_0.vtk",node,test_);
        }

      if(count == 1 )
        {
          vector<deque<pair<size_t,size_t> > > test_;
          test_.push_back(other_chain);
          dump_singularity_to_vtk("other_chain_1.vtk",node,test_);
        }
      ++count;
#endif
      size_t begin_ = -1;
      if(other_chain.front().first == begin_vertex)
        begin_ = other_chain.front().second;
      else
        if(other_chain.back().second == begin_vertex)
          begin_ = other_chain.back().first;

      deque<size_t> path;
      // conversion

      deque<size_t> black_chain_;
      black_chain_.push_back(black_chain[0].first);
      for(size_t k = 0; k < black_chain.size(); ++k){
          black_chain_.push_back(black_chain[k].second);
        }
      find_shortest_path_new(begin_,end_vertex,tet,node,black_chain_,path);

#undef visual
#if visual
      for(size_t i = 0; i < path.size() - 1; ++i)
        {
          black_line_other_path.push_back(path[i]);
          black_line_other_path.push_back(path[i+1]);
        }
#endif
      iterative_modify_face_type(fa,black_chain_,path,ortae,inner_face_jump_type,remove_faces);

      {// adjust the singularity chain
        deque<pair<size_t,size_t> > & other_chain = singularity_edges[other_chain_idx];
        if(other_chain.front().first == begin_vertex){
            other_chain.pop_front();
            for(size_t k = 0; k < path.size() - 1; ++k){
                other_chain.push_front(make_pair(path[k+1],path[k]));
              }
          }else if (other_chain.back().second == begin_vertex) {
            other_chain.pop_back();
            for(size_t k = 0; k < path.size() - 1; ++k){
                other_chain.push_back(make_pair(path[k],path[k+1]));
              }
          }

        // need to detect the other_chain, if parts of this chain is on surface, we should split it.
        // and remove those edges which are on surface
        vector<deque<pair<size_t,size_t> > > chain_after_split;
        bool is_other_chain_changed = false;
        for(size_t k = 0; k < other_chain.size(); ++k){
            const pair<size_t,size_t> & each_edge = other_chain[k];
            if(!is_outside_edge(each_edge.first,each_edge.second,outside_face)) {
                if(chain_after_split.empty()) {
                    deque<pair<size_t,size_t> > empty_deque;
                    chain_after_split.push_back(empty_deque);
                  }
                chain_after_split.back().push_back(each_edge);
              }else{
                is_other_chain_changed = true;
                if(!chain_after_split.empty() && !chain_after_split.back().empty()){
                    deque<pair<size_t,size_t> > empty_deque;
                    chain_after_split.push_back(empty_deque);
                  }
              }
          }
        if(is_other_chain_changed)
          {
            need_to_be_removed_chain[other_chain_idx] = true;
            if(!chain_after_split.empty())
              {
                if(chain_after_split.size() == 1 && !chain_after_split[0].empty()) //to ensure there is singularity edge
                  need_to_be_insert_chain.insert(need_to_be_insert_chain.end(),chain_after_split.begin(),chain_after_split.end());
              }
          }

        ////////////////////////////////////////////////////////////////////////////////////////////
        if(linked_chain_num > 2){ // need to merge the black line to 3rd chain
            size_t third_chain_idx = -1;
            for(set<size_t>::const_iterator sci = vertex_count[begin_vertex].begin();
                sci != vertex_count[begin_vertex].end(); ++sci){
                if(*sci != other_chain_idx && * sci != black_lines[t])
                  {
                    third_chain_idx = *sci;
                    break;
                  }
              }
            if(third_chain_idx == -1){
                cerr << "# strange, can not find 3rd chain, but this vertex has at least three chain linked." << endl;
                return __LINE__;
              }

            deque<pair<size_t,size_t> > & third_chain = singularity_edges[third_chain_idx];
            if(third_chain.front().first == begin_vertex){
                for(size_t k = 0; k < black_chain.size(); ++k){
                    third_chain.push_front(make_pair(black_chain[k].second,black_chain[k].first));
                  }
              }else if (third_chain.back().second == begin_vertex) {
                for(size_t k = 0; k < black_chain.size(); ++k){
                    third_chain.push_back(black_chain[k]);
                  }
              }

            // need to detect the third_chain, if parts of this chain is on surface, we should split it.
            chain_after_split.clear();
            bool is_third_chain_changed = false;
            for(size_t k = 0; k < third_chain.size(); ++k){
                const pair<size_t,size_t> & each_edge = third_chain[k];
                if(!is_outside_edge(each_edge.first,each_edge.second,outside_face)) {
                    if(chain_after_split.empty()) {
                        deque<pair<size_t,size_t> > empty_deque;
                        chain_after_split.push_back(empty_deque);
                      }
                    chain_after_split.back().push_back(each_edge);
                  }else{
                    is_third_chain_changed = true;
                    if(!chain_after_split.empty() && !chain_after_split.back().empty())
                      {
                        deque<pair<size_t,size_t> > empty_deque_;
                        chain_after_split.push_back(empty_deque_);
                      }
                  }
              }
            if(is_third_chain_changed){
                need_to_be_removed_chain[third_chain_idx] = true;
                if(!chain_after_split.empty())
                  {
                    if(chain_after_split.size() == 1 && !chain_after_split[0].empty())
                      need_to_be_insert_chain.insert(need_to_be_insert_chain.end(),chain_after_split.begin(),chain_after_split.end());
                  }
              }
            /////////////////////////////////////////

            need_to_be_removed_chain[black_lines[t]] = true;
          }
      }
    }// end for

#if visual
  ofstream ofs("path.vtk");
  ofstream ofs_face("faces.vtk");
  line2vtk(ofs,&node[0],node.size(2),&black_line_other_path[0],black_line_other_path.size()/2);
  tri2vtk(ofs_face,&node[0],node.size(2),&remove_faces[0],remove_faces.size()/3);
#endif

  // assemble the singularities edges
  vector<deque<pair<size_t,size_t> > > new_singularity_edges;
  new_singularity_edges.reserve(singularity_edges.size());
  for(size_t t = 0; t < need_to_be_removed_chain.size(); ++t){
      if(!need_to_be_removed_chain[t]){
          new_singularity_edges.push_back(singularity_edges[t]);
        }
    }
  new_singularity_edges.insert(new_singularity_edges.end(),need_to_be_insert_chain.begin(),need_to_be_insert_chain.end());
  swap(new_singularity_edges,singularity_edges);

#if 1
  for(size_t i = 0; i < singularity_edges.size(); ++i)
    {
      if(singularity_edges[i].size() == 0)
        cerr << "# [error] this singularity is empty." << endl;
    }
#endif

  return 0;
}

int relabel_zigzag_by_modify_face_jump_type(
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    vector<deque<pair<size_t,size_t> > > &singularity_edges,
    vector<deque<size_t> > &singularity_type,
    boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type)
{
  vector<bool> need_to_remove_chain(singularity_edges.size(),false);
  vector<deque<pair<size_t,size_t> > > need_to_insert_chain;

#define VISUAL
#ifdef VISUAL
  vector<size_t> zigzag_edges;
  vector<size_t> regular_edges;
#endif

  timer tr;
  tr.start();
  for(size_t t = 0; t < singularity_edges.size(); ++t){ // for each singularity chain
      deque<pair<size_t,size_t> > &single_chain = singularity_edges[t];
      deque<pair<size_t,size_t> > new_chain;
      new_chain.push_back(single_chain[0]);

      for(size_t i = 0; i < single_chain.size() - 1; ++i){
          const size_t face_idx = is_zigzag(new_chain.back(),single_chain[i+1],fa);
          if(face_idx != -1) { // is zigzag
#ifdef VISUAL
              {
                zigzag_edges.push_back(new_chain.back().first);
                zigzag_edges.push_back(new_chain.back().second);
                zigzag_edges.push_back(single_chain[i+1].first);
                zigzag_edges.push_back(single_chain[i+1].second);
              }
#endif
              remove_zigzag(new_chain.back(),single_chain[i+1],face_idx,ortae,fa,inner_face_jump_type,tet,node);
              pair<size_t,size_t> back_ = new_chain.back();
              new_chain.pop_back();
              new_chain.push_back(make_pair(back_.first,single_chain[i+1].second));
            }
          else
            new_chain.push_back(single_chain[i+1]);
        }
      if(new_chain.size() != single_chain.size()){ // detect zigzag, need to replace the chain

#if 1   // to handle the edges exist on surface
          vector<deque<pair<size_t,size_t> > > chain_after_split;
          bool is_chain_changed = false;
          for(size_t k = 0; k < new_chain.size(); ++k){
              const pair<size_t,size_t> & each_edge = new_chain[k];
              if(!is_outside_edge(each_edge.first,each_edge.second,outside_face)) {
                  if(chain_after_split.empty()) {
                      deque<pair<size_t,size_t> > empty_deque;
                      chain_after_split.push_back(empty_deque);
                    }
                  chain_after_split.back().push_back(each_edge);
                }else{
                  is_chain_changed = true;
                  if(!chain_after_split.empty() && !chain_after_split.back().empty()){
                      deque<pair<size_t,size_t> > empty_deque;
                      chain_after_split.push_back(empty_deque);
                    }
                }
            }
          if(is_chain_changed) {
              need_to_remove_chain[t] = true;
              if(!chain_after_split.empty())
                need_to_insert_chain.insert(need_to_insert_chain.end(),chain_after_split.begin(),chain_after_split.end());
            }
#endif
        }
#ifdef VISUAL
      {
        if(is_regular_type((singularity_type[t].front()))){
            for(size_t ei = 0; ei < single_chain.size(); ++ei){
                if(find(zigzag_edges.begin(), zigzag_edges.end(), single_chain[ei].first)
                   == zigzag_edges.end() ||
                   find(zigzag_edges.begin(), zigzag_edges.end(), single_chain[ei].second)
                   == zigzag_edges.end()){
                    regular_edges.push_back(single_chain[ei].first);
                    regular_edges.push_back(single_chain[ei].second);
                  }
              }
          }
      }
#endif
    }
  vector<deque<pair<size_t,size_t> > > new_singularities;
  new_singularities.reserve(singularity_edges.size());
  for(size_t t = 0;t < singularity_edges.size(); ++t){
      if(!need_to_remove_chain[t]) new_singularities.push_back(singularity_edges[t]);
    }

  new_singularities.insert(new_singularities.end(),need_to_insert_chain.begin(),
                           need_to_insert_chain.end());
  swap(new_singularities,singularity_edges);

  singularity_type.clear();
  singularity_type.resize(singularity_edges.size());
  for(size_t ci = 0; ci < singularity_edges.size(); ++ci){
      singularity_type[ci].resize(singularity_edges[ci].size());
      const deque<pair<size_t,size_t> > & one_chain = singularity_edges[ci];
      for(size_t ei = 0; ei < one_chain.size(); ++ei){
          singularity_type[ci][ei] =
              get_edge_type_with_part_face_type_map(ortae,one_chain[ei],
                                                    inner_face_jump_type);
        }
    }

  tr.finish();
  cerr << "# [zigzag time] cost " << tr.result()
       << "ms to remove inner zigzag" << endl;

#ifdef VISUAL
  {
    ofstream ofs("zigzag_edges.vtk");
    line2vtk(ofs, &node[0], node.size(2), &zigzag_edges[0], zigzag_edges.size()/2);
    cerr << "# [zigzag num] inner zigzag edges num " << zigzag_edges.size()/2 << endl;
    ofstream ofs2("regular_edges.vtk");
    line2vtk(ofs2, &node[0], node.size(2), &regular_edges[0], regular_edges.size()/2);
  }
#endif
  return 0;
}

int relabel_remove_near_surface_by_modify_face_jump_type(
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    vector<deque<pair<size_t,size_t > > > &singularity_edges,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type)
{
  ///////////////////////////////////////////////
  /// construct a graph with edge weighting  ////
  map<pair<size_t,size_t>,double> edge_weight_;
  for(size_t t = 0 ; t < outside_face.size(2); ++t){
      for(size_t j = 0; j < 3; ++j){
          pair<size_t,size_t> edge(outside_face(j,t),outside_face((j+1)%3,t));
          if(edge.first > edge.second) swap(edge.first,edge.second);
          edge_weight_[edge] = norm(node(colon(),edge.first) - node(colon(),edge.second)); // or use other weighting like distanct
        }
    }
  unique_ptr<vertex_connection<UNDIRECT> >  vc(vertex_connection<UNDIRECT>::create(edge_weight_));

  vector<bool> need_to_be_removed(singularity_edges.size(),false);
  for(size_t t = 0; t < singularity_edges.size(); ++t) {
      if(is_one_ring_near_surface_chain(singularity_edges[t],tet,outside_face)){

          remove_near_surface_chain(tet,node,outside_face,fa,ortae,singularity_edges[t],inner_face_jump_type,*vc);
          need_to_be_removed[t] = true;
        }
    }

  vector<deque<pair<size_t,size_t> > > temp_edges;
  temp_edges.reserve(singularity_edges.size());
  for(size_t t = 0; t < need_to_be_removed.size(); ++t){
      if(!need_to_be_removed[t]){
          temp_edges.push_back(singularity_edges[t]);
        }
    }
  swap(temp_edges,singularity_edges);
  return 0;
}

int split_fans_between_two_path(
    const std::deque<size_t> &black_chain,
    std::deque<size_t> &path,
    std::vector<size_t> &tet_array,
    std::vector<double> &node_array,
    std::vector<matrixd > & frame,
    jtf::mesh::face2tet_adjacent &fa,
    jtf::mesh::one_ring_tet_at_edge &ortae,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type)
{
  if(path.size() < 3) return 0;
  assert(path.back() == black_chain.back());
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  typedef deque<size_t>::const_iterator dci;
  dci black_begin = black_chain.begin();
  dci black_end = black_chain.end() - 1;
  dci other_begin = path.begin();

  deque<size_t> new_path;
  new_path.push_back(path.front());

  boost::unordered_map<pair<size_t,size_t>,vector<pair<size_t,size_t> > > split_edge_map;
  boost::unordered_map<size_t,size_t > old_tet_map_to_new;

  vector<size_t> next_point(2);
  vector<size_t> face_idx(2);

  while((*other_begin != *black_end)
        && (*black_begin != *black_end))
    {
      next_point[0] = *(black_begin + 1);
      next_point[1] = *(other_begin + 1);
      size_t t = 0;
      for(; t < 2; ++t){
          face_idx[t] =
              fa.get_face_idx(next_point[t],*black_begin,*other_begin);
          if(face_idx[t] == -1) continue;
          if(next_point[t] == *black_end) break; // reach the end
          pair<size_t,size_t> edge(next_point[t],
                                   (t == 0? *other_begin:*black_begin));
          split_tets_around_edge(
                tet_array, node_array, ortae, frame, edge, 0.5,
                split_edge_map,old_tet_map_to_new,inner_face_jump_type);
          new_path.push_back(node_array.size()/3-1); // add new_point
          break;
        }
      if(t == 2){
          cerr << "# [error] can not find a suitable face." << endl;
          return __LINE__;
        }else if(t == 0){
          ++black_begin;
        }else if(t == 1){
          ++other_begin;
        }
    }
  new_path.push_back(path.back());
  //assert(new_path.size() == path.size());
  swap(new_path, path);

  itr_matrix<size_t*> new_tet(4, tet_array.size()/4, &tet_array[0]);
  unique_ptr<jtf::mesh::face2tet_adjacent> fa_new(jtf::mesh::face2tet_adjacent::create(new_tet));
  fa = *fa_new;
  return 0;
}

int iterative_modify_face_type(
    const jtf::mesh::face2tet_adjacent &fa,
    const deque<size_t> & black_chain,
    const deque<size_t> &other_chain,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type,
    vector<size_t> &remove_faces)
{
  // this function assume that black_chain and other_chain has the same ending and different beginning
  // the two beginning's distance is one edge
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  typedef deque<size_t>::const_iterator dci;
  dci black_begin = black_chain.begin();
  dci black_end = black_chain.end() - 1;
  dci other_begin = other_chain.begin();
  dci next = other_begin;
  const size_t prev_removed_points = remove_faces.size();
  while((*other_begin != *black_end)
        && (*black_begin != *black_end))
    {
      next = other_begin;
      ++next;
      size_t face_idx = fa.get_face_idx(*next,*black_begin,*other_begin);
      if(face_idx != -1){ // it's a common face, need to remove the singularity edge <*black_begin,*other_begin>

          remove_faces.push_back(*next);
          remove_faces.push_back(*black_begin);
          remove_faces.push_back(*other_begin);

          const pair<size_t,size_t> & two_tets = fa.face2tet_[face_idx];
          oecit it = ortae.e2t_.find(make_pair(*black_begin,*other_begin));
          if(it == ortae.e2t_.end()) it = ortae.e2t_.find(make_pair(*other_begin,*black_begin));
          const vector<size_t> &around_tets = it->second;
          modify_face_jump_type_at_given_tets_edge(two_tets.first,two_tets.second,inner_face_jump_type,around_tets);
          //modify_face_jump_type_to_remove_singularity_edge(two_tets.first,two_tets.second,inner_face_jump_type,around_tets);
          other_begin = next;



#if 0 // check the around tets type
          matrixd rot = eye<double>(3);
          for(size_t t = 0; t < around_tets.size() - 1; ++t){
              rot = temp(rot * type_transition2(inner_face_jump_type[make_pair(around_tets[t],around_tets[t+1])]) );
            }
          cerr << "# rot = " << rot << endl;
#endif
        }
      else{
          next = black_begin;
          ++next;
          size_t face_idx = fa.get_face_idx(*next,*black_begin,*other_begin);
          if(face_idx != -1){ // it's a common face, need to remove the singularity edge <*black_begin,*other_begin>

              remove_faces.push_back(*next);
              remove_faces.push_back(*black_begin);
              remove_faces.push_back(*other_begin);

              const pair<size_t,size_t> & two_tets = fa.face2tet_[face_idx];
              oecit it = ortae.e2t_.find(make_pair(*black_begin,*other_begin));
              if(it == ortae.e2t_.end()) it = ortae.e2t_.find(make_pair(*other_begin,*black_begin));
              const vector<size_t> &around_tets = it->second;
              modify_face_jump_type_at_given_tets_edge(two_tets.first,two_tets.second,inner_face_jump_type,around_tets);
              black_begin = next;

#if 0 // check the around tets type
              matrixd rot = eye<double>(3);
              for(size_t t = 0; t < around_tets.size() - 1; ++t){
                  rot = temp(rot * type_transition2(inner_face_jump_type[make_pair(around_tets[t],around_tets[t+1])]) );
                }
              cerr << "# rot = " << rot << endl;
#endif
            }else{
              black_begin = next;
            }
        }
    }
  if(remove_faces.size() == prev_removed_points) {
      cerr << "# [error] has not remove this black edges" << endl;
      return __LINE__;
    }
  return 0;
}

//int remove_zigzag_new(const std::pair<size_t,size_t> & edge0,
//                      const std::pair<size_t,size_t> & edge1,
//                      const size_t face_idx,
//                      const one_ring_tet_at_edge &ortae,
//                      const jtf::mesh::face2tet_adjacent & fa,
//                      std::map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
//                      const matrixst &tet)
//{
//  const pair<size_t,size_t> & tet_pair = fa.face2tet_[face_idx];
//  if(fa.is_outside_face(tet_pair)){
//    cerr << "# [error] wrong tet pair" << endl;
//    return __LINE__;
//  }


//  return 0;
//}

// return -1 if is not zigzag, else return the common face idx
int remove_zigzag(
    const pair<size_t,size_t> &edge0,
    const pair<size_t,size_t> &edge1,
    const size_t face_idx,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::face2tet_adjacent &fa,
    boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const matrixst &tet,
    const matrixd &node)
{
  // this function assume the edge0 and edge1 is zigzag
  if((edge0.first == edge1.first) || (edge0.second == edge1.second))
    {
      cerr << "# this two edge has opposite order." << endl;
      return __LINE__;
    }

  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  oecit edge0_it = ortae.e2t_.find(edge0);
  oecit edge1_it = ortae.e2t_.find(edge1);
  vector<size_t> edge0_tet,edge1_tet;
  if(edge0_it == ortae.e2t_.end()) // this edge should be up_and_down
    {
      edge0_it = ortae.e2t_.find(make_pair(edge0.second,edge0.first));
      if( edge0_it == ortae.e2t_.end()) {
          cerr << "# strange can not find this edge." << endl;
          return __LINE__;
        }
      edge0_tet = edge0_it->second;
      reverse(edge0_tet.begin(),edge0_tet.end()); // reverse this around tets
    }else
    edge0_tet = edge0_it->second;

  if(edge1_it == ortae.e2t_.end()) // this edge should be up_and_down
    {
      edge1_it = ortae.e2t_.find(make_pair(edge1.second,edge1.first));
      if( edge1_it == ortae.e2t_.end()) {
          cerr << "# strange can not find this edge." << endl;
          return __LINE__;
        }
      edge1_tet = edge1_it->second;
      reverse(edge1_tet.begin(),edge1_tet.end()); // reverse this around tets
    }else
    edge1_tet = edge1_it->second;

  const pair<size_t,size_t> &two_tets = fa.face2tet_[face_idx];

#if 1 // check whether the both tets exist in each edges' arounding tets
  if((find(edge0_tet.begin(),edge0_tet.end(),two_tets.first) == edge0_tet.end())
     || (find(edge0_tet.begin(),edge0_tet.end(),two_tets.second) == edge0_tet.end())
     || (find(edge1_tet.begin(),edge1_tet.end(),two_tets.first) == edge1_tet.end())
     || (find(edge1_tet.begin(),edge1_tet.end(),two_tets.second) == edge1_tet.end()))
    {
      cerr << "# strange can not find both tets in these edges' arounding." << endl;
    }
#endif

  modify_face_jump_type_at_given_tets_edge(two_tets.first,two_tets.second,
                                           inner_face_jump_type,edge0_tet);

#if 0 // check the other edge, which may not be singulairty edge anymore
  matrixd rot = eye<double>(3);
  typedef map<pair<size_t,size_t>,size_t>::const_iterator mcit;
  vector<size_t> error_edges;

  for(size_t t = 0; t < edge1_tet.size() - 1; ++t){
      if((edge1_tet[t] == two_tets.first && (edge1_tet[t+1] == two_tets.second))
         || (edge1_tet[t] == two_tets.second && (edge1_tet[t+1] == two_tets.first))){ // match the two tets
          vector<size_t> reorder_tets;
          reorder_tets.reserve(edge1_tet.size());
          for(size_t i = t; i < edge1_tet.size() - 1; ++i) reorder_tets.push_back(edge1_tet[i]);
          for(size_t i = 0; i < t + 1; ++i) reorder_tets.push_back(edge1_tet[i]);
          for(size_t j = 0; j < reorder_tets.size() - 1; ++j){
              mcit it = inner_face_jump_type.find(make_pair(reorder_tets[j],reorder_tets[j+1]));
              if(it == inner_face_jump_type.end()) continue;
              rot = temp(rot * type_transition2(it->second));
            }
          if(fabs(norm(rot - eye<double>(3))) > 1e-8) {
              cerr << "# strange the zigzag modification is nor correct on the other edge." << endl;
              cerr << "# rot " << rot << endl;
              //return __LINE__;
              ofstream ofs("error_singularity_edges.vtk");
              error_edges.push_back(edge1.first);
              error_edges.push_back(edge1.second);
              line2vtk(ofs,&node[0],node.size(2),&error_edges[0],error_edges.size()/2);
            }
          break;
        }
    }
#endif

#if 0 // dump out singularity edges from modified face jump type
  vector<pair<size_t,size_t> > singularity_edges;
  extract_singulairty_edges(ortae,inner_face_jump_type,singularity_edges);
  ofstream ofs("recal_singularity_edges.vtk");
  vector<size_t> singularity_edges_;
  for(size_t t = 0; t < singularity_edges.size(); ++t){
      singularity_edges_.push_back(singularity_edges[t].first);
      singularity_edges_.push_back(singularity_edges[t].second);
    }
  line2vtk(ofs,&node[0],node.size(2),&singularity_edges_[0],singularity_edges_.size()/2);
#endif
  return 0;
}

//int relabel_face_type_to_remove_zigzag(
//    std::vector<std::deque<std::pair<size_t,size_t> > > &chain_list,
//    std::vector<std::deque<size_t> >  &singularities_type_,
//    const matrixst &tet,
//    const jtf::mesh::face2tet_adjacent & fa,
//    const matrixst &outside_face,
//    std::map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
//    matrixst &tet_rot_type,
//    const one_ring_tet_at_edge &ortae)
//{
//  map<pair<size_t,size_t>,size_t> singularity_edges;
//  for(size_t t = 0; t < chain_list.size(); ++t){
//    const deque<pair<size_t,size_t> > & chain = chain_list[t];
//    for(size_t j = 0; j < chain.size(); ++j){
//      const pair<size_t,size_t> & edge = chain[j];
//      if(edge.first > edge.second){
//        singularity_edges[make_pair(edge.second,edge.first)] = singularities_type_[t][j];
//      }else
//        singularity_edges[edge] = singularities_type_[t][j];
//    }
//  }
//  // detech zigzag singularity edges
//  typedef map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
//  for(size_t t = 0; t < fa.faces_.size(); ++t){
//    const vector<size_t> & face = fa.faces_[t];
//    vector<size_t> edge_idx;
//    vector<size_t> edge_type;
//    for(size_t i = 0; i < face.size(); ++i){
//      pair<size_t,size_t> edge(face[i],face[(i+1)%face.size()]);
//      if(edge.first > edge.second) swap(edge.second,edge.first);
//      mpscit mcit = singularity_edges.find(edge);
//      if(mcit != singularity_edges.end()){
//        edge_idx.push_back(i);
//        edge_type.push_back(mcit->second);
//      }
//    }
//    if(edge_idx.size() == 0 || // this face do not have singularity edges
//       edge_idx.size() == 1 || // this face has only one singularity edge
//       edge_idx.size() == 3) // this face has three singularity edges, which needs to be  removed
//      continue;

//    assert(edge_idx.size() == 2);
//    assert(edge_type.size() == 2);

//    if(edge_type.front() != edge_type.back()){
//      cerr << "# [error] two singularity edges is zigzag with different type: " << endl;
//      cerr << "#N ------ edge <" << face[edge_idx.front()] << ","
//           << face[(edge_idx.front()+1)%face.size()] << " type "
//           << edge_idx.front() << endl;
//      cerr << "#N ------ edge <" << face[edge_idx.back()] << ","
//           << face[(edge_idx.back()+1)%face.size()] << " type "
//           << edge_idx.back() << endl;
//      continue;

//    }
//    remove_zigzag_new(make_pair(face[edge_idx.front()],face[(edge_idx.front()+1)%face.size()]),
//                      make_pair(face[edge_idx.back()],face[(edge_idx.back()+1)%face.size()]),
//                      t,ortae,fa,inner_face_jump_type,tet);
//  }
//  return 0;
//}

int check_surface_vertex_to_make_face(const vector<size_t> &path,
                                      const pair<size_t,size_t> &one_edge,
                                      const jtf::mesh::face2tet_adjacent &fa,
                                      const size_t order = 0) // 0: from begin to end in path; 1: from end to begin
{
  if(order == 0){
      for(size_t t = 0; t < path.size(); ++t){
          size_t face_idx = fa.get_face_idx(one_edge.first,one_edge.second,path[t]);
          if(face_idx != -1) return path[t];
        }
    }else {
      for(size_t t = path.size() - 1; t != -1; --t){
          size_t face_idx = fa.get_face_idx(one_edge.first,one_edge.second,path[t]);
          if(face_idx != -1) return path[t];
        }
    }
  return -1;
}

int remove_near_surface_chain(
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    const deque<pair<size_t,size_t> > &chain,
    boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const vertex_connection<UNDIRECT> & vc)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;

  vector<size_t> path;
  if(vc.get_shortest_path(chain.front().first,chain.back().second,path))
    {
      cerr << "# [error] can not find a shortest path." << endl;
      return __LINE__;
    }

  ///////////////////////////////////////////////////////////////////
  matrixd center_of_surface_path = zeros<double>(3,1);
  for(size_t t = 0; t < path.size(); ++t)
    center_of_surface_path += node(colon(),path[t]);
  center_of_surface_path /= path.size();
  vector<size_t> one_ring_vertex;

  vector<pair<double,size_t> > distance_each_vertex;
  deque<pair<size_t,size_t> > modified_chain_;
  modified_chain_.resize(chain.size());
  copy(chain.begin(),chain.end(),modified_chain_.begin());

  vector<pair<size_t,size_t> > edges_on_surface;
  for(size_t t = 0; t < path.size() -1; ++t) {
      edges_on_surface.push_back(make_pair(path[t],path[t+1]));
      edges_on_surface.push_back(make_pair(path[t+1],path[t]));
    }


  //    while(!modified_chain_.empty())
  //    {
  //        cerr << "--- need_to_modified_chain_length: " << modified_chain_.size() << endl;
  //        const pair<size_t,size_t> & one_edge = modified_chain_.front();
  //        if(find(edges_on_surface.begin(),edges_on_surface.end(),one_edge) != edges_on_surface.end())
  //        {// this edge is already exist on surface, do not need to handle it
  //            modified_chain_.pop_front();
  //            continue;
  //        }

  //        distance_each_vertex.clear();
  //        find_one_ring_vertex_around_edge(tet,one_edge,ortae, one_ring_vertex);
  //        for(size_t t = 0; t < one_ring_vertex.size(); ++t){
  //            distance_each_vertex.push_back(make_pair(norm(node(colon(),one_ring_vertex[t]) - center_of_surface_path),t));
  //        }
  //        sort(distance_each_vertex.begin(),distance_each_vertex.end());
  //        size_t face_idx = fa.get_face_idx(one_edge.first,one_edge.second,one_ring_vertex[distance_each_vertex[0].second]);
  //        const pair<size_t,size_t> &tets = fa.face2tet_[face_idx];

  //        oecit it = ortae.e2t_.find(one_edge);
  //        if(it == ortae.e2t_.end()) it = ortae.e2t_.end();
  //        if(it == ortae.e2t_.end()) return __LINE__;

  //        modify_face_jump_type_at_given_tets_edge(tets.first,tets.second,inner_face_jump_type, it->second);
  //        modified_chain_.pop_front();
  //        modified_chain_.push_front(make_pair(one_edge.first,one_ring_vertex[distance_each_vertex[0].second]));
  //        modified_chain_.push_front(make_pair(one_edge.second,one_ring_vertex[distance_each_vertex[0].second]));
  //    }

#if 1 // check whether the chain edge is inner
  for(size_t t = 0; t < chain.size(); ++t){
      const pair<size_t,size_t> & each_edge_ = chain[t];
      if(is_outside_edge(each_edge_.first,each_edge_.second,outside_face))
        cerr << "# error this edge is on the surface." << endl;
    }
#endif

  static vector<size_t> changed_face;

  while(!modified_chain_.empty()){
      const pair<size_t,size_t>  one_edge = modified_chain_.front();

      if(find(edges_on_surface.begin(),edges_on_surface.end(),one_edge) != edges_on_surface.end())
        {// this edge is already exist on surface, do not need to handle it
          modified_chain_.pop_front();
          continue;
        }

      size_t vertex_ =  check_surface_vertex_to_make_face(path,one_edge,fa);
      if(vertex_ == -1){
          distance_each_vertex.clear();
          find_one_ring_vertex_around_edge(tet,one_edge,ortae, one_ring_vertex);
          for(size_t t = 0; t < one_ring_vertex.size(); ++t){
              distance_each_vertex.push_back(make_pair(norm(node(colon(),one_ring_vertex[t]) - center_of_surface_path),t));
            }
          sort(distance_each_vertex.begin(),distance_each_vertex.end());
          vertex_ = one_ring_vertex[distance_each_vertex[0].second];
        }
      size_t face_idx = fa.get_face_idx(one_edge.first,one_edge.second,vertex_);

      const pair<size_t,size_t> &tets = fa.face2tet_[face_idx];

      oecit it = ortae.e2t_.find(one_edge);
      if(it == ortae.e2t_.end()) it = ortae.e2t_.end();
      if(it == ortae.e2t_.end()) return __LINE__;

      modify_face_jump_type_at_given_tets_edge(tets.first,tets.second,inner_face_jump_type, it->second);
      modified_chain_.pop_front();
      modified_chain_.push_front(make_pair(one_edge.first,vertex_));
      modified_chain_.push_front(make_pair(one_edge.second,vertex_));

      changed_face.push_back(one_edge.first);
      changed_face.push_back(one_edge.second);
      changed_face.push_back(vertex_);
    }

  ofstream ofs("changed_face.vtk");
  tri2vtk(ofs,&node[0],node.size(2),&changed_face[0],changed_face.size()/3);
  //find_one_ring_vertex_around_edge(tet,edge,ortae,one_ring_vertex);
  return 0;
}

//! @brief split edge <first,second> at (1-lambda) * first + lambda * second
// WARNINGl: if you want to iteratorly subdivide the tets around edge, you should go along first to second
// because I used fa to get faces around point "edge_second" and I didn't update fa;
// only if you go along first--> second, can you get correct answer.
// after this function called, the node_array will add one point.
int split_tets_around_edge(
    vector<size_t> &tet_array,
    vector<double> &node_array,
    jtf::mesh::one_ring_tet_at_edge &ortae,
    vector<matrixd > &frame,
    const std::pair<size_t,size_t> &edge,
    const double lambda,
    boost::unordered_map<pair<size_t,size_t>,vector<pair<size_t,size_t> > > &split_edge_map,
    boost::unordered_map<size_t,size_t > & old_tet_map_to_new,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::iterator oeit;
  typedef boost::unordered_map<pair<size_t,size_t>,vector<size_t> >::iterator mvit;
  typedef boost::unordered_map<pair<size_t,size_t>,vector<pair<size_t,size_t> > >::iterator mvpit;
  //adjust_cross_edges_right_hand_order_of_tet(tet(colon(),tet_idx),node,black_edge,edge_need_to_sub);

  matrixd point_first(3,1),point_second(3,1);
  copy(node_array.begin() + 3 * edge.first,node_array.begin() + 3 * edge.first + 3, point_first.begin());
  copy(node_array.begin() + 3 * edge.second,node_array.begin() + 3 * edge.second + 3, point_second.begin());

  matrixd new_point = (1-lambda) * point_second + lambda * point_first;

  node_array.insert(node_array.end(),new_point.begin(),new_point.end());
  const size_t new_point_idx = node_array.size() / 3 -1;
#if 0 // debug
  if(edge.first == 82 && edge.second == 313){
      cerr << "# pause" << endl;
    }
#endif
  mvpit mvpit_ = split_edge_map.find(edge);
  if(mvpit_ == split_edge_map.end()) mvpit_ = split_edge_map.find(make_pair(edge.second,edge.first));
  if(mvpit_ == split_edge_map.end()) {
      split_edge_map[edge].push_back(make_pair(edge.first,new_point_idx));
      split_edge_map[edge].push_back(make_pair(new_point_idx,edge.second));
    }else{
      cerr << "# [error] can not split again edge " << edge.first << "-->" << edge.second << endl;
      return __LINE__;
    }

  oeit it = ortae.e2t_.find(make_pair(edge.first,edge.second));
  if(it == ortae.e2t_.end()){
      it = ortae.e2t_.find(make_pair(edge.second,edge.first));
      if(it == ortae.e2t_.end()) {
          cerr << "# [error] strange: can not find edge < " << edge.second << "," //edge_need_to_sub[0] << ","
               << edge.first << ">." << endl;
          return __LINE__;
        }else{
          vector<size_t> loop = it->second;
          reverse(loop.begin(),loop.end());
          ortae.e2t_.erase(it);
          ortae.e2t_[make_pair(edge.first,edge.second)] = loop;
          it = ortae.e2t_.find(make_pair(edge.first,edge.second));//edge_need_to_sub[1],edge_need_to_sub[0]));
        }
    }// after this it->first == make_pair(edge_need_to_sub[0],edge_need_to_sub[1])
  const vector<size_t> loop = it->second;

  set<size_t> around_tet_idx;

  for(size_t t = 0; t < loop.size(); ++t) if(loop[t] != -1) around_tet_idx.insert(loop[t]);

  //map<size_t,size_t> old_tet_map_to_new;
  map<size_t,size_t> new_tet_map_to_old;
  vector<pair<size_t,size_t> > around_edges; // store the original edges which are around "edge"
  {// preprocess, get useful information before tet has been changed
    itr_matrix<size_t *> tet(4,tet_array.size() / 4,&tet_array[0]);
    itr_matrix<double *> node(3,node_array.size() / 3,&node_array[0]);
    find_one_ring_edges_around_edge(tet,node,edge,ortae,around_edges);
  }
  ///// step 0: change the tet idx, and add new tets and frames
  ///// WARNING: be careful!!! after this step the vector and node will be increased which may result in memory realloc
  {
    for(set<size_t>::const_iterator scit = around_tet_idx.begin();
        scit != around_tet_idx.end(); ++scit){
        {// add new tet
          //matrixst new_tet = tet(colon(),*scit);
          vector<size_t> new_tet(4);
          copy(tet_array.begin() + 4 * (*scit),tet_array.begin() + 4 * (*scit) + 4,new_tet.begin());
          for(size_t t = 0; t < 4; ++t)
            //if(new_tet[t] == edge_need_to_sub[1]) {
            if(new_tet[t] == edge.first) {
                new_tet[t] = new_point_idx;
                break;
              }
          tet_array.insert(tet_array.end(),new_tet.begin(),new_tet.end());
          old_tet_map_to_new[*scit] = tet_array.size() / 4 - 1;
          new_tet_map_to_old[tet_array.size() / 4 - 1] = *scit;
        }
        {// modify original tet edge
          for(size_t t = 0; t < 4; ++t){
              //if(tet(t,*scit) == edge_need_to_sub[0]) {
              if(tet_array[4 * (*scit) + t] == edge.second){
                  tet_array[4 * (*scit) + t] = new_point_idx;
                  //          if(tet(t,*scit) == edge.second) {
                  //            tet(t,*scit) = new_point_idx;
                  break;}
            }
        }
        { // add frame for new tets
          const matrixd & old_frame_tet = frame[*scit];
          if(*scit > frame.size()){
              cerr << "# [error] copy frame error: frame size = "
                   << frame.size() << " "
                   << "tet_idx = " << *scit << endl;
            }
          if(norm(old_frame_tet) < 1e-8){
              cerr << "# [error] copy frame error. norm frame = 0."  << endl;
            }
          frame.push_back(old_frame_tet);
        }
      }// end for
  }

  itr_matrix<size_t *> tet(4,tet_array.size() / 4,&tet_array[0]);
  itr_matrix<double *> node(3,node_array.size() / 3,&node_array[0]);
  //unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet)); // TODO: need to speed up

  /////  step1: change all edges' arounding tet loop
  { // update the one_ring_neighbour_tet for inner edge
    {
      vector<size_t> new_loop = loop;
      for(size_t t = 0; t < new_loop.size(); ++t)
        if(new_loop[t] != -1) {
            boost::unordered_map<size_t,size_t>::const_iterator mcit = old_tet_map_to_new.find(loop[t]);
            if(mcit == old_tet_map_to_new.end()) {
                cerr << "# [error] can not find tet: " << loop[t] << " in old_tet_map_to_new." << endl;
                return __LINE__;
              }
            new_loop[t] = mcit->second;
          }

      ortae.e2t_[make_pair(edge.first,new_point_idx)] = loop;
      ortae.e2t_[make_pair(new_point_idx,edge.second)] = new_loop;
    }

    // modify parts tet loop around edge <edge_need_to_sub[0],p>,
    // while edge<edge_need_to_sub[1],p> do not need to be modified
    set<size_t> around_point_idx;
    for(size_t t = 0; t < around_edges.size(); ++t) {
        around_point_idx.insert(around_edges[t].first);
        around_point_idx.insert(around_edges[t].second);
      }
    //for(size_t t = 0; t < around_edges.size(); ++t){
    for(set<size_t>::const_iterator scit = around_point_idx.begin();
        scit != around_point_idx.end(); ++scit){
        oeit it = ortae.e2t_.find(make_pair(edge.second,*scit));
        if(it == ortae.e2t_.end())
          it = ortae.e2t_.find(make_pair(*scit,edge.second));
        if(it == ortae.e2t_.end()) return __LINE__;
        vector<size_t> &tet_loop = it->second;
        for(size_t i = 0; i < tet_loop.size(); ++i){
            if(tet_loop[i] == -1) continue;
            if(find(loop.begin(),loop.end(),tet_loop[i]) != loop.end()) tet_loop[i] = old_tet_map_to_new[tet_loop[i]];
          }
      }
    // modify the around edges' tets loop
    // since around_edge $edge_a$ is in right hand order of subdivided edge,
    // we just need to insert the new tet idx after the original tet idx in $edge_a$ 's right hand loop
    vector<size_t> new_tet_loop;
    for(size_t t = 0; t < around_edges.size(); ++t){
        oeit it = ortae.e2t_.find(around_edges[t]);
        if(it == ortae.e2t_.end()) {
            it = ortae.e2t_.find(make_pair(around_edges[t].second,around_edges[t].first));
            if(it == ortae.e2t_.end()) return __LINE__;
            vector<size_t> tet_loop = it->second;
            reverse(tet_loop.begin(),tet_loop.end());
            ortae.e2t_.erase(it);
            ortae.e2t_[around_edges[t]] = tet_loop;
            it = ortae.e2t_.find(around_edges[t]);
          }
        vector<size_t> &tet_loop = it->second;
        if(tet_loop.size() == 2 && (tet_loop.front() == -1 || tet_loop.back() == -1)){// the sharp edge will be regular surface edge after splitting
            const size_t inner_tet = (tet_loop.front() == -1)?tet_loop.back():tet_loop.front();
            tet_loop.resize(3);
            tet_loop.front() = -1;
            tet_loop.back() = -1;
            tet_loop[1] = inner_tet;
          }
        new_tet_loop.resize(tet_loop.size() + 1);
        for(size_t i = 0; i < tet_loop.size(); ++i){
            new_tet_loop[i] = tet_loop[i];
            if(tet_loop[i] == -1) continue;
            if(find(loop.begin(),loop.end(),tet_loop[i]) != loop.end()){
                new_tet_loop[i+1] = old_tet_map_to_new[tet_loop[i]];//new_loop[i];
                copy(tet_loop.begin() + i + 1,tet_loop.end(),new_tet_loop.begin() + i + 2);
                //assert(new_loop.front() == new_loop.back());
                break;
              }
          }
        swap(tet_loop,new_tet_loop);
      }

    // modfiy the inner edges's arounding tet loop
    boost::unordered_map<pair<size_t,size_t>,vector<size_t> > temp_loop_for_inner_edge;

    set<pair<size_t,size_t> > around_edges_set(around_edges.begin(),around_edges.end());
    for(set<size_t>::const_iterator scit = around_tet_idx.begin();
        scit != around_tet_idx.end(); ++scit){
        //cerr << tet(colon(),*scit) << endl;
        for(set<pair<size_t,size_t> >::iterator sit = around_edges_set.begin();
            sit != around_edges_set.end(); ++sit){
            if(find(tet(colon(),*scit).begin(),tet(colon(),*scit).end(),sit->first)
               != tet(colon(),*scit).end()
               &&find(tet(colon(),*scit).begin(),tet(colon(),*scit).end(),sit->second)
               != tet(colon(),*scit).end()){
                mvit itv = temp_loop_for_inner_edge.find(make_pair(new_point_idx,sit->second));
                if(itv == temp_loop_for_inner_edge.end()) {
                    vector<size_t> temp_tet_loop0(4,-1);
                    temp_tet_loop0[0] = old_tet_map_to_new[*scit];
                    temp_tet_loop0[1] = *scit;
                    temp_loop_for_inner_edge[make_pair(new_point_idx,sit->second)] = temp_tet_loop0;
                  }else{
                    vector<size_t> &temp_tet_loop0 = itv->second;
                    assert(temp_tet_loop0[0] == -1 && temp_tet_loop0[1] == -1); // not set
                    temp_tet_loop0[0] = old_tet_map_to_new[*scit];
                    temp_tet_loop0[1] = *scit;
                  }

                itv = temp_loop_for_inner_edge.find(make_pair(new_point_idx,sit->first));
                if(itv == temp_loop_for_inner_edge.end()) {
                    vector<size_t> temp_tet_loop1(4,-1);
                    temp_tet_loop1[2] = *scit;
                    temp_tet_loop1[3] = old_tet_map_to_new[*scit];
                    temp_loop_for_inner_edge[make_pair(new_point_idx,sit->first)] = temp_tet_loop1;
                  }else{
                    vector<size_t> &temp_tet_loop1 = itv->second;
                    assert(temp_tet_loop1[2] == -1 && temp_tet_loop1[3] == -1); // not set
                    temp_tet_loop1[2] = *scit;
                    temp_tet_loop1[3] = old_tet_map_to_new[*scit];
                  }
                around_edges_set.erase(sit);
                break;
              }
          }
      }
    for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator it =
        temp_loop_for_inner_edge.begin();
        it != temp_loop_for_inner_edge.end(); ++it){
        const vector<size_t> & vs = it->second;
        if(count(vs.begin(),vs.end(),-1) == 0){// inner not outside edge
            vector<size_t> temp_around_tet = vs;
            temp_around_tet.push_back(temp_around_tet.front());
            //        cerr << "# [info] around tet idx ";
            //        for(size_t k = 0; k < temp_around_tet.size(); ++k) cerr << temp_around_tet[k] << " ";
            ortae.e2t_[it->first] = temp_around_tet;
          }else{ // outside edge
            vector<size_t> temp_around_tet(vs.size());
            if(vs.front() == -1){
                copy(vs.begin() + 1,vs.end(),temp_around_tet.begin());
                temp_around_tet.back() = -1;
                ortae.e2t_[it->first] = temp_around_tet;
              }else{ // vs.back == -1
                temp_around_tet.front() = -1;
                copy(vs.begin(),vs.end() - 1,temp_around_tet.begin() + 1);
                ortae.e2t_[it->first] = temp_around_tet;
              }
          }
      }// end for
  }// end step 1

  ///// step2: change all faces jump type
  {
    if(loop.size() != 2){ // no need to handle sharp edge
        for(size_t t = 0; t < loop.size() - 1; ++t){
            if(loop[t] == -1) continue;
            boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mcit
                = inner_face_jump_type.find(make_pair(loop[t],loop[t+1]));
            if(mcit == inner_face_jump_type.end()) continue; // since the original face type is identity, no need to build a type

            inner_face_jump_type[make_pair(old_tet_map_to_new[loop[t]],
                old_tet_map_to_new[loop[t+1]])]
                = mcit->second;
            mcit = inner_face_jump_type.find(make_pair(loop[t+1],loop[t]));
            assert(mcit != inner_face_jump_type.end());
            inner_face_jump_type[make_pair(old_tet_map_to_new[loop[t+1]],
                old_tet_map_to_new[loop[t]])]
                = mcit->second;
          }

        // modify the faces around edge.second
        typedef  boost::unordered_map<pair<size_t,size_t>,size_t>::iterator mit;
        for(size_t t = 0; t < around_edges.size(); ++t){
            //const pair<size_t,size_t> &pair_tet_ = fa->query(around_edges[t].first,around_edges[t].second,edge.second);
            const pair<size_t,size_t> pair_tet = get_adj_tets(around_edges[t].first,around_edges[t].second,edge.second,&tet_array[0],tet_array.size()/4);
            //cerr << "# [comp] fa--pair_tet: " << pair_tet_.first << ";" << pair_tet_.second << endl;
            //cerr << "# [comp] new--pair_tet: " << pair_tet.first << ";" << pair_tet.second << endl;
            if(pair_tet.first == -1 || pair_tet.second == -1) continue;

            const size_t new_tet_idx = (find(tet(colon(),pair_tet.first).begin(),tet(colon(),pair_tet.first).end(),new_point_idx)
                                        != tet(colon(),pair_tet.first).end())?pair_tet.first:pair_tet.second;
            const size_t old_tet_idx = new_tet_map_to_old[new_tet_idx];
            const size_t adj_tet_idx = pair_tet.first == new_tet_idx ? pair_tet.second: pair_tet.first;
            mit it0 = inner_face_jump_type.find(make_pair(adj_tet_idx,old_tet_idx));
            mit it1 = inner_face_jump_type.find(make_pair(old_tet_idx,adj_tet_idx));
            // if original face type is identity, we do not need to set the new face, because it's still identity
            if(it0 == inner_face_jump_type.end() || it1 == inner_face_jump_type.end()) continue;
            if(new_tet_idx == pair_tet.first){
                inner_face_jump_type[pair_tet] = it1->second;
                inner_face_jump_type[make_pair(pair_tet.second,pair_tet.first)] = it0->second;
              }else{
                inner_face_jump_type[pair_tet] = it0->second;
                inner_face_jump_type[make_pair(pair_tet.second,pair_tet.first)] = it1->second;
              }

            inner_face_jump_type.erase(it0);
            inner_face_jump_type.erase(it1);
          }
      }
  }
  ortae.e2t_.erase(it); // delete original edge

#if 0 // debug
  if(edge.first == 82 && edge.second == 313){
      one_ring_tet_at_edge::e2tet_type::const_iterator oecit = ortae.e2t_.find(make_pair(313,82));
      if(oecit != ortae.e2t_.end()){
          const vector<size_t> & vec = oecit->second;
        }
    }
#endif
  return 0;
}



//! @brief this function is used to modify type of face<B,A,e>, while edge<B,A> and edge<D,C> is in each other's right hand order
//       A
//    / /| \
//  C--e-|--D
//    \ \| /
//       B
static int fix_face_type_with_stupid_method(
    const pair<size_t,size_t> &edge_BA,
    const pair<size_t,size_t> &edge_DC,
    const size_t &point_e_idx,
    const size_t &right_hand_face_type,
    const matrixst &tet,
    const matrixd &node,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oceit;

  //const size_t add_point_idx = node.size(2) /  - 1;
  oceit it = ortae.e2t_.find(make_pair(edge_DC.first,point_e_idx));
  if(it == ortae.e2t_.end()) it = ortae.e2t_.find(make_pair(point_e_idx,edge_DC.first));
  if(it == ortae.e2t_.end()) return __LINE__;

  const vector<size_t> & around_tet = it->second;
  size_t tet_idx = -1;
  for(size_t t = 0; t < around_tet.size(); ++t){
      if(around_tet[t] == -1) continue;
      if(find(tet(colon(),around_tet[t]).begin(),tet(colon(),around_tet[t]).end(),edge_BA.first)
         != tet(colon(),around_tet[t]).end()
         && find(tet(colon(),around_tet[t]).begin(),tet(colon(),around_tet[t]).end(),edge_BA.second)
         != tet(colon(),around_tet[t]).end()) {
          tet_idx = around_tet[t];break;
        }
    }
  if(tet_idx == -1) {
      cerr << "# [error] strange: can not find a tet contain black edge." << endl;
      return __LINE__;
    }

  oceit it_new = ortae.e2t_.find(make_pair(point_e_idx,edge_DC.second));
  if(it_new == ortae.e2t_.end()) it_new = ortae.e2t_.find(make_pair(edge_DC.second,point_e_idx));
  if(it_new == ortae.e2t_.end()) return __LINE__;
  const vector<size_t> & around_tet_new = it_new->second;
  size_t tet_idx_new = -1;
  for(size_t t = 0; t < around_tet_new.size(); ++t){
      if(around_tet_new[t] == -1) continue;
      if(find(tet(colon(),around_tet_new[t]).begin(),tet(colon(),around_tet_new[t]).end(),edge_BA.first)
         != tet(colon(),around_tet_new[t]).end()
         && find(tet(colon(),around_tet_new[t]).begin(),tet(colon(),around_tet_new[t]).end(),edge_BA.second)
         != tet(colon(),around_tet_new[t]).end()) {
          tet_idx_new = around_tet_new[t];break;
        }
    }
  if(tet_idx_new == -1) {
      cerr << "# [error] strange: can not find a tet contain black edge." << endl;
      return __LINE__;
    }

#if 0 // check the modification is correct, which means will reduce the compound axis num
  {
    oceit it_black = ortae.e2t_.find(make_pair(edge_BA.first,edge_BA.second));
    const vector<size_t> & tet_around_black = it_black->second;
    matrixd rot_black = eye<double>(3);
    for(size_t t = 0; t < tet_around_black.size() - 1; ++t){
        map<pair<size_t,size_t>,size_t>::const_iterator cit =
            inner_face_jump_type.find(make_pair(tet_around_black[t],tet_around_black[t+1]));
        if(cit == inner_face_jump_type.end()) continue;
        rot_black = temp(rot_black * type_transition2(cit->second));
      }
    cerr << "# [info] edge type = " << type_transition1(rot_black) << endl;
  }
#endif

  // right hand order around black edge is : tet_idx --> tet_idx_new
  inner_face_jump_type[make_pair(tet_idx,tet_idx_new)] = right_hand_face_type;
  inner_face_jump_type[make_pair(tet_idx_new,tet_idx)] = type_transition1(trans(type_transition2(right_hand_face_type)));

#if 0 // check the modification is correct, which means will reduce the compound axis num
  {
    oceit it_black = ortae.e2t_.find(make_pair(edge_BA.first,edge_BA.second));
    const vector<size_t> & tet_around_black = it_black->second;
    matrixd rot_black = eye<double>(3);
    for(size_t t = 0; t < tet_around_black.size() - 1; ++t){
        map<pair<size_t,size_t>,size_t>::const_iterator cit =
            inner_face_jump_type.find(make_pair(tet_around_black[t],tet_around_black[t+1]));
        if(cit == inner_face_jump_type.end()) continue;
        rot_black = temp(rot_black * type_transition2(cit->second));
      }
    cerr << "# [info] edge type = " << type_transition1(rot_black) << endl;
  }
#endif

  return 0;
}

int split_to_remove_black_line(
    vector<size_t> &tet_array,
    vector<double> &node_array,
    const size_t &tet_idx,
    const jtf::mesh::face2tet_adjacent &fa,
    jtf::mesh::one_ring_tet_at_edge & ortae,
    vector<matrixd > &frame,
    const std::pair<size_t,size_t> &black_edge,
    const matrixd &angle_zyx,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    boost::unordered_map<pair<size_t,size_t>,vector<pair<size_t,size_t> > > &split_edge_map,
    size_t black_line_type) // 1:compound of two axes, 0 means compound of three axes
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::iterator oeit;
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oceit;
  vector<size_t> edge_need_to_sub;

  {// find an edge in right hand order of black edges, which is need to be subdivied
    itr_matrix<size_t *> tet(4,tet_array.size() / 4,&tet_array[0]);
    itr_matrix<double *> node(3,node_array.size() / 3,&node_array[0]);

    adjust_cross_edges_right_hand_order_of_tet(tet(colon(),tet_idx),node,black_edge,edge_need_to_sub);
  }
  if(black_line_type == 1) // compound of two axes
    {
#if 0 //debug
      if((black_edge.first == 35941 && black_edge.second == 35958) ||
         (black_edge.first == 35958 && black_edge.second == 35941))
        {
          one_ring_tet_at_edge::e2tet_type::const_iterator oecit =
              ortae.e2t_.find(make_pair(35941,35958));
          if(oecit == ortae.e2t_.end()) oecit = ortae.e2t_.find(make_pair(35958,35941));
          if(oecit == ortae.e2t_.end()) {
              cerr << "# error." << endl;
            };

          const vector<size_t> &loop = oecit->second;
          vector<size_t> check_tets;
          //cerr << " tet pair " << loop[t] << " " << loop[t+1] << " type: " << endl;
          for(size_t t = 0; t < loop.size() - 1; ++t){
              cerr << inner_face_jump_type[make_pair(loop[t],loop[t+1])] << " ";
              check_tets.push_back(tet_array[4 * loop[t] + 0]);
              check_tets.push_back(tet_array[4 * loop[t] + 1]);
              check_tets.push_back(tet_array[4 * loop[t] + 2]);
              check_tets.push_back(tet_array[4 * loop[t] + 3]);
            }
          ofstream ofs("debug_0.vtk");

          tet2vtk(ofs,&node_array[0],node_array.size()/3,&check_tets[0],check_tets.size()/4);
          cerr << endl;
        }
#endif

      boost::unordered_map<size_t,size_t> split_tet_map;
      if(split_tets_around_edge(tet_array,node_array,ortae,frame,make_pair(edge_need_to_sub[1],edge_need_to_sub[0]),
                                1/2.0,split_edge_map,split_tet_map,inner_face_jump_type))
        cerr << "# [error] shit happens while splitting tets." << endl;

      // since the tet and node has increased, it may result in memory realloc, we should be very careful while using the memory address
      itr_matrix<size_t *> tet(4,tet_array.size() / 4,&tet_array[0]);
      itr_matrix<double *> node(3,node_array.size() / 3,&node_array[0]);

      size_t point_add_idx = node.size(2) - 1;
      matrixd modified_zyx = zeros<double>(3,1);
      size_t first_non_zero_angle = -1;
      for(size_t t = 0; t < 3; ++t){
          if(fabs(angle_zyx[t]) < 1e-8) continue;
          first_non_zero_angle = t;
        }
      if(first_non_zero_angle == -1) {
          cerr << "# [error] can not find first non zero angle for rot matrix." << endl;
          return __LINE__;
        }
      modified_zyx[first_non_zero_angle] = angle_zyx[first_non_zero_angle] * -1;
      matrixd rot = zeros<double>(3,3);
      convert_euler_angle_to_rot(&modified_zyx[0],&rot[0]);

      size_t rot_type = type_transition1(rot);

      fix_face_type_with_stupid_method(black_edge,make_pair(edge_need_to_sub[1],edge_need_to_sub[0]),
          point_add_idx,rot_type,tet,node,ortae,inner_face_jump_type);

#if 1 // ugly modification
      // since the insert edge will be of a strange type
      // for example jump faces around edge maybe 6,2,8, which is supposed to be 2,8,6,
      // different order may result in different type and
#endif

#if 1 // check the new type of face which contain black edge and the subdivied point, it should be modified to a regular singularity
      oceit it_black = ortae.e2t_.find(black_edge);
      if(it_black == ortae.e2t_.end()) it_black = ortae.e2t_.find(make_pair(black_edge.second,black_edge.first));
      if(it_black == ortae.e2t_.end()) {
          cerr << "# [error] can not find this black edge." << endl;
          return __LINE__;
        }
      const vector<size_t> & loop_around_black_edge = it_black->second;
      assert(loop_around_black_edge.front() == loop_around_black_edge.back());
      matrixd rot_ = eye<double>(3);
      for(size_t t = 0; t < loop_around_black_edge.size() - 1; ++t){
          boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator it =
              inner_face_jump_type.find(make_pair(loop_around_black_edge[t],loop_around_black_edge[t+1]));
          if(it == inner_face_jump_type.end()) continue;
          //cerr << "# [info] edge " << it->first.first << "," << it->first.second << " :" << it->second << endl;
          assert(it->second < 24);
          //cerr << "## right multi " << type_transition2(it->second);
          rot_ = temp(rot_ * type_transition2(it->second));
          //cerr << "## after right multi " << rot_;
        }
      assert(type_transition1(rot_) <= 9);
#endif


#if 0 //debug
      if((black_edge.first == 35941 && black_edge.second == 35958) ||
         (black_edge.first == 35958 && black_edge.second == 35941))
        {
          {
            one_ring_tet_at_edge::e2tet_type::const_iterator oecit =
                ortae.e2t_.find(make_pair(57881,35958));
            if(oecit == ortae.e2t_.end()) oecit = ortae.e2t_.find(make_pair(35958,57881));
            if(oecit == ortae.e2t_.end()) {
                cerr << "# error." << endl;
              };

            const vector<size_t> &loop = oecit->second;
            vector<size_t> check_tets;
            vector<size_t> jump_faces;
            size_t faces[3];
            itr_matrix<size_t *> tet(4,tet_array.size() / 4,&tet_array[0]);
            //cerr << " tet pair " << loop[t] << " " << loop[t+1] << " type: " << endl;
            for(size_t t = 0; t < loop.size() - 1; ++t){
                //cerr << inner_face_jump_type[make_pair(loop[t],loop[t+1])] << " ";
                check_tets.push_back(tet_array[4 * loop[t] + 0]);
                check_tets.push_back(tet_array[4 * loop[t] + 1]);
                check_tets.push_back(tet_array[4 * loop[t] + 2]);
                check_tets.push_back(tet_array[4 * loop[t] + 3]);
                map<pair<size_t,size_t>,size_t>::const_iterator mcit =
                    inner_face_jump_type.find(make_pair(loop[t],loop[t+1]));
                if(mcit == inner_face_jump_type.end()) continue;
                if(is_trivial_type(mcit->second)) continue;
                find_common_face(tet,loop[t],loop[t+1],&faces[0]);
                jump_faces.push_back(faces[0]);
                jump_faces.push_back(faces[1]);
                jump_faces.push_back(faces[2]);
              }
            cerr << endl;
            ofstream ofs("debug_1.vtk");
            ofstream ofs_face("debug_jump_face.vtk");
            tet2vtk(ofs,&node_array[0],node_array.size()/3,&check_tets[0],check_tets.size()/4);
            tri2vtk(ofs_face,&node_array[0],node_array.size()/3,&jump_faces[0],jump_faces.size()/3);
          }

          {
            one_ring_tet_at_edge::e2tet_type::const_iterator oecit =
                ortae.e2t_.find(make_pair(57881,35941));
            if(oecit == ortae.e2t_.end()) oecit = ortae.e2t_.find(make_pair(35941,57881));
            if(oecit == ortae.e2t_.end()) {
                cerr << "# error." << endl;
              };

            const vector<size_t> &loop = oecit->second;
            vector<size_t> check_tets;
            vector<size_t> jump_faces;
            size_t faces[3];
            itr_matrix<size_t *> tet(4,tet_array.size() / 4,&tet_array[0]);
            //cerr << " tet pair " << loop[t] << " " << loop[t+1] << " type: " << endl;
            for(size_t t = 0; t < loop.size() - 1; ++t){
                //cerr << inner_face_jump_type[make_pair(loop[t],loop[t+1])] << " ";
                check_tets.push_back(tet_array[4 * loop[t] + 0]);
                check_tets.push_back(tet_array[4 * loop[t] + 1]);
                check_tets.push_back(tet_array[4 * loop[t] + 2]);
                check_tets.push_back(tet_array[4 * loop[t] + 3]);
                map<pair<size_t,size_t>,size_t>::const_iterator mcit =
                    inner_face_jump_type.find(make_pair(loop[t],loop[t+1]));
                if(mcit == inner_face_jump_type.end()) continue;
                if(is_trivial_type(mcit->second)) continue;
                find_common_face(tet,loop[t],loop[t+1],&faces[0]);
                jump_faces.push_back(faces[0]);
                jump_faces.push_back(faces[1]);
                jump_faces.push_back(faces[2]);
              }
            cerr << endl;
            ofstream ofs("debug_2.vtk");
            ofstream ofs_face("debug_jump_face2.vtk");
            tet2vtk(ofs,&node_array[0],node_array.size()/3,&check_tets[0],check_tets.size()/4);
            tri2vtk(ofs_face,&node_array[0],node_array.size()/3,&jump_faces[0],jump_faces.size()/3);
          }
        }

#endif

    }else{ // compount of three axes
#if 1
      boost::unordered_map<size_t,size_t> split_tet_map;
      size_t first_sub_point = -1;
      size_t first_non_zero_angle = -1;
      size_t second_non_zero_angle = -1;
      {// split at first 1/3 point of edge "first --|------> second"
        if(split_tets_around_edge(tet_array,node_array,ortae,frame,make_pair(edge_need_to_sub[1],edge_need_to_sub[0]),
                                  2.0/3.0,split_edge_map,split_tet_map,inner_face_jump_type))
          cerr << "# [error] shit happens while splitting tets" << endl;

        itr_matrix<size_t *> tet(4,tet_array.size() / 4,&tet_array[0]);
        itr_matrix<double *> node(3,node_array.size() / 3,&node_array[0]);

        for(size_t t = 0; t < 3; ++t){
            if(fabs(angle_zyx[t]) < 1e-8) continue;
            first_non_zero_angle = t;
          }
        if(first_non_zero_angle == -1) {
            cerr << "# [error] can not find first non zero angle for rot matrix." << endl;
            return __LINE__;
          }

        first_sub_point = node.size(2) - 1;
        matrixd modified_zyx = zeros<double>(3,1);
        modified_zyx[first_non_zero_angle] = angle_zyx[first_non_zero_angle] * -1; // apply a rotation of Rx^{-1}
        matrixd rot = zeros<double>(3,3);
        convert_euler_angle_to_rot(&modified_zyx[0],&rot[0]);
        size_t rot_type = type_transition1(rot);

        // to modify the face introduced by splitting method
        fix_face_type_with_stupid_method(black_edge,make_pair(edge_need_to_sub[1],
                                         edge_need_to_sub[0]),
            first_sub_point,rot_type,tet,node,
            ortae,inner_face_jump_type);
      }
#endif
#if 1
      {// split at second 1/3 point of edge "first ------|--> second"
        if(split_tets_around_edge(tet_array,node_array,ortae,frame,make_pair(first_sub_point,edge_need_to_sub[0]),
                                  1.0/3.0,split_edge_map,split_tet_map,inner_face_jump_type))
          cerr << "# [error] shit happens while splitting tets" << endl;

        itr_matrix<size_t *> tet(4,tet_array.size() / 4,&tet_array[0]);
        itr_matrix<double *> node(3,node_array.size() / 3,&node_array[0]);

        for(size_t t = first_non_zero_angle; t < 3; ++t){
            if(fabs(angle_zyx[t]) < 1e-8) continue;
            second_non_zero_angle = t;
          }
        if(second_non_zero_angle == -1) {
            cerr << "# [error] can not find second non zero angle for rot matrix." << endl;
            return __LINE__;
          }

        size_t second_sub_point = node.size(2) - 1;
        matrixd modified_zyx = zeros<double>(3,1);
        modified_zyx[second_non_zero_angle] = angle_zyx[second_non_zero_angle] * -1; // apply a rotation of Ry^{-1}
        matrixd rot = zeros<double>(3,3);
        convert_euler_angle_to_rot(&modified_zyx[0],&rot[0]);
        size_t rot_type = type_transition1(rot);

        fix_face_type_with_stupid_method(black_edge,make_pair(edge_need_to_sub[1],edge_need_to_sub[0]),
            second_sub_point,rot_type,tet,node,ortae,inner_face_jump_type);
      }
#endif
    }
  return 0;
}

int cal_singularity_type_of_edge(
    const pair<size_t,size_t> & edge,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  oecit it = ortae.e2t_.find(edge);
  if(it == ortae.e2t_.end())
    it = ortae.e2t_.find(make_pair(edge.second,edge.first));
  const vector<size_t> & loop = it->second;
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mcit;

  matrixd rot = eye<double>(3);
  for(size_t t = 0; t < loop.size() - 1; ++t){
      mcit p = inner_face_jump_type.find(make_pair(loop[t],loop[t+1]));
      if(p == inner_face_jump_type.end()) continue;
      rot = temp(rot * type_transition2(p->second));
    }
  return type_transition1(rot);
}

static int apply_jump_on_cut_faces_for_each_group(
    const matrixst & tet,
    const vector<pair<size_t,size_t> >  &need_to_apply_jump_tet_pair,
    boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const deque<pair<size_t,size_t> > &chain,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge & ortae )
{
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mcit;

  const pair<size_t,size_t> & begin_edge = chain.front();
  pair<size_t,size_t> begin_edge_tet_pair(-1,-1);
  jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator it = ortae.e2t_.find(begin_edge);
  vector<size_t> loop;
  if(it != ortae.e2t_.end())
    loop = it->second;
  else{
      it = ortae.e2t_.find(make_pair(begin_edge.second,begin_edge.first));
      if(it == ortae.e2t_.end()) {
          cerr << "# [error] can not find this singularity edge in ortae." << endl;
          return __LINE__;
        }
      else{
          loop = it->second;
          reverse(loop.begin(),loop.end());
        }
    }
  // now the loop tets around begin_edge is in correct right hand format
  vector<size_t> face(3,-1);
  for(size_t t = 0; t < need_to_apply_jump_tet_pair.size(); ++t){
      const pair<size_t,size_t> & tet_pair = need_to_apply_jump_tet_pair[t];
      common_tet_face(&tet(0,tet_pair.first),
                      &tet(0,tet_pair.second),&face[0]);
      if(find(face.begin(),face.end(),begin_edge.first)!= face.end()
         && find(face.begin(),face.end(),begin_edge.second)!= face.end()){
          // begin_edge exist in this tet pair
          begin_edge_tet_pair = tet_pair;
          break;
        }
    }

  {// calculate the rotation jump from pair first to pair second
    if(begin_edge_tet_pair.first == -1 && begin_edge_tet_pair.second == -1) {
        cerr << "# [error] can not find singularity edge in cut faces." << endl;
        return 0;
      }

    size_t first_find_tet_idx = -1;
    for(size_t t = 0; t < loop.size() - 1; ++t){
        if(loop[t] == begin_edge_tet_pair.first||
           loop[t] == begin_edge_tet_pair.second){
            first_find_tet_idx = t;
            break;
          }
      }
    assert(first_find_tet_idx != -1);
    if(loop[first_find_tet_idx+1] == begin_edge_tet_pair.first||
       loop[first_find_tet_idx+1] == begin_edge_tet_pair.second){
        vector<size_t> reverse_loop;
        reverse_loop.insert(reverse_loop.end(),
                            loop.begin()+first_find_tet_idx,loop.end()-1);
        reverse_loop.insert(reverse_loop.end(),loop.begin(),loop.begin()+first_find_tet_idx);
        reverse_loop.push_back(reverse_loop.front());
        //swap(reverse_loop,loop);
        loop = reverse_loop;
      }else{
        if(first_find_tet_idx == 0){// loop order may be: firtst * *  * * second first
            vector<size_t> reverse_loop(loop.size(),0);
            copy(loop.begin(), loop.end() - 1,reverse_loop.begin() + 1);
            reverse_loop.front() = reverse_loop.back();
            //swap(reverse_loop,loop);
            loop = reverse_loop;
          }else if(loop[first_find_tet_idx-1] == begin_edge_tet_pair.first||
                   loop[first_find_tet_idx-1] == begin_edge_tet_pair.second)
          {
            vector<size_t> reverse_loop;
            reverse_loop.insert(reverse_loop.end(),
                                loop.begin()+first_find_tet_idx-1,loop.end()-1);
            reverse_loop.insert(reverse_loop.end(),loop.begin(),loop.begin()+first_find_tet_idx-1);
            reverse_loop.push_back(reverse_loop.front());
            //swap(reverse_loop,loop);
            loop = reverse_loop;
          }else{
            cerr << "# [error] adjacent tets order wrong." << endl;
            cerr << "# [error] adjacent tet pair: "
                 << begin_edge_tet_pair.first << " "
                 << begin_edge_tet_pair.second << endl;
            cerr << "# [error] ";
            for(size_t t = 0; t < loop.size(); ++t) cerr << loop[t] << " ";
            cerr << endl;
            return __LINE__;
          }
      }

#if 1
    {
      if((loop[0] != begin_edge_tet_pair.first && loop[0] != begin_edge_tet_pair.second) ||
         (loop[1] != begin_edge_tet_pair.first && loop[1] != begin_edge_tet_pair.second))
        {
          cerr << "# [error] wrong loop: " ;
          for(size_t t = 0; t < loop.size(); ++t) cerr << loop[t] << " ";
          cerr << endl;
          cerr << "# [error] tet pair = "
               << begin_edge_tet_pair.first << " "
               << begin_edge_tet_pair.second << endl;
        }
    }
#endif

    matrixd rot = eye<double>(3);

    for(size_t t = 0; t < loop.size()-1; ++t){
        boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator it =
            inner_face_jump_type.find(make_pair(loop[t],loop[t+1]));
        if(it == inner_face_jump_type.end()) continue;
        rot = temp(rot * type_transition2(it->second));
      }
    if(is_trivial_type(type_transition1(rot))){
        cerr << "# [info] this singularity edge is already changed to non-singularity." << endl;
        return 0;
      }

#if 1 // check the singularity type should be the same along each chain
    {
      vector<size_t> stvec;
      set<size_t> singularity_type;
      for(size_t t = 0; t < chain.size(); ++t){
          stvec.push_back(cal_singularity_type_of_edge(chain[t],ortae,
                                                       inner_face_jump_type));
          singularity_type.insert(stvec.back());
        }

      if(singularity_type.size() != 1){
          cerr << "# [error] singularity type along this chain is not the same." << endl;
          cerr << "# [error] singularity type is : " ;
          for(vector<size_t>::const_iterator scit = stvec.begin();
              scit != stvec.end(); ++scit) cerr << *scit << " ";
          cerr << endl;
        }
    }
#endif

    //assert(norm(rot) > 1e-6);
    matrixd rot_inv = trans(rot); //
    assert(norm(rot_inv - rot) > 1e-8);
    if(loop.front() == begin_edge_tet_pair.first){
        if(loop[1] != begin_edge_tet_pair.second){
            cerr << "# [error] wrong loop. " << __LINE__ << endl;
          }
        // need to left multiply rot^-1 to each cut tet pair
        for(size_t t = 0; t < need_to_apply_jump_tet_pair.size(); ++t){
            const pair<size_t,size_t> & tet_pair = need_to_apply_jump_tet_pair[t];
            mcit it = inner_face_jump_type.find(tet_pair);
            if(it == inner_face_jump_type.end()){
                inner_face_jump_type[tet_pair] = type_transition1(rot_inv);
                inner_face_jump_type[make_pair(tet_pair.second,tet_pair.first)]
                    = type_transition1(rot);
              }else{
                inner_face_jump_type[tet_pair]
                    = type_transition1(rot_inv * type_transition2(it->second));
                inner_face_jump_type[make_pair(tet_pair.second,tet_pair.first)]
                    = type_transition1(trans(type_transition2(inner_face_jump_type[tet_pair])));
              }
          }
      }else{ // need to left multiply rot^-1 to each cut tet reverse pair
        if(loop[0] != begin_edge_tet_pair.second &&
           loop[1] != begin_edge_tet_pair.first) {
            cerr << "# [error] wrong loop. " <<  __LINE__ << endl;
          }
        for(size_t t = 0; t < need_to_apply_jump_tet_pair.size(); ++t){
            const pair<size_t,size_t> & tet_pair = need_to_apply_jump_tet_pair[t];
            mcit it = inner_face_jump_type.find(tet_pair);
            if(it == inner_face_jump_type.end()){
                inner_face_jump_type[make_pair(tet_pair.second,tet_pair.first)]
                    = type_transition1(rot_inv);
                inner_face_jump_type[tet_pair] = type_transition1(rot);
              }else{
                inner_face_jump_type[make_pair(tet_pair.second,tet_pair.first)]
                    = type_transition1(rot_inv * trans(type_transition2(it->second)));
                //type_transition1(trans(type_transition2(inner_face_jump_type[tet_pair])));
                inner_face_jump_type[tet_pair]
                    = type_transition1( trans(type_transition2(
                                                inner_face_jump_type[make_pair(tet_pair.second,tet_pair.first)])));//type_transition2(it->second) * rot);
                //cerr << endl;
              }
          }
      }
  }// end calculate rotation jump

#if 1
  {// check whether the singularity edge is removed
    //matrixd rot = eye<double>(3);
    for(size_t t = 0; t < chain.size(); ++t){
        if(!is_trivial_type(cal_singularity_type_of_edge(chain[t],ortae,inner_face_jump_type)))
          cerr << "# [error] this singularity edge should be non-singularity." << endl;
      }

  }
#endif

  return 0;
}

bool is_edge_contained_in_group(const pair<size_t,size_t> & edge,
                                const vector<pair<size_t,size_t> > & tet_group,
                                const matrixst & tet)
{
  vector<size_t> face(3,-1);
  //set<pair<size_t,size_t> > all_edges;
  for(size_t t = 0; t < tet_group.size(); ++t){
      const pair<size_t,size_t> & tet_pair = tet_group[t];
      common_tet_face(&tet(0,tet_pair.first),
                      &tet(0,tet_pair.second),&face[0]);
      for(size_t i = 0; i < 3; ++i){
          pair<size_t,size_t> one_edge(face[i],face[(i+1)%3]);
          if((one_edge.first == edge.first && one_edge.second == edge.second)||
             (one_edge.first == edge.second && one_edge.second == edge.first))
            return true;
        }
    }
  return false;
}


// Warning: no used
//! @brief this function is used to calculate a roation jump on each cut face,
// and apply this jump for all cut faces
// since the input tet pairs are all in a special order, first of pair is in one side
// while second of pair is on the other side
int apply_jump_on_cut_faces(
    const matrixst & tet,
    const vector<vector<pair<size_t,size_t> > > &need_to_apply_jump_tet_pair,
    boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const deque<pair<size_t,size_t> > &chain,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge & ortae )
{
  //set<size_t> singularity_part_idx;
  vector<vector<size_t> > chain_hit_recorder(need_to_apply_jump_tet_pair.size());
  for(size_t t = 0; t < need_to_apply_jump_tet_pair.size(); ++t){
      const vector<pair<size_t,size_t> > & one_group = need_to_apply_jump_tet_pair[t];
      deque<pair<size_t,size_t> > chain_part;
      for(size_t i = 0; i < chain.size(); ++i){
          if(is_edge_contained_in_group(chain[i],one_group,tet)){
              chain_part.push_back(chain[i]);
              //singularity_part_idx.insert(i);
              chain_hit_recorder[t].push_back(i);
            }
        }
      if(!chain_part.empty()){
          if(apply_jump_on_cut_faces_for_each_group(tet,one_group,inner_face_jump_type,
                                                    chain_part,fa,ortae))
            return __LINE__;
        }
#if 1 // check
      if(chain_hit_recorder[t].size() == 1) continue;
      for(size_t i = 1; i < chain_hit_recorder[t].size(); ++i){
          if(chain_hit_recorder[t][i] - chain_hit_recorder[t][i-1] != 1){
              cerr << "# [error] this singularity chain is broken." << endl;
            }
        }
#endif
    }

  return 0;
}

int get_around_vertex(
    const pair<size_t,size_t> & edge,
    const matrixst & tet,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    vector<size_t> &around_vertex)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  oecit it = ortae.e2t_.find(edge);
  if(it == ortae.e2t_.end())
    it = ortae.e2t_.find(make_pair(edge.second,edge.first));

  if(it == ortae.e2t_.end()) {
      cerr << "# [error] can not find edge <" << edge.first
           << " " << edge.second << "> in ortae." << endl;
      return __LINE__;
    }

  const vector<size_t> & loop = it->second;
  set<size_t> vertex;
  for(size_t t = 0; t < loop.size(); ++t){
      if(loop[t] == -1) continue;
      for(size_t i = 0; i < 4; ++i){
          if(tet(i,loop[t]) != edge.first && tet(i,loop[t])!= edge.second)
            vertex.insert(tet(i,loop[t]));
        }
    }
  around_vertex.resize(vertex.size());
  copy(vertex.begin(),vertex.end(),around_vertex.begin());

  return 0;
}

static int add_face(
    const size_t v0, const size_t v1, const size_t v2,
    const matrixst & tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    vector<pair<size_t,size_t> > &path_face_with_tets)
{
  if(fa.get_face_idx(v0,v1,v2) == -1){
      cerr << "# [error] three vertex can not make a face." << endl;
      return 2;
    }
  const pair<size_t,size_t> & tet_pair = fa.query(v0,v1,v2);
  const matrixd edge_loop_order = cross(node(colon(),v1) - node(colon(),v0),
                                        node(colon(),v2) - node(colon(),v1));
  if(fa.is_outside_face(tet_pair)){
      cerr << "# [info] this face is outside surface." << endl;
      return 3;
    }

  if(find(path_face_with_tets.begin(),path_face_with_tets.end(),
          tet_pair) != path_face_with_tets.end()
     || find(path_face_with_tets.begin(),path_face_with_tets.end(),
             make_pair(tet_pair.second,tet_pair.first)) != path_face_with_tets.end()){
      return 4;
    }

  const matrixd face_center =
      (node(colon(),v0) + node(colon(),v1) + node(colon(),v2))/3.0;

  const size_t other_vertex[2] =
  {accumulate(tet(colon(),tet_pair.first).begin(),tet(colon(),tet_pair.first).end(),0)
   - (v0+v1+v2),
   accumulate(tet(colon(),tet_pair.second).begin(),tet(colon(),tet_pair.second).end(),0)
   - (v0+v1+v2),};

  const matrixd dir = face_center - node(colon(),other_vertex[0]);
  if(dot(dir,edge_loop_order) > 0)
    path_face_with_tets.push_back(tet_pair);
  else
    path_face_with_tets.push_back(make_pair(tet_pair.second,tet_pair.first));

  return 0;
}

static int add_face_by_modify_type(
    const size_t v0, const size_t v1, const size_t v2,
    const matrixst & tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type,
    vector<pair<size_t,size_t> > &path_visited_faces)
{
  if(fa.get_face_idx(v0,v1,v2) == -1){
      cerr << "# [error] three vertex can not make a face." << endl;
      return 2;
    }

  const pair<size_t,size_t> & tet_pair = fa.query(v0,v1,v2);
  const matrixd edge_loop_order = cross(node(colon(),v1) - node(colon(),v0),
                                        node(colon(),v2) - node(colon(),v1));
  if(fa.is_outside_face(tet_pair)){
      cerr << "# [info] this face is outside surface." << endl;
      return 3;
    }

  const matrixd face_center =
      (node(colon(),v0) + node(colon(),v1) + node(colon(),v2))/3.0;

  const size_t other_vertex[2] =
  {accumulate(tet(colon(),tet_pair.first).begin(),tet(colon(),tet_pair.first).end(),0)
   - (v0+v1+v2),
   accumulate(tet(colon(),tet_pair.second).begin(),tet(colon(),tet_pair.second).end(),0)
   - (v0+v1+v2),};

  const matrixd dir = face_center - node(colon(),other_vertex[0]);
  pair<size_t,size_t> order_tet;
  if(dot(dir,edge_loop_order) > 0)
    order_tet = tet_pair;
  else{
      order_tet.first = tet_pair.second;
      order_tet.second = tet_pair.first;
    }

  if(find(path_visited_faces.begin(),path_visited_faces.end(),order_tet)
     != path_visited_faces.end()){
      return 4;
    }else
    path_visited_faces.push_back(order_tet);

  deque<pair<size_t,size_t> > chain;
  if(ea.get_edge_idx(v0,v1) == -1)
    chain.push_back(make_pair(v0,v1));
  else if(ea.get_edge_idx(v1,v2) == -1)
    chain.push_back(make_pair(v1,v2));
  else if(ea.get_edge_idx(v2,v0) == -1)
    chain.push_back(make_pair(v2,v0));
  else{
      return 5;
      cerr << "# [error] three edges are all on surface." << endl;
    }
  vector<pair<size_t,size_t> > need_to_apply_jump_tet_pair;
  need_to_apply_jump_tet_pair.push_back(order_tet);
  apply_jump_on_cut_faces_for_each_group(tet,need_to_apply_jump_tet_pair,
                                         inner_face_jump_type,chain,fa,ortae);
  return 0;
}

int modify_face_type_to_remove_near_surface_with_distance(
    const matrixst & tet,
    const matrixd & node,
    const matrixd & dist,
    const matrixst & prev_node,
    const jtf::mesh::edge2cell_adjacent & ea,
    const deque<pair<size_t,size_t> > & chain){

  // preprocess
  deque<list<pair<size_t,size_t> > > split_chains;
  {
    list<pair<size_t,size_t> > chain_list;
    for(deque<pair<size_t,size_t> >::const_iterator it = chain.begin();
        it != chain.end(); ++it){
        chain_list.push_back(*it);
      }
    split_chains.push_back(chain_list);
  }

  typedef list<pair<size_t,size_t> >::iterator lit;
  while(!split_chains.empty()){
      list<pair<size_t,size_t> > & chain_segment = split_chains.front();

      for(lit it = chain_segment.begin(); it != chain_segment.end();){
          if(ea.get_edge_idx(it->first,it->second) != -1) // surface edge
            {
              list<pair<size_t,size_t> > one_segment;
              one_segment.insert(one_segment.end(),chain_segment.begin(),it);
              if(!one_segment.empty()) split_chains.push_back(one_segment);
              chain_segment.erase(it++);
            }else
            ++it;
        }
    }
  return 0;
}

//! @brief for each near surface singularity chain, get a surface path which can link
// inner singularity chain into a loop, then find a patch whose boundary is such a loop
// modify face jump type of each face in such patch to remove the inner singularity
// WARNING: This function is not robust, for it used geometry direction to indicate which
// face should be handled. If singularity twist a lot, this process possibly failed.
// If can not get such a patch, return non-zero
int modify_face_type_to_remove_near_surface(
    const matrixst &tet,
    const matrixd &node,
    const matrixst outside_face,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::face2tet_adjacent &fa,
    const vertex_connection<UNDIRECT> &vc,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const deque<pair<size_t,size_t> > &chain,
    boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type)
{
#define visualization
#ifdef visualization
  static size_t iter = 0;
  ostringstream oos;
  oos << iter++;
  const string iter_str = oos.str();
#endif
  assert(find(outside_face.begin(),outside_face.end(),chain.front().first) != outside_face.end());
  assert(find(outside_face.begin(),outside_face.end(),chain.back().second) != outside_face.end());


  vector<size_t> cut_faces;
  vector<pair<size_t,size_t> > path_visited_faces; // store the adjacent tets idx
  list<pair<size_t,size_t> > near_surface_singularity_list_;
  for(size_t t = 0; t < chain.size(); ++t)
    near_surface_singularity_list_.push_back(chain[t]);

  vector<size_t> around_vertex;

  deque<list<pair<size_t,size_t> > > nssl_vec;
  nssl_vec.push_back(near_surface_singularity_list_);

  // during the calculation, this singularity will be splitted
  while(!nssl_vec.empty()){
      list<pair<size_t,size_t> > & near_surface_singularity_list = nssl_vec.front();
      if(near_surface_singularity_list.empty()) {
          nssl_vec.pop_front();
          continue;
        }
      if(near_surface_singularity_list.front().first ==
         near_surface_singularity_list.back().second)
        {
          cerr << "# [error] generate an inner suface loop." << endl;
          cerr << "# [errro] Sorry, I can not handle this singularity currently." << endl;
          return __LINE__;
        }
      ////////////////////////////////////////////////////////////
      ///// calculate the surface path
      vector<size_t> path;
      vc.get_shortest_path(near_surface_singularity_list.front().first,
                           near_surface_singularity_list.back().second,path);
      assert(path.front() == near_surface_singularity_list.front().first &&
             path.back() == near_surface_singularity_list.back().second);

      deque<size_t> variable_path(path.size());
      copy(path.begin(),path.end(),variable_path.begin());

      matrixd sign_of_loop =
          cross(node(colon(),path[1]) - node(colon(),path[0]),
          node(colon(),chain.front().second) - node(colon(),chain.front().first));
      ////////////////////////////////////////////////////////////////

      while(!near_surface_singularity_list.empty()){
          list<pair<size_t,size_t> > temp_near_surface_singularity_list
              = near_surface_singularity_list;

          cerr << "# [info] ***** left near surface singularity edge num: "
               << near_surface_singularity_list.size() << endl;
          //if inner singularity is only one edge, we should shrink the surface edge and add inner edges

          assert(near_surface_singularity_list.front().first == variable_path.front() &&
                 near_surface_singularity_list.back().second == variable_path.back());

          // update the loop direction
          if(variable_path.size() > 1){
              sign_of_loop =
                  cross(node(colon(),variable_path[1]) - node(colon(),variable_path[0]),
                  node(colon(),near_surface_singularity_list.front().second)
                  - node(colon(),near_surface_singularity_list.front().first));
            }

          bool is_chain_splitted = false;
          // need to split near surface singularity if it reach surface
          for(list<pair<size_t,size_t> >::iterator lit = near_surface_singularity_list.begin();
              lit != near_surface_singularity_list.end();){
              if(ea.get_edge_idx(lit->first,lit->second) != -1){// surface edge
                  list<pair<size_t,size_t> > split_list;
                  for(list<pair<size_t,size_t> >::const_iterator it0 = near_surface_singularity_list.begin();
                      it0 != lit;  ++it0){
                      split_list.push_back(*it0);
                    }
                  if(!split_list.empty())
                    nssl_vec.push_back(split_list);
                  near_surface_singularity_list.erase(near_surface_singularity_list.begin(),
                                                      ++lit);
                  if(near_surface_singularity_list.empty()) nssl_vec.pop_front();
                  is_chain_splitted = true;
                }else  ++lit;
            }

          if(is_chain_splitted) break;

          if(near_surface_singularity_list.size() == 1){
              if(fa.get_face_idx(variable_path[2],variable_path[1],variable_path[0])
                 != -1){// these three vertex are of one face
                  int rtn = add_face_by_modify_type(variable_path[1],variable_path[0],variable_path[2],tet,node,
                      fa,ea,ortae,inner_face_jump_type,path_visited_faces);
                  if(rtn == 3) { // this face is outside face  or visited
                      // need to remove variable_path[1]
                      const size_t temp = variable_path.front();
                      variable_path.pop_front();
                      variable_path.pop_front();
                      variable_path.push_front(temp);
                      continue; // this face is not valid
                    }
#ifdef visualization
                  cut_faces.push_back(variable_path[1]);
                  cut_faces.push_back(variable_path[0]);
                  cut_faces.push_back(variable_path[2]);
#endif
                  near_surface_singularity_list.insert(
                        near_surface_singularity_list.begin(),
                        make_pair(variable_path[2], variable_path[0]));

                  variable_path.pop_front();
                  variable_path.pop_front();
                  if(variable_path.size() > 2)
                    sign_of_loop =
                        cross(node(colon(),variable_path[1]) - node(colon(),variable_path[0]),
                        node(colon(),chain.front().second) - node(colon(),chain.front().first));
                }else{ // need to find a vertex around <v0,v1> which is close to v2
                  assert(variable_path.size() > 1);
                  get_around_vertex(make_pair(near_surface_singularity_list.front().first,
                                              near_surface_singularity_list.front().second),tet,ortae,
                                    around_vertex);
                  //        get_around_vertex(make_pair(variable_path[0],variable_path[1]),
                  //                          tet,ortae,around_vertex);
                  //        double min_length = std::numeric_limits<double>::max();
                  //        size_t idx_ = -1;
                  vector<pair<double,size_t> > length_idx;
                  for(size_t t = 0; t < around_vertex.size(); ++t){
                      const double length_ = norm(node(colon(),around_vertex[t])
                                                  - node(colon(),variable_path[1]))
                          + norm(node(colon(),around_vertex[t])
                                 - node(colon(),variable_path[0]))
                          + norm(node(colon(),around_vertex[t])
                                 - node(colon(),near_surface_singularity_list.front().second));
                      length_idx.push_back(make_pair(length_,t));
                      //          if(length_ < min_length){
                      //            min_length = length_;
                      //            idx_ = t;
                      //          }
                    }// end for
                  size_t idx_ = -1;
                  //if(idx_ == -1) cerr << "# [error] can not find a minimal length." << endl;
                  sort(length_idx.begin(),length_idx.end());
                  size_t vertex_idx = 0;
                  for(vertex_idx = 0; vertex_idx < length_idx.size(); ++vertex_idx)
                    {
                      if(add_face_by_modify_type(near_surface_singularity_list.front().first,
                                                 near_surface_singularity_list.front().second,
                                                 around_vertex[length_idx[vertex_idx].second],
                                                 tet,node,fa,ea,ortae,inner_face_jump_type,path_visited_faces)){
                          // this face is already inserted
                          cerr << "# [info] this face has already been visited."
                               << " face: " << near_surface_singularity_list.front().first << " "
                               << near_surface_singularity_list.front().second << " "
                               << around_vertex[length_idx[vertex_idx].second] << endl;
                          continue;
                        }
#ifdef visualization
                      cut_faces.push_back(near_surface_singularity_list.front().first);
                      cut_faces.push_back(near_surface_singularity_list.front().second);
                      cut_faces.push_back(around_vertex[length_idx[vertex_idx].second]);
#endif
                      break;
                    }

                  if(vertex_idx == length_idx.size()){
                      cerr << "# [error] strange can not find a valid face." << endl;
                      cerr << "# [error] Sorry. I can not handle this singularity." << endl;
                      return __LINE__;
                    }else idx_ = length_idx[vertex_idx].second;

                  const pair<size_t,size_t> edge_need_to_insert_prev(near_surface_singularity_list.front().first,
                                                                     around_vertex[idx_]);
                  const pair<size_t,size_t> edge_need_to_insert_next(around_vertex[idx_],
                                                                     near_surface_singularity_list.front().second);

                  near_surface_singularity_list.clear();
                  near_surface_singularity_list.insert(near_surface_singularity_list.begin(),
                                                       edge_need_to_insert_prev);
                  near_surface_singularity_list.insert(near_surface_singularity_list.end(),
                                                       edge_need_to_insert_next);

                  //        variable_path.pop_front();

                  //        list<pair<size_t,size_t> >::iterator lit= near_surface_singularity_list.begin();
                  //        near_surface_singularity_list.insert(lit,edge_need_to_insert_prev);
                  //        near_surface_singularity_list.insert(lit,edge_need_to_insert_next);
                  continue;
                }
            }
          // inner singularity edge length > 1, we can shrink go along it
          for(list<pair<size_t,size_t> >::iterator lit = near_surface_singularity_list.begin();
              lit != near_surface_singularity_list.end();){

              if(lit != near_surface_singularity_list.begin()){
                  list<pair<size_t,size_t> >::iterator prev = --lit;
                  ++lit;
                  assert(prev->second == lit->first);
                  if(prev->first == lit->second){ // if this two edge has already degenerate to a point, remove them
                      near_surface_singularity_list.erase(prev);
                      near_surface_singularity_list.erase(lit++);
                      continue;
                    }
                  const double path_order =
                      dot(sign_of_loop,cross(
                            node(colon(),prev->first) - node(colon(),prev->second),
                            node(colon(),lit->second) - node(colon(),lit->first)));
                  if(path_order < 0) {
                      ++lit;
                      continue; // wrong order
                    }

                  if(fa.get_face_idx(prev->first,prev->second,lit->second) !=-1){ // is in one face
                      const pair<size_t,size_t> edge_need_to_insert(prev->first,lit->second);
                      int rtn = add_face_by_modify_type(prev->first,prev->second,lit->second,
                                                        tet,node,fa,ea,ortae,inner_face_jump_type,path_visited_faces);
                      if(rtn == 2){
                          return __LINE__; // three edges can not make a face
                        }
                      if(rtn != 0 ){ // TODO: this face is already outside face, needs to split it
                          ++lit;
                          continue;
                        }
#ifdef visualization
                      cut_faces.push_back(prev->first);
                      cut_faces.push_back(prev->second);
                      cut_faces.push_back(lit->second);
#endif
                      near_surface_singularity_list.erase(prev);
                      near_surface_singularity_list.erase(lit++);
                      //if(ea.get_edge_idx(edge_need_to_insert.first,edge_need_to_insert.second) == -1) //if this edge is not on surface
                      near_surface_singularity_list.insert(lit,edge_need_to_insert);
                      // else{// TODO: this needs to split variable path or other things
                      //  if(edge_need_to_insert.first == variable_path)
                      //  variable_path.push_back(edge_need_to_insert.first);
                      //   break;
                      // }
                      continue;
                    }else{ // two edge is not in one face, find a close face to pre vertex
                      // inner near surface singularity chain become a calabash shape
                      if((near_surface_singularity_list.front().first ==
                          near_surface_singularity_list.back().second) &&
                         (near_surface_singularity_list.front().second ==
                          near_surface_singularity_list.back().first) )
                        {
                          if(find(outside_face.begin(),
                                  outside_face.end(),
                                  near_surface_singularity_list.front().second)
                             != outside_face.end()){
                              near_surface_singularity_list.pop_front();
                              near_surface_singularity_list.pop_back();
                              continue;
                            }else{
                              cerr << "# [error] near surface becomes an inner loop." << endl;
                              cerr << "# [error] Sorry. I can not handle this singularity." << endl;
                              return __LINE__;
                            }
                        }


                      get_around_vertex(*lit,tet,ortae,around_vertex);
                      double min_length = std::numeric_limits<double>::max();
                      size_t idx_ = -1;
                      // to find a vertex which is close to these vertex
                      for(size_t t = 0; t < around_vertex.size(); ++t){
                          const double length_ = norm(node(colon(),around_vertex[t])
                                                      - node(colon(),prev->first))
                              + norm(node(colon(),around_vertex[t])
                                     - node(colon(),prev->second))
                              + norm(node(colon(),around_vertex[t])
                                     -node(colon(),lit->second));
                          if(length_ < min_length){
                              min_length = length_;
                              idx_ = t;
                            }
                        }// end for
                      if(idx_ == -1) cerr << "# [error] can not find a minimal length." << endl;
                      int rtn = add_face_by_modify_type(lit->first,lit->second,around_vertex[idx_],tet,node,fa,ea,ortae,
                                                        inner_face_jump_type,path_visited_faces);
                      if(rtn == 2){
                          return __LINE__;
                        }
                      if(rtn != 0){ //TODO: this face is already outside face, needs to split it
                          ++lit;
                          continue;
                        }
#ifdef visualization
                      cut_faces.push_back(lit->first);
                      cut_faces.push_back(lit->second);
                      cut_faces.push_back(around_vertex[idx_]);
#endif
                      const pair<size_t,size_t> edge_need_to_insert_prev(lit->first,around_vertex[idx_]);
                      const pair<size_t,size_t> edge_need_to_insert_next(around_vertex[idx_],lit->second);
                      near_surface_singularity_list.erase(lit++);
                      near_surface_singularity_list.insert(lit,edge_need_to_insert_prev);
                      near_surface_singularity_list.insert(lit,edge_need_to_insert_next);
                      continue;
                    }
                }// end if current edge is not on the edge of singularity
              else{
                  assert(lit == near_surface_singularity_list.begin());
                  //        get_around_vertex(*lit,tet,ortae,around_vertex);
                  // this face is adjacent surface
                  if(variable_path.size() > 1){
                      if(fa.get_face_idx(variable_path[1],lit->first,lit->second) != -1)
                        {
                          const pair<size_t,size_t> edge_need_to_insert(variable_path[1],lit->second);
                          if(add_face_by_modify_type(lit->first,lit->second,variable_path[1],
                                                     tet,node,fa,ea,ortae,inner_face_jump_type,
                                                     path_visited_faces))
                            continue;
#ifdef visualization
                          cut_faces.push_back(lit->first);
                          cut_faces.push_back(lit->second);
                          cut_faces.push_back(variable_path[1]);
#endif
                          const size_t check_vertex = lit->second;
                          variable_path.pop_front();
                          near_surface_singularity_list.erase(lit++);
                          if(ea.get_edge_idx(edge_need_to_insert.first,edge_need_to_insert.second) == -1){
                              near_surface_singularity_list.insert(lit,edge_need_to_insert);
                            }
                          else{
                              if(find(variable_path.begin(),variable_path.end(),check_vertex)
                                 != variable_path.end())
                                {
                                  while(variable_path.front() != check_vertex) variable_path.pop_front();
                                }else {
                                  variable_path.push_front(check_vertex);
                                }
                            }
                          assert(variable_path[0] == near_surface_singularity_list.front().first);
                          if(variable_path.size() > 2){ // update sign
                              sign_of_loop =
                                  cross(node(colon(),variable_path[1]) - node(colon(),variable_path[0]),
                                  node(colon(),chain.front().second) - node(colon(),chain.front().first));
                            }
                          continue;
                        }else ++lit;
                    }else ++lit;
                } // end else
            }
          // if this singulariy has been degenerate to a point
          if(near_surface_singularity_list.size() == 2 &&
             near_surface_singularity_list.front().first ==
             near_surface_singularity_list.back().second){
              near_surface_singularity_list.clear();
              continue;
            }

          if(near_surface_singularity_list.empty()) continue;

          if(temp_near_surface_singularity_list
             == near_surface_singularity_list){
              cerr << "# [error] Sorry. I can not handle this singularity." << endl;
              return __LINE__;
            }
        }
    }
#ifdef visualization
  ofstream os_cut_faces((iter_str + "_cut_faces.vtk").c_str());
  tri2vtk(os_cut_faces,&node[0],node.size(2),&cut_faces[0],cut_faces.size()/3);
#endif

  return 0;
}

//! @brief find a face patch which connect inner near surface singularity
// and surface edges
int find_cut_faces_along_singularity(const matrixst &tet,
                                     const matrixd &node,
                                     const jtf::mesh::edge2cell_adjacent &ea,
                                     const jtf::mesh::face2tet_adjacent &fa,
                                     const vertex_connection<UNDIRECT> & vc,
                                     const matrixst &outside_face,
                                     vector<vector<pair<size_t,size_t> > > &need_to_apply_jump_tet_pair,
                                     const deque<pair<size_t,size_t> > &chain,
                                     const jtf::mesh::one_ring_tet_at_edge &ortae)
{
#define visualization
#ifdef visualization
  static size_t iter = 0;
  ostringstream oos;
  oos << iter++;
  const string iter_str = oos.str();
#endif
  assert(find(outside_face.begin(),outside_face.end(),chain.front().first) != outside_face.end());
  assert(find(outside_face.begin(),outside_face.end(),chain.back().second) != outside_face.end());



  //#ifdef visualization
  //  ofstream ofs_path((iter_str + "_surface_path.vtk").c_str());
  //  vector<size_t> surface_edge;
  //  for(size_t t = 0; t < path.size()-1; ++t) {
  //    surface_edge.push_back(path[t]);
  //    surface_edge.push_back(path[t+1]);
  //  }
  //  line2vtk(ofs_path,&node[0],node.size(2),&surface_edge[0],surface_edge.size()/2);
  //#endif

  list<pair<size_t,size_t> > near_surface_singularity_list;
  for(size_t t = 0; t < chain.size(); ++t)
    near_surface_singularity_list.push_back(chain[t]);

  vector<pair<size_t,size_t> > path_face_with_tets; // pair first is one group while second is the other
  vector<size_t> around_vertex;

  //  deque<list<pair<size_t,size_t> > > nssl_vec;
  //  nssl_vec.push_back(near_surface_singularity_list_);

  // during the calculation, this singularity will be splitted
  //  while(!nssl_vec.empty()){
  //    list<pair<size_t,size_t> > & near_surface_singularity_list = nssl_vec.front();

  ////////////////////////////////////////////////////////////
  ///// calculate the surface path
  vector<size_t> path;
  vc.get_shortest_path(chain.front().first,chain.back().second,path);
  assert(path.front() == chain.front().first && path.back() == chain.back().second);

  deque<size_t> variable_path(path.size());
  copy(path.begin(),path.end(),variable_path.begin());

  matrixd sign_of_loop =
      cross(node(colon(),path[1]) - node(colon(),path[0]),
      node(colon(),chain.front().second) - node(colon(),chain.front().first));
  ////////////////////////////////////////////////////////////////

  while(!near_surface_singularity_list.empty()){
      list<pair<size_t,size_t> > temp_near_surface_singularity_list
          = near_surface_singularity_list;

      cerr << "# [info] ***** left near surface singularity edge num: "
           << near_surface_singularity_list.size() << endl;
      //if inner singularity is only one edge, we should shrink the surface edge and add inner edges

      assert(near_surface_singularity_list.front().first == variable_path.front() &&
             near_surface_singularity_list.back().second == variable_path.back());

      //      bool is_chain_splitted = false;
      for(list<pair<size_t,size_t> >::iterator lit = near_surface_singularity_list.begin();
          lit != near_surface_singularity_list.end();){
          if(ea.get_edge_idx(lit->first,lit->second) != -1){// surface edge
              if(lit == near_surface_singularity_list.begin()){
                  variable_path.pop_front();
                  near_surface_singularity_list.erase(lit++);
                  if(near_surface_singularity_list.empty()) break;
                  if(variable_path.front() != near_surface_singularity_list.front().first){
                      // if path first is not the same with near surface singularity first,
                      // need to find a path to link them
                      vector<size_t> relink_path;
                      vc.get_shortest_path(variable_path.front(),
                                           near_surface_singularity_list.front().first,
                                           relink_path);
                      assert(relink_path.size() > 1);
                      for(size_t t = 1; t < relink_path.size(); ++t)
                        variable_path.push_front(relink_path[t]);
                    }
                }else if(*lit == near_surface_singularity_list.back()){
                  if(lit->first != variable_path[variable_path.size()-2]) // this surface edge is not contained in variable path
                    variable_path.push_back(lit->first);
                  else
                    variable_path.pop_back();
                  near_surface_singularity_list.erase(lit++);
                } else{
                  //            list<pair<size_t,size_t> > split_list;
                  //            for(list<pair<size_t,size_t> >::const_iterator it0 = near_surface_singularity_list.begin();
                  //                it0 != lit;  ++it0){
                  //              split_list.push_back(*it0);
                  //            }
                  //            nssl_vec.push_back(split_list);
                  //            near_surface_singularity_list.erase(near_surface_singularity_list.begin(),
                  //                                                ++lit);
                  //            is_chain_splitted = true;
                  cerr << "# [error] singularity chain touch surface too early, needs to split." << endl;
                  cerr << "# [error] Sorry, I can not handle this currently." << endl;
                  return __LINE__;
                  //            cerr << "# [error] strange shrink of near surface singularity." << endl;
                  //            cerr << "# [info] Sorry, I can not handle this singularity." << endl;
                  //            return __LINE__;
                }
            }else  ++lit;
        }

      //      if(is_chain_splitted) break;

      if(near_surface_singularity_list.size() == 1){
          if(fa.get_face_idx(variable_path[2],variable_path[1],variable_path[0])
             != -1){// these three vertex are of one face
              int rtn = add_face(variable_path[1],variable_path[0],variable_path[2],tet,node,
                  fa,ortae,path_face_with_tets);
              if(rtn == 3) { // this face is outside face  or visited
                  // need to remove variable_path[1]
                  const size_t temp = variable_path.front();
                  variable_path.pop_front();
                  variable_path.pop_front();
                  variable_path.push_front(temp);
                  continue; // this face is not valid
                }

              near_surface_singularity_list.insert(
                    near_surface_singularity_list.begin(),
                    make_pair(variable_path[2], variable_path[0]));

              variable_path.pop_front();
              variable_path.pop_front();
              if(variable_path.size() > 2)
                sign_of_loop =
                    cross(node(colon(),variable_path[1]) - node(colon(),variable_path[0]),
                    node(colon(),chain.front().second) - node(colon(),chain.front().first));
            }else{ // need to find a vertex around <v0,v1> which is close to v2
              assert(variable_path.size() > 1);
              get_around_vertex(make_pair(near_surface_singularity_list.front().first,
                                          near_surface_singularity_list.front().second),tet,ortae,
                                around_vertex);
              //        get_around_vertex(make_pair(variable_path[0],variable_path[1]),
              //                          tet,ortae,around_vertex);
              //        double min_length = std::numeric_limits<double>::max();
              //        size_t idx_ = -1;
              vector<pair<double,size_t> > length_idx;
              for(size_t t = 0; t < around_vertex.size(); ++t){
                  const double length_ = norm(node(colon(),around_vertex[t])
                                              - node(colon(),variable_path[1]));
                  length_idx.push_back(make_pair(length_,t));
                  //          if(length_ < min_length){
                  //            min_length = length_;
                  //            idx_ = t;
                  //          }
                }// end for
              size_t idx_ = -1;
              //if(idx_ == -1) cerr << "# [error] can not find a minimal length." << endl;
              sort(length_idx.begin(),length_idx.end());
              size_t vertex_idx = 0;
              for(vertex_idx = 0; vertex_idx < length_idx.size(); ++vertex_idx)
                {
                  if(add_face(near_surface_singularity_list.front().first,
                              near_surface_singularity_list.front().second,
                              around_vertex[length_idx[vertex_idx].second],
                              tet,node,fa,ortae,path_face_with_tets)){
                      // this face is already inserted
                      cerr << "# [info] this face has already been visited."
                           << " face: " << near_surface_singularity_list.front().first << " "
                           << near_surface_singularity_list.front().second << " "
                           << around_vertex[length_idx[vertex_idx].second] << endl;
                      continue;
                    }

                  break;
                }

              if(vertex_idx == length_idx.size()){
                  cerr << "# [error] strange can not find a valid face." << endl;
                  cerr << "# [error] Sorry. I can not handle this singularity." << endl;
                  return __LINE__;
                }else idx_ = length_idx[vertex_idx].second;

              const pair<size_t,size_t> edge_need_to_insert_prev(near_surface_singularity_list.front().first,
                                                                 around_vertex[idx_]);
              const pair<size_t,size_t> edge_need_to_insert_next(around_vertex[idx_],
                                                                 near_surface_singularity_list.front().second);

              near_surface_singularity_list.clear();
              near_surface_singularity_list.insert(near_surface_singularity_list.begin(),
                                                   edge_need_to_insert_prev);
              near_surface_singularity_list.insert(near_surface_singularity_list.end(),
                                                   edge_need_to_insert_next);

              //        variable_path.pop_front();

              //        list<pair<size_t,size_t> >::iterator lit= near_surface_singularity_list.begin();
              //        near_surface_singularity_list.insert(lit,edge_need_to_insert_prev);
              //        near_surface_singularity_list.insert(lit,edge_need_to_insert_next);
              continue;
            }
        }
      // inner singularity edge length > 1, we can shrink go along it
      for(list<pair<size_t,size_t> >::iterator lit = near_surface_singularity_list.begin();
          lit != near_surface_singularity_list.end();){

          if(lit != near_surface_singularity_list.begin()){
              list<pair<size_t,size_t> >::iterator prev = --lit;
              ++lit;
              assert(prev->second == lit->first);
              if(prev->first == lit->second){ // if this two edge has already degenerate to a point, remove them
                  near_surface_singularity_list.erase(prev);
                  near_surface_singularity_list.erase(lit++);
                  continue;
                }
              const double path_order =
                  dot(sign_of_loop,cross(
                        node(colon(),prev->first) - node(colon(),prev->second),
                        node(colon(),lit->second) - node(colon(),lit->first)));
              if(path_order < 0) {
                  ++lit;
                  continue; // wrong order
                }

              if(fa.get_face_idx(prev->first,prev->second,lit->second) !=-1){ // is in one face
                  const pair<size_t,size_t> edge_need_to_insert(prev->first,lit->second);
                  int rtn = add_face(prev->first,prev->second,lit->second,tet,node,fa,ortae,path_face_with_tets);
                  if(rtn == 2){
                      return __LINE__; // three edges can not make a face
                    }
                  if(rtn != 0 ){ // TODO: this face is already outside face, needs to split it
                      ++lit;
                      continue;
                    }

                  near_surface_singularity_list.erase(prev);
                  near_surface_singularity_list.erase(lit++);
                  //if(ea.get_edge_idx(edge_need_to_insert.first,edge_need_to_insert.second) == -1) //if this edge is not on surface
                  near_surface_singularity_list.insert(lit,edge_need_to_insert);
                  // else{// TODO: this needs to split variable path or other things
                  //  if(edge_need_to_insert.first == variable_path)
                  //  variable_path.push_back(edge_need_to_insert.first);
                  //   break;
                  // }
                  continue;
                }else{ // two edge is not in one face, find a close face to pre vertex
                  get_around_vertex(*lit,tet,ortae,around_vertex);
                  double min_length = std::numeric_limits<double>::max();
                  size_t idx_ = -1;
                  for(size_t t = 0; t < around_vertex.size(); ++t){
                      const double length_ = norm(node(colon(),around_vertex[t])
                                                  - node(colon(),prev->first));
                      if(length_ < min_length){
                          min_length = length_;
                          idx_ = t;
                        }
                    }// end for
                  if(idx_ == -1) cerr << "# [error] can not find a minimal length." << endl;
                  int rtn = add_face(lit->first,lit->second,around_vertex[idx_],tet,node,fa,ortae,path_face_with_tets);
                  if(rtn == 2){
                      return __LINE__;
                    }
                  if(rtn != 0){ //TODO: this face is already outside face, needs to split it
                      ++lit;
                      continue;
                    }

                  const pair<size_t,size_t> edge_need_to_insert_prev(lit->first,around_vertex[idx_]);
                  const pair<size_t,size_t> edge_need_to_insert_next(around_vertex[idx_],lit->second);
                  near_surface_singularity_list.erase(lit++);
                  near_surface_singularity_list.insert(lit,edge_need_to_insert_prev);
                  near_surface_singularity_list.insert(lit,edge_need_to_insert_next);
                  continue;
                }
            }// end if current edge is not on the edge of singularity
          else{
              assert(lit == near_surface_singularity_list.begin());
              //        get_around_vertex(*lit,tet,ortae,around_vertex);
              // this face is adjacent surface
              if(variable_path.size() > 1){
                  if(fa.get_face_idx(variable_path[1],lit->first,lit->second) != -1)
                    {
                      const pair<size_t,size_t> edge_need_to_insert(variable_path[1],lit->second);
                      if(add_face(lit->first,lit->second,variable_path[1],
                                  tet,node,fa,ortae,path_face_with_tets))
                        continue;
                      const size_t check_vertex = lit->second;
                      variable_path.pop_front();
                      near_surface_singularity_list.erase(lit++);
                      if(ea.get_edge_idx(edge_need_to_insert.first,edge_need_to_insert.second) == -1){
                          near_surface_singularity_list.insert(lit,edge_need_to_insert);
                        }
                      else{
                          if(find(variable_path.begin(),variable_path.end(),check_vertex)
                             != variable_path.end())
                            {
                              while(variable_path.front() != check_vertex) variable_path.pop_front();
                            }else {
                              variable_path.push_front(check_vertex);
                            }
                        }
                      assert(variable_path[0] == near_surface_singularity_list.front().first);
                      if(variable_path.size() > 2){ // update sign
                          sign_of_loop =
                              cross(node(colon(),variable_path[1]) - node(colon(),variable_path[0]),
                              node(colon(),chain.front().second) - node(colon(),chain.front().first));
                        }
                      continue;
                    }else ++lit;
                }else ++lit;
            } // end else
        }
      // if this singulariy has been degenerate to a point
      if(near_surface_singularity_list.size() == 2 &&
         near_surface_singularity_list.front().first ==
         near_surface_singularity_list.back().second){
          near_surface_singularity_list.clear();
          continue;
        }

      if(near_surface_singularity_list.empty()) continue;

      if(temp_near_surface_singularity_list
         == near_surface_singularity_list){

          // last try, since the geometry stander is not accurate
          if(near_surface_singularity_list.size() == 2
             && ea.get_edge_idx(near_surface_singularity_list.front().first,
                                near_surface_singularity_list.back().second)
             != -1) // this singularity is a simple U, and left one face to cover it
            {
              int rtn = add_face(near_surface_singularity_list.front().first,
                                 near_surface_singularity_list.front().second,
                                 near_surface_singularity_list.back().second,
                                 tet,node,fa,ortae,path_face_with_tets);
              near_surface_singularity_list.clear();
            }else{
              cerr << "# [error] Sorry. I can not handle this singularity." << endl;
              return __LINE__;
            }
        }
    }

  need_to_apply_jump_tet_pair.push_back(path_face_with_tets);
  //nssl_vec.pop_front();
#ifdef visualization
  matrixst faces(3,path_face_with_tets.size());
  matrixst left_tet(4,path_face_with_tets.size()),right_tet(4,path_face_with_tets.size());
  for(size_t t = 0; t < path_face_with_tets.size(); ++t){
      common_tet_face(&tet(0,path_face_with_tets[t].first),
                      &tet(0,path_face_with_tets[t].second),
                      &faces(0,t));
      left_tet(colon(),t) = tet(colon(),path_face_with_tets[t].first);
      right_tet(colon(),t) = tet(colon(),path_face_with_tets[t].second);
    }
  ofstream ofs_f((iter_str + "_cut_faces_new.vtk").c_str());
  ofstream ofs_left((iter_str + "_cut_left_tet_new.vtk").c_str());
  ofstream ofs_right((iter_str + "_cut_right_tet_new.vtk").c_str());
  tri2vtk(ofs_f,&node[0],node.size(2),&faces[0],faces.size(2));
  tet2vtk(ofs_left,&node[0],node.size(2),&left_tet[0],left_tet.size(2));
  tet2vtk(ofs_right,&node[0],node.size(2),&right_tet[0],right_tet.size(2));
#endif

  return 0;
}

int split_to_remove_zigzag(
    std::vector<size_t> &tet_array,
    std::vector<double> &node_array,
    const std::vector<std::deque<std::pair<size_t,size_t > > > &chain_list,
    jtf::mesh::one_ring_tet_at_edge &ortae,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    std::vector<matrixd > &frame,
    vector<size_t> & tet_rot_type)
{
  map<pair<size_t,size_t>,size_t> edge_idx_map;
  deque<list<pair<size_t,size_t> > > singularitiy_chain_list; // use list to support dynamic splitting edges
  {// initialize
    for(size_t t = 0; t < chain_list.size(); ++t){
        const deque<pair<size_t,size_t> > & singularity_chain = chain_list[t];
        list<pair<size_t,size_t> > singularity_list;
        for(size_t i = 0;i < singularity_chain.size(); ++i){
            singularity_list.push_back(singularity_chain[i]);
            edge_idx_map[singularity_chain[i]] = t;
          }
        singularitiy_chain_list.push_back(singularity_list);
      }
  }

  typedef list<pair<size_t,size_t> >::iterator lpit;
  typedef boost::unordered_map<pair<size_t,size_t>,vector<pair<size_t,size_t > > >::const_iterator mvpcit;
  boost::unordered_map<pair<size_t,size_t>,vector<pair<size_t,size_t > > > split_edge_map;
  boost::unordered_map<size_t,size_t> split_tet_map;
  while(!singularitiy_chain_list.empty()){
      cerr << "# [info] left " << singularitiy_chain_list.size() << " to be checked for removing zigzag." << endl;
      list<pair<size_t,size_t> > single_chain = singularitiy_chain_list.front();
      singularitiy_chain_list.pop_front();
      lpit lpit_next = single_chain.begin();
      lpit lpit_pre = single_chain.begin();
      while(lpit_next != single_chain.end()){
          lpit_pre = lpit_next;
          mvpcit mvpcit_ = split_edge_map.find(*lpit_pre);
          if(mvpcit_ != split_edge_map.end()){// this edge has beed splitted
              const vector<pair<size_t,size_t> > &split_edges = mvpcit_->second;
              assert(split_edges.size() == 2);
              assert(split_edges.front().second == split_edges.back().first);
              assert(split_edges.front().first == lpit_pre->first
                     && split_edges.back().second == lpit_pre->second);
              lpit lpit_temp = lpit_next++; // lpit_temp == lpit_pre
              *lpit_pre = split_edges.front(); // modify the edge which is need to be splitted to part0 of split edge
              lpit_next = single_chain.insert(lpit_next,split_edges.back());
              lpit_next = ++lpit_temp; // set lpit_next back to lpit_pre+1
            }else{
              mvpcit_ = split_edge_map.find(make_pair(lpit_pre->second,lpit_pre->first));
              if(mvpcit_ != split_edge_map.end()){
                  const vector<pair<size_t,size_t> > &split_edges = mvpcit_->second;
                  assert(split_edges.size() == 2);
                  assert(split_edges.front().second == split_edges.back().first);
                  assert(split_edges.front().first == lpit_pre->second
                         && split_edges.back().second == lpit_pre->first);
                  lpit lpit_temp = lpit_next++; // lpit_temp == lpit_pre
                  *lpit_pre = make_pair(split_edges.back().second,split_edges.back().first); // modify the edge which is need to be splitted to part0 of split edge
                  lpit_next = single_chain.insert(lpit_next,make_pair(split_edges.front().second,split_edges.front().second));
                  lpit_next = ++lpit_temp; // set lpit_next back to lpit_pre+1
                }else{
                  ++lpit_next;
                }
            }

          if(lpit_next == single_chain.end()) break;
          assert(lpit_pre->second == lpit_next->first);

          if(is_zigzag_edge(*lpit_pre,*lpit_next,&tet_array[0],tet_array.size()/4)){
              pair<size_t,size_t> split_edge(lpit_pre->first,lpit_next->second);
              if(split_tets_around_edge(tet_array,node_array,ortae,frame,split_edge,0.5,
                                        split_edge_map,split_tet_map,
                                        inner_face_jump_type))
                cerr << "# [error] shit happens while splitting tets." << endl;

#if 0 // debug
              one_ring_tet_at_edge::e2tet_type::const_iterator oecit = ortae.e2t_.find(split_edge);
              assert(oecit == ortae.e2t_.end());
              oecit = ortae.e2t_.find(make_pair(split_edge.second,split_edge.first));
              assert(oecit == ortae.e2t_.end());
              mvpcit mvpcit_split = split_edge_map.find(split_edge);
              const vector<pair<size_t,size_t> > & split_edges = mvpcit_split->second;
              one_ring_tet_at_edge::e2tet_type::const_iterator oecit0 = ortae.e2t_.find(split_edges.front());
              one_ring_tet_at_edge::e2tet_type::const_iterator oecit1 = ortae.e2t_.find(make_pair(split_edges.front().second,
                                                                                                  split_edges.front().first));
              assert(oecit0 != ortae.e2t_.end() || oecit1 != ortae.e2t_.end());
              oecit0 = ortae.e2t_.find(split_edges.back());
              oecit1 = ortae.e2t_.find(make_pair(split_edges.back().second, split_edges.back().first));
              assert(oecit0 != ortae.e2t_.end() || oecit1 != ortae.e2t_.end());
#endif
            }
        }
    }

  //  assert(tet_rot_type.size() + split_tet_map.size() == tet_array.size()/4);
  //  while(tet_rot_type.size() < tet_array.size()/4){
  //    tet_rot_type.push_back(-1);
  //  }
  //  for(map<size_t,size_t>::const_iterator mcit = split_tet_map.begin();
  //      mcit != split_tet_map.end(); ++mcit){
  //    tet_rot_type[mcit->second] = tet_rot_type[mcit->first];
  //  }
  return 0;
}

static size_t get_potential_center_node_idx(
    map<pair<size_t,size_t>,size_t> & potential_node,
    pair<size_t,size_t>  edge,
    const size_t start_idx)
{
  typedef map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
  if(edge.first > edge.second)
    swap(edge.first,edge.second);

  mpscit mit = potential_node.find(edge);
  if(mit != potential_node.end())
    return mit->second;
  const size_t idx = potential_node.size() + start_idx;
  potential_node.insert(make_pair(edge,idx));
  return idx;
}

static int get_potential_center_node_with_edge(
    const map<pair<size_t,size_t>,size_t> & potential_node,
    const size_t idx,
    pair<size_t,size_t> & edge)
{
  for(map<pair<size_t,size_t>,size_t>::const_iterator mcit = potential_node.begin();
      mcit != potential_node.end(); ++mcit){
      if(mcit->second != idx) continue;
      edge = mcit->first;
      return 0;
    }
  return -1;
}

static int orient_tet_pair_group_along_path(
    const vector<size_t> &tet_array,
    const vector<double> &node_array,
    vector<pair<size_t,size_t> > &tet_pair,
    const vector<size_t> & fan_point,
    const size_t corresponding_point)
{
  //assert(fan_point.size() == tet_pair.size() + 1);
  itr_matrix<const double*> node(3,node_array.size()/3,&node_array[0]);
  itr_matrix<const size_t*> tet(4,tet_array.size()/4,&tet_array[0]);

  matrixd edge0 = zeros<double>(3,1),
      edge1 = zeros<double>(3,1),
      direction = zeros<double>(3,1);
  matrixd into_tet0_direction = zeros<double>(3,1),
      into_tet1_direction = zeros<double>(3,1);
  size_t face[3];

  for(size_t t = 0; t < tet_pair.size(); ++t){
      pair<size_t,size_t> & one_tet_pair = tet_pair[t];
      jtf::mesh::find_common_face(tet(colon(),one_tet_pair.first),tet(colon(),one_tet_pair.second),&face[0]);
      into_tet0_direction =
          (  node(colon(),tet(0,one_tet_pair.first))
             + node(colon(),tet(1,one_tet_pair.first))
             + node(colon(),tet(2,one_tet_pair.first))
             + node(colon(),tet(3,one_tet_pair.first)))/4.0
          - ( node(colon(),face[0])
          + node(colon(),face[1])
          + node(colon(),face[2]))/3.0;
      into_tet0_direction /= norm(into_tet0_direction);
      into_tet1_direction =
          (  node(colon(),tet(0,one_tet_pair.second))
             + node(colon(),tet(1,one_tet_pair.second))
             + node(colon(),tet(2,one_tet_pair.second))
             + node(colon(),tet(3,one_tet_pair.second)))/4.0
          - ( node(colon(),face[0])
          + node(colon(),face[1])
          + node(colon(),face[2]))/3.0;
      into_tet1_direction /= norm(into_tet1_direction);

      edge0 = node(colon(),fan_point[t+1]) - node(colon(),fan_point[t]);
      edge1 = node(colon(),corresponding_point) - node(colon(),fan_point[t+1]);
      direction = cross(edge0,edge1);
      direction /= norm(direction);

      // common face to t1 is more like the direction
      if(dot(direction,into_tet1_direction)
         > dot(direction,into_tet0_direction)){
          swap(one_tet_pair.first,one_tet_pair.second);
        }
    }
  return 0;
}

static size_t calc_compensate_jump_type(
    const map<pair<size_t,size_t>,size_t> & inner_face_jump_type,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    const pair<size_t,size_t> &tet_pair,
    const pair<size_t,size_t> &regular_edge)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  typedef map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
  oecit ocit = ortae.e2t_.find(regular_edge);
  if(ocit == ortae.e2t_.end())
    ocit = ortae.e2t_.find(make_pair(regular_edge.second,regular_edge.first));
  if(ocit == ortae.e2t_.end()){
      cerr << "# [error] can not find regular edge <" << regular_edge.first
           << "," << regular_edge.second << ">."<< endl;
      return -1;
    }

  const vector<size_t> &loop = ocit->second;
  if(find(loop.begin(),loop.end(),-1) != loop.end()){
      cerr << "# [error] this edge is on surface." << endl;
      return -1;
    }

  assert(loop.front() == loop.back());
  deque<size_t> loop_reorder(loop.size()-1);//(loop.begin(),loop.end());
  copy(loop.begin(),loop.end()-1,loop_reorder.begin());

  const size_t first_idx =
      static_cast<size_t>(find(loop_reorder.begin(),
                               loop_reorder.end(),
                               tet_pair.first)
                          - loop_reorder.begin());
  const size_t second_idx =
      static_cast<size_t>(find(loop_reorder.begin(),
                               loop_reorder.end(),
                               tet_pair.second)
                          - loop_reorder.begin());

  if(fabs(first_idx - second_idx) != 1 && fabs(second_idx - first_idx) != 1){
      if(first_idx != 0 && second_idx != 0){
          cerr << "# [error] tet pair should be adjacent in this tet loop." << endl;
          cerr << "# tet pair: " << tet_pair.first << " "
               << tet_pair.second << endl;
          for(size_t t = 0; t < loop.size(); ++t) cerr << loop[t] << " ";
          return -1;
        }
      loop_reorder.push_back(loop_reorder.front());
      loop_reorder.pop_front();
    }else{
      while(loop_reorder.back() != tet_pair.first &&
            loop_reorder.back() != tet_pair.second){
          loop_reorder.push_front(loop_reorder.back());
          loop_reorder.pop_back();
        }
    }

  assert(loop_reorder.back() == tet_pair.first ||
         loop_reorder.back() == tet_pair.second);

  assert(loop_reorder[loop_reorder.size()-2] == tet_pair.first ||
      loop_reorder[loop_reorder.size()-2] == tet_pair.second);

  matrixd rot = eye<double>(3);
  for(size_t t = 0; t < loop_reorder.size(); ++t){
      mpscit mit =
          inner_face_jump_type.find(make_pair(
                                      loop_reorder[t],
                                      loop_reorder[(t+1)%loop_reorder.size()]));
      if(mit == inner_face_jump_type.end())   continue;
      rot = temp(rot * type_transition2(mit->second));
    }
  const size_t edge_type = type_transition1(rot);
  assert(is_regular_type(edge_type)); // regular singularity edge
  pair<size_t,size_t> tet_pair_in_loop(loop_reorder[loop_reorder.size()-2],
      loop_reorder.back());
  if(tet_pair_in_loop == tet_pair)
    return type_transition1(trans(type_transition2(edge_type)));
  else
    return edge_type;
}

static int dump_tets_around_point(const char * filename,
                                  const vector<size_t> & tet_array,
                                  const vector<double> & node_array,
                                  const size_t point_idx)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      cerr << "# [error] can not open file." << endl;
      return __LINE__;
    }

  itr_matrix<const size_t*> tet(4,tet_array.size()/4,&tet_array[0]);
  itr_matrix<const double*> node(3,node_array.size()/3,&node_array[0]);
  vector<size_t> selected_tets;
  for(size_t t = 0; t < tet.size(2); ++t){
      if(find(tet(colon(),t).begin(),tet(colon(),t).end(),point_idx)
         != tet(colon(),t).end()){
          selected_tets.insert(selected_tets.end(),tet(colon(),t).begin(),
                               tet(colon(),t).end());
        }
    }
  tet2vtk(ofs,&node_array[0],node_array.size()/3,&selected_tets[0],selected_tets.size()/4);

  return 0;
}



//! @brief this function is used to create triagle fan associate black_edge.second and regular_edge.first
//                A  regular_edge.first
//               | \
//               |  \
//               B'..B black_edge.first(regular_edge.second)
//                \  |
//                 \ |
//                   C black_edge.second
// we want to remove <A,B>,<B,C> by modify the create triangle fan, and the new singularity will be
// <A.B'> <B',C>
int split_tet_to_remove_black_edge_per_edge(
    vector<size_t> & tet_array,
    vector<double> & node_array,
    vector<matrixd > & frame,
    const pair<size_t,size_t> & black_edge,
    const pair<size_t,size_t> & regular_edge,
    jtf::mesh::one_ring_tet_at_edge & ortae,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type,
    pair<size_t,size_t> & new_regular_edge,
    vector<size_t> & tet_rot_global)
{
  assert(regular_edge.second == black_edge.first);

#if 1 // debug
  static size_t count = 0;
  ostringstream os ;
  os << count++;
#endif

  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  oecit it_black = ortae.e2t_.find(black_edge);
  if(it_black == ortae.e2t_.end())
    it_black = ortae.e2t_.find(make_pair(black_edge.second,black_edge.first));
  if(it_black == ortae.e2t_.end())
    {
      cerr << "# [error] can not find black edge in one_ring_tet_at_edge." << endl;
      return __LINE__;
    }

  oecit it_regular = ortae.e2t_.find(regular_edge);
  if(it_regular == ortae.e2t_.end())
    it_regular = ortae.e2t_.find(make_pair(regular_edge.second,regular_edge.first));
  if(it_regular == ortae.e2t_.end())
    {
      cerr << "# [error] can not find regular  edge in one_ring_tet_at_edge." << endl;
      return __LINE__;
    }

  const vector<size_t> & black_loop = it_black->second;
  const vector<size_t> & regular_loop = it_regular->second;

  if(find(black_loop.begin(),black_loop.end(),-1) != black_loop.end() ||
     find(regular_loop.begin(),regular_loop.end(),-1) != regular_loop.end()){
      cerr << "# [error] black edge or regular edge is surface edge." << endl;
      return __LINE__;
    }

#if 1 // this method will introduce error if black edge and regular edge are on the same tet, split it will be helpful
  itr_matrix<size_t*> tet_(4,tet_array.size()/4,&tet_array[0]);
  {
    boost::unordered_map<pair<size_t,size_t>,vector<pair<size_t,size_t> > > split_edge_map;
    boost::unordered_map<size_t,size_t> split_tet_map;
    unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet_,"topology"));
    if(fa->get_face_idx(black_edge.first,black_edge.second,regular_edge.first)
       != -1){ // exit in one tet
        split_tets_around_edge(tet_array,node_array,ortae,frame,make_pair(regular_edge.first,black_edge.second),0.5,
                               split_edge_map,split_tet_map,inner_face_jump_type);

      }
  }
  itr_matrix<size_t*> tet(4,tet_array.size()/4,&tet_array[0]);
  itr_matrix<double*> node(3,node_array.size()/3,&node_array[0]);
#endif

  map<pair<size_t,size_t>,size_t> potential_node; // record the center node with adjacent two vertex and its index

  //  const size_t regular_start_idx = 0;
  //  const size_t black_end_idx = 1;

  //set<size_t> original_point;
  vector<size_t> original_point;
  {// add original around points
    set<size_t> original_point_set;
    for(size_t t = 0; t < tet.size(2); ++t){
        if(find(tet(colon(),t).begin(),
                tet(colon(),t).end(),
                black_edge.first)
           != tet(colon(),t).end()){
            for(size_t i = 0; i < tet.size(1); ++i){
                if(tet(i,t) != black_edge.first){
                    original_point_set.insert(tet(i,t));
                  }
              }
          }
      }

    original_point.resize(original_point_set.size());
    copy(original_point_set.begin(),original_point_set.end(),original_point.begin());
  }

  const vector<size_t>::const_iterator regular_start_ptr =
      find(original_point.begin(),original_point.end(),regular_edge.first);
  const vector<size_t>::const_iterator black_end_ptr =
      find(original_point.begin(),original_point.end(),black_edge.second);
  if(regular_start_ptr == original_point.end() ||
     black_end_ptr == original_point.end()){
      cerr << "# [error] can not find original_start or black end point." << endl;
      return __LINE__;
    }

  const size_t regular_start_idx = static_cast<size_t>(regular_start_ptr - original_point.begin());
  const size_t black_end_idx = static_cast<size_t>(black_end_ptr - original_point.begin());


  set<size_t> visited_tet;
  map<pair<size_t,size_t>,double> edge_with_weight;
  const size_t origina_point_num = original_point.size();
  {
    // add around black edge center node
    set<size_t> visited_tet;
    for(size_t t = 0; t < black_loop.size() - 1; ++t){
        visited_tet.insert(black_loop[t]);
        vector<size_t> edge_;
        for(size_t i = 0; i < tet.size(1); ++i){
            if(tet(i,black_loop[t]) != black_edge.second &&
               tet(i,black_loop[t]) != black_edge.first){
                edge_.push_back(tet(i,black_loop[t]));
              }
          }
        assert(edge_.size() == 2);
        const size_t new_node_idx = get_potential_center_node_idx(potential_node,
                                                                  make_pair(edge_.front(),edge_.back()),
                                                                  origina_point_num);

        edge_with_weight[make_pair(new_node_idx,black_end_idx)] =
            norm((node(colon(),edge_.front()) + node(colon(),edge_.back()))/2.0
                 - node(colon(),black_edge.second));// the weight can be real distance
      }

    // add around regular edge center node
    for(size_t t = 0; t < regular_loop.size() - 1; ++t){
        visited_tet.insert(regular_loop[t]);
        vector<size_t> edge_;
        for(size_t i = 0; i < tet.size(1); ++i){
            if(tet(i,regular_loop[t]) != regular_edge.second &&
               tet(i,regular_loop[t]) != regular_edge.first){
                edge_.push_back(tet(i,regular_loop[t]));
              }
          }
        assert(edge_.size() == 2);
        const size_t new_node_idx = get_potential_center_node_idx(potential_node,
                                                                  make_pair(edge_.front(),edge_.back()),
                                                                  origina_point_num);

        edge_with_weight[make_pair(new_node_idx,regular_start_idx)] =
            norm((node(colon(),edge_.front()) + node(colon(),edge_.back()))/2.0
                 - node(colon(),regular_edge.first));// the weight can be real distance
      }

    // add around black_edge.first center node
    vector<size_t> one_point_link_to_black_tets;
    for(size_t t = 0; t < tet.size(2); ++t){
        if(find(tet(colon(),t).begin(),tet(colon(),t).end(),black_edge.first)
           != tet(colon(),t).end()){
            if(visited_tet.find(t) == visited_tet.end())
              one_point_link_to_black_tets.push_back(t);
          }
      }

#if 0 // debug
    {
      vector<size_t> adjacnet_tet(visited_tet.size() + one_point_link_to_black_tets.size());
      copy(one_point_link_to_black_tets.begin(),
           one_point_link_to_black_tets.end(),
           adjacnet_tet.begin());
      copy(visited_tet.begin(),visited_tet.end(),adjacnet_tet.begin() + one_point_link_to_black_tets.size());
      ofstream ofs((os.str() + "adjacent_tets.vtk").c_str());
      vector<size_t> adjacent_tet_;//(adjacnet_tet.size()*4);
      for(size_t t = 0; t < adjacnet_tet.size(); ++t){
          adjacent_tet_.push_back(tet_array[4*adjacnet_tet[t]]);
          adjacent_tet_.push_back(tet_array[4*adjacnet_tet[t]+1]);
          adjacent_tet_.push_back(tet_array[4*adjacnet_tet[t]+2]);
          adjacent_tet_.push_back(tet_array[4*adjacnet_tet[t]+3]);
        }
      tet2vtk(ofs,&node_array[0],node_array.size()/3,&adjacent_tet_[0],adjacent_tet_.size()/4);
    }
#endif

    vector<size_t> face(3);
    size_t three_node_idx[3];
    matrixd three_node(3,3);
    for(size_t t = 0; t < one_point_link_to_black_tets.size(); ++t){
        for(size_t j = 0; j < tet.size(1); ++j){
            face[0] = tet(j,one_point_link_to_black_tets[t]);
            face[1] = tet((j+1)%tet.size(1),one_point_link_to_black_tets[t]);
            face[2] = tet((j+2)%tet.size(1),one_point_link_to_black_tets[t]);

            if(find(face.begin(),face.end(),black_edge.first) != face.end()) continue;

            for(size_t k = 0; k < 3; ++k){
                three_node_idx[k] = get_potential_center_node_idx(potential_node,
                                                                  make_pair(face[k],face[(k+1)%3]),
                    origina_point_num);
                three_node(colon(),k) = (node(colon(),face[k]) + node(colon(),face[(k+1)%3]))/2.0;
                const vector<size_t>::const_iterator  anti_point_ptr =
                    find(original_point.begin(),original_point.end(),face[(k+2)%3]);
                if(anti_point_ptr == original_point.end()) {
                    cerr << "# [error] can not find anti point" << endl;
                    return __LINE__;
                  }
                edge_with_weight[make_pair(three_node_idx[k],
                                           static_cast<size_t>(anti_point_ptr-original_point.begin()))] =
                    norm(three_node(colon(),k) - node(colon(),
                                                      static_cast<size_t>(anti_point_ptr-original_point.begin())));
              }

            //        edge_with_weight[make_pair(three_node_idx[0],three_node_idx[1])] =
            //            norm(three_node(colon(),0) - three_node(colon(),1));
            //        edge_with_weight[make_pair(three_node_idx[1],three_node_idx[2])] =
            //            norm(three_node(colon(),1) - three_node(colon(),2));
            //        edge_with_weight[make_pair(three_node_idx[2],three_node_idx[0])] =
            //            norm(three_node(colon(),2) - three_node(colon(),0));
          }
      }
  }

  {
    unique_ptr<vertex_connection<UNDIRECT> > vc(vertex_connection<UNDIRECT>::create(edge_with_weight));

    vector<size_t> path;
    vc->get_shortest_path(regular_start_idx,black_end_idx,path);

    cerr << " path " << endl;
    for(size_t t = 0; t < path.size(); ++t)
      cerr << path[t] << " " ;
    cerr << endl;

    if(path.empty()){
        cerr << "# [error]  can not find a path from regular_edge.first to black_edge.second" << endl;
        return __LINE__;
      }
    pair<size_t,size_t> edge;
    vector<size_t> fan_point;
    boost::unordered_map<pair<size_t,size_t>,vector<pair<size_t,size_t> > > split_edge_map;
    boost::unordered_map<size_t,size_t> split_tet_map;
    //fan_point.push_back(regular_edge.first);
    for(size_t t = 0; t < path.size(); ++t){
        if(path[t] < original_point.size()){
            fan_point.push_back(original_point[path[t]]);
            continue;
          }
        if(get_potential_center_node_with_edge(potential_node,path[t],edge))
          return __LINE__;
        split_tets_around_edge(tet_array,node_array,ortae,frame,edge,0.5,
                               split_edge_map,split_tet_map,
                               inner_face_jump_type);
        fan_point.push_back(split_edge_map[edge].front().second);
        assert(fan_point.back() == node_array.size()/3 -1);
      }

    //    map<pair<size_t,size_t>,vector<pair<size_t,size_t> > > split_edge_map;
    //    map<size_t,size_t> split_tet_map;
    //    for(size_t t = 1; t < path.size() - 1; ++t){
    //      if(get_potential_center_node_with_edge(potential_node,path[t],edge))
    //        return __LINE__;

    //      //split_tets_around_edge(tet_array,node_array,ortae,frame,edge,0.5,split_edge_map,inner_face_jump_type);
    //      split_tets_around_edge(tet_array,node_array,ortae,frame,edge,0.5,
    //                             split_edge_map,split_tet_map,
    //                             inner_face_jump_type);
    //      fan_point.push_back(split_edge_map[edge].front().second);
    //      assert(fan_point.back() == node_array.size()/3 -1);
    //    }
    //    fan_point.push_back(black_edge.second);



    //    while(tet_rot_global.size() != tet_array.size()/4){
    //      tet_rot_global.push_back(-1);
    //    }
    //    for(map<size_t,size_t>::const_iterator mcit = split_tet_map.begin();
    //        mcit != split_tet_map.end(); ++mcit){
    //      tet_rot_global[mcit->second] = tet_rot_global[mcit->first];
    //    }

    //    assert(find(tet_rot_global.begin(),tet_rot_global.end(),-1) == tet_rot_global.end());

#if 1//
    {
      itr_matrix<size_t*>  tet(4,tet_array.size()/4,&tet_array[0]);
#if 0 // debug
      {
        size_t face[3];
        for(map<size_t,size_t>::const_iterator mcit = split_tet_map.begin();
            mcit != split_tet_map.end(); ++mcit){
            find_common_face(tet,mcit->first,mcit->second,&face[0]);
          }
      }
      dump_tets_around_point((os.str() + "tet_around_black_first.vtk").c_str(),tet_array,node_array,black_edge.first);
#endif


      typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
      typedef boost::unordered_map<pair<size_t,size_t>,size_t>::iterator mpsit;

      unique_ptr<jtf::mesh::face2tet_adjacent> fa(
            jtf::mesh::face2tet_adjacent::create(tet,"topology"));
      vector<size_t> new_fan;
      vector<pair<size_t,size_t> > tet_pair;
      for(size_t t = 0; t < fan_point.size()-1; ++t){
          const size_t face_idx = fa->get_face_idx(black_edge.first,fan_point[t],fan_point[t+1]);
          if(face_idx == -1){
              cerr << "# [error] strange it's not a valid face "
                   << black_edge.first << " "
                   << fan_point[t] << " "
                   << fan_point[t+1] << endl;
              return __LINE__;
            }
          tet_pair.push_back(fa->face2tet_[face_idx]);
          new_fan.push_back(black_edge.first);
          new_fan.push_back(fan_point[t]);
          new_fan.push_back(fan_point[t+1]);
        }
#define visual
#ifdef visual
      ofstream ofs((os.str() + "new_fan.vtk").c_str());
      tri2vtk(ofs,&node_array[0],node_array.size()/3,&new_fan[0],new_fan.size()/3);
      ofstream ofs_temp((os.str() + "new_tet.vtk").c_str());
      tet2vtk(ofs_temp,&node_array[0],node_array.size()/3,&tet_array[0],tet_array.size()/4);
#endif

      //      orient_tet_pair_group_along_path(tet_array,node_array,tet_pair,fan_point,
      //                                       black_edge.first);
      ////      if(tet_pair.front().first == 301324 && tet_pair.front().second == 301331)
      ////        cerr << "# pause" << endl;
      //      const size_t comp_jump_type =
      //          calc_compensate_jump_type(inner_face_jump_type,ortae,
      //                                    tet_pair.front(),regular_edge);
      //      if(comp_jump_type == -1) return __LINE__;

      //      for(size_t t = 0; t < tet_pair.size(); ++t){
      //        mpsit mit = inner_face_jump_type.find(tet_pair[t]);
      //        if(mit == inner_face_jump_type.end()){
      //          inner_face_jump_type[tet_pair[t]] = comp_jump_type;
      //          inner_face_jump_type[make_pair(tet_pair[t].second,
      //                                         tet_pair[t].first)]
      //              = type_transition1(trans(type_transition2(comp_jump_type)));
      //        }else{
      //          mit->second = type_transition1(type_transition2(mit->second)
      //                                         * type_transition2(comp_jump_type));
      //          inner_face_jump_type[make_pair(tet_pair[t].second,
      //                                         tet_pair[t].first)]
      //              = type_transition1(trans(type_transition2(mit->second)));
      //        }
      //      }
      assert(fan_point.size() == tet_pair.size() + 1);
      for(size_t t = 0; t < tet_pair.size(); ++t){
          pair<size_t,size_t> edge(fan_point[t],black_edge.first);
          modidy_edge_type_by_change_face_type(edge,tet_pair[t],tet,ortae,inner_face_jump_type);
        }
      new_regular_edge.second = fan_point.back();
      new_regular_edge.first = fan_point[fan_point.size()-2];

#if 1 // check the black edge
      {
        oecit oit = ortae.e2t_.find(black_edge);
        if(oit == ortae.e2t_.end())
          oit = ortae.e2t_.find(make_pair(black_edge.second,
                                          black_edge.first));
        if(oit == ortae.e2t_.end()){
            cerr << "# [error] can not find this black edge. " << endl;
            return __LINE__;
          }

        const vector<size_t> & black_loop = oit->second;
        matrixd rot = eye<double>(3);
        for(size_t t = 0; t < black_loop.size()-1; ++t){
            mpsit mit = inner_face_jump_type.find(make_pair(black_loop[t],
                                                            black_loop[t+1]));
            if(mit == inner_face_jump_type.end()) continue;
            rot = temp(rot * type_transition2(mit->second));
          }
        //        rot = temp(type_transition2(tet_rot_global[black_loop.front()]) * rot);
        //        rot = temp(rot * trans(type_transition2(tet_rot_global[black_loop.front()])));
        if(!is_regular_type(type_transition1(rot))){
            cerr << "# [error] modified black edge' type is still balck "
                 << type_transition1(rot) << endl;
            return __LINE__;
          }
      }
#endif
    }
#endif

  }
  return 0;
}

static bool is_edge_type_compatible_with_black(const size_t t1,
                                               const size_t t2)
{
  return true;
  assert(is_black_line_new(t1) || is_black_line_new(t2));
  size_t regular_type = t1<t2?t1:t2;
  size_t black_type = t1+t2 - regular_type;
  const matrixd black_type_mat = type_transition2(black_type);
  const matrixd rev_regular_type = trans(type_transition2(regular_type));
  if(is_regular_type(type_transition1(black_type_mat * rev_regular_type)) ||
     is_regular_type(type_transition1(rev_regular_type * black_type_mat))){
      return true;
    }
  return false;
}

int reverse_singularity_with_type(deque<pair<size_t,size_t> > & chain,
                                  deque<size_t> & chain_type)
{
  if(chain.size() != chain_type.size()) {
      cerr << "# [error] chain num and chain type is not compatiable." << endl;
      return __LINE__;
    }

  reverse(chain.begin(),chain.end());
  reverse(chain_type.begin(),chain_type.end());
  for(size_t t = 0; t < chain.size(); ++t){
      pair<size_t,size_t> & edge = chain[t];
      swap(edge.first,edge.second);
    }
  return 0;
}

int reverse_singularity(deque<pair<size_t,size_t> > & chain)
{
  reverse(chain.begin(),chain.end());
  for(size_t t = 0; t < chain.size(); ++t){
      pair<size_t,size_t> & edge = chain[t];
      swap(edge.first,edge.second);
    }
  return 0;
}

struct linked_info
{
  size_t black_chain_idx;
  size_t other_end_point_of_black_chain;
  vector<size_t> regular_chain;
};

int relabel_singularity_chain_by_splitting(
    jtf::tet_mesh &tm,
    zjucad::matrix::matrix<matrixd > &frame,
    std::vector<std::deque<std::pair<size_t,size_t> > > &chain_list,
    std::vector<std::deque<size_t> > &singularities_type_,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    boost::property_tree::ptree &pt,
    std::vector<size_t> &tet_rot_type)
{
  vector<size_t> tet_array(tm.tetmesh_.mesh_.size());
  copy(tm.tetmesh_.mesh_.begin(),tm.tetmesh_.mesh_.end(),tet_array.begin());
  vector<double> node_array(tm.tetmesh_.node_.size());
  copy(tm.tetmesh_.node_.begin(),tm.tetmesh_.node_.end(),node_array.begin());
  vector<matrixd > frame_array;
  for(size_t t = 0; t < frame.size(); ++t) {
      frame_array.push_back(frame[t]);
    }

  // remove balck_edge
  map<size_t,linked_info>  point_info;
  typedef map<size_t,linked_info>::const_iterator  mslcit;
  const linked_info init_linked_info={size_t(-1),size_t(-1)};
  for(size_t t = 0; t < chain_list.size(); ++t){
      const deque<pair<size_t,size_t> >& chain = chain_list[t];
      //      if(chain.front().first == 703 ||
      //         chain.back().second == 703){
      //        cerr << "pause" << endl;
      //      }
      mslcit li_front_ptr = point_info.find(chain.front().first);
      if(li_front_ptr == point_info.end()){
          point_info.insert(make_pair(chain.front().first,init_linked_info));
        }
      mslcit li_back_ptr = point_info.find(chain.back().second);
      if(li_back_ptr == point_info.end()){
          point_info.insert(make_pair(chain.back().second,init_linked_info));
        }

      linked_info &li_front = point_info[chain.front().first];
      if(is_black_line_new(singularities_type_[t].front())){
          li_front.black_chain_idx = t;
          li_front.other_end_point_of_black_chain = chain.back().second;
        }else {
          //li_front.black_chain_idx = -1;
          li_front.regular_chain.push_back(t);
        }
      linked_info &li_back = point_info[chain.back().second];
      if(is_black_line_new(singularities_type_[t].back())){
          li_back.black_chain_idx = t;
          li_back.other_end_point_of_black_chain = chain.front().first;
        }else {
          //li_back.black_chain_idx = -1;
          li_back.regular_chain.push_back(t);
        }
    }


  for(map<size_t,linked_info>::iterator mslit = point_info.begin();
      mslit != point_info.end();){
      if(mslit->second.black_chain_idx == -1)
        point_info.erase(mslit++); // remove those point which do not link black edges
      else
        ++mslit;
    }
  // each black edge just needs to handle once
  for(map<size_t,linked_info>::iterator mslit = point_info.begin();
      mslit != point_info.end();){
      const linked_info & li = mslit->second;
      map<size_t,linked_info>::iterator mit = point_info.find(li.other_end_point_of_black_chain);
      if(mit == point_info.end()) {
          ++mslit;
        }else{
          if(li.regular_chain.size() >  mit->second.regular_chain.size()) {
              point_info.erase(mit);
              ++mslit;
            }else{
              point_info.erase(mslit++);
            }
        }
    }

  cerr << "# [info] left " << point_info.size() << "  black chain to handle." << endl;
  for(map<size_t,linked_info>::const_iterator mcit = point_info.begin();
      mcit != point_info.end(); ++mcit){
      const linked_info & li = mcit->second;
      deque<pair<size_t,size_t> > &black_chain = chain_list[li.black_chain_idx];
      deque<size_t> & chain_type = singularities_type_[li.black_chain_idx];

      if(mcit->first == black_chain.back().second)
        reverse_singularity_with_type(black_chain,chain_type);
      assert(mcit->first == black_chain.front().first);
      size_t regular_chain = -1;
      for(size_t t = 0; t < li.regular_chain.size(); ++t){
          if(is_edge_type_compatible_with_black(chain_type.front(),singularities_type_[li.regular_chain[t]].front()))
            {
              regular_chain = li.regular_chain[t];
              break;
            }
        }
      if(regular_chain == -1){
          cerr << "# [error] can not find a chain compatible with black edge." << endl;
          return __LINE__;
        }

      pair<size_t,size_t> regular_edge;
      if(chain_list[regular_chain].front().first == mcit->first){
          regular_edge.first = chain_list[regular_chain].front().second;
          regular_edge.second = chain_list[regular_chain].front().first;
        }else{
          assert(chain_list[regular_chain].back().second == mcit->first);
          regular_edge = chain_list[regular_chain].back();
        }
      pair<size_t,size_t> new_regular_edge;
      for(size_t t = 0; t < black_chain.size(); ++t){
          split_tet_to_remove_black_edge_per_edge(tet_array,node_array,frame_array,black_chain[t],
                                                  regular_edge,tm.ortae_,inner_face_jump_type,new_regular_edge,tet_rot_type);
          regular_edge = new_regular_edge;
        }
      cerr << "# [info] finish one chain." << endl;
      //break;
    }


  itr_matrix<size_t*> tet(4,tet_array.size()/4,&tet_array[0]);
  itr_matrix<double*> node(3,node_array.size()/3,&node_array[0]);
  orient_tet_raw(&node[0],node.size(2),&tet[0], tet.size(2));

  jtf::tet_mesh tm_modify(tet, node);
  singularity_extractor se(tm_modify);
  vector<deque<pair<size_t,size_t> > > singularity_edges;
  std::vector<deque<size_t> > singularities_type;
  se.extract(inner_face_jump_type, singularity_edges, singularities_type);

  dump_singularity_chain_to_vtk_2("singulatiy_after_remove_black.vtk",node,chain_list,singularities_type_);


#if 1
  { // remove zigzag
#if 1 // original
    split_to_remove_zigzag(tet_array,node_array,chain_list,tm_modify.ortae_,
                           inner_face_jump_type,frame_array,tet_rot_type);
    cerr << "# [info] remove zigzag step 2" << endl;
    split_to_remove_zigzag(tet_array,node_array,chain_list,tm_modify.ortae_,
                           inner_face_jump_type,frame_array,tet_rot_type);
#endif
    //    ortae.find_singularities_global_align_with_type(tet_rot_type,inner_face_jump_type,outside_face,chain_list,singularities_type_);
    //    relabel_face_type_to_remove_zigzag(chain_list,singularities_type_,tet_,*fa_new,
    //                                       outside_face,inner_face_jump_type,
    //                                       tet_rot_type,ortae);
  }
#endif

  { // post-process
    cerr << "# [info] post process" << endl;
    tet.resize(4,tet_array.size() / 4);
    node.resize(3,node_array.size() / 3);
    copy(tet_array.begin(),tet_array.end(),tet.begin());
    copy(node_array.begin(),node_array.end(),node.begin());
    orient_tet_raw(&node[0], node.size(2), &tet[0], tet.size(2));
    tm_modify.load(tet, node);

    frame.resize(frame_array.size());
    for(size_t t = 0; t < frame_array.size(); ++t){
        frame[t].resize(3,3);
        frame[t] = frame_array[t];
      }
  }


#if 1 // visuallization
  {
#if 1
    {
      tet_rot_type.resize(tet.size(2));
      calc_tet_rotation(tm_modify,inner_face_jump_type,tet_rot_type);
      itr_matrix<size_t*> tet_rot_type_(1,tet_rot_type.size(), &tet_rot_type[0]) ;
      find_singularities_global_align_with_type(tm_modify.ortae_, tet_rot_type_,inner_face_jump_type,
                                                tm_modify.outside_face_,chain_list,singularities_type_);
      dump_singularity_chain_to_vtk_2("singulatiy_after_near_surface.vtk",node,chain_list,singularities_type_);
    }
#endif
    ofstream ofs("new-tet.vtk");
    tet2vtk(ofs,&node[0],node.size(2),&tet[0],tet.size(2));

    set<vector<size_t> > face_set;
    vector<size_t> shared_face(3);
    for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mcit = inner_face_jump_type.begin();
        mcit != inner_face_jump_type.end(); ++mcit){
        if(is_trivial_type(mcit->second))  continue;
        const pair<size_t,size_t> & tet_pair = mcit->first;

        if(jtf::tetmesh::get_shared_face(tet_pair,tet,&shared_face[0])) {
            cerr << "# [error] can not find shared face for tet "
                 << tet_pair.first << " " << tet_pair.second << " : type " << mcit->second << endl;
            continue;
          }
        sort(shared_face.begin(),shared_face.end());
        face_set.insert(shared_face);
      }
    vector<size_t> face_array;
    face_array.reserve(3 * face_set.size());
    for(set<vector<size_t> >::const_iterator scit = face_set.begin();
        scit != face_set.end(); ++scit){
        const vector<size_t> &face = *scit;
        face_array.insert(face_array.end(),face.begin(),face.end());
      }
    ofstream ofs_tri("inner_face.vtk");
    tri2vtk(ofs_tri,&node[0],node.size(2),&face_array[0],face_set.size());
  }
#endif

  return 0;
}
