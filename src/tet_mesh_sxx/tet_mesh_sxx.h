#ifndef TET_MESH_SXX_H
#define TET_MESH_SXX_H

#include <set>
#include <map>
#include <list>
#include <vector>
#include <iterator>
#include <iostream>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>
#include <boost/tuple/tuple.hpp>
#include <zjucad/matrix/matrix.h>
#include "../tetmesh/tetmesh.h"
//#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/mesh.h>

using namespace std;

namespace sxx{
  template <typename T>
  class my_hash
  {
  public:
    size_t operator()(const typename list<T>::iterator &it) const
    {
      boost::hash<T> t_hash;
      return t_hash(*it);
    }
    size_t operator()(const typename list<T>::const_iterator &it) const
    {
      boost::hash<T> t_hash;
      return t_hash(*it);
    }
  };

  template <typename T>
  class cmp
  {
  public:
    bool operator()(const pair<int, T> &p1, const pair<int, T> &p2) const
    {
      if(p1.first < p2.first)
        return true;
      else
        return false;
    }
  };

  template <typename It>
  int signed_distance(const It &it1, const It &it2, const It &it3)
  {
    It it;
    int d = 0;
    for(it = it1; it != it3; ++it, ++d)
      if(it == it2)
        return d;
    d = 0;
    for(it = it2; it != it1; ++it, --d);
    return d;
  }

  typedef vector<size_t> tet;
  typedef vector<size_t> face;
  typedef pair<size_t, size_t> edge;
  typedef list<tet>::iterator tet_list_it;
  typedef list<edge>::iterator edge_list_it;
  typedef list<face>::iterator face_list_it;
  typedef boost::unordered_map<tet_list_it, list<edge_list_it> ,my_hash<tet> > tet2edge_type;
  typedef boost::unordered_map<edge_list_it, list<tet_list_it> ,my_hash<edge> > edge2tet_type;
  typedef boost::unordered_map<face_list_it, list<tet_list_it> ,my_hash<face> > face2tet_type;
  typedef boost::unordered_map<tet_list_it, list<face_list_it> ,my_hash<tet> > tet2face_type;
  typedef boost::unordered_map<size_t, list<tet_list_it> > vertex2tet_type;
  typedef boost::unordered_map<edge, edge_list_it> edge2it_type;
  typedef boost::unordered_map<face, const face_list_it> face2it_type;
  typedef boost::unordered_map<face, tet_list_it> tet2it_type;
  typedef std::set<pair<int, tet_list_it>, cmp<tet_list_it> > tet_it_set_type;
  typedef std::set<pair<int, edge_list_it>, cmp<edge_list_it> > edge_it_set_type;
  typedef std::set<pair<int, face_list_it>, cmp<face_list_it> > face_it_set_type;
  typedef boost::unordered_map<set<size_t>, size_t> tet2original_type;
  typedef boost::unordered_map<face, face> face2orginal_type;

  inline size_t get_max(const size_t id1, const size_t id2)
  {
    return id2 < id1 ? id1 : id2;
  }

  inline size_t get_min(const size_t id1, const size_t id2)
  {
    return id1 < id2 ? id1 : id2;
  }

  class tet_mesh
  {
  public:
    tet_mesh();
    int create_tetmesh(const char *path);
    int create_tetmesh(const zjucad::matrix::matrix<double> &node,
                       const zjucad::matrix::matrix<size_t> &tet_matrix);
    int clear(){
      tet_list.clear();
      edge_list.clear();
      face_list.clear();
      vertex_vec.clear();
      tet2edge.clear();
      edge2tet.clear();
      face2tet.clear();
      tet2face.clear();
      vertex2tet.clear();
      edge2it.clear();
      face2it.clear();
      tet2it.clear();
      //    current_tet_it = nullptr;
      //    current_edge_it = nullptr;
      //    current_face_it = nullptr;

      tet2original.clear();
      face2orginal.clear();
    }
    int write_tetmesh_to_matrix(zjucad::matrix::matrix<double> &node,
                                zjucad::matrix::matrix<size_t> &tet_matrix,
                                bool remove_extra_node_flag = true) const;
    size_t split_edge(const edge &e, const bool debug_info = false);
    int split_tet(const vector<size_t> &v_id);
    int split_tet(const size_t & p0, const size_t & p1,
                  const size_t & p2, const size_t & p3, const bool debug_info = false);

    int del_tet(const size_t & p0, const size_t & p1,
                const size_t & p2, const size_t & p3,
                const bool debug_info = false);

    int split_tet(const size_t * tet_p, const bool debug_info = false);
    int split_face(const size_t id0, const size_t id1, const size_t id2,
                   const bool debug_info = false);
    int split_face(const size_t * face_p, const bool debug_info = false);
    int collapse_edge(const edge &e, const bool debug_info = true);

    int write_tetmesh_to_file(const char *path1);
    int test_topology_operation();

    int get_tet2orginal_index(const zjucad::matrix::matrix<size_t> &tet_matrix,
                              zjucad::matrix::matrix<size_t> &tet_index_map) const;
    const tet2original_type &get_tet2orginal() const
    {
      return tet2original;
    }

    face2orginal_type &get_face2orginal()
    {
      vector<size_t> null_vector(3);
      face2orginal_type::iterator face2orginal_it, face2orginal_it1;
      for(face2orginal_it = face2orginal.begin(); face2orginal_it != face2orginal.end(); )
        {
          face2orginal_it1 = face2orginal_it++;
          if(face2orginal_it1->second == null_vector)
            face2orginal.erase(face2orginal_it1);
        }
      return face2orginal;
    }


    int output_all_edges(const char * path);

    ~tet_mesh();

  private:
    int tet_split(tet_list_it &tet_it,
                  tet_it_set_type *del_tet_set,
                  edge_it_set_type *del_edge_set,
                  face_it_set_type *del_face_set);
    //  int tet_split_at_edges(const tet_list_it &tet_it,
    //                         const tet_it_set_type *del_tet_set,
    //                         const edge_it_set_type *del_edge_set,
    //                         const face_it_set_type *del_face_set);
    int face_split( face_list_it &face_it,
                    tet_it_set_type *del_tet_set,
                    edge_it_set_type *del_edge_set,
                    face_it_set_type *del_face_set);
    int face_split_at_edge( face_list_it &face_it,
                            tet_it_set_type *del_tet_set,
                            edge_it_set_type *del_edge_set,
                            face_it_set_type *del_face_set);
    int create_tet(const vector<size_t> &tet_veretx);
    int create_tet(const set<size_t> &tet_set);
    int delete_tet(const tet_list_it &tet_it,
                   tet_it_set_type *del_tet_set,
                   edge_it_set_type *del_edge_set,
                   face_it_set_type *del_face_set);
    int add_info_ver2tet(const size_t id, const tet_list_it &it);
    int add_info_edge2tet(const size_t id1, const size_t id2, const tet_list_it &it);
    int add_info_face2tet(const size_t id1, const size_t id2,
                          const size_t id3, const tet_list_it &it);
    int del_info_ver2tet(const size_t id, const tet_list_it &it);
    int del_info_edge2tet(const edge_list_it &edge_it, const tet_list_it &tet_it,
                          edge_it_set_type *del_edge_set);
    int del_info_face2tet(const face_list_it &face_it, const tet_list_it &tet_it,
                          face_it_set_type *del_face_set);
    int sort_face_index(const size_t id1, const size_t id2, const size_t id3,
                        face &f);
    int delete_all_tet();
    int delete_info(const tet_it_set_type *del_tet_set, const edge_it_set_type *del_edge_set,
                    const face_it_set_type *del_face_set);

    int edge_split( edge_list_it &edge_it,
                    tet_it_set_type *del_tet_set,
                    edge_it_set_type *del_edge_set,
                    face_it_set_type *del_face_set,
                    vector<vector<size_t> > *face_adj_tet = NULL);
    int edge_collapse( edge_list_it &edge_it,
                       tet_it_set_type *del_tet_set,
                       edge_it_set_type *del_edge_set,
                       face_it_set_type *del_face_set, const bool is_debug = true);
    int edge_to_face( edge_list_it &edge_it, vector<size_t> &v_id,
                      tet_it_set_type *del_tet_set,
                      edge_it_set_type *del_edge_set,
                      face_it_set_type *del_face_set);
    int face_to_edge( face_list_it &face_it, vector<size_t> &v_id,
                      tet_it_set_type *del_tet_set,
                      edge_it_set_type *del_edge_set,
                      face_it_set_type *del_face_set);
    int edge_to_edge( edge_list_it &edge_it, const vector<edge> &swap_edge,
                      size_t selected_edge_index,
                      tet_it_set_type *del_tet_set,
                      edge_it_set_type *del_edge_set,
                      face_it_set_type *del_face_set);
    edge get_opposite_edge(const edge_list_it &edge_it, const tet_list_it &tet_it);
    int  get_edges_vertex(const vector<edge> &edge_vec, vector<size_t> &v_id);
    size_t get_id_except_face(const face_list_it &face_it, const tet_list_it &tet_it);
    int get_face_except_point(const tet &temp_tet, const size_t id, face &f);

    bool can_face_to_edge(const face_list_it &face_it, vector<size_t> &v_id);
    bool can_edge_to_face(const edge_list_it &edge_it, vector<size_t> &v_id);
    bool can_edge_to_edge(const edge_list_it &edge_it, vector<edge> &swap_edge);

    bool intersect_triangle(const size_t id1, const size_t id2, const face &f);

    int create_tet_bin(const char *path1, const char *path2);
    template <typename T>
    int get_next_it(const std::set<pair<int, T>, cmp<T> > &del_set, T &it);

    double compute_area(const face &f);

    double compute_quality(const tet &temp_tet);

    int dump_edge_adj_tet(const char *path, const edge &e);
    int dfs(const size_t index, vector<bool> &is_visited,
            const zjucad::matrix::matrix<size_t> &face,
            const jtf::mesh::edge2cell_adjacent &edge_adj);
    bool test();

    bool test(const zjucad::matrix::matrix<size_t> &tet_ver);

    template<typename T>
    bool is_ele_in_vec(const T &ele, const vector<T> &vec, size_t &index)
    {
      for(size_t i = 0; i < vec.size(); ++i)
        {
          if(ele == vec[i])
            {
              index = i;
              return true;
            }
        }
      index = -1;
      return false;
    }

    list<tet> tet_list;
    list<edge> edge_list;
    list<face> face_list;
    vector<zjucad::matrix::matrix<double> > vertex_vec;
    tet2edge_type tet2edge;
    edge2tet_type edge2tet;
    face2tet_type  face2tet;
    tet2face_type tet2face;
    vertex2tet_type vertex2tet;
    edge2it_type edge2it;
    face2it_type face2it;
    tet2it_type tet2it;
    tet_list_it current_tet_it;
    edge_list_it current_edge_it;
    face_list_it current_face_it;

    tet2original_type tet2original;
    face2orginal_type face2orginal;
    size_t initial_tet_size;


    size_t del_face_index;

  };
}
/* class tet_mesh_operation */
/* { */
/*  public: */
/*   int edge_split(const edge_list_it &edge_it); */
/*   int edge_collapse(const edge_list_it &edge_it); */
/*   int edge_to_face(const edge_list_it &edge_it); */
/*   int face_to_edge(const face_list_it &face_it); */
/*   int edge_to_edge(const edge_list_it &edge_it); */
/*   edge get_opposite_edge(const tet_list_it &tet_it); */
/* }; */

#endif 
