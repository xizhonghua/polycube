#ifndef TYPE_PATCH_GRAPH_H
#define TYPE_PATCH_GRAPH_H

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <vector>

#include "../tetmesh/tetmesh.h"
//#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/mesh.h>
#include "../common/def.h"

class type_patch_graph
{
public:
  typedef boost::unordered_map<size_t, boost::unordered_set<size_t> > linking_type;

  type_patch_graph(boost::unordered_map<size_t,size_t> & surface_type)
    : surface_type_(surface_type){}

  //! @brief build a graph
  int build_graph(
      const matrixst & tet,
      const matrixst & cut_tet,
      const matrixd & node,
      const matrixst & cut_tet2tet,
      const matrixst & outside_face,
      const matrixst & outside_face_idx,
      const jtf::mesh::face2tet_adjacent & fa,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const jtf::mesh::edge2cell_adjacent & ea,
      const std::vector<std::pair<size_t,size_t> > & g_unknown_face_pair);

  //! @brief clean all graph.
  int clean_graph(){ patches_.clear(); patch_linking_info_.clear();}

  //! @brief merge group g0 to g1
  int merge_group(const size_t & g0, const size_t & g1);

  //! @brief del greoup g
  int del_group(const size_t & g);

  //! @brief get one ring adjacent groups of g
  //         if g is not exist, throw an exception
  const boost::unordered_set<size_t> &
  get_linking_of_group(const size_t & g) const;

  //! @brief get group faces, if g is not exist, throw an exception
  const boost::unordered_set<size_t> &
  get_group_faces(const size_t &g) const { return patches_.at(g); }
  
  //! @brief get group type, if g is not exist, return -1, use "check "
  //         to check whether all faces of this group has the same type.
  size_t get_group_type(const size_t & g,
                        const bool check = false)const;

  //! @brief get the whole linking info
  const linking_type& get_linking_info() const {return patch_linking_info_;}

  //! @brief get the patch degree, if g is not exist, return -1
  size_t get_patch_degree(const size_t & g) const;

private:
  boost::unordered_map<size_t, size_t> & surface_type_;
  std::vector<boost::unordered_set<size_t> > patches_;
  linking_type patch_linking_info_;
};

#endif
