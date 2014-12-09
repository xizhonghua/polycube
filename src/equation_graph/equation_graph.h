#ifndef JTF_EQUATION_GRAPH_H
#define JTF_EQUATION_GRAPH_H

#include <deque>
#include <set>
#include <stack>
#include <memory>
#include <tuple>
#include <boost/unordered_set.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/dynamic_bitset.hpp>

#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/util/vertex_connection.h>
#include <hjlib/sparse/sparse.h>
#include <jtflib/algorithm/gauss_elimination.h>
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"
#include "../hex_param/find_singularities.h"
#include <jtflib/algorithm/gauss_elimination.h>
#include "../common/transition_type.h"

class step_state;
class state_each_step;
class edge_with_type;

template <typename T>
class group : public boost::unordered_set<T>
{
public:
  typedef typename boost::unordered_set<T>::const_iterator const_iterator;
public:

  void operator << (const T & a){
    this->insert(a);
  }
  void operator << (const group<T> & a){
    this->insert(a.begin(), a.end());
  }
  void operator >> (group<T> &a)const{
    a.insert(this->begin(), this->end());
  }
  void merge(group<T> &a){
    if(&a == this)
      return ;
    this->insert(a.begin(), a.end());
    a.clear();
  }
  friend std::ostream& operator << (std::ostream &output,
                                    const group<T> &g)
  {
    output << "# [info] group has: ";
    for(group<T>::const_iterator it = g.begin(); it != g.end(); ++it)
      output << *it << " ";
    output << std::endl;
    return output;
  }
};

///
/// @brief eliminate_gaps, eliminate gaps for those gap = 0 and rot  = I faces
/// @param cut_tet
/// @param uncut_tet
/// @param cut_tet2tet
/// @param inner_type
/// @param cut_tet_new
/// @return
///
int eliminate_gaps(
    const zjucad::matrix::matrix<size_t> &cut_tet,
    const zjucad::matrix::matrix<size_t> &uncut_tet,
    const zjucad::matrix::matrix<size_t> &cut_tet2tet,
    zjucad::matrix::matrix<double> & cut_node,
    const zjucad::matrix::matrix<double> & uncut_node,
    const boost::unordered_map<std::pair<size_t, size_t>, size_t> &inner_type,
    zjucad::matrix::matrix<size_t> &cut_tet_new);

class transition_elimination
{
public:

  /// @breief this function create transition elimination function,
  ///         it assumes all inner singularity are exposed outside
  /// @param tetmesh  input tetmesh
  /// @param cut_mesh input cut mesh
  /// @param cut_tet2tet input cut tet to tet mapping
  /// @param tet_node input mesh node
  /// @param inner_face_type input inner face jump type
  /// @param surface_type input surface restricted type, surface index is of original fa
  /// @return transition_elimination *
  static transition_elimination * create(
      const zjucad::matrix::matrix<size_t> & tetmesh,
      const zjucad::matrix::matrix<size_t> & cut_mesh,
      const zjucad::matrix::matrix<size_t> & cut_tet2tet,
      const zjucad::matrix::matrix<double> & tet_node,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type,
      const boost::unordered_map<size_t,size_t> & surface_type,
      const bool is_restricted_type);

  ///
  /// @brief create_with_gap, create transition with gaps
  /// @param tetmesh
  /// @param cut_mesh
  /// @param cut_tet2tet
  /// @param tet_node
  /// @param inner_face_type
  /// @param surface_type
  /// @param is_restricted_type
  /// @return
  ///
  static transition_elimination * create_with_gap(
      const zjucad::matrix::matrix<size_t> & tetmesh,
      const zjucad::matrix::matrix<size_t> & cut_mesh,
      const zjucad::matrix::matrix<size_t> & cut_tet2tet,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type,
      const boost::unordered_map<size_t,size_t> & surface_type,
      const bool is_restricted_type,
      const zjucad::matrix::matrix<double> &node);

  ///
  /// @brief export node_group
  /// @return
  ///
  const std::vector<group<size_t> > out() const {
    return node_group_;
  }

  ///
  /// \brief get_integer_group_idx, this function defines which variable
  ///        should be integer in post hexmeshing
  /// \return
  ///
  const std::vector<size_t> get_integer_group_idx() {
    if(integer_groups_.size())
      return integer_groups_;

    std::set<size_t> integer_groups_set;
    for(const auto & v_idx : integer_variable_){
        integer_groups_set.insert(node2group_[v_idx]);
      }
    integer_groups_.resize(integer_groups_set.size());
    std::copy(integer_groups_set.begin(), integer_groups_set.end(),
              integer_groups_.begin());
    return integer_groups_;
  }

  ///
  /// \brief get_equation, output variable index and coeff in each equation
  /// \param eqn_idx
  /// \param eqn_coeff
  ///
  void get_equation(std::vector<std::vector<size_t> > & eqn_idx,
                    std::vector<std::vector<double> > & eqn_coeff)const;

  int get_trivial_cut_face_pair(
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_type,
      const zjucad::matrix::matrix<size_t> & cut_tet,
      const zjucad::matrix::matrix<size_t> & cut_tet2tet,
      std::vector<zjucad::matrix::matrix<size_t> > & trivial_cut_face_pairs) const ;

private:
  int init(const zjucad::matrix::matrix<size_t> & tetmesh,
           const zjucad::matrix::matrix<size_t> & cut_mesh,
           const zjucad::matrix::matrix<size_t> & cut_tet2tet,
           const zjucad::matrix::matrix<double> & tet_node,
           const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type,
           const boost::unordered_map<size_t,size_t> & surface_type,
           const bool is_restricted_type);


  int init_with_gap(
      const zjucad::matrix::matrix<size_t> & tetmesh,
      const zjucad::matrix::matrix<size_t> & cut_mesh,
      const zjucad::matrix::matrix<size_t> & cut_tet2tet,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type,
      const boost::unordered_map<size_t,size_t> & surface_type,
      const bool is_restricted_type,
      const zjucad::matrix::matrix<double> & node);

private:

  int get_gaps_stupid(
      const zjucad::matrix::matrix<size_t> & tetmesh,
      const zjucad::matrix::matrix<size_t> & cut_mesh,
      const zjucad::matrix::matrix<size_t> & cut_tet2tet,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type,
      const zjucad::matrix::matrix<double> & node);

  int get_gaps_elimination(
      const zjucad::matrix::matrix<size_t> &tetmesh,
      const zjucad::matrix::matrix<size_t> &cut_mesh,
      const zjucad::matrix::matrix<size_t> &cut_tet2tet,
      const boost::unordered_map<std::pair<size_t, size_t>, size_t> &inner_face_type,
      const zjucad::matrix::matrix<double> & node);

  ///
  /// @brief construct_basic_info
  /// @param tetmesh
  /// @param cut_tet
  /// @param cut_tet2tet
  /// @param tet_node input tetmesh node
  /// @param node variants for gauss eliminator
  /// @param node_flag variant_flag for gauss eliminator
  /// @param gnode gap variants for gauss eliminator
  /// @param gnode_flag gap variant_flag for gauss eliminator
  /// @return 0 if ok or non-zeros
  ///
  int construct_basic_info(const zjucad::matrix::matrix<size_t> & tetmesh,
                           const zjucad::matrix::matrix<size_t> & cut_tet,
                           const zjucad::matrix::matrix<size_t> & cut_tet2tet,
                           const zjucad::matrix::matrix<double> & tet_node,
                           std::vector<double> & node,
                           boost::dynamic_bitset<> & node_flag,
                           std::vector<double> & gnode,
                           boost::dynamic_bitset<> & gnode_flag);

  ///
  /// @brief add_inner_transition for each inner face, equations are built based
  ///        on gauss eliminator.
  /// @param tetmesh
  /// @param cut_mesh
  /// @param cut_tet2tet
  /// @param inner_face_type
  /// @return 0 if ok, or non-zeros
  ///
  int add_inner_transition(
      const zjucad::matrix::matrix<size_t> & tetmesh,
      const zjucad::matrix::matrix<size_t> & cut_mesh,
      const zjucad::matrix::matrix<size_t> & cut_tet2tet,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type);

  ///
  /// @brief add inner transition based on edge type
  ///        For an edge (p_i,p_j): (I-type) (f_t(p_i)-f_t(p_j)) = 0
  /// @param tetmesh
  /// @param cut_mesh
  /// @param cut_tet2tet
  /// @param inner_face_type
  /// @return 0 if ok, non-zeros if error.
  ///
  int add_inner_transition_edge(
      const zjucad::matrix::matrix<size_t> &tetmesh,
      const zjucad::matrix::matrix<size_t> &cut_mesh,
      const zjucad::matrix::matrix<size_t> &cut_tet2tet,
      const boost::unordered_map<std::pair<size_t, size_t>, size_t> &inner_face_type);

  ///
  /// @brief add_inner_transition_edge2
  /// @param tetmesh
  /// @param cut_mesh
  /// @param cut_tet2tet
  /// @param inner_face_type
  /// @return
  ///
  int add_inner_transition_edge2(
      const zjucad::matrix::matrix<size_t> &tetmesh,
      const zjucad::matrix::matrix<size_t> &cut_mesh,
      const zjucad::matrix::matrix<size_t> &cut_tet2tet,
      const boost::unordered_map<std::pair<size_t, size_t>, size_t> &inner_face_type);


  ///
  /// @brief add_inner_transition_gap based on gap equation:
  ///        for a cut face:
  ///        f_t(p) = Pi_st*f_s(p) + g_st
  /// @param tetmesh
  /// @param cut_mesh
  /// @param cut_tet2tet
  /// @param inner_face_type
  /// @return 0 if ok, or non-zeros
  ///
  int add_inner_transition_gap(
      const zjucad::matrix::matrix<size_t> &tetmesh,
      const zjucad::matrix::matrix<size_t> &cut_mesh,
      const zjucad::matrix::matrix<size_t> &cut_tet2tet,
      const boost::unordered_map<std::pair<size_t, size_t>, size_t> &inner_face_type);

  ///
  /// @brief add_surface_transition for each outside_face of original tetmesh.
  /// @param tetmesh
  /// @param cut_mesh
  /// @param cut_tet2tet
  /// @param surface_type
  /// @param is_restricted_type if true, surface_type contains [0/1/2], or its type contains [0,23]
  /// @return 0 if ok, or non-zeros
  ///
  int add_surface_transition(
      const zjucad::matrix::matrix<size_t> & tetmesh,
      const zjucad::matrix::matrix<size_t> & cut_mesh,
      const zjucad::matrix::matrix<size_t> & cut_tet2tet,
      const boost::unordered_map<size_t,size_t> & surface_type,
      const bool is_restricted_type);

  ///
  /// @brief merge_node from to "to"
  /// @param from
  /// @param to
  ///
  void merge_node(const size_t & from, const size_t &to);

  ///
  /// @brief analysis equations from gauss elimination, all equations
  ///        (after gauss elimination) containing only two variants whose
  ///        coefficients are 1 and -1 should be grouped.
  ///
  void collect_group_variants_from_equations();

  //////////////////////////////////////////////////////////////////////////////
  ///                   NODE WITH GAPS
  //////////////////////////////////////////////////////////////////////////////
  /// \brief construct_basic_info_with_gaps
  /// \param tetmesh
  /// \param cut_mesh
  /// \param cut_tet2tet
  /// \param tet_node
  /// \param inner_face_type
  /// \return
  ///
  int construct_basic_info_with_gaps(
      const zjucad::matrix::matrix<size_t> & tetmesh,
      const zjucad::matrix::matrix<size_t> & cut_mesh,
      const zjucad::matrix::matrix<size_t> & cut_tet2tet,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type,
      const zjucad::matrix::matrix<double> & node);

  ///
  /// \brief cluster_transition_patch, based on walking on cut faces, ang regard
  /// each neigubour face with the same transition type as belongs to the same patch
  ///
  void cluster_transition_patch(
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
      const zjucad::matrix::matrix<double> * node_ptr = 0);

  ///
  /// \brief cluster_transition_patch_each, this function regard each jump face as
  /// a patch, too slow!!!
  /// \param inner_type
  /// \param node_ptr
  ///
  void cluster_transition_patch_each(
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_type,
      const zjucad::matrix::matrix<double> * node_ptr);

  ///
  /// \brief cluster_transition_patch_by_elimination, this function first eliminate
  /// gap variants, and group them, based on these gap groups, it will build patches
  /// \param inner_type
  /// \param node_ptr
  ///
  void cluster_transition_patch_by_elimination(
      const zjucad::matrix::matrix<size_t> & tetmesh,
      const zjucad::matrix::matrix<size_t> & cut_mesh,
      const zjucad::matrix::matrix<size_t> & cut_tet2tet,
      const zjucad::matrix::matrix<double> & tet_node,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_type);

  size_t get_face_rot_type(
      const zjucad::matrix::matrix<size_t> & one_face_0,
      const zjucad::matrix::matrix<size_t> & one_face_1,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_type);

  ///
  /// \brief build_transitions
  /// \param uncut_mesh
  /// \param cut_mesh
  /// \param cut_tet2tet
  /// \param inner_face_type
  /// \param surface_type
  ///
  void build_transitions(
      const zjucad::matrix::matrix<size_t> &uncut_mesh,
      const zjucad::matrix::matrix<size_t> &cut_mesh,
      const zjucad::matrix::matrix<size_t> &cut_tet2tet,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_type,
      const boost::unordered_map<size_t,size_t> & surface_type,
      bool is_restricted_type);

  ////
  /// \brief reorder_nodes
  /// \param cut_tet
  /// \param M_new2old
  /// \param M_old2new
  /// \param P_out_end_idx
  ///
  void reorder_nodes(
      const zjucad::matrix::matrix<size_t> & cut_tet,
      zjucad::matrix::matrix<size_t> & M_new2old, // recorded mapping from new index to new one
      zjucad::matrix::matrix<size_t> & M_old2new,
      size_t &P_out_end_idx);

  ///
  /// \brief assemble_inner_transition_and_surface_type
  /// \param ge
  /// \param cut_face_pairs
  /// \param cut_face_patches
  /// \param cut_face_patch_rot_type
  /// \param surface_type
  /// \param is_restricted_type
  ///
  void assemble_inner_transition_and_surface_type(
      const zjucad::matrix::matrix<size_t> & tetmesh,
      const zjucad::matrix::matrix<size_t> & cut_mesh,
      const zjucad::matrix::matrix<size_t> & cut_tet2tet,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_type,
      const boost::unordered_map<size_t,size_t> & surface_type,
      const bool is_restricted_type,
      const size_t original_node_number);


  int add_surface_transition_with_gaps(
      const zjucad::matrix::matrix<size_t> & tetmesh,
      const zjucad::matrix::matrix<size_t> & cut_mesh,
      const zjucad::matrix::matrix<size_t> & cut_tet2tet,
      const boost::unordered_map<size_t,size_t> & surface_type,
      const bool is_restricted_type);

  ///
  /// @brief add_inner_transition_gap2, original variables with gap variables
  /// @param tetmesh
  /// @param cut_mesh
  /// @param cut_tet2tet
  /// @param original_node_number
  ///
  void add_inner_transition_gap2(
      const zjucad::matrix::matrix<size_t> & tetmesh,
      const zjucad::matrix::matrix<size_t> & cut_mesh,
      const zjucad::matrix::matrix<size_t> & cut_tet2tet,
      const size_t original_node_number);

  ///
  /// \brief add_inner_transition_gap3, add equations about pure gap variables
  /// \param tetmesh
  /// \param cut_mesh
  /// \param cut_tet2tet
  /// \param original_node_number
  ///
  void add_inner_transition_gap3(
      const zjucad::matrix::matrix<size_t> & tetmesh,
      const zjucad::matrix::matrix<size_t> & cut_mesh,
      const zjucad::matrix::matrix<size_t> & cut_tet2tet,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t > & inner_face_type,
      const size_t original_node_number);

  ///
  /// \brief compress_equation
  /// \param M_old2new
  /// \param P_out_end_index
  /// \param eqns
  ///
  void compress_equation(
      const zjucad::matrix::matrix<size_t> &M_old2new,
      const size_t P_out_end_index,
      std::vector<jtf::algorithm::equation<double> > & all_eqn,
      size_t & independent_variable_number);

public:
  ///
  /// \brief get_node_comp
  /// \param eqn
  ///
  void get_node_comp(
      hj::sparse::csc<double,int32_t> & MT) const;

  ///
  /// \brief get_zero_idx
  /// \param zero_set
  ///
  void get_zero_idx(std::set<size_t> & zero_set) const;

  ///
  void get_ordered_cut_face_patches(
      zjucad::matrix::matrix<size_t> & cut_faces,
      std::vector<std::vector<std::pair<size_t,size_t> > > & patches,
      std::vector<size_t> & rot_type) const;
  //////////////////////////////////////////////////////////////////////////////
private:
  class basic_info{
  public:
    ///
    /// @brief get gap node defined of tet pair, each node is a float defined
    ///        of its uvw variants index.
    /// @param tet_pair input tet_pair
    /// @param gnode output gnode
    /// @return 0 if ok, 1 if reversed order, other if error
    ///
    int get_gnode_idx(const std::pair<size_t,size_t> & tet_pair,
                      zjucad::matrix::matrix<double> & gnode_idx)const {
      gnode_idx.resize(3,1);
      int reversed = 0;
      auto it = tet_pair2g_idx.find(tet_pair);
      if(it == tet_pair2g_idx.end()){
          it = tet_pair2g_idx.find(std::make_pair(tet_pair.second, tet_pair.first));
          reversed = 1;
        }

      if(it == tet_pair2g_idx.end()) {// this face is glued, do not need to calculated
          //std::cerr << "# [error] can not find tet_pair_to_gnode." << std::endl;
          return __LINE__;
        }

      for(size_t di = 0; di < 3; ++di)
        gnode_idx[di] = 3 * it->second + di;
      return reversed;
    }

    ///
    /// \brief get_gnode_idx_of_patch, return gnode_idx of tet_pair
    /// \param tet_pair
    /// \param tet_pair2cut_patch_idx
    /// \param gnode_idx
    /// \return
    ///
    int get_gnode_idx_of_patch(const std::pair<size_t, size_t> &tet_pair,
                               const std::map<std::pair<size_t,size_t>,size_t> & tet_pair2cut_patch_idx,
                               zjucad::matrix::matrix<double> &gnode_idx) const{
      int reversed = 0;
      auto it = tet_pair2cut_patch_idx.find(tet_pair);
      if(it == tet_pair2cut_patch_idx.end()){
          it = tet_pair2cut_patch_idx.find(std::make_pair(tet_pair.second, tet_pair.first));
          reversed = 1;
        }
      if(it == tet_pair2cut_patch_idx.end()){
          return __LINE__;
        }
      if(gnode_idx.size() != 3)
        gnode_idx.resize(3,1);
      gnode_idx = 3 * it->second * zjucad::matrix::ones<double>(3,1);
      gnode_idx[1] += 1;
      gnode_idx[2] += 2;
      return reversed;
    }

    ///
    /// @brief get_fnode_idx
    /// @param idx input point idx
    /// @param node_idx output point variant
    ///
    void get_fnode_idx(const size_t idx,
                       zjucad::matrix::matrix<double> & node_idx)const
    {
      node_idx.resize(3,1);
      for(size_t di = 0; di < 3; ++di)
        node_idx[di] = 3 * idx + di;
    }
    ///
    /// @brief get_inner_face_type get inner rotation type based on tet pair
    /// @param tet_pair
    /// @param inner_face_type
    /// @return
    ///
    size_t get_inner_face_type(
        const std::pair<size_t,size_t> & tet_pair,
        const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type)const
    {
      auto it = inner_face_type.find(tet_pair);
      if(it == inner_face_type.end()) return TRIVIAL_TYPE;
      return it->second;
    }

    ///
    /// @brief get_tet_pair_from_cut_face
    /// @param cut_face0
    /// @param cut_face1
    /// @return
    ///
    const std::pair<size_t,size_t>  get_tet_pair_from_cut_face(
        const std::vector<size_t> & cut_face0,
        const std::vector<size_t> & cut_face1)const
    {
      const size_t face_idx_0 = fa_cut->get_face_idx(&cut_face0[0]);
      const size_t face_idx_1 = fa_cut->get_face_idx(&cut_face1[0]);

      if(face_idx_0 == -1 || face_idx_1 == -1)
        return std::make_pair(-1,-1);

      const std::pair<size_t,size_t> & tet_pair_0 =
          fa_cut->face2tet_[face_idx_0];
      const std::pair<size_t,size_t> & tet_pair_1 =
          fa_cut->face2tet_[face_idx_1];

      return std::make_pair(
            tet_pair_0.first==-1?tet_pair_0.second:tet_pair_0.first,
            tet_pair_1.first==-1?tet_pair_1.second:tet_pair_1.first);
    }

  public:
    std::shared_ptr<jtf::mesh::face2tet_adjacent> fa,fa_cut;
    zjucad::matrix::matrix<size_t> outside_face_cut;
    jtf::mesh::one_ring_tet_at_edge ortae;

    zjucad::matrix::matrix<size_t> cut_face_pairs;
    std::shared_ptr<jtf::mesh::edge2cell_adjacent> ea_cut_face_pair;
    std::vector<std::vector<std::pair<size_t,size_t> > > cut_face_patches_;
    std::vector<size_t> cut_face_patch_rot_type_;

    typedef std::vector<size_t> one_face;
    std::map<one_face, std::vector<one_face> > orig_face2_cut_faces;
    std::map<std::pair<size_t,size_t>, size_t> tet_pair2g_idx;
    std::map<one_face, size_t> orig_face2g_idx;
  };

  std::unique_ptr<jtf::algorithm::gauss_eliminator<double> > ge;
  basic_info bi_;
  std::vector<size_t> node2group_;
  std::vector<group<size_t> > node_group_;

  std::unique_ptr<jtf::algorithm::gauss_eliminator<double> > ge_gap;
  std::vector<double> gnodes_;
  boost::dynamic_bitset<> gnode_flag_;

  std::set<size_t> integer_variable_;
  std::vector<size_t> integer_groups_;

  //  // only used in init_with_gap
  //  std::vector<std::vector<std::pair<size_t,double> > > node_comp_eqn;
  //  size_t independent_variable_num;

  hj::sparse::csc<double,int32_t>  node_comp_mapping;
};

class singularity_graph
{
public:
  static singularity_graph* create(
      const jtf::mesh::one_ring_tet_at_edge & ortae_cut,
      const jtf::mesh::one_ring_tet_at_edge & ortae,
      const matrixst & cut_tet2tet,
      const matrixst & outside_face_cut,
      const matrixst & cut_face_pair,
      const matrixst & outside_face_idx_cut,
      const boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & jump_face_to_rot_idx,
      const bool no_surface = false);
  ~singularity_graph(){}

  /**
   * @brief This function is used to insert transistion on a given face
   *
   * @param face        input face
   * @param cut_tet     input cut_tet mesh
   * @param cut_tet2tet input cut_tet mesh to original mesh index mapping
   * @param fa_cut      inputjtf::mesh::face2tet_adjacent
   * @param ortae_original  input original ortae
   * @param rot_type    input rot_type
   * @param rot_type_map_idx  input rot_type_map_idx tet_pair to rot_type_idx
   *                          which is indexed in rot_type
   * @param mode        mode = 0: do not return non-zero if meet compound edge
   *                    mode = 1: return non-zeros if meet compound edge
   * @return int  return 0 if works well, or return non-zero
   */
  int insert_transition(
      const std::pair<size_t,size_t> & face,
      const matrixst & cut_tet,
      const matrixst & cut_tet2tet,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const jtf::mesh::one_ring_tet_at_edge & ortae_original,
      const std::vector<size_t> & rot_type,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_pair_to_rot_idx,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair2rot_idx,
      const size_t mode = 0);

  int insert_transition_group(
      const std::vector<std::pair<size_t,size_t> > & face_groups,
      const matrixst & cut_tet,
      const matrixst & cut_tet2tet,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const jtf::mesh::one_ring_tet_at_edge & ortae_original,
      const std::vector<size_t> & rot_type,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_pair_to_rot_idx,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair2rot_idx,
      const size_t mode = 0);


  /**
   * @brief This function is used to insert transistion on a given surface face
   *
   * @param rot_type    input rot_type
   * @param face_idx_to_rot_tye_idx  input face index to rotation type index
   *                                 which is indexed in rot_type
   * @param tet_pair2rot_type_idx_map input tet pair to rotation type map
   * @param outside_face_idx_cut     input outside_face_idx of cut tet mesh
   * @param fa_cut      inputjtf::mesh::face2tet_adjacent of cut tet mesh
   * @param e2a_cut     input edge2adjacent of cut tet mesh
   * @param cut_tet2tet input cut tet to tet mapping
   * @param orig_ortae       input original one ring tet at edge
   * @param cut_tet2tet       input cut tet 2 orig tet index
   * @param outside_face_normal_cut input outside face normal
   * @param orig_ortae      input original one_ring_tet_at_edge
   * @param face_idx    surface triangle index
   * @param step_stamp  used to record the step order
   * @return int  return 0 if works well, or return non-zero
   */
  int insert_surface_transition(
      const size_t & rot_type,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const matrixst & cut_tet2tet,
      const size_t & face_idx,
      const std::vector<size_t> & rot_type_vec,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair2rot_idx);

  int insert_surface_transition_group(
      const size_t & rot_type,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const matrixst & cut_tet2tet,
      const std::vector<size_t> & face_group,
      const std::vector<size_t> & rot_type_vec,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair2rot_idx);
  /**
   * @brief This function is used to check whether the configuration is valid
   *
   * @param cut_tet2tet  input cut_tet2tet info
   * @param mode         mode = 0; check all loops
   *                     mode = 1: return false if meet one loop
   * @return int  return 0 if works well, or return non-zero
   */
  bool is_valid_with_info(
      const matrixst & cut_tet2tet,
      const std::vector<size_t> & rot_type,
      const jtf::mesh::one_ring_tet_at_edge & ortae,
      const boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
      const size_t mode = 0 );

  //  bool is_valid_with_info_accelerate(
  //      const matrixst & cut_tet2tet,
  //      const std::vector<size_t> & rot_type,
  //      const jtf::mesh::one_ring_tet_at_edge & ortae,
  //      const boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
  //      const size_t mode = 0 );

  /**
   * @brief convert node idx to point idx in cut mesh
   *
   * @param node_idx  input node_idx
   * @return pair<size_t,size_t>  return point_idx and u/v/w, if can not find
   *                              this node, return <-1,-1>
   */
  std::pair<size_t,size_t> get_point_idx_from_node_idx(
      const size_t &node_idx) const;


  /**
   * @brief convert point idx with uvw type to node idx
   *
   * @param point_idx  input node_idx
   * @param type       u/v/w of this point
   * @return size_t    if find this point return node_idx, or return -1
   */
  size_t get_node_idx_from_point_idx( const size_t & point_idx,
                                      const size_t & type)const;

  size_t get_gnode_idx_from_rot_idx( const size_t & rot_idx,
                                     const size_t & type)const;

  int get_step_state(step_state & ss)const;
  int set_step_state(step_state & ss);

  int save_state(const std::string & file_name )const;
  int load_state(const std::string & file_name );

  int clear_state_mem();
  int save_state_mem(state_each_step & ss)const;
  int load_state_mem(const state_each_step & ss);

  std::pair<size_t,size_t> get_rot_idx_from_gnode_idx(const size_t & g_idx)const;
  const std::pair<size_t,size_t>& get_jump_face_from_gnode_idx(
      const size_t gnode_idx)const;
  /**
   * @brief This function is used to add singularity constraints, which contains
   *        two equal constraints and one non-equal constraint
   *
   * @param se  input singularity edge
   * @return int  return 0 if works well, or return non-zero
   */
  int add_singularity(
      const std::vector<size_t> &tet_loop_vec,
      const matrixst & cut_tet,
      const matrixst & cut_tet2tet,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & rot_type_map_idx,
      const std::vector<size_t> &rot_type,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const edge_with_type & se);

  /**
   * @brief This function is used to combine two different fnode into one group
   * @param fnode_idx_0 fnode_idx_combined_to
   * @param fnode_idx_1 fnode_idx_combined_from
   * @return int  return 0 if works well, or return non-zero
   */
  int combine_fnode(const size_t & fnode_idx_0,
                    const size_t & fnode_idx_1);

  /**
   * @brief This function is used to add trivial edge constraint, which is
   *        about transition, and as the result the fnode will be grouped
   *        according to gnode.
   *
   * @param tet_loop_vec  tet loop around this edge
   * @param cut_tet       cut tet mesh
   * @param cut_tet2tet   cut tet vertex to original tet mesh mapping
   * @param rot_type      rot type of each jump face pair current
   * @param fa_cut        face to tets adjacent of cut tet mesh
   * @param rot_type_map_idx  a map which record the tet pair to index of rot type
   * @return int  return 0 if works well, or return non-zero
   */
  int add_trivial_edge(
      const std::vector<size_t> &tet_loop_vec,
      const matrixst & cut_tet,
      const matrixst & cut_tet2tet,
      const std::vector<size_t> &rot_type,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair2rot_idx);


  /**
   * @brief This function is used to update fnode according to gnode
   *        information
   *
   * @param face_pair  input face_pair where transition is calculated
   * @param gnode_idx  input gnode idx which should be
   * @return int  return 0 if works well, or return non-zero
   */
  int connect_fnode_according_to_gnode_at_modified_face_pair(
      const matrixst & cut_tet2tet,
      const std::vector<size_t> & rot_type,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & rot_type_map_idx,
      const jtf::mesh::face2tet_adjacent &fa_cut,
      const std::pair<size_t,size_t> & face_pair,
      const size_t & gnode_idx);


  int check_unready_edges(
      const jtf::mesh::one_ring_tet_at_edge & ortae_original,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const matrixst & cut_tet,
      const matrixst & cut_tet2tet,
      const std::vector<size_t> & rot_type,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pari2rot_idx,
      const size_t mode = 0);

  bool is_group_info_broken()const;

  int convert_g_eq_to_f_eq(
      const std::vector<size_t> & rot_type,
      const matrixst & cut_tet2tet,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & jump_face_to_rot_idx,
      std::shared_ptr<jtf::algorithm::gauss_eliminator<double> > &fe,
      std::vector<double> &fnode_vec,
      boost::dynamic_bitset<> & fnode_flag)const;
  // there are 3 * surface point node, store the group info
  // which stands for the equal constraints
  std::vector<size_t> fnode_;
  std::vector<double> gnode_; // there are 3 * jump face node
  boost::dynamic_bitset<> gnode_flag_;
  boost::unordered_map<size_t,size_t> outside_point_idx_;
  std::vector<size_t> outside_points_;
  std::vector<std::pair<size_t,size_t> > jump_faces_; // jump faces of cut mesh
  //! @brief this uset stores the non-euqal constraints which are represented
  //  as unordered edges, to check whether there is a loop, need to build a graph
  // of group, stack is used to store the modification stamp for each edge
  //boost::unordered_map<std::pair<size_t,size_t>,std::stack<size_t> > unordered_edge_;
  //boost::unordered_map<std::pair<size_t,size_t> > unordered_edge_;

  //boost::unordered_set<std::pair<size_t,size_t> > unordered_edge_;
  //  boost::unordered_set<size_t> unordered_edge_idx_;
  //  boost::unordered_set<size_t> possible_edge_idx_;

  //  std::vector<std::pair<size_t,size_t> > edges_;
  //  boost::unordered_map<std::pair<size_t,size_t>,size_t> edge2idx_;

  std::list<std::tuple<size_t,size_t,size_t> > unready_edges_;

  //boost::unordered_map<std::pair<size_t,size_t>,std::stack<size_t> > ordered_edge_;
  //! During type assignment, whether an edge is a restricted edge can not be determined
  // accoring to current types, those type which is given later may make it become
  // singularity later, we gather them in possible edges
  // bool stands for inner near surface edge(true), or others (false)
  boost::unordered_map<std::pair<size_t,size_t>,bool> possible_edge_;
  std::vector<group<size_t> > groups_; // group of fnode_

  boost::unordered_map<size_t,size_t> rot_idx2g_idx_;
  std::vector<size_t> gnode2rot_idx_;
  boost::unordered_map<std::pair<size_t,size_t>, size_t> jump_face_to_gnode_idx_;
  //boost::unordered_map<size_t, std::pair<size_t,size_t> > gnode_idx_to_jump_face_;
  std::shared_ptr<jtf::algorithm::gauss_eliminator<double> > ge_ptr;

  // to speed up, use point_around_faces_with_rot_type_idx instead of surface_idx
  boost::unordered_map<size_t, boost::unordered_set<size_t> > point_around_faces_with_rot_idx;
  std::vector<std::deque<std::pair<size_t,size_t> > > chains_;

  boost::unordered_set<std::pair<size_t,size_t> > graph_edges;

private:
  int add_singularity_equation_for_gnode(
      const std::vector<size_t> & tet_loop_vec,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & rot_type_map_idx,
      const std::vector<size_t> & rot_type,
      const matrixst & cut_tet,
      const matrixst & cut_tet2tet,
      const jtf::mesh::face2tet_adjacent &fa_cut,
      const size_t & singularity_axis_type,
      jtf::algorithm::equation<double> &eq);

  int add_trivial_edge_equation_for_gnode(
      const std::vector<size_t> & tet_loop_vec,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & rot_type_map_idx,
      const std::vector<size_t> & rot_type,
      const matrixst & cut_tet,
      const matrixst & cut_tet2tet,
      const jtf::mesh::face2tet_adjacent &fa_cut,
      std::vector<jtf::algorithm::equation<double> > &eqs);

  size_t extract_free_axis_of_edge(const size_t &p0, const size_t & p1);

  int extract_directed_chain_from_node_edges(
      const boost::unordered_set<std::pair<size_t,size_t> > & edges,
      std::vector<std::deque<std::pair<size_t,size_t> > > & chains,
      const matrixst &  cut_tet2tet,
      const std::vector<size_t> & rot_type,
      const jtf::mesh::one_ring_tet_at_edge & ortae,
      const boost::unordered_map<size_t,size_t>& surface_face2rot_idx,
      const size_t mode = 0);

  //! @brief return -1 : shit happen
  //         return -2 : this edge is not ready
  //         return 0-23: edge type
  int cal_singularity_type_at_given_tet(
      const std::vector<size_t> & tet_loop,
      const size_t &begin_tet,
      const std::vector<size_t> &rot_type,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & rot_type_map_idx,
      std::vector<size_t> * reordered_tet_loop_ptr = 0);


  int init(
      const jtf::mesh::one_ring_tet_at_edge & ortae_cut,
      const jtf::mesh::one_ring_tet_at_edge & ortae,
      const matrixst & cut_tet2tet,
      const matrixst & outside_face_cut,
      const matrixst & cut_face_pair,
      const matrixst & outside_face_idx_cut,
      const boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & jump_face_to_rot_idx,
      const bool no_surface = false);

  singularity_graph(){}

public: // debug info
  std::vector<std::deque<std::pair<size_t,size_t> > > singularity_edges_of_cut_tet_;
  std::vector<std::deque<size_t> > singularity_type_of_cut_tet_;

  std::vector<std::deque<std::pair<size_t,size_t> > > compound_edges_;
  std::vector<std::deque<std::pair<size_t,size_t> > > loop_edges_;
  //std::vector<size_t> jump_faces_;
  // singularity_graph& operator = (const singularity_graph & sg);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

class equation_graph
{
public:
  equation_graph(const size_t & varaint_number)
    :dim_(3){
    assert(varaint_number % dim_ == 0);
    reset(varaint_number);
  }
  equation_graph():dim_(3){}

  void reset(const size_t &variant_number){
    groups_.clear();
    v2g_idx_.clear();
    edges_.clear();
    groups_.resize(variant_number);
    v2g_idx_.resize(variant_number);
    for(size_t  i = 0; i < variant_number; ++i){
        groups_[i] << i;
        v2g_idx_[i] = i;
      }
    node_ = 0;
    //point_idx_ = 0;
  }
  // based operation

  class virtual_point
  {
  public:
    std::vector<size_t> orig_point;
    std::vector<size_t> uvw_group; // uvw with group_idx
  };


  ///
  /// @brief get_variant_idx
  /// @param vp input virtual point
  /// @param di input variant di
  /// @return variant idx
  ///
  size_t get_variant_idx(const virtual_point& vp, const size_t &di) const
  {
    return vp.uvw_group[di];
  }

  ///
  /// @brief get_variant_idx
  /// @param real_point_idx
  /// @param di
  /// @return
  ///
  size_t get_variant_idx(const size_t& real_point_idx, const size_t &di) const
  {
    return  dim_ * real_point_idx + di;
  }

  enum class graph_state {ABSOLUTE_FAIL, COMMON_FAIL, SUCCEED};

public:
  //! @brief this function is used to build equaltiy constraints and inequality constraints
  //! @param node_group input node_group
  //! @param edges input edges to be extracted inequaltiy constraints
  void assemble_equation_graph_generally(
      const std::vector<std::vector<size_t> > & node_group,
      const std::vector<std::pair<size_t,size_t> > &edges);

  ///
  /// \brief assemble_equation_graph_generally
  /// \param node_group_dim
  /// \param edges edges should not contain duplicates
  /// \return if find duplicate edges in virtual points, then degeneration happens, and return false, or return true
  void assemble_equation_graph_generally(
      const std::vector<std::vector<std::vector<size_t> > > & node_group_dim,
      const std::vector<std::pair<size_t,size_t> > & edges);

  //! @brief merge group which contains "from" to one containing "to"
  //! @param from
  //! @param to    merge "from" to "to"
  void add_equal_constraint(const size_t &from,const size_t &to);

  //! @brief add variant edges, I will check the virtual points associated to each edge,
  //!       if two edges have the same virtual points, this means they are cutted open,
  //!       and can be glued
  //! @param variant_edge input edge with variant index
  void add_inequal_constraint(const std::tuple<size_t,size_t,size_t> & edge);

  //! @brief check whether the graph is valid or not
  bool check_valid(bool show_info = false);

  //! @brief bind node and point_idx
  void bind_node(const zjucad::matrix::matrix<double> *node){
    node_ = node;
  }

  //! @brief test whether node_ and point_idx_ is assigned
  //! @param return true if bind well, or false
  bool has_node()const{return (node_);}

public: // util tool

  //! @desc: add all edges which has only one free variant,
  //!        and degenerated(zero-free variants) edges(this is for a complete check)
  //!        Warning: this function require equal_constraint to be built, input
  //!        edges should not contain duplicates
  //! @param edges: all edges of this mesh
  //! @return if find duplicate edges with virtual points, then return false, or return true.
  void auto_add_inequal_cons_based_on_equal_cons(
      const std::vector<std::pair<size_t,size_t> > &edges);

protected: // new scheme
  typedef std::deque<std::tuple<size_t,size_t, size_t> > path_type;

  //! @brief extract edges to path according to edge degree
  //! @param path output pathes
  void edge2path(std::vector<path_type> & path);

  //! @brief collect_direction_variant based on report
  //! @param dir_path output_path with direction
  //! @param dir_v direction variant
  void collect_direction_variant(std::vector<path_type> & dir_path,
                                 std::vector<group<int> > & dir_v)const;



  //! @brief unify edge order along each axis
  //! @param directional_path unified path
  //! @param direction_v grouped direction variant
  void unify_axis_order(const std::vector<size_t> & path_idx_of_one_axis,
                        std::vector<path_type> & directional_path,
                        std::vector<group<int> > & direction_v) const;

  //! @brief find next order assignment
  //! @param order  next generated order
  //! @param orig_path input orig path
  //! @param output path output path adjusted by order
  //! @param dir_v_idx  input order idx to grouped directional variant
  //! @return return true if find next assigment, false if all assignment
  bool next_assignment(boost::dynamic_bitset<> & order,
                       const std::vector<path_type> & orig_path,
                       std::vector<path_type> & output_path,
                       const std::vector<group<int> > & dir_v_idx)const;

  //! @brief check direction edge graph based on input direction path, this step
  //!        only check absolutely error
  //! @param dir_path a given path which defines the direction
  //! @param show_info show information if needed
  //! @return if graph is ok return true, otherwise return false
  graph_state check_directional_edge_graph_step0(const std::vector<path_type> & dir_path,
                                                 bool show_info = false);

  //! @brief check direction edge graph based on input direction path, this step
  //!        only check possible error
  //! @param dir_path a given path which defines the direction
  //! @param show_info show information if needed
  //! @return if graph is ok return true, otherwise return false
  graph_state check_directional_edge_graph_step1(const std::vector<path_type> & dir_path,
                                                 bool show_info = false);

  //! @brief convert all points to virtual point in case of duplicated points
  //!        introduced by arbitrary cutting
  //! @param variant2group input variant2group
  //! @param vp output virtual points
  void convert2virtual_point(const std::vector<size_t> & variant2group,
                             std::vector<virtual_point> & vp);

  //! @brief check whether input dir_path contains self loop
  //! @param dir_path input directional path
  //! @param show_info show information if needed
  //! @return return true if find self loop, or false
  bool check_self_loop(const std::vector<path_type> & dir_path,
                       bool show_info = false)const;

  //! @brief check whether two type path have the same endings
  //! @param dir_path input directional path
  //! @param show_info show information if needed
  //! @return return true if find duplicate type path, or false
  bool check_duplicate_type_path(const std::vector<path_type> & dir_path,
                                 bool show_info = false)const;

  //! @brief check whether the direction assignment is ok, depends on report def 4.3
  //! @param dir_path input directional path
  //! @param show_info show information if needed or false
  //! @return return true if path direction is ok, or return false
  bool check_path_direction(const std::vector<path_type> & dir_path,
                            bool show_info = false) const;

  //! @brief check whether two edges intersection based on given direction assignment
  //! @param dir_path input directional path
  //! @param show_info  show information if needed or false
  //! @return return true if find intersection or false
  bool check_edge_intersection(const std::vector<path_type> & dir_path,
                               bool show_info = false) const;

  //! @brief check whether there is a loop in graph
  //! @param dir_path input directional path
  //! @parma show_info show information if needed or false
  //! @return return true if find loop or false
  bool check_loop(const std::vector<path_type> & dir_path,
                  bool show_info = false);

  //! @brief this function is used to detech two path intersection
  //! @param path_idx_0 one path group with path index
  //! @param path_idx_1 one path group with path index
  //! @param dir_path   input directional pathes
  //! @return if find intersection, return the path_idx, or return <-1,-1>
  std::pair<size_t,size_t> find_path_intersection(
      const std::vector<size_t> & path_idx_0,
      const std::vector<size_t> & path_idx_1,
      const std::shared_ptr<vertex_connection<DIRECT> >  vc0,
      const std::shared_ptr<vertex_connection<DIRECT> >  vc1,
      const std::vector<path_type> &dir_path)const;

private:
  const size_t dim_; // dimension of variant. 3: uvw
  std::vector<group<size_t> > groups_;
  std::vector<size_t> v2g_idx_; // variant to group index

  std::set<std::tuple<size_t,size_t,size_t> > edges_; // edge with
  const zjucad::matrix::matrix<double> *node_;

  std::vector<virtual_point> vp_;
  std::vector<size_t> p2vp_;

  // this information is only used for debugging
  boost::unordered_map<std::pair<size_t,size_t>, std::vector<std::pair<size_t,size_t> > >
  vpedge2rawedge_;
};

#endif //JTF_EQUATION_GRAPH_H
