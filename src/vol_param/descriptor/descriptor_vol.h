#ifndef DESCRIPTOR_VOL_H
#define DESCRIPTOR_VOL_H

#include <memory>
#include <jtflib/mesh/mesh.h>
#include <jtflib/algorithm/gauss_elimination.h>
#include <jtflib/optimizer/opt.h>
#include "def.h"
#include "descriptor_base.h"
#include "../../common/util.h"

class descriptor_vol : public descriptor_base
{
public:
  virtual ~descriptor_vol(){}
public:
  ///
  /// @brief init all related basic information
  /// @param mesh
  /// @param node
  /// @return
  ///
  virtual int init(
      const zjucad::matrix::matrix<size_t> & mesh,
      const zjucad::matrix::matrix<double> & node,
      boost::property_tree::ptree &pt);

  ///
  /// @brief set_objective, set objective of problem
  /// @param mesh input mesh
  /// @param node input mesh nodes
  /// @param pt  input configuration
  ///
  virtual void set_objective(
      const zjucad::matrix::matrix<size_t> & mesh,
      const zjucad::matrix::matrix<double> & node,
      boost::property_tree::ptree &pt);

  ///
  /// @brief set_constraint, set set_constraint of problem
  /// @param mesh input mesh
  /// @param node input mesh node
  /// @param pt input configuration
  ///
  virtual void set_constraint(
      const zjucad::matrix::matrix<size_t> & mesh,
      const zjucad::matrix::matrix<double> & node,
      boost::property_tree::ptree &pt);

  ///
  /// @brief add_eqn_constraint
  /// @param fc
  /// @param hard_or_soft
  /// @param eqn_type
  ///
  virtual void add_eqn_constraint(jtf_func_ptr fc,int hard_or_soft,const char * eqn_type =0);

  ///
  /// @brief get_objective
  /// @return shared_ptr
  ///
  virtual jtf_func_ptr
  get_objective() const {return obj_;}

  ///
  /// @brief get_objective
  /// @param obj_type input obj_type,
  /// @return jtf_func_cons_ptr
  ///
  virtual jtf_func_ptr
  get_objective(const std::string obj_type) const;

  ///
  /// @brief get_constraint
  /// @return shared_ptr
  ///
  virtual const std::vector<jtf_func_ptr> &
  get_eqn_constraint()const {return eqn_cons_;}

  virtual std::vector<jtf_func_ptr> & get_eqn_constraint() {
    return eqn_cons_;}

  virtual const std::pair<size_t,size_t>
  get_eqn_constraint(const std::string type) const {
    const auto & it = eqn_cons_type_.find(type);
    if(it == eqn_cons_type_.end()) return std::make_pair(-1,-1);
    return it->second;
  }
  ///
  /// @brief get_constraint
  /// @return shared_ptr
  ///
  virtual const std::vector<jtf_func_ptr> &
  get_ineqn_constraint()const {return ineqn_cons_;}

  ///
  /// @brief get_constraint_type, each item contains <type, range_begin, range_end>
  ///        here [range_begin, range_end) is the container order.
  /// @return
  virtual const std::map<std::string, std::pair<size_t,size_t> > &
  get_eqn_constraint_type() const {return eqn_cons_type_;}

  ///
  /// @brief has_node_mapping
  /// @return
  ///
  virtual bool has_node_mapping()const{
    return (bi_.NM.ZT.size(1) > 0 && bi_.NM.ZT.size(2) > 0?true:false);
  }

  ///
  /// @brief get_node_mapping
  /// @return
  ///
  virtual const node_mapping & get_node_mapping() const{
    return bi_.NM;
  }

  ///
  /// @brief local_stiffening
  ///
  void frame_local_stiffening(const zjucad::matrix::matrix<size_t> & mesh,
                              const zjucad::matrix::matrix<double> & orig_node,
                              const zjucad::matrix::matrix<double> & real_node);

  ///
  /// @brief reset      reset all objective and constraints
  /// @param uncut_tet
  /// @param uncut_node
  /// @param cut_tet
  /// @param cut_node
  /// @param ptree
  ///
  void reset(const zjucad::matrix::matrix<size_t> & uncut_tet,
             const zjucad::matrix::matrix<double> & uncut_node,
             const zjucad::matrix::matrix<size_t> & cut_tet,
             const zjucad::matrix::matrix<double> & cut_node,
             boost::property_tree::ptree & pt);

  ///
  /// \brief recover_gap_node, when node mapping is required, gap
  ///        node should be recovered from tet nodes
  /// @brief node input node should be with gap nodes
  void recover_gap_node(zjucad::matrix::matrix<double> & node) const;

  ///
  /// \brief update_Z_q update Z and q according to new gauss elimination result
  /// \return 0 : Z,q has been changed, 1: nothing changed, others error
  int update_Z_q(const std::pair<size_t,double> & fix_integer);

  jtf::opt::callbacks * get_callback() const
  {
    if(cb_.get()) return cb_.get();
    return 0;
  }
private:
  jtf_func_ptr obj_;
  std::vector<jtf_func_ptr> obj_vec_;
  std::map<std::string, jtf_func_ptr> obj_type_;

  std::vector<jtf_func_ptr> eqn_cons_;
  //std::vector<std::tuple<std::string,size_t,size_t> > eqn_cons_type_;
  std::map<std::string, std::pair<size_t,size_t> > eqn_cons_type_;

  std::vector<jtf_func_ptr> ineqn_cons_;
  std::vector<std::tuple<std::string,size_t,size_t> > ineqn_cons_type_;

  std::shared_ptr<jtf::opt::callbacks> cb_;
  std::vector<std::pair<jtf_func *, double> > wf_;

  class basic_infor
  {
  public:
    ///
    /// @brief load_node_mapping This function load node mapping
    /// @param node_mapping_file
    /// @param node_mapping
    /// @param variable_number
    /// @return
    ///
    int load_node_mapping(const char* node_mapping_file,
                          zjucad::matrix::matrix<size_t> & node_mapping,
                          const size_t variable_number);

    ///
    /// \brief load_node_comp_eqn, load node compression equation, it contains
    ///        x equation and u independent variables
    /// \param node_mapping_file
    /// \param node_mapping it store node_mappingT
    /// \return
    ///
    int load_node_comp_eqn(const char * node_mapping_file,
                           hj::sparse::csc<double,int32_t> & node_mappingT);

    int load_zero_index(const char * zero_idx_file,
                        std::set<size_t> & zero_idx_set);
    ///
    /// @brief get_independent_node
    /// @param node_mapping
    /// @param node
    /// @param independent_node
    ///
    void get_independent_node(const zjucad::matrix::matrix<size_t> & node_mapping,
                              const zjucad::matrix::matrix<double> & node,
                              zjucad::matrix::matrix<double> & independent_node)const;

  public:
    std::unique_ptr<jtf::mesh::face2tet_adjacent> fa;
    zjucad::matrix::matrix<size_t> outside_face;
    zjucad::matrix::matrix<size_t> outside_face_idx;

    zjucad::matrix::matrix<size_t> outside_face_uncut;

    zjucad::matrix::matrix<size_t> uncut_tet;
    zjucad::matrix::matrix<double> uncut_node;

    zjucad::matrix::matrix<size_t> cut_tet;
    zjucad::matrix::matrix<double> cut_node;

    zjucad::matrix::matrix<size_t> cut_tet2tet;
    std::unique_ptr<jtf::mesh::face2tet_adjacent> fa_uncut;

    node_mapping NM;
    std::set<size_t> zero_index_; // some gaps will be zero after elimination

    zjucad::matrix::matrix<zjucad::matrix::matrix<double> > frames;
    zjucad::matrix::matrix<double> vol_weight;
    std::vector<hj_func_cons_ptr> vol_func_;

    zjucad::matrix::matrix<size_t> orig_face_in_cut; // store cut face which belong to original surface.
    std::unique_ptr<jtf::mesh::edge2cell_adjacent> ea_orig_in_cut;
//    std::vector<std::vector<size_t> > variable_idx; // equations
//    std::vector<std::vector<double> > coefficient;

    boost::unordered_map<std::pair<size_t,size_t>,size_t> inner_type;
    boost::unordered_map<size_t,size_t> surface_type;

    zjucad::matrix::matrix<size_t> gap_faces;
    std::vector<std::vector<std::pair<size_t,size_t> > > cut_patches;
    std::vector<size_t> rot_type;

    std::vector<std::vector<size_t> > node_group;
    std::vector<bool> integer_group_flag ;
  };

public:
  SIMPLEX sim_;
  basic_infor bi_;
};
#endif // DESCRIPTOR_VOL_H
