#ifndef DESCRIPTOR_BASE_H
#define DESCRIPTOR_BASE_H

#include <zjucad/matrix/matrix.h>
#include <zjucad/ptree/ptree.h>
#include <map>
#include <hjlib/sparse/sparse.h>
#include "def.h"
#include "../common/util.h"

class descriptor_base{
public:
  descriptor_base(){}
  virtual ~descriptor_base(){}

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
      boost::property_tree::ptree &pt) {
    throw std::logic_error(std::string( __PRETTY_FUNCTION__) + " empty func.");}

  ///
  /// @brief set_objective, set objective of problem
  /// @param mesh input mesh
  /// @param node input mesh nodes
  /// @param pt  input configuration
  ///
  virtual void set_objective(
      const zjucad::matrix::matrix<size_t> & mesh,
      const zjucad::matrix::matrix<double> & node,
      boost::property_tree::ptree &pt){
    throw std::logic_error(std::string( __PRETTY_FUNCTION__) + " empty func.");}

  ///
  /// @brief set_constraint, set set_constraint of problem
  /// @param mesh input mesh
  /// @param node input mesh node
  /// @param pt input configuration
  ///
  virtual void set_constraint(
      const zjucad::matrix::matrix<size_t> & mesh,
      const zjucad::matrix::matrix<double> & node,
      boost::property_tree::ptree &pt){
    throw std::logic_error(std::string( __PRETTY_FUNCTION__) + " empty func.");}

  ///
  /// @brief add_eqn_constraint
  /// @param fc
  /// @param hard_or_sof, 0: soft, non-zeros: hard
  /// @param eqn_type
  ///
  virtual void add_eqn_constraint(jtf_func_ptr fc, int hard_or_sof, const char * eqn_type =0){
    throw std::logic_error(std::string( __PRETTY_FUNCTION__) + " empty func.");
  }

  virtual void add_ineqn_constraint(jtf_func_ptr fc){
    throw std::logic_error(std::string( __PRETTY_FUNCTION__) + " empty func.");
  }

  ///
  /// @brief get_objective
  /// @return shared_ptr
  ///
  virtual jtf_func_ptr
  get_objective() const {
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " empty func.");}

  ///
  /// @brief get_objective
  /// @param obj_type input obj_type
  /// @return jtf_func_cons_ptr
  ///
  virtual jtf_func_ptr
  get_objective(const std::string obj_type) const {
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " empty func.");}

  ///
  /// @brief get_constraint
  /// @return shared_ptr
  ///
  virtual const std::vector<jtf_func_ptr> &
  get_eqn_constraint()const {
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " empty func.");}

  virtual std::vector<jtf_func_ptr> & get_eqn_constraint() {
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " empty func.");}
  ///
  /// @brief get_constraint
  /// @return shared_ptr
  ///
  virtual const std::pair<size_t,size_t>
  get_eqn_constraint(const std::string type)const {
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " empty func.");}

  virtual const std::vector<jtf_func_ptr> &
  get_ineqn_constraint()const {
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " empty func.");}

  ///
  /// @brief get_constraint_type, each item contains <type, range_begin, range_end>
  ///        here [range_begin, range_end) is the container order.
  /// @return
  virtual const std::map<std::string, std::pair<size_t,size_t> > &
  get_eqn_constraint_type() const {
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " empty func.");}

  ///
  /// @brief get_constraint_type, each item contains <type, range_begin, range_end>
  ///        here [range_begin, range_end) is the container order.
  /// @return
  virtual const std::vector<std::tuple<std::string, size_t, size_t> > &
  get_ineqn_constraint_type()const{
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " empty func.");}

  ///
  /// @brief check whether this descriptor has node mapping
  /// @return ture/false
  ///
  virtual bool has_node_mapping() const{
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " empty func.");
  }

  ///
  /// @brief get_node_mapping
  /// @return
  ///
  virtual const node_mapping &
  get_node_mapping()const{
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " empty func.");
  }


};

#endif // DESCRIPTOR_BASE_H
