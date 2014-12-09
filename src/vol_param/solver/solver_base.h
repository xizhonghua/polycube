#ifndef SOLVER_BASE_H
#define SOLVER_BASE_H

#include <zjucad/matrix/matrix.h>
#include <zjucad/ptree/ptree.h>
#include <jtflib/optimizer/opt.h>

class descriptor_base;

class solver_base
{
public :
  solver_base(){}
  virtual ~solver_base(){}

  virtual std::string get_solver_type()const {
    throw std::logic_error(std::string( __PRETTY_FUNCTION__) + " empty func.");
  }
public:
  virtual int solve(zjucad::matrix::matrix<double> & init_node,
                    const std::shared_ptr<descriptor_base> &,
                    boost::property_tree::ptree &)  {
    throw std::logic_error(std::string( __PRETTY_FUNCTION__) + " empty func.");}
  virtual int set_callback(jtf::opt::callbacks * cb){
    throw std::logic_error(std::string( __PRETTY_FUNCTION__) + " empty func.");
  }
  virtual jtf::opt::callbacks* get_callback(){
    throw std::logic_error(std::string( __PRETTY_FUNCTION__) + " empty func.");
  }
};

#endif // SOLVER_BASE_H
