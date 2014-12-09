#ifndef SOLVER_JTF_H
#define SOLVER_JTF_H

#include "../descriptor/descriptor_base.h"
#include "solver_base.h"

class solver_jtf : public solver_base
{
public:
  solver_jtf():cb_(nullptr){}
  virtual ~solver_jtf(){}
  virtual std::string get_solver_type()const {return "jtf";}
public:
  virtual int solve(zjucad::matrix::matrix<double> &init_node,
                    const std::shared_ptr<descriptor_base> & desc,
                    boost::property_tree::ptree & pt) ;
  virtual int set_callback(jtf::opt::callbacks * cb){
    cb_ = cb;
  }
  virtual jtf::opt::callbacks* get_callback(){
    return cb_;
  }
private:
  jtf::opt::callbacks * cb_;
};

#endif
