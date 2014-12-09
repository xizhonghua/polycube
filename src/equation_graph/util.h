#ifndef EQUATION_GRAPH_UTIL_H
#define EQUATION_GRAPH_UTIL_H


#include <deque>
#include <tuple>
#include <list>
#include <jtflib/mesh/mesh.h>
#include <boost/dynamic_bitset.hpp>
#include <zjucad/matrix/matrix_expression.h>

#include <jtflib/algorithm/gauss_elimination.h>

class step_state{
public:
  std::list<jtf::algorithm::equation<double> > es_;
  boost::dynamic_bitset<> gnode_flag_;
  std::vector<double> gnode_;
  std::vector<size_t> fnode_;
  std::list<std::tuple<size_t,size_t,size_t> >  unready_edges_;
};

class state_each_step{
public:
  std::vector<double> nodes_;
  boost::dynamic_bitset<> node_flag_;
  std::vector<jtf::algorithm::equation<double> > eq_vec_;
  std::vector<size_t> fnode_;
  std::list<std::tuple<size_t,size_t,size_t> > unready_edges_;
};

class edge_with_type
{
public:
  edge_with_type(){}
  edge_with_type(const size_t & point0_, const size_t & point1_,
                 const size_t & type_)
  {
    type() = type_;
    point0() = point0_;
    point1() = point1_;
  }
  size_t& type(){ return std::get<2>(ewt);}
  const size_t& type()const{ return std::get<2>(ewt);}
  size_t& point0(){ return std::get<0>(ewt); }
  const size_t& point0()const{ return std::get<0>(ewt); }
  size_t& point1(){ return std::get<1>(ewt); }
  const size_t& point1()const{ return std::get<1>(ewt); }
private:
  std::tuple<size_t,size_t,size_t> ewt;
};

template <typename T1, typename T2>
bool matrix_equal(const zjucad::matrix::matrix_expression<T1> & A,
                  const zjucad::matrix::matrix_expression<T2> & B)
{
  if(A().size(1) != B().size(1) || A().size(2) != B().size(2))
     return false;
  double a = 0;
  for(size_t i = 0; i < A().size(); ++i)
    a += (A()[i] - B()[i]) * (A()[i] - B()[i]);
  if(fabs(a) < 1e-8)
    return true;
  return false;
}

template <typename T>
void reverse_container(std::pair<T,T> & v){
  std::swap(v.first,v.second);
}

template <typename T1, typename T2>
void reverse_container(std::tuple<T1,T1,T2> & v){
  std::swap(std::get<0>(v),std::get<1>(v));
}


template <typename ITERATOR >
void reverse_container(ITERATOR first, ITERATOR last)
{
  reverse(first, last);
  for(ITERATOR it = first; it != last; ++it)
    reverse_container(*it);
}

template <typename T>
void reverse_container(std::vector<T> & v){
  reverse_container(v.begin(), v.end());
}

template <typename T>
void reverse_container(std::deque<T> & v){
  reverse_container(v.begin(), v.end());
}

template <typename T1, typename T2>
void assert_iterator(const T1 & container,  T2 & it,
                     const std::string & msg)
{
  if(it == container.end())
    throw std::out_of_range(msg);
}

class equation_graph;

//! @brief this function is used to assemble equation graph with input node group
//! and surface type, which contains six types [0,5].
//! @param tm input tetmesh
//! @param node_group variants group
//! @param surface type input surface type [0,5]: {+u,-u,+v,-v,+w,-w}
//! @param eg output equation graph
void assemble_equation_graph_generally(
    const jtf::mesh::meshes & tm,
    const std::vector<std::vector<size_t> > & node_group,
    const std::map<size_t,int> & surface_type,// this surfac type contains signs: 0~5
    equation_graph &eg);

#endif // UTIL_H
