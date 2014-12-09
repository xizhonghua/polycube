#ifndef DESCRIPTOR_DEF_H
#define DESCRIPTOR_DEF_H

#include <jtflib/function/function.h>
#include <hjlib/function/function.h>

typedef jtf::function::functionN1_t<double,int32_t> jtf_func;
typedef std::shared_ptr<const jtf_func> jtf_func_cons_ptr;
typedef std::shared_ptr<jtf_func> jtf_func_ptr;

typedef hj::function::function_t<double,int32_t> hj_func;
typedef std::shared_ptr<const hj_func> hj_func_cons_ptr;
typedef std::shared_ptr<hj_func> hj_func_ptr;

enum SIMPLEX{TRI, TET};
#endif
