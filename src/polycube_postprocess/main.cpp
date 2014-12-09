#include <boost/property_tree/ptree.hpp>
#include <zjucad/ptree/ptree.h>

#include <string>
#include <iostream>
#include "../mesh_func/tri-normal.h"
using namespace std;
using boost::property_tree::ptree;

int main(int argc, char *argv[])
{
  ptree pt;

  try {
     zjucad::read_cmdline(argc,argv,pt);

     pt.put("prog.desc", "[polycube_improve,update_surface_node, simplify_polycube_struct]");
#define CALL_SUB_PROG(prog)         \
     int prog(ptree &pt);                         \
     if(pt.get<string>("prog.value") == #prog)    \
     return prog(pt);
    CALL_SUB_PROG(polycube_improve);
    CALL_SUB_PROG(update_surface_node);
  } catch(std::exception& e) {
    cerr << endl;
    cerr << "[error] "<<e.what() << endl;
    cerr << "Usage:" << endl;
     zjucad::show_usage_info(std::cerr, pt);
  }
  return 0;
}
