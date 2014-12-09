#include <boost/property_tree/ptree.hpp>
#include <zjucad/ptree/ptree.h>

#include <string>
#include <iostream>

using namespace std;
using boost::property_tree::ptree;

int main(int argc, char *argv[])
{
  ptree pt;

  try {
     zjucad::read_cmdline(argc,argv,pt);

#define CALL_SUB_PROG(prog)         \
     int prog(ptree &pt);                         \
     if(pt.get<string>("prog.value") == #prog)    \
     return prog(pt);

    CALL_SUB_PROG(polycube_with_type);
    CALL_SUB_PROG(polycube3);
    CALL_SUB_PROG(polycube3_control);
    CALL_SUB_PROG(polycube_obj);
    CALL_SUB_PROG(polycube_fix_normal);
    CALL_SUB_PROG(polycube_generate);
    //CALL_SUB_PROG(polycube_local);

  } catch(std::exception& e) {
    cerr << endl;
    cerr << "[error] "<<e.what() << endl;
    cerr << "Usage:" << endl;
     zjucad::show_usage_info(std::cerr, pt);
  }
  return 0;
}
