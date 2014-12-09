#include <string>
#include <iostream>
#include <zjucad/ptree/ptree.h>

using namespace std;
using boost::property_tree::ptree;

#define CALL_SUB_PROG(prog)						\
  int prog(ptree &pt);						  \
  if(pt.get<string>("prog.value") == #prog)				\
  return prog(pt);

static void pt_description(ptree &pt)
{
  pt.put("prog.desc", "sub-prog: param3d");
}


int main(int argc, char *argv[])
{
  ptree pt;

  try{
    zjucad::read_cmdline(argc,argv,pt);
    pt_description(pt);

    CALL_SUB_PROG(param3d);
    CALL_SUB_PROG(param3d_integer);
   // CALL_SUB_PROG(test);

  }catch(std::exception& e){
    cerr << endl;
    cerr << "[error] "<<e.what() << endl;
    cerr << "Usage:" << endl;
    zjucad::show_usage_info(std::cerr, pt);
  }
  return 0;
}
