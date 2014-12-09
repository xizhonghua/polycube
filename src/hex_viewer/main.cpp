#include <zjucad/ptree/ptree.h>
#include <zjucad/matrix/matrix.h>
#include <string>
#include <iostream>

using namespace std;
using namespace zjucad::matrix;
using boost::property_tree::ptree;

int main(int argc, char *argv[])
{
  ptree pt;

  try {
    zjucad::read_cmdline(argc,argv,pt);
    pt.put("prog.desc","<frame, fram_inner, frame_polycube, frame2vtk, cut_tet,cut_tet_inner, int_param, find_singularities,find_singularities2, find_singularities3,int_pts, init_zyz, init_zyz_inner, init_zyz_polycube, label_polycube, polycube_param, remove_black_lines, find_singularities_type_consistency, map_hex_to_tet, cut_hex_to_tet>");
    pt.put("tet.desc","file name to the input tet model");

#define CALL_SUB_PROG(prog)						\
  int prog(ptree &pt);						  \
  if(pt.get<string>("prog.value") == #prog)				\
  return prog(pt);

    CALL_SUB_PROG(view_tet);

  }catch(std::exception &e){
    cerr << endl;
    cerr << "[error] " << e.what() << endl;
    cerr << "Usage: " << endl;
    zjucad::show_usage_info(std::cerr, pt);
  }
  return 0;
}
