#include "hex_process.h"
#include <string>
#include <iostream>
using namespace std;
using namespace boost::property_tree;
using namespace zjucad::matrix;

int extrapolate_hexmesh_from_surface(const matrixst &hex,
                                     const matrixd &node,
                                     matrixst &hex_new,
                                     matrixd &node_new,
                                     boost::property_tree::ptree &pt)
{
    pt.put("ext_alg.desc","extrapolate methods name");
    string ext_alg = pt.get<string>("ext_alg.value");
//    if(ext_alg == "linear"){

//    }
    std::cerr << "# [error] this function is not finished." << std::endl;
    return 0;
}
