#include "../common/def.h"
#include <jtflib/mesh/io.h>
#include <iostream>
using namespace zjucad::matrix;
using namespace std;

int count_separated_patches(int argc, char * argv[])
{
    if(argc != 3){
        cerr <<  "# [usage] count_separated_patches quad fl" << endl;
        return __LINE__;
    }

//    matrix<size_t> mesh;
//    matrix<double> node;

//    if(jtf::mesh::load_obj(mesh, node, 4)){
//        return __LINE__;
//    }




    return 0;
}
