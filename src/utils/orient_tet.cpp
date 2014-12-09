#include <iostream>
#include <jtflib/mesh/io.h>

#include "../common/def.h"
#include "../tetmesh/tetmesh.h"


using namespace std;

int orient_tet(int argc, char * argv[])
{
    if(argc != 3){
        cerr << "# [usage] orient_tet input_tet output_tet" << endl;
    }
    matrixst tet;
    matrixd node;

    if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &node, &tet)){
        cerr << "# [error] can not open tet file." << endl;
        return __LINE__;
    }

    orient_tet(node, tet);

    if(jtf::mesh::tet_mesh_write_to_zjumat(argv[2], &node, &tet)){
        cerr << "# [error] can not save out tet." << endl;
        return __LINE__;
    }
    return 0;
}
