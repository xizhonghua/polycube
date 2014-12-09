#include <iostream>

#include <jtflib/mesh/io.h>

#include "../common/util.h"
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"
#include "../hexmesh/io.h"

using namespace std;

int remove_extra_nodes(int argc, char * argv[])
{
    if(argc != 4){
        cerr << "# [usage] remove_extra_nodes tet/hex input_tet output_tet." << endl;
        return __LINE__;
    }

    jtf::mesh::meshes tm;
    const string tet_or_hex = argv[1];
    if(tet_or_hex == "tet")
        if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &tm.node_, &tm.mesh_))
            return __LINE__;
    if(tet_or_hex == "hex")
        if(jtf::hexmesh::hex_mesh_read_from_wyz(argv[2], tm.mesh_, tm.node_, 1))
            return __LINE__;

    remove_extra_node(tm.mesh_, tm.node_);

    if(tet_or_hex == "tet")
        if(jtf::mesh::tet_mesh_write_to_zjumat(argv[3], &tm.node_, &tm.mesh_))
            return __LINE__;
    if(tet_or_hex == "hex")
        if(jtf::hexmesh::hex_mesh_write_to_wyz(argv[3], tm.mesh_, tm.node_))
            return __LINE__;

    return 0;
}

