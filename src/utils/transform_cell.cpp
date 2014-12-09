#include "../tetmesh/hex_io.h"
#include "../numeric/util.h"
#include "../hexmesh/io.h"
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>

#include <zjucad/matrix/matrix.h>
#include <iostream>

using namespace std;
using namespace zjucad::matrix;

matrixd rotate_x_angle(const double &angle)
{
    matrixd rot = zeros<double>(3,3);
    rot[0] = 1.0;
    rot(1,1) = rot(2,2) = cos(angle/180.0 * My_PI());
    rot(1,2) = -1* sin(angle/180.0 * My_PI());
    rot(2,1) = sin(angle/180.0 * My_PI());
    return rot;
}

matrixd rotate_y_angle(const double &angle)
{
    matrixd rot = zeros<double>(3,3);
    rot(1,1) = 1.0;
    rot(0,0) = rot(2,2) = cos(angle/180.0 * My_PI());
    rot(0,2) = sin(angle/180.0 * My_PI());
    rot(2,0) = -1*sin(angle/180.0 * My_PI());
    return rot;
}

matrixd rotate_z_angle(const double &angle)
{
    matrixd rot = zeros<double>(3,3);
    rot(2,2) = 1.0;
    rot(0,0) = rot(1,1) = cos(angle/180.0 * My_PI());
    rot(0,1) = -1* sin(angle/180.0 * My_PI());
    rot(1,0) = sin(angle/180.0 * My_PI());
    return rot;
}

int transform_cell(int argc, char * argv[])
{
    if(argc != 9){
        //cerr << "# [usage] recenter_tet input_tet angle_x angle_y angle_z output_tet" << endl;
        cerr << "# [usage] transform_cell tet/hex/tri/quad cell_input "
             << "angle_x angle_y angle_z scale [recenter:y/n] output." << endl;
        return __LINE__;
    }

    const string cell_type = argv[1];

    matrixd node;
    matrixst cell;

    if(cell_type == "tet"){
        if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &node, &cell))
            return __LINE__;
    }else if(cell_type == "hex"){
        if(jtf::hexmesh::hex_mesh_read_from_wyz(argv[2], cell, node, 1))
            return __LINE__;
    }else if(cell_type == "tri"){
        if(jtf::mesh::load_obj(argv[2], cell, node))
            return __LINE__;
    }else if (cell_type == "quad"){
        if(jtf::mesh::load_obj(argv[2], cell, node))
            return __LINE__;
    }else {
        cerr << "# [error] unsupported type." << endl;
        return __LINE__;
    }


    assert(node.size(1) == 3);
    matrixd center = node * ones<double>(node.size(2),1);
    center /= node.size(2);

    cerr << "# [info] center " << center << endl;

    for(size_t pi = 0; pi < node.size(2); ++pi){
        node(colon(),pi) -= center;
    }

    double rot_xyz[] = {atof(argv[3]), atof(argv[4]), atof(argv[5])};

    node = temp( rotate_x_angle(rot_xyz[0])
                 * rotate_y_angle(rot_xyz[1])
                 * rotate_y_angle(rot_xyz[2]) * node);

    double scale_  = atof(argv[6]);
    node *= scale_;

    const string recenter_flag = argv[7];
    if(recenter_flag == "n" || recenter_flag == "N"){
      for(size_t pi = 0; pi < node.size(2); ++pi)
        node(colon(),pi) += center;
    }

    if(cell_type == "tet"){
        if(jtf::mesh::tet_mesh_write_to_zjumat(argv[8], &node, &cell)){
            cerr << "# [error] can not write to tetmesh." << endl;
            return __LINE__;
        }
    }else if(cell_type == "hex"){
        if(jtf::hexmesh::hex_mesh_write_to_wyz(argv[8], cell, node))
            return __LINE__;
    }else if(cell_type == "tri" || cell_type == "quad"){
        if(jtf::mesh::save_obj(argv[8], cell, node))
            return __LINE__;
    }else{
        cerr << "# [error] strange, can not be here." << endl;
        return __LINE__;
    }

    return 0;
}
