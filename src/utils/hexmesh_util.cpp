#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include "../common/vtk.h"
#include "../common/def.h"
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>

using namespace std;
using namespace zjucad::matrix;

int read_hexmesh(const char * file, matrixd &node, matrixst &hexmesh)
{
    ifstream ifs(file);
    if(ifs.fail()){
        cerr << "# read hexmesh fail." << endl;
        return __LINE__;
    }
    string temp;
    size_t vertex_num;
    size_t hex_num;
    ifs >> temp >> vertex_num;
    ifs >> temp >> hex_num;
    node.resize(3,vertex_num);
    hexmesh.resize(8,hex_num);
    for(size_t t = 0; t < vertex_num * 3; ++t){
        ifs >> node[t];
    }
    for(size_t t = 0; t < hex_num * 8; ++t){
        ifs >> hexmesh[t];
    }
    return 0;
}

int pull_hexmesh(int argc, char *argv[]){
    if(argc != 4){
        cerr << "input_hexmesh_file  output_file pull_delta[0,1]" << endl;
        return __LINE__;
    }

    matrixd node;
    matrixst hexmesh;
    if(read_hexmesh(argv[1],node,hexmesh)){
        return __LINE__;
    }
    double pull_delta = atof(argv[3]);
    assert(node.size(2));
    matrixd min_position = node(colon(),0);
    matrixd max_position = node(colon(),0);
    {
        for(size_t t = 1; t < node.size(2); ++t){
            for(size_t i = 0; i < 3; ++i){
                if(node(i,t) < min_position[i]) min_position[i] = node(i,t);
                if(node(i,t) > max_position[i]) max_position[i] = node(i,t);
            }
        }
    }
    const double radius = norm(max_position - min_position) / 2.0;
    const matrixd center = (max_position + min_position)/2.0;
    matrixd new_node(3, 8 * hexmesh.size(2));
    for(size_t t = 0; t < hexmesh.size(2); ++t){
        matrixd center_of_each_hex = zeros<double>(3,1);
        for(size_t j = 0; j < 8; ++j){
            center_of_each_hex += node(colon(),hexmesh(j,t));
        }
        center_of_each_hex /= 8.0;

        double step[3]; // assign step for x,y,z axis
        for(size_t k = 0; k < 3; ++k)
            step[k] = pull_delta * fabs(max_position[k] - min_position[k])
                    * static_cast<size_t>(100 * fabs(center_of_each_hex[k] - center[k])
                                          / fabs(max_position[k] - min_position[k])) / 100.0;

        double sign[3]; // to get the sign for x,y,z axis
        matrixd one_axis = zeros<double>(3,1);
        for(size_t j = 0; j < 8; ++j){

            for(size_t i = 0; i < 3; ++i){
                one_axis = zeros<double>(3);
                one_axis[i] = max_position[i] - min_position[i];
                sign[i] = 1.0;
                if(dot(center_of_each_hex - center, one_axis) < -1e-8)
                    sign[i] = -1;

                new_node(i, t * 8 + j) = node(i,hexmesh(j,t)) + step[i] * sign[i];// * (center_of_each_hex[i] - center[i]);
            }
            hexmesh(j,t) = t * 8 + j;
        }
    }

    {// to vtk
        string input_hex_str = argv[1];
        input_hex_str += ".pull_" ;
        input_hex_str += argv[3];
        input_hex_str += ".vtk";
        ofstream ofs(input_hex_str.c_str());
        hex2vtk(ofs,&new_node[0],new_node.size(2),&hexmesh[0],hexmesh.size(2));
    }
    return 0;
}
