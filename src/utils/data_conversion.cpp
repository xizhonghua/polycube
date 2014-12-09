#include "../common/IO.h"
#include "../tetmesh/hex_io.h"
#include "../common/zyz.h"

#include <iostream>

using namespace std;
using namespace zjucad::matrix;

int load_yf_tet(const char* filename,
                matrix<size_t> & tet,
                matrix<double> & node)
{
    ifstream ifs(filename);
    if(ifs.fail()){
        cerr << "# [error] can not open yf_tet file." << endl;
        return __LINE__;
    }

    string temp;
    size_t cell_num,node_num;

    ifs >> node_num >> temp; // node_num "vertices"
    ifs >> cell_num >> temp; // cell_num "cells"

    node = zeros<double>(3, node_num);
    tet = zeros<size_t>(4, cell_num);

    for(size_t vi = 0; vi < node_num; ++vi){
        ifs >> node(0,vi) >> node(1,vi) >> node(2,vi);
    }

    for(size_t ti = 0; ti < cell_num; ++ti){
        ifs >> temp >> tet(0,ti) >> tet(1,ti) >> tet(2,ti) >> tet(3,ti);
    }

    orient_tet(node, tet);
    return 0;
}

int load_yf_frame(const char * filename,
                  matrix<matrix<double> > & frame)
{
    ifstream ifs(filename);
    if(ifs.fail()){
        cerr << "# [error] can not open yf_frame." << endl;
        return __LINE__;
    }

    size_t tet_num = 0;
    ifs >> tet_num;

    frame.resize(tet_num,1);
    for(size_t ti = 0; ti < tet_num; ++ti){
        frame[ti].resize(3,3);
        for(size_t i = 0; i < 9; ++i){
            ifs >> frame[ti][i];
        }
    }

    return 0;
}

int tet_yf2jtf(const char * input_tet,
               const char * output_tet)
{
    matrix<size_t> tet;
    matrix<double> node;

    if(load_yf_tet(input_tet, tet, node))
        return __LINE__;


    if(jtf::mesh::tet_mesh_write_to_zjumat(output_tet, &node, &tet))
        return __LINE__;

    return 0;
}

int fld_yf2jtf(const char * input_fld,
               const char * output_fld)
{
    matrix<matrix<double> > frame_matrix;
    if(load_yf_frame(input_fld, frame_matrix))
        return __LINE__;

    cerr << "# [info] frame number " << frame_matrix.size() << endl;
    matrix<double> zyz = zeros<double>(3, frame_matrix.size());
    for(size_t ti = 0; ti < frame_matrix.size(); ++ti){
        rotation_matrix_2_zyz_angle(&frame_matrix[ti][0], &zyz(0,ti), 0);
    }

    if(jtf::mesh::write_matrix(output_fld, zyz))
        return __LINE__;

    return 0;
}

int yf2jtf(int argc, char * argv[])
{
    if(argc != 4){
        cerr << "# [usage] yf2jtf tetfld input_tet[input_fld] output_tet[output_zyz]." << endl;
        return __LINE__;
    }

    const string opt = argv[1];
    if(opt == "tet"){
        return tet_yf2jtf(argv[2], argv[3]);
    }else if(opt == "fld"){
        return fld_yf2jtf(argv[2], argv[3]);
    }else{
        cerr << "# [error] unsupported option." << endl;
        return __LINE__;
    }

    cerr << "# [info] success." << endl;
    return 0;
}
