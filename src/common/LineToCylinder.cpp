#include <cmath>
#include <iostream>
#include <cassert>
#include <fstream>

#include <zjucad/matrix/io.h>
#include "LineToCylinder.h"
#include "Cosserat.h"
#include "../numeric/util.h"
using namespace zjucad::matrix;

typedef matrixd matrixd;
typedef zjucad::matrix::matrix<int> matrixi;
const int segment_num = 16;

int circle_global_coord(int _num,double radius,matrixd& circle_point)
{
    double unit_angle = 2 * My_PI() / _num;
    circle_point.resize(2,_num);
    //	circle_point(0, colon()) = cos(2*PI*colon(0, _num-1)/_num);
    for(size_t i = 0; i < _num; ++i){
        circle_point(0,i) = radius * cos(unit_angle * i);
        circle_point(1,i) = radius * sin(unit_angle * i);
    }
    return 0;
}

int face_of_cylinder(int _num,int line_point_num,zjucad::matrix::matrix<matrixi>& face_index){
    const int quad_face_num = (line_point_num - 1) * _num;
    const int tri_face_num = 2 * _num;
    face_index.resize(quad_face_num + tri_face_num,1);
    for(size_t i = 0; i < (line_point_num - 1); ++i){
        for(size_t j = 0; j < _num; ++j){
            face_index[i * _num + j].resize(4,1);
            face_index[i * _num + j][0] = i * _num + j;
            face_index[i * _num + j][1] = i * _num + (j + 1) % _num;
            face_index[i * _num + j][2] = (i + 1) * _num + (j + 1) % _num;
            face_index[i * _num + j][3] = (i + 1) * _num + j;
        }
    }
    for(size_t i = 0; i < 2; ++i){
        for(size_t j = 0; j < _num; ++j){
            face_index[quad_face_num + i * _num + j].resize(3,1);
            face_index[quad_face_num + i * _num + j][0] = i * _num * (line_point_num - 1) + j;
            face_index[quad_face_num + i * _num + j][1] = i * _num * (line_point_num - 1) + (j + 1) % _num;
            face_index[quad_face_num + i * _num + j][2] = line_point_num * _num +i;
        }
    }
    return 0;
}

int line_to_cylinder(const matrixd& line_point,matrixd& point_coord,zjucad::matrix::matrix<matrixi>& face_index,double radius)
{
    zjucad::matrix::matrix<matrixd> Fn(line_point.size(2));
    matrixd circle_coord;
    point_coord.resize(3,segment_num * line_point.size(2) + 2);
    fill(Fn.begin(), Fn.end(), zeros<double>(3, 3));
    position_to_Bishop_frame_at_node(line_point,Fn);
    circle_global_coord(segment_num,radius,circle_coord);
    for(size_t i = 0; i < Fn.size(1); ++i){
        point_coord(colon(),colon(i * segment_num,(i + 1) * segment_num - 1)) = Fn[i](colon(),colon(1,2)) * circle_coord + line_point(colon(),i) * ones<double>(1,segment_num);
    }
    point_coord(colon(),line_point.size(2) * segment_num) = line_point(colon(),0);
    point_coord(colon(),line_point.size(2) * segment_num + 1) = line_point(colon(),line_point.size(2) - 1);
    face_of_cylinder(segment_num,line_point.size(2),face_index);
    return 0;
}


int convertToCylinder(
        const std::vector<std::vector<line_type> >& global_streamline,
        std::vector<std::vector<matrixd > >& point_coord,
        std::vector<std::vector<zjucad::matrix::matrix<zjucad::matrix::matrix<int> > > >& face_index,
        double radius)
{
    point_coord.resize(global_streamline.size());
    face_index.resize(global_streamline.size());
    for(size_t i = 0;i < global_streamline.size(); ++i){
        point_coord[i].resize(global_streamline[i].size());
        face_index[i].resize(global_streamline[i].size());
        for(size_t j = 0; j < global_streamline[i].size(); ++j){
            //line
            matrixd line_point_i(3, global_streamline[i][j].size());
            //printf("%d\n", global_streamline[i][j].size());
            //std::cerr << "# [info] lines segments size "
            //     << global_streamline[i][j].size() << std::endl;
            for(size_t k = 0; k < global_streamline[i][j].size(); ++k){
                for(size_t t = 0 ; t < 3; ++t){
                    line_point_i(t, k) = global_streamline[i][j][k][t];
                }
            }
            if(line_point_i.size(2) > 1){
                //std::cerr << std::endl;
                line_to_cylinder(line_point_i, point_coord[i][j], face_index[i][j], radius);
            }
        }
    }
    return 0;
}

int saveCylinder(
        const std::string& filename,
        const std::vector<std::vector<matrixd > >& point_coord,
        const std::vector<std::vector<zjucad::matrix::matrix<zjucad::matrix::matrix<int> > > >& face_index)
{
    std::ofstream ofs(filename.c_str());
    if(!ofs)
        return 1;
    for(size_t i = 0; i < point_coord.size(); ++i) {
        for(size_t j = 0; j < point_coord[i].size(); ++j) {
            for(size_t k = 0; k < point_coord[i][j].size(2); ++k) {
                ofs << "v";
                for(size_t t = 0; t < point_coord[i][j].size(1); ++t) {
                    ofs << " " << point_coord[i][j](t, k);
                }
                ofs << std::endl;
            }
        }
    }
    int point_num = 0;
    for(size_t i = 0; i < point_coord.size(); ++i) {
        for(size_t j = 0; j < point_coord[i].size(); ++j) {
            //ofs << "g group" << i << j % (point_coord[i].size() / 2) << std::endl;
            ofs << "g group" << i << std::endl;
            for(size_t k = 0; k < face_index[i][j].size(1); ++k) {
                if(face_index[i][j][k].size(1) == 3){
                    ofs << "f";
                    for(size_t t = 0; t < face_index[i][j][k].size(1); ++t) {
                        ofs << " " << face_index[i][j][k][t] + 1 + point_num;
                    }
                    ofs << std::endl;
                }else{
                    assert(face_index[i][j][k].size(1) == 4);
                    for(size_t kk = 0; kk < 2; ++kk){
                        size_t det = (kk==0?0:2);
                        ofs << "f";
                        for(size_t t = 0; t < 3; ++t) {
                            ofs << " " << face_index[i][j][k][(t+det)%4] + 1 + point_num;
                        }
                        ofs << std::endl;
                    }
                }
            }
            point_num += point_coord[i][j].size(2);
        }
    }
    return 0;

}
