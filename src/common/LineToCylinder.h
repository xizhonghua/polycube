#ifndef LINETOCYLINDER_H
#define LINETOCYLINDER_H

#include <zjucad/matrix/matrix.h>
#include "../common/def.h"
#include <vector>
//const int segment_num = 16;
typedef std::vector<matrixd > line_type;

int line_to_cylinder(const matrixd& line_point,
					 matrixd& point_coord,
					 zjucad::matrix::matrix<zjucad::matrix::matrix<int> >& face_index,
					 double radius);

int convertToCylinder(
    const std::vector<std::vector<line_type> >& line_segments,
    std::vector<std::vector<matrixd > >& point_coord,
    std::vector<std::vector<zjucad::matrix::matrix<zjucad::matrix::matrix<int> > > >& face_index,
    double radius);

int saveCylinder(
    const std::string& filename,
    const std::vector<std::vector<matrixd > >& point_coord,
    const std::vector<std::vector<zjucad::matrix::matrix<zjucad::matrix::matrix<int> > > >& face_index);

#endif
