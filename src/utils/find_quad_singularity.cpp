#include <iostream>
#include <fstream>
#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/io.h>

#include "../common/visualize_tool.h"


using namespace std;
using namespace zjucad::matrix;

int dump_quad_face_singularity_point(
    const matrix<size_t> & outside_face,
    const matrix<double> & node,
    const double radius,
    const char * sphere_file,
    const string singularity_point_str_pref);

int find_quad_singularity(int argc, char * argv[])
{
  if(argc != 4 && argc != 5 && argc != 6){
    cerr << "# [usage] find_quad_singularity quad_obj output_singularity.obj input_unit_sphere [radius] [feature_line] " << endl;
    return __LINE__;
  }

  matrix<size_t> quad;
  matrix<double> nodes;

  if(jtf::mesh::load_obj(argv[1], quad, nodes))
    return __LINE__;

  cerr << "# [info] read quad seccess." << endl;

  double radius  = 0.002;
  if(argc > 4)
    radius = atof(argv[4]);

  string singularity_point_str = argv[2];
  singularity_point_str += "_point_";
  dump_quad_face_singularity_point(quad, nodes, radius, argv[3], singularity_point_str);

  if(argc > 5){
    string feature_line_str = argv[2];
    feature_line_str += "_line.obj";
    vector<deque<pair<size_t,size_t> > > lines;
    jtf::mesh::load_feature_line(argv[5], lines);
    cerr << "# [info] read feature line seccess." << endl;
    dump_singularity_to_cylinder(feature_line_str.c_str(), nodes, lines,radius);
  }
  cerr << "# [info] seccess." << endl;
  return 0;
}
