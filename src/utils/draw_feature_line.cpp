#include <vector>
#include <deque>

#include "../common/visualize_tool.h"
#include "../common/def.h"
#include <jtflib/mesh/io.h>

using namespace std;

int draw_feature_line(int argc, char * argv[])
{
    if(argc != 5){
        cerr << "# [usage] draw_feature_line obj fl radius output_fl" << endl;
        return __LINE__;
    }
    matrixst mesh;
    matrixd node;
    if(jtf::mesh::load_obj(argv[1], mesh, node)){
        return __LINE__;
    }

    vector<deque<pair<size_t,size_t> > > feature_lines;

    if(jtf::mesh::load_feature_line(argv[2], feature_lines))
        return __LINE__;

    dump_singularity_to_cylinder(argv[4], node, feature_lines, atof(argv[3]));
    return 0;
}
