#include <numeric>

#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/mesh.h>

#include "../common/def.h"
#include "../common/util.h"

using namespace std;
using namespace zjucad::matrix;

int remove_degenerated_tri(int argc, char * argv[])
{
    if(argc != 4){
        cerr << "# [usage] remove_degenerated_tet tri_obj output_obj avg_threshold." << endl;
        return __LINE__;
    }

    matrix<size_t> tri;
    matrix<double> node;

    if(jtf::mesh::load_obj(argv[1], tri, node)){
        return __LINE__;
    }

    matrix<double> face_area(tri.size(2),1);
    for(size_t fi = 0; fi < tri.size(2); ++fi){
        face_area[fi] = jtf::mesh::cal_face_area(tri(colon(),fi), node);
    }

    const double avg_area = std::accumulate(face_area.begin(), face_area.end(),0.0)/face_area.size();
    cerr << "# [info] avg_area " << avg_area << endl;
    vector<bool> is_degenerated(tri.size(2), false);

    const double avg_threshold = atof(argv[3]);
    for(size_t fi = 0; fi < tri.size(2); ++fi){
        if(face_area[fi] < avg_area/avg_threshold){
            is_degenerated[fi] = true;
        }
    }

    matrix<size_t> left_tri(3, count(is_degenerated.begin(), is_degenerated.end(),false));
    for(size_t fi = 0, i = 0; fi < tri.size(2); ++fi){
        if(!is_degenerated[fi]){
            left_tri(colon(),i++) = tri(colon(),fi);
        }
    }

    remove_extra_node(left_tri, node);

    if(jtf::mesh::save_obj(argv[2], left_tri,node)){
        return __LINE__;
    }

    cerr << "# [info] sucess." << endl;
    return 0;
}
