#include "../common/def.h"
#include "../tetmesh/hex_io.h"
#include <fstream>
#include <jtflib/mesh/io.h>

#include "../hex_param/io.h"

using namespace std;

int fit_position(int argc, char *argv[])
{
    if(argc != 3){
        cerr << "# [usage] fit_position tet group_file" << endl;
        return __LINE__;
    }

    matrixst tet;
    matrixd node;
    if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &node, &tet))
        return __LINE__;

    vector<std::vector<size_t> > groups;
    if(load_group_file(argv[2], groups))
        return __LINE__;
    ofstream ofs("fit_position");
    for(size_t gi = 0; gi < groups.size(); ++gi){
        double avg = 0;
        const vector<size_t> & one_group = groups[gi];
        for(size_t ti = 0; ti < one_group.size(); ++ti){
            avg += node[one_group[ti]];
        }
        avg /= one_group.size();
        for(size_t ti = 0; ti < one_group.size(); ++ti)
            ofs << one_group[ti] << " " << avg << endl;
    }

    return 0;
}
