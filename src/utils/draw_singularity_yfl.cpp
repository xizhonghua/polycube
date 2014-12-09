#include <fstream>
#include <vector>

#include "../common/def.h"
#include "../tetmesh/hex_io.h"
#include "../common/visualize_tool.h"
#include "../common/util.h"

using namespace std;
using namespace zjucad::matrix;

int load_singularity_list(const char * filename,
                          vector<pair<size_t,size_t> > & singularity_list)
{
    ifstream ifs(filename);

    if(ifs.fail()){
        cerr << "# [error] can not open singularity list file." << endl;
        return __LINE__;
    }

    singularity_list.clear();

    size_t p0,p1;
    while(!ifs.eof()){
        ifs >> p0 >> p1;
        if(ifs.eof()) break;
        singularity_list.push_back(make_pair(p0,p1));
    }

    return 0;
}


int draw_singularity_yfl(int argc, char * argv[])
{
    if( argc != 3){
        cerr << "# [usage] draw_singularity_yfl tet singularity_list." << endl;
        return __LINE__;
    }

    matrixst tet;
    matrixd node;

    if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &node, &tet))
        return __LINE__;

    vector<pair<size_t,size_t> > singularity_list;
    if(load_singularity_list(argv[2], singularity_list))
        return __LINE__;

    vector<deque<pair<size_t,size_t> > > singularity_chains;
    jtf::util::extract_chain_from_edges(singularity_list, singularity_chains);

    const double radius = 0.003;
    dump_singularity_to_cylinder("singularity.obj", node, singularity_chains, radius);

    return 0;
}
