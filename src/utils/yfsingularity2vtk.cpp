#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
#include "../common/visualize_tool.h"
#include "../common/util.h"

#include <fstream>

using namespace std;

int load_yfsingularity(const char * singularity_file,
                       vector<pair<size_t,size_t> > & singularity)
{
    ifstream ifs(singularity_file);
    if(ifs.fail()){
        cerr << "# [error] can not open singularity file." << endl;
        return __LINE__;
    }

    singularity.clear();
    size_t p0,p1;
    while(!ifs.eof()){
        ifs >> p0 >> p1;
        if(ifs.eof()) break;
        singularity.push_back(make_pair(p0,p1));
    }

    return 0;
}

int yfsingularity2vtk(int argc, char * argv[])
{
    if(argc != 4){
        cerr << "# [usage] yfsingularity2vtk tet singularity vtk." << endl;
        return __LINE__;
    }

    matrixst tet;
    matrixd node;

    if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &node, &tet))
        return __LINE__;

    vector<pair<size_t,size_t> > singularity;
    if(load_yfsingularity(argv[2],  singularity))
        return __LINE__;

    vector<deque<pair<size_t,size_t> > > singularity_chain;
    jtf::util::extract_chain_from_edges(singularity, singularity_chain);

    dump_singularity_to_vtk(argv[3],node, singularity_chain);
    cerr << "# [info] success." << endl;
    return 0;
}
