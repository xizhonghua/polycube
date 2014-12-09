#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>

using namespace std;

int add_noise(int argc, char * argv[])
{
    if(argc != 4){
        cerr << "# [usage] add_noise input_obj output_obj weight" << endl;
        return __LINE__;
    }

    jtf::mesh::meshes trim;
    if(jtf::mesh::load_obj(argv[1], trim.mesh_, trim.node_)){
        cerr << "# [error] can not load obj." << endl;
        return __LINE__;
    }

    const double w = atof(argv[3]);

    trim.node_ += zjucad::matrix::rand<double>(trim.node_.size(1), trim.node_.size(2)) * w;

    if(jtf::mesh::save_obj(argv[2], trim.mesh_, trim.node_)){
        cerr << "# [error] can not save obj." << endl;
        return __LINE__;
    }
    return 0;
}
