#include "../common/IO.h"
#include <jtflib/mesh/io.h>
#include "../common/def.h"
#include <sstream>

using namespace std;
using namespace zjucad::matrix;

static int load_surface_patch(
        const char * surface_patch_file,
        std::vector<matrixst> &patch_faces)
{
    ifstream ifs_patch(surface_patch_file);
    if(ifs_patch.fail()){
        cerr << "# [error] can not open patch surface file." << endl;
        return __LINE__;
    }

    size_t patch_num = 0;
    ifs_patch >> patch_num;

    cerr << "# [info] patch number " << patch_num << endl;

    if(patch_num == 0) return 0;

    patch_faces.clear();

    patch_faces.resize(patch_num);
    size_t face_num = 0, face_idx, trash_1 ,trash_2;

    for(size_t pi = 0; pi < patch_num; ++pi){
        ifs_patch >> face_num >> trash_1;
        patch_faces[pi].resize(face_num,1);

        for(size_t fi = 0; fi < face_num; ++fi){
            ifs_patch >> patch_faces[pi][fi];
        }
        for(size_t li = 0; li < trash_1; ++li)
            ifs_patch >> trash_2;
    }
    return 0;
}


int dump_surface_patch(int argc, char * argv[])
{
    if(argc != 3){
        cerr << "# [usage] dump_surface_patch obj patch output_patch_obj" << endl;
        return __LINE__;
    }

    matrixst mesh;
    matrixd node;

    if(jtf::mesh::load_obj(argv[1], mesh, node)){
        return __LINE__;
    }

    vector<matrixst> patch_faces;
    if(load_surface_patch(argv[2],  patch_faces))
        return __LINE__;

    matrixst one_patch;
    string obj_name;
    for(size_t pi = 0; pi < patch_faces.size(); ++pi){
        ostringstream os ;
        os << "surface_patch_" << pi << ".obj";
        one_patch.resize(3, patch_faces[pi].size());
        for(size_t fi = 0; fi < one_patch.size(2); ++fi){
            one_patch(colon(),fi) = mesh(colon(),patch_faces[pi][fi]);
        }
        obj_name = os.str();
        if(jtf::mesh::save_obj(obj_name.c_str(), one_patch,node))
            return __LINE__;
    }
    return 0;
}
