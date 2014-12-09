#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include "../common/IO.h"
#include "../tetmesh/tetmesh.h"
#include "../common/zyz.h"

#include "../common/transition_type.h"
#include "../common/transition.h"
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>

using namespace std;
using namespace zjucad::matrix;

int dump_jump_type(int argc,char *argv[])
{
    if(argc != 4){
        cerr << "# tet zyz output_jump_type." << endl;
        return __LINE__;
    }
    jtf::mesh::meshes tm;
    matrixd zyz;

    if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1],&tm.node_, &tm.mesh_)) {
        cerr << "# read tet fail." << endl;
        return __LINE__;
    }

    ifstream ifs(argv[2],ifstream::binary);
    if(!ifs.fail())
    {
        jtf::mesh::read_matrix(ifs,zyz);
        if(!zyz.size()) {
            cerr << "#read zyz fail." << endl;
            return __LINE__;
        }else if(zyz.size(2) != tm.mesh_.size(2)){
            cerr << "# error zyz format." << endl;
            return __LINE__;
        }
    }else{
        cerr << "# open zyz file fail." << endl;
        return __LINE__;
    }

    matrix<matrixd > frame_inner(zyz.size(2));
    for(size_t t = 0; t < zyz.size(2); ++t){
        frame_inner[t].resize(3,3);
        zyz_angle_2_rotation_matrix1(&zyz(0,t),&frame_inner[t](0,0));
    }

    unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
    map<pair<size_t,size_t>,size_t> tet_jump_type;
    matrixd rot(3,3);
    for(size_t t = 0;t < fa->face2tet_.size(); ++t){
        const pair<size_t,size_t> & tet_pair = fa->face2tet_[t];
        if(tet_pair.first == -1 || tet_pair.second == -1) continue;
         get_best_alignment(&frame_inner[tet_pair.first](0,0),
                            &frame_inner[tet_pair.second](0,0),
                            &rot[0]);
         if(norm(rot - eye<double>(3)) < 1e-8) continue;
         const size_t type_A2B = type_transition1(rot);
         const size_t type_B2A = type_transition1(trans(rot)); // the inverse matrix of rot from B-->A
         assert(  type_A2B != 25 && type_B2A != 25 // not invalid and not identity
               && type_A2B != 9 && type_B2A != 9);
         tet_jump_type[make_pair(tet_pair.first,tet_pair.second)] = type_A2B;
         tet_jump_type[make_pair(tet_pair.second,tet_pair.first)] = type_B2A;
    }

    {// dump out the tet_pair_type
        ofstream ofs(argv[3]);
        if(ofs.fail()){
            cerr << "# error. can not open file." << endl;
            return __LINE__;
        }
        typedef map<pair<size_t,size_t>,size_t>::const_iterator mci;
        for(mci ci = tet_jump_type.begin(); ci != tet_jump_type.end(); ++ci){
            ofs << ci->first.first << " " << ci->first.second << " " << ci->second << endl;
        }

    }
    return 0;
}
