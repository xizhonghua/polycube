#include "edge2tri_test.h"
#include "../trimesh/trimesh_io.h"
#include "../trimesh/trimesh.h"
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>
#include <iostream>
#include <string>
#include <memory>
using namespace std;
using namespace jtf::io;
using namespace jtf::trimesh;
void edge2tri_test::setUp()
{
    //pt.put("obj.value","../../da2@jtf/fandisk/301k/tet/fandisk-301k.tet.obj");
    pt.put("obj.value","../../dat@jtf/cube/408b/tet/cube-408b.tet.obj");
}

void edge2tri_test::tearDown(){}

void edge2tri_test::test_edge2tri_adjacent()
{
    string obj_file = pt.get<string>("obj.value");
    trimesh trm;
    if(load_from_obj(obj_file.c_str(),trm))
        return ;
    auto_ptr<edge2tri_adjacent> ea(edge2tri_adjacent::create(trm.tri));
    matrixst edge_idx;
    matrixst edge;
    get_boundary_edge_idx(*ea,edge_idx);
    get_boundary_edge(*ea,edge);
    cout << edge_idx << endl;
    cout << edge << endl;
}
