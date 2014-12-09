#include "find_singularities2_test.h"
#include "../tetmesh/tetmesh.h"
#include "../hex_ui/prog.h"
#include "../common/IO.h"

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <string>

using namespace std;
using namespace zjucad::matrix;

void find_singularities2_test::tearDown()
{}

void find_singularities2_test::setUp()
{
    pt.put("prog.value","find_singularities2_inner");
    pt.put("tet.value","../../dat@jtf/sphere/80k/tet/sphere-80k.tet.obj");
    pt.put("zyz.value","../../dat@jtf/sphere/80k/test/sphere-80k.test.zyz");
    pt.put("package.value","alglib");
    pt.put("alg.value","lbfgs");
    pt.put("lbfgs-len.value",7);
    pt.put("iter.value",1000);
    pt.put("align_w.value",1);
    pt.put("fix_w.value",1);
    frame_inner(pt);
}

void find_singularities2_test::test_find_singularities2()
{
    tetmesh tm;
    matrixst tri;
    pt.put("tet.desc","tet file");
    CPPUNIT_ASSERT_MESSAGE("# read tet fail.", !tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm.node, &tm.tet, &tri));

    auto_ptr<face2tet_adjacent> fa(face2tet_adjacent::create(tm.tet));

    matrixd tet_zyz;

    ifstream ifs(pt.get<string>("zyz.value").c_str(), ifstream::binary);
    CPPUNIT_ASSERT_MESSAGE("# not a zyz file.", !read_matrix(ifs, tet_zyz));
}
