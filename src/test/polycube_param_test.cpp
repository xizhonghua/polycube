#include "polycube_param_test.h"
#include <string.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>

#include <omp.h>

#include "../common/zyz.h"
#include "../common/IO.h"
#include "../common/util.h"
#include "../common/vtk.h"
#include "../common/check_jac.h"
#include "../common/transition.h"

#include "../spherical_harmonics/rot_cubic_f_SH.h"

#include <hjlib/sparse/sparse.h>
//#include <hjlib/arg_opts/arg_opts.h>
#include <hjlib/function/function.h>
//#include <hjlib/optimizer/optimizer.h>
#include <zjucad/optimizer/optimizer.h>

#include <hjlib/math/blas_lapack.h>

#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/lapack.h>

#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"

#include "../hex_frame_opt/hex_frame_opt.h"

#include "../hex_param/hex_param.h"
#include "../hex_param/common.h"

#include <boost/property_tree/ptree.hpp>
#include <zjucad/ptree/ptree.h>

#include <minpack.h>


using namespace std;
using namespace zjucad::matrix;
using namespace hj::function;

void polycube_param_test::setUp()
{
    pt.put("tet.value","../../dat@jtf/sphere/5k/tet/sphere-5k.tet");
    pt.put("prog.value","polycube_param");
    //pt.put("init.value","../../dat@jtf/sphere/5k/polycube_zyz/sphere-5k-r390-11-29-14-22@1-1e3-1000.inner.init-polycube.with_surface_normal.zyz");
    pt.put("param_algin_fram_w.value",1);
    pt.put("frame_align_normal_w.value",1000);
    pt.put("param_planar_w.value",1);
    pt.put("package.value","alglib");
    pt.put("alg.value","lbfgs");
    pt.put("lbfgs-len.value",7);
    pt.put("iter.value",1000);
    pt.put("comp_zyz.value","../../dat@jtf/sphere/5k/polycube_zyz/sphere-5k-r390-11-29-15-05@1-1000-1e3-1000.inner-polycube.zyz");
}

void polycube_param_test::tearDown()
{}

void polycube_param_test::test_frame()
{
    tetmesh tm;
    matrixd zyz;
    sym_frame_opt sys;
    CPPUNIT_ASSERT(!tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm.node, &tm.tet));

    pt.put("init.desc","need inner frame zyz file");
    if(zjucad::has("init.value",pt))
    {
        ifstream ifs(pt.get<string>("init.value").c_str(),ifstream::binary);
        CPPUNIT_ASSERT_MESSAGE("#open init zyz fail.",!ifs.fail());
        read_matrix(ifs,zyz);
        CPPUNIT_ASSERT_MESSAGE("#read zyz fail.", zyz.size());
        CPPUNIT_ASSERT_EQUAL_MESSAGE("# error zyz format",tm.tet.size(2),zyz.size(2));
    }
    else
        zyz = zeros<double>(3,tm.tet.size(2)) + 3.1415926/4;

    // weight[0]:param align frame
    // weight[1]:surface normal alignment
    // weigth[2]:planar alignment
    const double weight[3] = {
            pt.get<double>("param_algin_fram_w.value"),
            pt.get<double>("frame_align_normal_w.value"),
            pt.get<double>("param_planar_w.value")
    };/*add default value 1e6*/
    matrixd stiff = ones<double>(tm.tet.size(2), 1); // do not know how to use the stiff

    if(zjucad::has("stiff.value",pt)){
        ifstream ifs(pt.get<string>("stiff.value").c_str(), ifstream::binary);
        CPPUNIT_ASSERT_MESSAGE("#open stiff fail.",!ifs.fail());
        read_matrix(ifs, stiff);
    }

    matrixd frame_uvw = zeros<double>(3, tm.node.size(2) + tm.tet.size(2));
    frame_uvw(colon(),colon(0,tm.tet.size(2)-1)) = zyz;

    frame_uvw(colon(),colon(tm.tet.size(2), tm.tet.size(2) + tm.node.size(2)-1)) = tm.node;
    //matrixd uvw_frame = zyz;
    polycube_param(pt, sys, tm.tet, tm.node, frame_uvw, weight, stiff);

   // const matrixd & frame = uvw_frame(colon(),colon(0, tm.tet.size(2) -1));

    matrixd comp_zyz;
    ifstream ifs(pt.get<string>("comp_zyz.value").c_str(),ifstream::binary);
    CPPUNIT_ASSERT_MESSAGE("#open comp_zyz fail.",!ifs.fail());
    read_matrix(ifs,comp_zyz);
    CPPUNIT_ASSERT_MESSAGE("#read zyz fail.", comp_zyz.size());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("# error zyz format",tm.tet.size(2),comp_zyz.size(2));

    const matrixd & frame = frame_uvw(colon(),colon(0,tm.tet.size(2)-1));

    CPPUNIT_ASSERT(fabs(norm(comp_zyz - frame)) < 1e-8);
}
