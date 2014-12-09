#include "init_cut_inner_test.h"
#include "../common/IO.h"
#include "../tetmesh/tetmesh.h"
#include "../hex_ui/function_term.h"

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>


using namespace std;
using namespace zjucad::matrix;
using namespace hj::function;

void init_cut_inner_test::setUp()
{
  pt.put("prog.value","init_cut_inner");
  pt.put("tet.value","../../dat@jtf/sphere/80k/tet/sphere-80k.tet");
  pt.put("obj.value","../../dat@jtf/sphere/80k/tet/sphere-80k.tet.obj");
  pt.put("fix_w.value",1e3);
  //pt.put("obj.value","../../dat@jtf/sphere/80k/tet/sphere-80k.tet.obj");
}

void init_cut_inner_test::test_init_cut_inner()
{
  pt.put("align_w.value",0.01);
  pt.put("output.value","../../dat@jtf/sphere/80k/test/sphere-80k@0.01-1e3-7.init.with_surface.zyz");
  init_zyz_inner(pt);

  pt.put("align_w.value",1e6);
  pt.put("output.value","../../dat@jtf/sphere/80k/test/sphere-80k@1e6-1e3-7.init.with_surface.zyz");
  init_zyz_inner(pt);

  matrixd zyz0;
  matrixd zyz1;
  ifstream ifs0,ifs1;
  ifs0.open("../../dat@jtf/sphere/80k/test/sphere-80k@0.01-1e3-7.init.with_surface.zyz");
  //if(ifs.fail())
  CPPUNIT_ASSERT (!ifs0.fail());
  read_matrix(ifs0, zyz0);
  CPPUNIT_ASSERT( zyz0.size()!= 0);


  ifs1.open("../../dat@jtf/sphere/80k/test/sphere-80k@1e6-1e3-7.init.with_surface.zyz");
  CPPUNIT_ASSERT (!ifs1.fail());
  read_matrix(ifs1, zyz1);
  CPPUNIT_ASSERT( zyz1.size()!= 0);

  CPPUNIT_ASSERT_EQUAL( zyz0.size(), zyz1.size() );
  CPPUNIT_ASSERT_EQUAL( zyz0.size(1), zyz1.size(1) );
  CPPUNIT_ASSERT_EQUAL( zyz0.size(2), zyz1.size(2) );
  vector<double> difference(zyz0.size(2));
  for(size_t ti = 0; ti < zyz0.size(2); ++ti)
  {
        difference[ti] = fabs(norm(zyz0(colon(), ti) - zyz1(colon(), ti)));
  }
  //copy(difference.begin(),difference.end(),ostream_iterator<size_t>(cerr, " "));
  CPPUNIT_ASSERT_MESSAGE("strange: the difference is zero!" ,   *max_element(difference.begin(),difference.end()) > 1e-8);
    cerr << "# max difference " << *max_element(difference.begin(),difference.end()) << endl;
}

void init_cut_inner_test::test_init_value_lsq()
{
    tetmesh tm;
    int flag = tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm.node, &tm.tet);
    CPPUNIT_ASSERT(flag == 0);

    pt.put("align_w.value",1e6);
    pt.put("package.value","hj");
    pt.put("iter.value",1000);
    auto_ptr<face2tet_adjacent> fa(face2tet_adjacent::create(tm.tet));

    std::vector<double> volume_of_tets(tm.tet.size(2));
    for(size_t ti = 0; ti < volume_of_tets.size(); ++ti)
            volume_of_tets[ti] = fabs(tet_volume(tm.node(colon(), tm.tet(colon(), ti))));

    const double volume_of_object = accumulate(volume_of_tets.begin(), volume_of_tets.end(), 0.0);

    cerr << "# begin to assumble the function" << endl;
    cerr << "# align_w = " << pt.get<double>("align_w.value") << endl;
    vector<boost::shared_ptr<function> > funcs0,funcs1;
    for(size_t i = 0; i < fa->faces_.size(); ++i) {// for each element, i.e. inner face
        const size_t &vi = fa->face2tet_[i].first;
        const size_t &vj = fa->face2tet_[i].second;

        if(fa->is_outside_face(fa->face2tet_[i]))
        {
            const size_t v = (vi == -1)?vj:vi;
            const double weight = sqrt(volume_of_tets[v]/(4 * volume_of_object) * pt.get<double>("align_w.value"));
            const vector<size_t> &vert_idx = fa->faces_[i];
            const matrixd edge[2] = {
                    tm.node(colon(), vert_idx[1]) - tm.node(colon(), vert_idx[0]),
                    tm.node(colon(), vert_idx[2]) - tm.node(colon(), vert_idx[0]),
            };
            matrixd normal = cross(edge[0], edge[1]);
            const double len = norm(normal);
            //cerr << "# i " << i << "normal = " << len << endl;
            if(len > 1e-8) {
                    normal /= len;
            }
            else {
                    cerr << "# degenerated surface triangle at node: " << i << endl;
                    copy(vert_idx.begin(), vert_idx.end(), ostream_iterator<size_t>(cerr, " "));
                    cerr << edge[0] << edge[1] << endl;
            }
            funcs0.push_back(boost::shared_ptr<function>(
                                                    new align_sh_inner(tm.tet.size(2), v, &normal[0], weight)));
        }
        else
        {
            funcs0.push_back(boost::shared_ptr<function>(
                                new inner_smooth(tm.tet,vi,vj,*fa,
                                                 sqrt((volume_of_tets[vi] + volume_of_tets[vj] )/(4*volume_of_object)))));
            //funcs1.push_back(funcs0.back());
        }
    }

    std::auto_ptr<hj::function::function> all_f0;
    all_f0.reset(new_catenated_function<double, int32_t>(funcs0));
    matrixd residual0(all_f0->dim_of_f());
    matrixd sh0 = zeros<double>(9,tm.tet.size(2)) ;
    cerr << "# solve full functions" << endl;
    cerr << "# functions num = " << funcs0.size() << endl;
    cerr << "# resideual0 size = " << residual0.size() << endl;
    zjucad::optimize(*all_f0, sh0, residual0, pt);
//    cerr << "# solve only align functions" << endl;
//    zjucad::optimize(*all_f1, sh1, residual1, pt);
//    for(size_t shi = 0; shi < sh0.size(2); ++shi)
//    {
//        cerr << "# ";
//        for(size_t t = 0; t < 9; ++t)
//            cerr << sh0(t,shi) << " " ;
//        cerr << endl;
//    }
    double sh_diff = norm(sh0);
    double res_diff = norm(residual0);
    CPPUNIT_ASSERT( fabs(sh_diff) > 1e-8);
    CPPUNIT_ASSERT( fabs(res_diff) > 1e-8);
}

void init_cut_inner_test::test_lsq_value_and_lbfgs_result()
{
    pt.put("prog.value","init_cut_inner");
    pt.put("tet.value","../../dat@jtf/fandisk/301k/tet/fandisk-301k.tet");
    pt.put("package.value","alglib");
    pt.put("alg.value","lbfgs");
    pt.put("lbfgs-len.value",7);
    pt.put("iter.value",100);
    pt.put("align_w.value",10);

    pt.put("output.value","../../dat@jtf/fandisk/301k/test/test_init_fandisk_zyz_inner.zyz");
    init_zyz_inner(pt);
    pt.put("init.value","../../dat@jtf/fandisk/301k/test/test_init_fandisk_zyz_inner.zyz");
    pt.put("zyz.value","../../dat@jtf/fandisk/301k/test/test_fandisk_zyz.zyz");
    frame_inner(pt);

    matrixd zyz_init;
    matrixd zyz_result;
    ifstream init_io,final_io;
    init_io.open(pt.get<string>("output.value").c_str());
    final_io.open(pt.get<string>("zyz.value").c_str());
    CPPUNIT_ASSERT_MESSAGE("can not read ../../dat@jtf/fandisk/301k/test/test_init_fandisk_zyz_inner.zyz!" ,  !init_io.fail() );
    CPPUNIT_ASSERT_MESSAGE("can not read ../../dat@jtf/fandisk/301k/test/test_fandisk_zyz.zyz!" ,  !final_io.fail() );

    read_matrix(init_io,zyz_init);
    read_matrix(final_io,zyz_result);
    CPPUNIT_ASSERT_EQUAL(zyz_init.size(), zyz_result.size());
    CPPUNIT_ASSERT_EQUAL(zyz_init.size(1), zyz_result.size(1));
    CPPUNIT_ASSERT_EQUAL(zyz_init.size(2), zyz_result.size(2));
    matrixd difference = zyz_init - zyz_result;
    CPPUNIT_ASSERT(norm(difference) > 1e-8);
    for(size_t t = 0; t < difference.size(2); ++t)
    {
        CPPUNIT_ASSERT(norm(difference(colon(), t)) > 1e-8);
    }

}
void init_cut_inner_test::tearDown(){}
