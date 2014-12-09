//#include "frame_test.h"
//#include "find_singularities_inner_test.h"
//#include "init_cut_inner_test.h"
//#include "weight_test.h"
//#include "edge2tri_test.h"
//#include "polycube_param_test.h"
//#include "vertex_connection.h"
//#include "endianness_test.h"
//#include "face2hex_adjacent.h"
//#include "rot2euler_angle_test.h"
//#include "split_for_black_edge_test.h"
#include "gauss_elimination_test.h"
//#include "feature_line_align_test.h"
//#include "hex2obj_test.h"
//#include "extrac_loop_from_undirected_edges_test.h"
//#include "gnuplot_test.h"
//#include "test_my_mesh.h"
//CPPUNIT_TEST_SUITE_REGISTRATION(frame_test); // 自动注册测试包
//CPPUNIT_TEST_SUITE_REGISTRATION(find_singularities_test);
//CPPUNIT_TEST_SUITE_REGISTRATION(init_cut_inner_test);
//CPPUNIT_TEST_SUITE_REGISTRATION(weight_test);
//CPPUNIT_TEST_SUITE_REGISTRATION(polycube_param_test);
//CPPUNIT_TEST_SUITE_REGISTRATION(vertex_connection_test);
//CPPUNIT_TEST_SUITE_REGISTRATION(gnuplot_test);
//CPPUNIT_TEST_SUITE_REGISTRATION(mymesh_test);
CPPUNIT_TEST_SUITE_REGISTRATION(gauss_eliminator_test);
//CPPUNIT_TEST_SUITE_REGISTRATION(feature_line_align_test);
//CPPUNIT_TEST_SUITE_REGISTRATION(hex2obj_test);
//CPPUNIT_TEST_SUITE_REGISTRATION(extract_loop_from_undirected_edges_test);
//CPPUNIT_TEST_SUITE_REGISTRATION(Endianness_test);
//CPPUNIT_TEST_SUITE_REGISTRATION(face2hex_adjacent_test);
//CPPUNIT_TEST_SUITE_REGISTRATION(rot2euler_angle_test);
//CPPUNIT_TEST_SUITE_REGISTRATION(split_for_black_edge_test);

//CPPUNIT_TEST_SUITE_REGISTRATION(edge2tri_test);

int main()
{
    CppUnit::TestResult r; 
    CppUnit::TestResultCollector rc;
    r.addListener(&rc); // 准备好结果收集器 

    CppUnit::TestRunner runner; // 定义执行实体
    runner.addTest(CppUnit::TestFactoryRegistry::getRegistry().makeTest());
    runner.run(r); // 运行测试

    CppUnit::TextOutputter o(&rc, std::cout);
    o.write(); // 将结果输出

    return rc.wasSuccessful() ? 0 : -1;
}
