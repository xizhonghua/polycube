#ifndef FACE2HEX_ADJACENT_H
#define FACE2HEX_ADJACENT_H


#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../hexmesh/hexmesh.h"
#include "../hexmesh/io.h"
//#include <boost/property_tree/ptree.hpp>
//using boost::property_tree::ptree;

class face2hex_adjacent_test: public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(face2hex_adjacent_test);
    CPPUNIT_TEST(get_outside_face_test);
//    CPPUNIT_TEST(get_outside_face_test);
    CPPUNIT_TEST_SUITE_END();
public:
    void setUp();
    void tearDown();
    //void hex_mesh_read_from_wyz_test();
    void get_outside_face_test();
    //void test_edge2tri_adjacent();
private:
    jtf::hexmesh::hexmesh hm;
};


#endif // FACE2HEX_ADJACENT_H
