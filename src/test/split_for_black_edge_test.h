#ifndef SPLIT_FOR_BLACK_EDGE_TEST_H
#define SPLIT_FOR_BLACK_EDGE_TEST_H
#include <memory>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

#include <boost/unordered_map.hpp>
#include <zjucad/matrix/matrix.h>
#include "../hex_param/singularity_adjustment.h"
#include "../hex_param/find_singularities.h"

class split_for_black_edge_test : public CppUnit::TestFixture
{
                CPPUNIT_TEST_SUITE(split_for_black_edge_test);
                CPPUNIT_TEST(split_for_black_edge_test_with_trivial_model);
                CPPUNIT_TEST_SUITE_END();
public:
                void setUp();
                void tearDown(){}
                void split_for_black_edge_test_with_trivial_model();
private:
  matrixst tet_;
  matrixd node_;
  boost::unordered_map<std::pair<size_t,size_t>,size_t> inner_face_jump_type;
  jtf::tetmesh::one_ring_tet_at_edge ortae;
  std::auto_ptr<face2tet_adjacent> fa;
  boost::property_tree::ptree pt;
};

#endif // SPLIT_FOR_BLACK_EDGE_TEST_H
