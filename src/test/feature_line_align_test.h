#ifndef FEATURE_LINE_ALIGN_TEST_H
#define FEATURE_LINE_ALIGN_TEST_H

//#include "../common/gauss_elimination.h"
#include "../trimesh/util.h"
#include "../quadmesh/util.h"
#include "../hex_process/hex_process.h"
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

class feature_line_align_test: public CppUnit::TestFixture
{
  //CPPUNIT_TEST_SUITE(gauss_eliminator_test);
  CPPUNIT_TEST_SUITE(feature_line_align_test);
  //CPPUNIT_TEST(gauss_elimination);
  CPPUNIT_TEST(feature_line_align);
  CPPUNIT_TEST_SUITE_END();
public:
    void setUp();
    void tearDown();
    //    void gauss_elimination();
    void feature_line_align();
private:
   // jtf::algorithm::gauss_eliminator ge;
};


#endif // GAUSS_ELIMINATION_TEST_H
