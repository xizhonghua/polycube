#ifndef FIND_SINGULARITIES_INNER_TEST_H
#define FIND_SINGULARITIES_INNER_TEST_H
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <boost/property_tree/ptree.hpp>

#include "../hex_ui/prog.h"
#include "../hex_param/hex_param.h"

class find_singularities_test: public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(find_singularities_test);
    CPPUNIT_TEST(test_find_singularities2_inner);
    CPPUNIT_TEST_SUITE_END();
public:
    void setUp();
    void tearDown();
    void test_find_singularities2_inner();
public:
    ptree pt;
};

#endif // FIND_SINGULARITIES2_INNER_TEST_H
