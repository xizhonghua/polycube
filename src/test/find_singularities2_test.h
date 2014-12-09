#ifndef FIND_SINGULARITIES2_TEST_H
#define FIND_SINGULARITIES2_TEST_H

#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

class find_singularities2_test: public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(find_singularities2_test);
    CPPUNIT_TEST(test_find_singularities2);
    CPPUNIT_TEST_SUITE_END();
public:
    void setUp();
    void tearDown();
    void test_find_singularities2();
private:
    ptree pt;
};

#endif // FIND_SINGULARITIES2_TEST_H
