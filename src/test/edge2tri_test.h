#ifndef EDGE2TRI_TEST_H
#define EDGE2TRI_TEST_H

#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

class edge2tri_test: public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(edge2tri_test);
    CPPUNIT_TEST(test_edge2tri_adjacent);
    CPPUNIT_TEST_SUITE_END();
public:
    void setUp();
    void tearDown();
    void test_edge2tri_adjacent();
private:
    ptree pt;
};
#endif // EDGE2TRI_TEST_H
