#ifndef EXTRAC_LOOP_FROM_UNDIRECTED_EDGES_TEST_H
#define EXTRAC_LOOP_FROM_UNDIRECTED_EDGES_TEST_H

#include "../common/extract_loop_from_undirected_edges.h"
#include <memory>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>


class extract_loop_from_undirected_edges_test : public CppUnit::TestFixture
{
        CPPUNIT_TEST_SUITE(extract_loop_from_undirected_edges_test);
        CPPUNIT_TEST(extract_loops);
        CPPUNIT_TEST_SUITE_END();
 public:
    void setUp(){}
    void tearDown(){}
    void extract_loops();
};

#endif // EXTRAC_LOOP_FROM_UNDIRECTED_EDGES_TEST_H
