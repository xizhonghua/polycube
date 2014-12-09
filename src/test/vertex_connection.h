#ifndef VERTEX_CONNECTION_H
#define VERTEX_CONNECTION_H

#include "../common/vertex_connection.h"
#include <memory>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>


class vertex_connection_test : public CppUnit::TestFixture
{
        CPPUNIT_TEST_SUITE(vertex_connection_test);
        CPPUNIT_TEST(get_shortest_path_test);
        CPPUNIT_TEST_SUITE_END();
 public:
    void setUp();
    void tearDown();
    void get_shortest_path_test();
private:
    std::auto_ptr<vertex_connection<DIRECT> > vc;
};

#endif // VERTEX_CONNECTION_H
