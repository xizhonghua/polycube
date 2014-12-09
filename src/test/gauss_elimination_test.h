#ifndef GAUSS_ELIMINATION_TEST_H
#define GAUSS_ELIMINATION_TEST_H

#include "../common/gauss_elimination.h"

#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

class gauss_eliminator_test: public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(gauss_eliminator_test);
    CPPUNIT_TEST(gauss_elimination1);
//    CPPUNIT_TEST(gauss_elimination2);
//    CPPUNIT_TEST(gauss_elimination8);
//    CPPUNIT_TEST(gauss_elimination4);
    CPPUNIT_TEST_SUITE_END();
public:
    void setUp();
    void tearDown();
    void gauss_elimination();
    void gauss_elimination2();
    void gauss_elimination8();
    void gauss_elimination4();
    void gauss_elimination1();
    void gauss_elimination3();
private:
   // jtf::algorithm::gauss_eliminator ge;
};


#endif // GAUSS_ELIMINATION_TEST_H
