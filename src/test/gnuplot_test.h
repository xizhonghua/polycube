#ifndef _TEST_GNUPLOT_TEST_H_
#define _TEST_GNUPLOT_TEST_H_

#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

class gnuplot_test: public CppUnit::TestFixture
{
  //CPPUNIT_TEST_SUITE(gauss_eliminator_test);
  CPPUNIT_TEST_SUITE(gnuplot_test);
  //CPPUNIT_TEST(gauss_elimination);
  CPPUNIT_TEST(gnuplot_test_0);
  CPPUNIT_TEST_SUITE_END();
public:
    void setUp(){}
    void tearDown(){}
    void gnuplot_test_0();
    virtual ~gnuplot_test(){}
private:
};


#endif // GNUPLOT
