#ifndef POLYCUBE_PARAM_TEST_H
#define POLYCUBE_PARAM_TEST_H

#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <boost/property_tree/ptree.hpp>
//#include "../hex_ui/prog.h"
using boost::property_tree::ptree;

class polycube_param_test : public CppUnit::TestFixture
{
        CPPUNIT_TEST_SUITE(polycube_param_test);
        CPPUNIT_TEST(test_frame);
        //CPPUNIT_TEST(weight_lbfgs_zyz_scale_test);
        //CPPUNIT_TEST(weight_lbfgs_zyz_subdivision_test);
        CPPUNIT_TEST_SUITE_END();
 public:
    void setUp();
    void tearDown();
    void test_frame();
private:
    ptree pt;

};

#endif // POLYCUBE_PARAM_TEST_H
