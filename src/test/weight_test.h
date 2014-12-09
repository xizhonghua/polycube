#ifndef WEIGHT_TEST_H
#define WEIGHT_TEST_H

#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <boost/property_tree/ptree.hpp>
//#include "../hex_ui/prog.h"
using boost::property_tree::ptree;

class weight_test : public CppUnit::TestFixture
{
        CPPUNIT_TEST_SUITE(weight_test);
        CPPUNIT_TEST(weight_init_zyz_scale_test);
        CPPUNIT_TEST(weight_lbfgs_zyz_scale_test);
        //CPPUNIT_TEST(weight_lbfgs_zyz_subdivision_test);
        CPPUNIT_TEST_SUITE_END();
 public:
    void setUp();
    void tearDown();
    void weight_init_zyz_scale_test();
    void weight_lbfgs_zyz_scale_test();
    //void weight_init_zyz_subdivision_test();
    void weight_lbfgs_zyz_subdivision_test();
private:
    ptree pt;
};

#endif // WEIGHT_TEST_H
