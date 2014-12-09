#ifndef INIT_CUT_INNER_TEST_H
#define INIT_CUT_INNER_TEST_H

#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <boost/property_tree/ptree.hpp>
#include "../hex_ui/prog.h"
using boost::property_tree::ptree;

class init_cut_inner_test : public CppUnit::TestFixture
{
                CPPUNIT_TEST_SUITE(init_cut_inner_test);
                //CPPUNIT_TEST(test_init_cut_inner);
                //CPPUNIT_TEST(test_init_value_lsq);
                CPPUNIT_TEST(test_lsq_value_and_lbfgs_result);
                CPPUNIT_TEST_SUITE_END();
public:
                void setUp();
                void tearDown();
                void test_init_cut_inner();
                void test_init_value_lsq();
                void test_lsq_value_and_lbfgs_result(); // the lsq result shoudl be different with the lbfgs result
private:

                ptree pt;
};

#endif // INIT_CUT_INNER_TEST_H
