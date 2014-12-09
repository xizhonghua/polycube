#ifndef ROT2EULER_ANGLE_TEST_H
#define ROT2EULER_ANGLE_TEST_H

#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../numeric/util.h"


class rot2euler_angle_test : public CppUnit::TestFixture
{
        CPPUNIT_TEST_SUITE(rot2euler_angle_test);
        //CPPUNIT_TEST(convert_rot_to_euler_zyx_test);
        CPPUNIT_TEST(delete_axis_test);
        CPPUNIT_TEST_SUITE_END();
public:
        void setUp(){}
        void tearDown(){}
        void convert_rot_to_euler_zyx_test();
        void delete_axis_test();
};


#endif // ROT2EULER_ANGLE_TEST_H
