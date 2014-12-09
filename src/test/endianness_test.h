#ifndef ENDIANNESS_TEST_H
#define ENDIANNESS_TEST_H

#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../common/byte_config.h"

class Endianness_test : public CppUnit::TestFixture
{
                CPPUNIT_TEST_SUITE(Endianness_test);
                CPPUNIT_TEST(swap_each_word_test);
                CPPUNIT_TEST(swap_endian_test);
                CPPUNIT_TEST_SUITE_END();
public:
                void setUp(){}
                void tearDown(){}
                void swap_each_word_test();
                void swap_endian_test();
};


#endif // ENDIANNESS_TEST_H
