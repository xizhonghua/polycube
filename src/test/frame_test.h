#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextOutputter.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <boost/property_tree/ptree.hpp>
#include "../hex_ui/prog.h"
using boost::property_tree::ptree;


class frame_test : public CppUnit::TestFixture
{
		CPPUNIT_TEST_SUITE(frame_test);
		CPPUNIT_TEST(test_alglib_lbfgs);
		CPPUNIT_TEST_SUITE_END();
public:
		void setUp();
		void tearDown();
		void test_alglib_lbfgs();
private:

		ptree pt;
};
