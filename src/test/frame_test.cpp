#include "frame_test.h"
#include <boost/date_time/posix_time/posix_time.hpp>  
#define BOOST_DATE_TIME_SOURCE
#include <string>

void frame_test::setUp()
{
  pt.put("prog.value","frame_inner");
  pt.put("tet.value","../../dat/sphere/20k/tet/sphere-20k.tet");
}

void frame_test::test_alglib_lbfgs()
{
  std::string strTime = boost::posix_time::to_iso_string(boost::posix_time::second_clock::local_time());  
  
  pt.put("package.value","alglib");
  pt.put("align_w.value",100);
  pt.put("fix_w.value",1e3);
  pt.put("alg.value","lbfgs");
  pt.put("lbfgs-len.value",7);
  pt.put("iter.value",100);
  std::string output_name = "../../dat/sphere/20k/test/sphere-20k_alglib_lbfgs_len7_iter100_";
  output_name += strTime;
  pt.put("zyz.value", output_name + ".test.zyz");
  frame(pt);
}
void frame_test::tearDown(){}
