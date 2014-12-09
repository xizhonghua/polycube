#include "find_singularities_inner_test.h"
#include <boost/date_time/posix_time/posix_time.hpp>
#define BOOST_DATE_TIME_SOURCE

void find_singularities_test::setUp()
{
    //std::string strTime = boost::posix_time::to_iso_string(boost::posix_time::second_clock::local_time());
    pt.put("prog.value","frame_inner");
    pt.put("tet.value","../../dat/sphere/20k/tet/sphere-20k.tet");
//    std::string output_name = "../../dat/sphere/20k/test/sphere-20k@10-1e3-1000_";
//    output_name += strTime;

    pt.put("zyz.value","../../dat/sphere/20k/test/sphere-20k@1000-1e3-1000.test.zyz");
    pt.put("align_w.value","1000");
    pt.put("fix_w.value","1e3");
    pt.put("iter.value", "1000");
    pt.put("package.value","alglib");
    pt.put("alg.value","lbfgs");
    pt.put("lbfgs-len.value","7");
    frame_inner(pt);
}

void find_singularities_test::tearDown()
{}

void find_singularities_test::test_find_singularities2_inner()
{
    pt.put("prog.value","find_singularities2_inner");
    find_singularities2_inner(pt);
}
