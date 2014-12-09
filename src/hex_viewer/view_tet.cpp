#include <zjucad/ptree/ptree.h>
#include <jtflib/igl/readtet.h>
#include <jtflib/igl/viewer_tet.h>
#include <eigen3/Eigen/Core>

using namespace std;
using boost::property_tree::ptree;

int view_tet(ptree &pt)
{
  igl::Viewer_tet viewer;
  viewer.launch(pt.get<string>("input/tet.value").c_str());
  return 0;
}

