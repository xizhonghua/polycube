#include <zjucad/ptree/ptree.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>

#include <jtflib/mesh/util.h>
#include <string>

#include "../tetmesh/tetmesh.h"
#include "../common/zyz.h"
#include "../hex_frame_opt/optimize_frame_field_L1.h"
using namespace std;
using namespace jtf::mesh;
using boost::property_tree::ptree;
using namespace zjucad::matrix;


int frame_inner_L1(ptree &pt)
{
  pt.put("input/tet.desc", "input tetmesh.");
  pt.put("input/init_zyz.desc", "input init zyz.");
  pt.put("output/zyz.desc" , "output zyz");

  jtf::tet_mesh tm(pt.get<string>("input/tet.value").c_str());

  matrix<double> zyz = zeros<double>(3,tm.tetmesh_.mesh_.size(2));
  if(zjucad::has("input/init_zyz.value",pt)){
      if(read_matrix(pt.get<string>("input/init_zyz.value").c_str(), zyz))
        return __LINE__;
      if(zyz.size(2) != tm.tetmesh_.mesh_.size(2)){
          cerr << "# [error] wrong zyz size." << endl;
          return __LINE__;
        }
      cerr << "# [info] use init_zyz." << endl;
    }

  matrix<double> frame_init(9, zyz.size(2));
  for(size_t i = 0; i < tm.tetmesh_.mesh_.size(2); ++i){
      zyz_angle_2_rotation_matrix1(&zyz(0,i), &frame_init(0,i));
    }


  optimize_frame_field_L1(tm, frame_init, pt);

  for(size_t i = 0; i < tm.tetmesh_.mesh_.size(2); ++i){
      rotation_matrix_2_zyz_angle(&frame_init(0,i), &zyz(0,i),0);
    }

  if(write_matrix(pt.get<string>("output/zyz.value").c_str(), zyz))
    return __LINE__;
  return 0;
}

int frame_inner_L1_euler(ptree &pt)
{
  pt.put("input/tet.desc", "input tetmesh.");
  pt.put("input/init_zyz.desc", "input init zyz.");
  pt.put("output/zyz.desc" , "output zyz");

  jtf::tet_mesh tm(pt.get<string>("input/tet.value").c_str());

  matrix<double> zyz = zeros<double>(3,tm.tetmesh_.mesh_.size(2));
  if(zjucad::has("input/init_zyz.value",pt)){
      if(read_matrix(pt.get<string>("input/init_zyz.value").c_str(), zyz))
        return __LINE__;
      if(zyz.size(2) != tm.tetmesh_.mesh_.size(2)){
          cerr << "# [error] wrong zyz size." << endl;
          return __LINE__;
        }
    }

  optimize_frame_field_L1_euler(tm, zyz, pt);

  if(write_matrix(pt.get<string>("output/zyz.value").c_str(), zyz))
    return __LINE__;
  return 0;
}
