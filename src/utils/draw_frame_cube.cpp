#include "../tetmesh/tetmesh.h"
#include "../common/zyz.h"
#include "../common/util.h"

using namespace std;
using namespace zjucad::matrix;

int draw_frame_cube(int argc, char * argv[])
{
  if(argc != 4 && argc != 5){
      cerr << "# [usage] draw_frame_cube tet zyz cube.obj [scale]" << endl;
      return __LINE__;
    }

  jtf::tet_mesh tm(argv[1]);
  matrix<double> zyz;
  if(jtf::mesh::read_matrix(argv[2], zyz)){
      cerr << "# [error] can not open zyz file." << endl;
      return __LINE__;
    }
  if(zyz.size(2) != tm.tetmesh_.mesh_.size(2)){
      cerr << "# [error] zyz is not compatible with tet." << endl;
      return __LINE__;
    }

  matrix<double> cube_node; matrix<size_t> cube_mesh;
  {
    if(jtf::mesh::load_obj(argv[3], cube_mesh, cube_node)) {
        cerr << "# [error] can not open cube obj." << endl;
        return __LINE__;
      }

    const double len = calc_bounding_sphere_size(cube_node) ;
    cube_node /= len;
  }

  matrix<matrix<double> > rot_m(zyz.size(2),1);
  for(size_t ti = 0; ti < zyz.size(2); ++ti){
      rot_m[ti].resize(3,3);
      zyz_angle_2_rotation_matrix1(&zyz(0,ti), &rot_m[ti][0]);
    }

  const double len = jtf::mesh::cal_average_edge(tm.tetmesh_.mesh_, tm.tetmesh_.node_);
  cube_node *= len;

  if(argc == 5){
      const double scale = atof(argv[4]);
      cube_node *= scale;
    }
  matrix<double> new_node(3, cube_node.size(2) * tm.tetmesh_.mesh_.size(2));
  matrix<size_t> new_face(cube_mesh.size(1), cube_mesh.size(2) * tm.tetmesh_.mesh_.size(2));
  matrix<double> deformed_cube;
  matrix<size_t> deformed_face;
  matrix<double> center(3,1);
  for(size_t ti = 0; ti < tm.tetmesh_.mesh_.size(2); ++ti){
      center *= 0;
      deformed_cube = rot_m[ti] * cube_node;
      for(size_t pi = 0; pi < tm.tetmesh_.mesh_.size(1); ++pi){
          center += tm.tetmesh_.node_(colon(), tm.tetmesh_.mesh_(pi,ti));
        }
      center /= 4;
      deformed_cube += center*ones<double>(1,cube_node.size(2));
      deformed_face = cube_mesh;
      deformed_face += ti * cube_node.size(2);
      new_node(colon(), colon(ti*cube_node.size(2),(ti+1)*cube_node.size(2)-1)) = deformed_cube;
      new_face(colon(), colon(ti*cube_mesh.size(2),(ti+1)*cube_mesh.size(2)-1)) = deformed_face;
    }

  jtf::mesh::save_obj("frame_cube.obj", new_face, new_node);
  cerr << "# [info] success save to frame_cube.obj" << endl;
  return 0;
}
