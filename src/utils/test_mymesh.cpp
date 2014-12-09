#include "../tet_mesh_sxx/tet_mesh_sxx.h"
//#include "../tet_mesh_sxx/edge_split_in_mid.h"
#include "../tetmesh/tetmesh.h"
#include <ctime>

using namespace std;
using namespace sxx;

int test_mymesh(int argc, char **argv)
{
  if(argc != 3)
    {
      cerr << "please input the tet mesh file" << endl;
      return 1;
    }
  tet_mesh t_mesh;
  time_t t_begin, t_end;
  double duration;
  t_begin = clock();
  t_mesh.create_tetmesh(argv[1]);
  t_mesh.test_topology_operation();
  t_mesh.write_tetmesh_to_file(argv[2]);
  t_end = clock();
  duration = static_cast<double>(t_end - t_begin) / CLOCKS_PER_SEC;
  cout << "run time: " << duration << endl; 	
  //t_mesh.create_tet_bin(argv[1], argv[2]);
  //t_mesh.delete_all_tet();	
  // zjucad::matrix::matrix<double> node;
  // zjucad::matrix::matrix<size_t> tet_matrix;
  // jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &node, &tet_matrix);
  // sxx::create_dynamic_mesh(t_mesh, node, tet_matrix);
  // pair<size_t, size_t> e(27, 411);
  // boost::tuple<size_t, size_t, size_t> tup;
  // sxx::split_edge_in_mid(t_mesh, e, tup);
  // cout << tup.get<0>() <<endl;

  return 0;
}		
