#include "tet_mesh_sxx.h"
#include <boost/tuple/tuple.hpp>
#include <zjucad/matrix/matrix.h>

using namespace std;
using namespace sxx;

namespace sxx
{
  int create_dynamic_mesh(tet_mesh &my_mesh, 
			  const zjucad::matrix::matrix<double> &node,
			  const zjucad::matrix::matrix<size_t> &tet_matrix)
  {
    my_mesh.create_tetmesh(node, tet_matrix);
    return 0;
  }

  //the first element in tup is the index of the middle point of the edge
  // and the 2nd and 3rd are the index of the edge
  size_t split_edge_in_mid(tet_mesh &my_mesh,
                           const pair<size_t, size_t> &e)
  {
    size_t mid_index;
    mid_index = my_mesh.split_edge(e);
    return mid_index;
  }

  int recover_dynamic_mesh_to_mat(const tet_mesh &my_mesh,
				  zjucad::matrix::matrix<double> &node,
				  zjucad::matrix::matrix<size_t> &tet_matrix)
  {
    my_mesh.write_tetmesh_to_matrix(node, tet_matrix);
    return 0;
  }

  int get_tet2orginal_index(const tet_mesh &my_mesh,
			    const zjucad::matrix::matrix<size_t> &tet_matrix,
			    zjucad::matrix::matrix<size_t> &tet_index_map)
  {
    my_mesh.get_tet2orginal_index(tet_matrix, tet_index_map);
  }

  int get_face2orginal(tet_mesh &my_mesh,
                       boost::unordered_map<vector<size_t>, vector<size_t> > &f2omap)
  {
    f2omap = my_mesh.get_face2orginal();
    return 0;
  }

}
