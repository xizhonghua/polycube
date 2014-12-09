#include <vector>

#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/io.h>
#include <zjucad/matrix/matrix.h>

#include "../tetmesh/tetmesh.h"

using namespace std;
using namespace zjucad::matrix;

int dfs(const size_t index, vector<bool> &is_visited, 
	const zjucad::matrix::matrix<size_t> &face,
	const jtf::mesh::edge2cell_adjacent_general &edge_adj)
{
  size_t face_index, edge_index;
  vector<size_t> edge2tri;
  if(!is_visited[index])
    {
      is_visited[index] = true;
      for(size_t i = 0; i < 3; ++i)
	for(size_t j = i+1; j < 3; ++j)
	  {
	    edge_index = edge_adj.get_edge_idx(face(i, index), face(j, index));
	    edge2tri = edge_adj.edge2cell_[edge_index];
	    if(edge2tri.size() != 2)
	      continue;
	    face_index = edge2tri[0] + edge2tri[1] - index;
	    dfs(face_index, is_visited, face, edge_adj);
	  }
    }
  return 0;
}

bool test(const char *path)
{
  zjucad::matrix::matrix<double> node;
  zjucad::matrix::matrix<size_t> tet, face;
  jtf::mesh::tet_mesh_read_from_zjumat(path, &node, &tet);

  double volume = 0.0;
  for(size_t i = 0; i < tet.size(2); ++i)
    {
      cout << jtf::mesh::cal_tet_vol(node(colon(), tet(colon(), i))) << " ";
      volume += jtf::mesh::cal_tet_vol(node(colon(), tet(colon(), i)));
    }
  cout << "volume: " << volume << endl;

  static jtf::mesh::face2tet_adjacent *face_adj = face_adj->create(tet);
  get_outside_face(*face_adj, face);
  static jtf::mesh::edge2cell_adjacent_general *edge_adj = edge_adj->create(face);
  vector<bool> is_visited(face.size(2), false);
  cout << is_visited.size() << endl;
  dfs(0, is_visited, face, *edge_adj);

  for(size_t i = 0; i < is_visited.size(); ++i)
    if(!is_visited[i])
      return false;
  return true;
}


int test_topology(int argc, char **argv)
{
  if(argc != 2)
    {
      cerr << "please input the tet mesh file" << endl;
      return 1;
    }
  if(test(argv[1]))
    cout << "OK" << endl;
  else
    cout << "not ok" <<endl;
  return 0;
}
