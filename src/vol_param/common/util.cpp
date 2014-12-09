#include "util.h"
#include <jtflib/mesh/mesh.h>

using namespace std;
using namespace zjucad::matrix;

void smooth_mesh_surface(
    const zjucad::matrix::matrix<size_t> &tri_faces,
    zjucad::matrix::matrix<double> &node,
    const size_t iter)
{
  // fix boundary
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(tri_faces));
  matrix<size_t> boundary_edges;
  jtf::mesh::get_boundary_edge(*ea, boundary_edges);
  set<size_t> boundary_points(boundary_edges.begin(), boundary_edges.end());
  unique_ptr<jtf::mesh::one_ring_point_at_point>  orpap(jtf::mesh::one_ring_point_at_point::create(tri_faces));
  matrix<double> avg_node = zeros<double>(3,1);
  for(size_t i = 0 ;i < iter; ++i){
      for(const auto & one_p2p: orpap->p2p_){
          avg_node *= 0;
          if(boundary_points.find( one_p2p.first) == boundary_points.end()){
              const vector<size_t> & other_points = one_p2p.second;
              for(size_t pi = 0; pi < other_points.size(); ++pi){
                  avg_node += node(colon(), other_points[pi]);
                }
              avg_node /= other_points.size();
            }
          node(colon(), one_p2p.first) = avg_node;
        }
    }
}
