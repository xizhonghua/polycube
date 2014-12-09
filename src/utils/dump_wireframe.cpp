#include <fstream>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <zjucad/matrix/matrix.h>
#include "../common/visualize_tool.h"

using namespace std;
using namespace jtf::mesh;
using namespace zjucad::matrix;

int dump_wireframe(int argc, char *argv[])
{
  if(argc != 4){
      cerr << "# [error] dump_wireframe input_obj output_obj radius" << std::endl;
      return __LINE__;
    }

  double radius = atof(argv[3]);
  matrix<size_t> surface;
  matrix<double> node;

  if(jtf::mesh::load_obj(argv[1], surface, node))
    return __LINE__;

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(surface));
  if(!ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << std::endl;
      return __LINE__;
    }

  vector<deque<pair<size_t,size_t> > > edges_chains;
  deque<pair<size_t,size_t> > one_chain;
  for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
      one_chain.clear();
      one_chain.push_back(ea->edges_[ei]);
      edges_chains.push_back(one_chain);
    }
  dump_singularity_to_cylinder(argv[2], node, edges_chains, radius);
  return 0;
}
