#include <iostream>
#include <fstream>
#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <jtflib/util/util.h>

#include "../common/def.h"
#include "../quadmesh/util.h"
#include "../common/util.h"
#include "../common/vtk.h"
#include "../numeric/util.h"

#include "../common/visualize_tool.h"
#include "../common/def.h"
#include <jtflib/mesh/io.h>


using namespace zjucad::matrix;
using namespace jtf::quadmesh;
using namespace std;

int dump_surface_feature(const char * file,
                      const matrixd & node,
                      const vector<deque<pair<size_t,size_t> > > & feature)
{
  ofstream ofs(file);
  if(ofs.fail()){
      cerr << "# [error] can not open file." << endl;
      return __LINE__;
    }
  ofs << feature.size() << endl;
  for(size_t fi = 0; fi < feature.size(); ++fi){
      const deque<pair<size_t,size_t> > & one_chain = feature[fi];
      ofs << one_chain.size() + 1 << endl;
      for(size_t i = 0; i < one_chain.size(); ++i){
          ofs << one_chain[i].first << " ";
        }
      ofs << one_chain.back().second << " " << endl;
    }
  return 0;
}


int dump_surface_feature(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] dump_surface_feature mesh min_angle_threshold[degree] output_feature " << endl;
      return __LINE__;
    }

  matrixst mesh;
  matrixd node;
  if(jtf::mesh::load_obj(argv[1], mesh, node)){
      cerr << "# [error] can not open quad file." << endl;
      return __LINE__;
    }

  const double min_angle_threshold = atof(argv[2]);
  matrixst feature_lines;
  jtf::mesh::extract_mesh_feature_line(mesh, node, feature_lines, cos(min_angle_threshold*My_PI()/180.0));

  vector<pair<size_t,size_t> > edges;
  for(size_t fi = 0; fi < feature_lines.size(2); ++fi){
      edges.push_back(make_pair(feature_lines(0,fi), feature_lines(1,fi)));
    }

  cerr << "edges " << edges.size() << endl;
  vector<deque<pair<size_t,size_t> > > chains;
  jtf::util::extract_chain_from_edges(edges, chains);
  dump_surface_feature(argv[3], node, chains);

  size_t edge_num = 0;
  for(size_t ci = 0; ci < chains.size(); ++ci){
      for(const auto & one_edge : chains[ci]){
          ++edge_num;
        }
    }
  cerr << "edges " << edge_num << endl;

  ofstream ofs("feature.vtk");
  line2vtk(ofs, &node[0], node.size(2), &feature_lines[0], feature_lines.size(2));

  {
    //debug
     vector<deque<pair<size_t,size_t> > > feature_lines_vec(feature_lines.size(2));
    for(size_t ei = 0 ;ei < feature_lines.size(2); ++ei){
        feature_lines_vec[ei].push_back(make_pair(feature_lines(0,ei),feature_lines(1,ei)));
      }
    dump_singularity_to_cylinder("feature_line_cylinder.obj", node, feature_lines_vec, 0.001);
  }
  return 0;
}
