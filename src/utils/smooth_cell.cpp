#include <iostream>
#include <fstream>


#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>

#include <zjucad/matrix/matrix.h>

#include "../common/vtk.h"
using namespace std;
using namespace jtf::mesh;
using namespace zjucad::matrix;

int smooth_cell(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] smooth_cell input_obj output_obj iterations" << endl;
      return __LINE__;
    }

  meshes input_obj, output_obj;
  if(jtf::mesh::load_obj(argv[1], input_obj.mesh_, input_obj.node_))
    return __LINE__;

  std::shared_ptr<edge2cell_adjacent> ea(edge2cell_adjacent::create(input_obj.mesh_));
  if(!ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }

  matrix<size_t> edges;
  jtf::mesh::get_boundary_edge(*ea, edges);

  set<size_t> boundary_points(edges.begin(), edges.end());

  std::shared_ptr<one_ring_point_at_point> orpap(
        one_ring_point_at_point::create(input_obj.mesh_));
  if(!orpap.get()){
      cerr << "# [error] can not build one_ring_point_at_point." << endl;
      return __LINE__;
    }

//  {
//    vector<size_t> fl;
//    ifstream ifs("compress.tri.obj.fl");
//    size_t num = 0;
//    ifs >> num >> num;
//    fl.resize(num);
//    for(size_t i = 0; i < num; ++i)
//      ifs >> fl[i];

//    matrix<double> delta = output_obj.node_(colon(),fl.back())
//        - output_obj.node_(colon(),fl.front());
//    delta /= num-1;
//    for(size_t i = 1; i < num-1; ++i){
//        output_obj.node_(colon(),fl[i]) = output_obj.node_(colon(),fl.front())
//            + i * delta;
//      }
//  }

  output_obj = input_obj;
  size_t iteritons = atoi(argv[3]);
  for(size_t j = 0; j < iteritons; ++j)
    for(one_ring_point_at_point::p2p_type::const_iterator cit = orpap->p2p_.begin();
        cit != orpap->p2p_.end(); ++cit){
        const vector<size_t> & linked_points = cit->second;
        if(boundary_points.find(cit->first) != boundary_points.end()) continue;
        output_obj.node_(colon(),cit->first) *= 0;
        for(size_t i = 0; i < linked_points.size(); ++i)
          output_obj.node_(colon(), cit->first) += output_obj.node_(colon(), linked_points[i]);
        output_obj.node_(colon(),cit->first) /= linked_points.size();
      }


  if(jtf::mesh::save_obj(argv[2], output_obj.mesh_, output_obj.node_))
    return __LINE__;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// just test

int load_size_field(const char * file,
                    matrix<double> & size_field)
{
  ifstream ifs(file);
  if(ifs.fail()){
      cerr << "# [error] can not open size_field" << endl;
      return __LINE__;
    }
  size_t num = 0;
  ifs >> num;
  if(size_field.size() != num)
    size_field.resize(num,1);

  for(size_t i = 0; i < num; ++i)
    ifs >> size_field[i];
  return 0;
}

void find_closet_point(const matrix<double> & p0,
                       const matrix<double> & p_set,
                       size_t &idx)
{
  static vector<double> dis(p_set.size(2));
  for(size_t pi = 0; pi < p_set.size(2); ++pi) {
      dis[pi] = norm(p0 - p_set(colon(),pi));
    }
  vector<double>::const_iterator min_it = min_element(dis.begin(), dis.end());
  idx = static_cast<size_t>(min_it - dis.begin());
}

int tri_size2quad(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] tri_size2quad tri size quad" << endl;
      return __LINE__;
    }
  jtf::mesh::meshes tri, quad;
  if(jtf::mesh::load_obj(argv[1], tri.mesh_, tri.node_))
    return __LINE__;
  if(jtf::mesh::load_obj(argv[3], quad.mesh_, quad.node_))
    return __LINE__;
  matrix<double> size_field = zeros<double>(tri.node_.size(2),1);

  if(load_size_field(argv[2], size_field))
    return __LINE__;
  matrix<double> quad_size_field = zeros<double>(quad.node_.size(2),1);
  size_t idx = 0;
  for(size_t pi = 0; pi < quad.node_.size(2); ++pi){
      find_closet_point(quad.node_(colon(),pi), tri.node_, idx);
      quad_size_field[pi] = size_field[idx];
      if(quad_size_field[pi] < 0) {
          cerr << "quad node " << pi << " tri node " << idx << " tri size " << size_field[idx] << endl;
        }
    }

  ofstream ofs("quad_size_field.vtk");
  quad2vtk(ofs, &quad.node_[0], quad.node_.size(2), &quad.mesh_[0], quad.mesh_.size(2));
  point_data(ofs, &quad_size_field[0], quad_size_field.size(), "size");

  ofstream ofs_tri("tri_size_field.vtk");
  tri2vtk(ofs_tri, &tri.node_[0], tri.node_.size(2), &tri.mesh_[0], tri.mesh_.size(2));
  point_data(ofs_tri, &size_field[0], size_field.size(), "tri_size");

  return 0;
}
