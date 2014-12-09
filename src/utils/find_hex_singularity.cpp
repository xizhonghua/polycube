#include <iostream>
#include <zjucad/matrix/matrix.h>
#include <vector>
#include <deque>

#include <jtflib/mesh/io.h>
#include "../hexmesh/io.h"
#include "../hexmesh/util.h"
#include <jtflib/mesh/mesh.h>
#include "../hexmesh/hexmesh.h"
#include "../common/visualize_tool.h"
#include "../common/util.h"

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
using namespace std;
using namespace zjucad::matrix;

int dump_surface_singularity_point(
    const vector<vector<size_t> > & singularity_point,
    const matrix<double> & node,
    const double radius,
    const char * sphere_file,
    const string &singularity_point_str);

int dump_quad_face_singularity_point(
    const matrix<size_t> & outside_face,
    const matrix<double> & node,
    const double radius,
    const char * sphere_file,
    const string singularity_point_str_pref)
{
  boost::unordered_map<size_t, boost::unordered_set<size_t> > p2f;
  for(size_t fi = 0; fi < outside_face.size(2); ++fi){
      for(size_t pi = 0; pi < outside_face.size(1); ++pi){
          p2f[outside_face(pi,fi)].insert(fi);
        }
    }

  std::shared_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent" << endl;
      return __LINE__;
    }
  matrix<size_t> edges;
  jtf::mesh::get_boundary_edge(*ea, edges);

  set<size_t> boundary_points(edges.begin(), edges.end());

  vector<vector<size_t> > singularity_point(2); // < 4; > 4
  for(boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
      cit = p2f.begin(); cit != p2f.end(); ++cit){
      if(cit->second.size() < 4 &&
         boundary_points.find(cit->first) == boundary_points.end()){
          singularity_point[0].push_back(cit->first);
        }else if(cit->second.size() > 4){
          singularity_point[1].push_back(cit->first);
        }
    }

  dump_surface_singularity_point(singularity_point, node, radius,
                                 sphere_file, singularity_point_str_pref);
}

int dump_surface_singularity_point(
    const vector<vector<size_t> > & singularity_point,
    const matrix<double> & node,
    const double radius,
    const char * sphere_file,
    const string &singularity_point_str)
{
  assert(singularity_point.size() == 2);

  jtf::mesh::meshes trm;
  if(jtf::mesh::load_obj(sphere_file, trm.mesh_, trm.node_))
    return __LINE__;

  matrix<double> bb(3,2);
  calc_bounding_box(trm.node_, &bb[0]);
  const matrix<double> center = (bb(colon(),0) + bb(colon(),1))/2.0;
  for(size_t pi = 0; pi < trm.node_.size(2); ++pi)
    trm.node_(colon(),pi) -= center;
  trm.node_ *= 2 * radius * calc_bounding_sphere_size(node);

  cerr << "# [info] positive singularity point " << singularity_point[0].size() << endl;
  cerr << "# [info] negative singularity point " << singularity_point[1].size() << endl;

  const string singularity_point_positive = singularity_point_str + "_positive.obj";
  const string singularity_point_negative = singularity_point_str + "_negative.obj";

  jtf::mesh::meshes trm_pos, trm_neg;
  trm_pos.node_.resize(3, trm.node_.size(2) * singularity_point[0].size());
  trm_pos.mesh_.resize(3, trm.mesh_.size(2) * singularity_point[0].size());
  for(size_t pi = 0; pi < singularity_point[0].size(); ++pi){
      trm_pos.node_(colon(), colon(pi * trm.node_.size(2), (pi+1)* trm.node_.size(2)-1))
          = trm.node_;
      for(size_t pp = pi * trm.node_.size(2); pp <(pi+1)* trm.node_.size(2); ++pp )
        trm_pos.node_(colon(), pp) += node(colon(), singularity_point[0][pi]);

      trm_pos.mesh_(colon(), colon(pi* trm.mesh_.size(2), (pi+1)* trm.mesh_.size(2)-1))
          = trm.mesh_ ;
      for(size_t pp = pi * trm.mesh_.size(2); pp <(pi+1)* trm.mesh_.size(2); ++pp)
        trm_pos.mesh_(colon(), pp) += pi * trm.node_.size(2);
    }

  trm_neg.node_.resize(3, trm.node_.size(2) * singularity_point[1].size());
  trm_neg.mesh_.resize(3, trm.mesh_.size(2) * singularity_point[1].size());
  for(size_t pi = 0; pi < singularity_point[1].size(); ++pi){
      trm_neg.node_(colon(), colon(pi * trm.node_.size(2), (pi+1)* trm.node_.size(2)-1))
          = trm.node_ ;

      for(size_t pp = pi * trm.node_.size(2); pp <(pi+1)* trm.node_.size(2); ++pp )
        trm_neg.node_(colon(), pp) += node(colon(), singularity_point[1][pi]);

      trm_neg.mesh_(colon(), colon(pi* trm.mesh_.size(2), (pi+1)* trm.mesh_.size(2) - 1))
          = trm.mesh_ ;
      for(size_t pp = pi * trm.mesh_.size(2); pp <(pi+1)* trm.mesh_.size(2); ++pp)
        trm_neg.mesh_(colon(), pp) += pi * trm.node_.size(2);
    }

  jtf::mesh::save_obj(singularity_point_positive.c_str(), trm_pos.mesh_, trm_pos.node_);
  jtf::mesh::save_obj(singularity_point_negative.c_str(), trm_neg.mesh_, trm_neg.node_);
  return 0;
}

int find_hex_singularity(int argc, char *argv[])
{
  if(argc != 5 && argc != 6){
      cerr << "# [usage] find_hex_singularity hex hex_format  singularity_obj input_unit_sphere_obj"
           << "[radius_of_cylinder]" << endl;
      return __LINE__;
    }

  jtf::hex_mesh hm(argv[1]);
  jtf::hexmesh::hex_singularity_extractor hse(hm);
  matrix<size_t> surface_singularity = hse.get_surface_singularity_edges();
  matrix<size_t> inner_singularity = hse.get_inner_singularity_edges();

  vector<pair<size_t,size_t> >  singularity_inner, singularity_surface;
  for(size_t ei = 0; ei < surface_singularity.size(2); ++ei)
    singularity_surface.push_back(make_pair(surface_singularity(0,ei),surface_singularity(1,ei)));
  for(size_t ei = 0; ei < inner_singularity.size(2); ++ei)
    singularity_inner.push_back(make_pair(inner_singularity(0,ei),inner_singularity(1,ei)));

  vector<deque<pair<size_t,size_t> > > singularity_inner_vec;
  vector<deque<pair<size_t,size_t> > > singularity_surface_vec;

  jtf::util::extract_chain_from_edges(singularity_inner, singularity_inner_vec);
  jtf::util::extract_chain_from_edges(singularity_surface, singularity_surface_vec);

  string inner_singularity_str = "inner_";
  inner_singularity_str += argv[3];
  string surface_singularity_str = "surface_";
  surface_singularity_str += argv[3];

  double radius = 0.002;
  if(argc == 6){
      radius = atof(argv[5]);
    }

  if(!singularity_inner.empty())
    dump_singularity_to_cylinder(
          inner_singularity_str.c_str(), hm.hexmesh_.node_, singularity_inner_vec,radius);
  if(!singularity_surface.empty())
    dump_singularity_to_cylinder(
          surface_singularity_str.c_str(), hm.hexmesh_.node_, singularity_surface_vec,radius);

  string singularity_point_str = argv[3] ;
  singularity_point_str += "_point";

  dump_quad_face_singularity_point(hm.outside_face_, hm.hexmesh_.node_, radius, argv[4], singularity_point_str);
  return 0;
}
