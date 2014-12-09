#include <boost/property_tree/ptree.hpp>

#include <string>
#include <fstream>

//#include "../mesh_func/tet-vol.h"
#include <zjucad/optimizer/optimizer.h>
#include <zjucad/ptree/ptree.h>

#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"
#include "../common/util.h"
#include "../common/IO.h"
#include <jtflib/mesh/mesh.h>
#include "util.h"
#include "io.h"

using namespace std;
using boost::property_tree::ptree;
using namespace zjucad::matrix;

int polycube_improve(ptree &pt)
{
  jtf::mesh::meshes tm,polycube_tm;

  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm.node_, &tm.mesh_))
    return __LINE__;

  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("polycube_tet.value").c_str(),
                               &polycube_tm.node_, &polycube_tm.mesh_))
    return __LINE__;

  if(tm.mesh_.size() != polycube_tm.mesh_.size()){
    cerr << "# [error] orig_tet is not compatible with polycube tet." << endl;
    return __LINE__;
  }

  cerr << "# read in tet " << tm.node_.size(2) << " " << tm.mesh_.size(2) << endl;

  pt.put("iter_b.desc", "smooth iteration on boundary, default is 5");
  pt.put("iter_p.desc", "smooth iteration in patch, default is 5");
  pt.put("iter_v.desc", "smooth iteration in volume, default is 5");
  pt.put("output.desc", "output tet");
  pt.put("surface_type.desc", "surface type file.");

  boost::unordered_map<size_t,size_t> surface_type;

  if(load_surface_restricted_type_static(
       pt.get<string>("surface_type.value").c_str(), surface_type))
    return __LINE__;

  const size_t iter_b = pt.get<size_t>("iter_b.value",5);
  const size_t iter_p = pt.get<size_t>("iter_p.value",5);
  const size_t iter_v = pt.get<size_t>("iter_v.value",5);

  vector<matrix<size_t> > surface_patches;
  boost::unordered_set<pair<size_t,size_t> > patch_boundary;

  convert_surface_type_to_surface_patches(
        tm.mesh_, surface_type, surface_patches, patch_boundary);

  matrixd node = polycube_tm.node_;
  { // smooth boundary
    smooth_boundary(patch_boundary, tm.node_, node, iter_b);
    //    vector<size_t> flipped_patchs;
    //    check_flipped_surface(surface_patches, node, flipped_patchs);
    //    dump_flipped_surface

    smooth_patch(surface_patches, patch_boundary, node, iter_p);
    smooth_volume(tm.mesh_, surface_patches, node, iter_v);

    //remove_flipped_tet(tm.mesh_, tm.node_, polycube_tm.mesh_, polycube_tm.node_);

    cerr << "# [info] smooth difference " << norm(node - tm.node_) << endl;

    ofstream ofs("after_smooth.vtk");
    tet2vtk(ofs, &node[0], node.size(2), &tm.mesh_[0], tm.mesh_.size(2));
  }

  jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("output.value").c_str(), &node, &tm.mesh_);

  cerr << "success." << endl;
  return 0;
}
