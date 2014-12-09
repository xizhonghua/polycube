#include <vector>
#include <fstream>
#include <numeric>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/property_tree/ptree.hpp>

#include <jtflib/util/vertex_connection.h>
#include <jtflib/mesh/io.h>

#include "../common/vtk.h"
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
#include "../common/util.h"
#include "io.h"
#include "util.h"

using namespace std;
using namespace boost;
using namespace zjucad::matrix;
using namespace boost::property_tree;

class quad_patch
{
public:
  std::vector<size_t> face_idx_;
  std::vector<size_t> boundary_chain_;
  //std::vector<deque<pair<size_t,size_t> > > & chains_;
};

int cal_node_avg_variance(const matrix<double> & node,
                          matrix<double> & avg,
                          matrix<double> & variance)
{
  avg = zeros<double>(node.size(1),1);
  variance = zeros<double>(node.size(1),1);
  for(size_t di = 0; di < node.size(1); ++di)
    avg[di] = std::accumulate(node(di,colon()).begin(), node(di,colon()).end(),0.0);
  avg /= node.size(2);
  for(size_t di = 0; di < 3; ++di){
    for(size_t p = 0; p < node.size(2); ++p){
      variance[di] += (node(di, p) - avg[di]) * (node(di, p) - avg[di]);
    }
  }
  variance /= node.size(2);

  return 0;
}

int analysis_patch_to_extract_planes(
    const vector<matrix<size_t> > &surface_patches,
    const matrix<double> &node,
    vector<vector<double> > &uvw_planes,
    vector<size_t> & patch_major_dir)
{
  uvw_planes.clear();
  uvw_planes.resize(3);

  patch_major_dir.clear();

  matrix<double> patch_node;
  matrix<double> avg = zeros<double>(3,1), variance = zeros<double>(3,1);

  patch_major_dir.resize(surface_patches.size());
  for(size_t pi = 0; pi < surface_patches.size(); ++pi){
    const matrix<size_t> & one_patch = surface_patches[pi];
    set<size_t> patch_node_set(one_patch.begin(), one_patch.end());
    matrix<size_t> patch_node_mat(patch_node_set.size(),1);
    copy(patch_node_set.begin(), patch_node_set.end(), patch_node_mat.begin());
    patch_node = node(colon(), patch_node_mat);

    // calculate avg and variance
    cal_node_avg_variance(patch_node, avg, variance);

    // find the minimal variance, and record the crosspoding avg as a plane.
    const size_t idx = min_element(variance.begin(), variance.end()) - variance.begin();
    uvw_planes[idx].push_back(avg[idx]);
    patch_major_dir[pi] = idx;
  }

  for(size_t i = 0; i < 3; ++i){
    sort(uvw_planes[i].begin(), uvw_planes[i].end());
  }

  return 0;
}

int split_each_patch(const vector<vector<double> > &uvw_planes,
                     const matrix<double> & node,
                     const matrix<size_t> &patch_faces,
                     const size_t & major_direction,
                     const double tolerance,
                     boost::unordered_set<pair<size_t,size_t> > & boundary_edges)
{

  // assume the uvw_plane must be ordered
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(patch_faces));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  map<pair<size_t,size_t>, double> edge_length_map;


  boost::unordered_set<size_t> boundary_points;
  for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
    const pair<size_t,size_t> & one_edge = ea->edges_[ei];
    const pair<size_t,size_t> & tri_pari = ea->edge2cell_[ei];
    edge_length_map[one_edge] = norm(node(colon(), one_edge.first) -
                                     node(colon(), one_edge.second));
    if(ea->is_boundary_edge(tri_pari)){
      boundary_points.insert(one_edge.first);
      boundary_points.insert(one_edge.second);
    }
  }

  unique_ptr<vertex_connection<UNDIRECT> > vc(
        vertex_connection<UNDIRECT>::create(edge_length_map));
  if(!vc.get()){
    cerr << "# [error] can not build vertex_connection." << endl;
    return __LINE__;
  }

  matrix<size_t> boundary_point_mat(boundary_points.size(),1);
  copy(boundary_points.begin(), boundary_points.end(), boundary_point_mat.begin());

  matrix<double> boundary_node = node(colon(), boundary_point_mat);

  matrix<double> bb = zeros<double>(3,2);

  calc_bounding_box(boundary_node, &bb[0]);

  cerr << "# bb " << bb << endl;

  vector<vector<double> > uvw_plane_space(3);
  for(size_t i = 0; i < 2; ++i){ // except the major direction, find the closet point
    const size_t di = (major_direction + i + 1) % node.size(1);
    size_t pi = 0;
    for(; pi < uvw_planes[di].size(); ++pi){
      if(uvw_planes[di][pi] < bb(di,0) && fabs(uvw_planes[di][pi] - bb(di,0)) > tolerance){
        if(uvw_plane_space[di].empty())
          uvw_plane_space[di].push_back(uvw_planes[di][pi]);
        else
          uvw_plane_space[di].back() = uvw_planes[di][pi];
      }else{
        if(!uvw_plane_space[di].empty() && fabs(uvw_planes[di][pi] - bb(di,0)) < tolerance)
          uvw_plane_space[di].back() = uvw_planes[di][pi];
        else
          uvw_plane_space[di].push_back(uvw_planes[di][pi]);
        break;
      }
    }
    if(pi == uvw_planes[di].size()){
      cerr << "# [error] patch bounding box does not match any uvw planes." << endl;
      return __LINE__; // to0 large tolerance?
    }

    for(; pi < uvw_planes[di].size(); ++pi){
      if(fabs(uvw_planes[di][pi] - bb(di,1)) < tolerance ||
         uvw_planes[di][pi] < bb(di,1)){
        if(fabs(uvw_plane_space[di].back() - uvw_planes[di][pi]) > tolerance)
        uvw_plane_space[di].push_back(uvw_planes[di][pi]);
      }else{
        break;
      }
    }
  }

  // find the space of plances, then we should find the points which most near
  // the cross points of planes
  for(size_t i = 0; i < 2;  ++i){
    // for each plane, find the nearest
    const size_t di = (major_direction + i + 1) % node.size(1);
    const size_t other_di = (1 + 2 - major_direction - di) % node.size(1);

    matrix<double> node_on_plane = zeros<double>(2,1);
    matrix<double> node_on_boundary = zeros<double>(2,1);

    for(size_t dii = 0; dii < uvw_plane_space[di].size(); ++dii){
      vector<size_t> points_need_to_find_short_path;
      for(size_t odii = 0; odii < uvw_plane_space[other_di].size(); ++odii){
        node_on_plane[0] = uvw_plane_space[di][dii];
        node_on_plane[1] = uvw_plane_space[other_di][odii];

        vector<pair<double, size_t> > length2plane_point;
        for(boost::unordered_set<size_t>::const_iterator cit = boundary_points.begin();
            cit != boundary_points.end(); ++cit){
          node_on_boundary[0] = node(di, *cit);
          node_on_boundary[1] = node(other_di, *cit);
          length2plane_point.push_back(
                make_pair(norm(node_on_plane - node_on_boundary), *cit));
        }
        sort(length2plane_point.begin(), length2plane_point.end());
        points_need_to_find_short_path.push_back(length2plane_point.front().second);
      }

      vector<size_t> path;
      for(size_t p = 0; p < points_need_to_find_short_path.size()-1; ++p){
        int rtn = vc->get_shortest_path(points_need_to_find_short_path[p],
                                        points_need_to_find_short_path[p+1], path);
        if(rtn == 0){
          for(size_t pi = 0; pi < path.size()-1; ++pi){
            if(path[pi] > path[pi+1])
              boundary_edges.insert(make_pair(path[pi+1], path[pi]));
            else
              boundary_edges.insert(make_pair(path[pi], path[pi+1]));
          }
        }
      }
    }
  }

  return 0;
}

//! @WARNING!!! This function takes an assumption that input tet must be polycube tet
int split_patch_to_quad_patch(ptree & pt)
{
  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("polycube_tet.value").c_str(), &tm.node_, &tm.mesh_))
    return __LINE__;

  boost::unordered_map<size_t,size_t> surface_type;
  if(load_surface_restricted_type_static(pt.get<string>("surface_type.value").c_str(), surface_type))
    return __LINE__;

  vector<matrix<size_t> > surface_patches;
  boost::unordered_set<pair<size_t,size_t> > patch_boundary;
  convert_surface_type_to_surface_patches(tm.mesh_, surface_type, surface_patches,
                                          patch_boundary);
  vector<quad_patch> split_patches;
  vector<size_t> patch_major_dir;
  vector<vector<double> > uvw_planes;
  analysis_patch_to_extract_planes(surface_patches, tm.node_, uvw_planes, patch_major_dir);

  {// debug
    ofstream ofs_iso_lines("iso_lines.vtk");
    vector<double> iso_lines;
    vector<size_t> lines;
    for(size_t x = 0; x < uvw_planes[0].size(); ++x)
      for(size_t y = 0; y < uvw_planes[1].size(); ++y){
        for(size_t z = 0; z < uvw_planes[2].size(); ++z){
          iso_lines.push_back(uvw_planes[0][x]);
          iso_lines.push_back(uvw_planes[1][y]);
          iso_lines.push_back(uvw_planes[2][z]);
          if(z != uvw_planes[2].size()-1){
            lines.push_back(iso_lines.size()/3-1);
            lines.push_back(iso_lines.size()/3);
          }
        }
      }

    for(size_t x = 0; x < uvw_planes[0].size(); ++x)
      for(size_t z = 0; z < uvw_planes[2].size(); ++z){
        for(size_t y = 0; y < uvw_planes[1].size(); ++y){
          iso_lines.push_back(uvw_planes[0][x]);
          iso_lines.push_back(uvw_planes[1][y]);
          iso_lines.push_back(uvw_planes[2][z]);
          if(y != uvw_planes[1].size()-1){
            lines.push_back(iso_lines.size()/3-1);
            lines.push_back(iso_lines.size()/3);
          }
        }
      }

    for(size_t z = 0; z < uvw_planes[2].size(); ++z){
      for(size_t y = 0; y < uvw_planes[1].size(); ++y){
        for(size_t x = 0; x < uvw_planes[0].size(); ++x){
          iso_lines.push_back(uvw_planes[0][x]);
          iso_lines.push_back(uvw_planes[1][y]);
          iso_lines.push_back(uvw_planes[2][z]);
          if(x != uvw_planes[0].size()-1){
            lines.push_back(iso_lines.size()/3-1);
            lines.push_back(iso_lines.size()/3);
          }
        }
      }
    }


    line2vtk(ofs_iso_lines, &iso_lines[0], iso_lines.size()/3, &lines[0], lines.size()/2);
  }
  boost::unordered_set<pair<size_t,size_t> > boundary_edges;

  double tolerance = 0.0;

  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit =
      patch_boundary.begin(); cit != patch_boundary.end(); ++cit){
    tolerance += norm(tm.node_(colon(), cit->first) - tm.node_(colon(), cit->second));
  }
  tolerance /= patch_boundary.size();
  tolerance /= 10.0;

  for(size_t pi = 0; pi < surface_patches.size(); ++pi){
    split_each_patch(uvw_planes, tm.node_, surface_patches[pi], patch_major_dir[pi], tolerance, boundary_edges);
  }

  {// visual
    ofstream ofs("boundary.vtk");
    vector<size_t> boundary_edges_vec;
    for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit
        = boundary_edges.begin(); cit != boundary_edges.end(); ++cit){
      boundary_edges_vec.push_back(cit->first);
      boundary_edges_vec.push_back(cit->second);
    }
    line2vtk(ofs, &tm.node_[0], tm.node_.size(2), &boundary_edges_vec[0], boundary_edges_vec.size()/2);

    vector<size_t> surface_patch_faces;
    vector<size_t> surface_patch_faces_type;
    for(size_t pi = 0; pi < surface_patches.size(); ++pi){
      const matrix<size_t> & one_patch = surface_patches[pi];
      surface_patch_faces.insert(surface_patch_faces.end(), one_patch.begin(), one_patch.end());
      for(size_t fi = 0; fi < one_patch.size(2); ++fi)
        surface_patch_faces_type.push_back(pi);
    }

    ofstream ofs_face("patch.vtk");
    tri2vtk(ofs_face, &tm.node_[0], tm.node_.size(2), &surface_patch_faces[0], surface_patch_faces.size()/3);
    cell_data(ofs_face, &surface_patch_faces_type[0], surface_patch_faces_type.size(), "patch");

  }

  return 0;
}
