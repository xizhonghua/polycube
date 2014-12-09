#include <boost/property_tree/ptree.hpp>

#include <string>
#include <fstream>
#include <numeric>
#include <hjlib/function/func_aux.h>
#include <hjlib/math/polar.h>

#include <zjucad/optimizer/optimizer.h>
#include <zjucad/ptree/ptree.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/util/util.h>
#include <jtflib/function/func_aux.h>
#include <jtflib/util/vertex_connection.h>
#include <jtflib/optimizer/optimizer.h>
#include "io.h"
#include "util.h"
#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"
#include "../common/util.h"
#include "../common/IO.h"
#include "../tetmesh/tetmesh.h"

#include "io.h"

#include "../mesh_func/tri-area.h"

#include "../mesh_func/tri-normal.h"

#include "util.h"

double lap_w, fix_w_3;

using namespace std;
using boost::property_tree::ptree;
using namespace hj::function;
using namespace zjucad::matrix;

double calc_path_length(
    const matrix<double> & node,
    vector<size_t> & path)
{
  double len  = 0;
  for(size_t pi = 0; pi < path.size()-1; ++pi){
    len += norm(node(colon(), path[pi]) - node(colon(),path[pi+1]));
  }
  return len;
}

double calc_path_length(
    const matrix<double> & node,
    deque<size_t> & path)
{
  double len  = 0;
  for(size_t pi = 0; pi < path.size()-1; ++pi){
    len += norm(node(colon(), path[pi]) - node(colon(),path[pi+1]));
  }
  return len;
}

bool is_large_distortion_chain(const matrix<double> & node,
                               const deque<pair<size_t,size_t> > & chain)
{
  if(chain.size() == 1) return false;
  const double PI=3.1415926;
  // this function assume each chain should be close to a straight line
  vector<double> cos_angle(chain.size()-1, 0);
  matrix<double> e(3,2);
  for(size_t ei = 0; ei < chain.size()-1; ++ei){
    const pair<size_t,size_t> & current_edge = chain[ei];
    const pair<size_t,size_t> & next_edge = chain[ei+1];
    e(colon(),0) = node(colon(), current_edge.second) -
                   node(colon(), current_edge.first);
    e(colon(),1) = node(colon(), next_edge.second) -
                   node(colon(), next_edge.first);
    double len_0 = norm(e(colon(),0));
    if(len_0 < 1e-6) len_0 = 1;
    double len_1 = norm(e(colon(),1));
    if(len_1 < 1e-6)  len_1 = 1;
    cos_angle[ei] = dot(e(colon(),0) , e(colon(),1)) / (len_0 * len_1);
  }

  double max_cos = *max_element(cos_angle.begin(), cos_angle.end());
  double avg_cos = std::accumulate(cos_angle.begin(), cos_angle.end(), 0.0) / cos_angle.size();

  if(max_cos < cos(PI/3.0) || avg_cos < cos(PI/12.0))
    return true;
  return false;
}

int retype_each_patch(
    boost::unordered_map<size_t,size_t> &surface_type,
    const vector<vector<size_t> > &patches)
{
  for(size_t pi = 0; pi < patches.size(); ++pi){
    const vector<size_t> & one_patch = patches[pi];
    boost::unordered_map<size_t,
        vector<boost::unordered_map<size_t,size_t>::iterator > > type2it;
    for(size_t fi = 0; fi < one_patch.size(); ++fi){
      boost::unordered_map<size_t,size_t>::iterator it =
          surface_type.find(one_patch[fi]);
      if(it == surface_type.end()){
        cerr << "# [error] strange can not find surface type of "
             << one_patch[fi] << endl;
      }
      type2it[it->second].push_back(it);
    }
    size_t max = 0;
    size_t type = -1;
    for(boost::unordered_map<size_t,
        vector<boost::unordered_map<size_t,size_t>::iterator> >::const_iterator
        cit = type2it.begin(); cit != type2it.end(); ++cit){
      if(max < cit->second.size())   {
        type = cit->first;
        max = cit->second.size();
      }
    }

    for(boost::unordered_map<size_t,
        vector<boost::unordered_map<size_t,size_t>::iterator> >::const_iterator
        cit = type2it.begin(); cit != type2it.end(); ++cit){
      if(cit->first != type){
        for(size_t i = 0; i < cit->second.size(); ++i){
          cit->second[i]->second = type;
        }
      }
    }
  }
  return 0;
}

int update_surface_type_according_to_boundary(
    boost::unordered_map<size_t,size_t> &surface_type,
    const std::vector<vector<size_t> > &modified_boundary,
    const matrix<size_t> & outside_face,
    const matrix<size_t> & outside_face_idx)
{
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  boost::unordered_set<pair<size_t,size_t> > boundary_edges;
  {
    for(size_t bi = 0; bi < modified_boundary.size(); ++bi){
      const vector<size_t> & one_boundary = modified_boundary[bi];
      for(size_t pi = 0; pi < one_boundary.size()-1 ; ++pi){
        pair<size_t,size_t> one_edge(one_boundary[pi], one_boundary[pi+1]);
        if(one_edge.first > one_edge.second)
          swap(one_edge.first,one_edge.second);
        boundary_edges.insert(one_edge);
      }
    }
  }

  // breadth first travelling to spit the surface patch
  vector<vector<size_t> > patches;
  get_face_patches_according_to_boundary(outside_face,*ea,
                                         boundary_edges,patches);
  {// transfer the face index to real one
    for(size_t pi = 0; pi < patches.size(); ++pi){
      vector<size_t> & one_patch = patches[pi];
      for(size_t fi = 0; fi < one_patch.size(); ++fi){
        one_patch[fi] = outside_face_idx[one_patch[fi]];
      }
    }
  }

  retype_each_patch(surface_type, patches);
  //vector<bool> is_face_visited(outside_face.size(2),false);
  //std::stack<size_t> face_stack;


  return 0;
}

int optimize_patch_boundary(
    const matrix<size_t> & tet,
    const matrix<double> & node,
    const vector<matrixst> & surface_patches,
    const boost::unordered_set<pair<size_t,size_t> > & patch_boundary,
    boost::unordered_map<size_t,size_t> & surface_type)
{
  vector<deque<pair<size_t,size_t> > > chains;
  {
    vector<pair<size_t,size_t> > possible_edges(patch_boundary.begin(),
                                                patch_boundary.end());
    jtf::util::extract_chain_from_edges(possible_edges, chains);
  }

  vector<vector<size_t> > chain_link_patch(chains.size());

  for(size_t ci = 0; ci < chains.size(); ++ci){
    const size_t & p0 = chains[ci].front().first;
    const size_t & p1 = chains[ci].back().second;
    for(size_t pi = 0; pi < surface_patches.size(); ++pi){
      if(find(surface_patches[pi].begin(), surface_patches[pi].end(), p0) !=
         surface_patches[pi].end() &&
         find(surface_patches[pi].begin(), surface_patches[pi].end(), p1) !=
         surface_patches[pi].end()){
        chain_link_patch[ci].push_back(pi);
      }
    }
  }


  // collect all linking patch to each chain
  {// check
    for(size_t ci = 0; ci < chain_link_patch.size(); ++ci){
      if(chain_link_patch[ci].size() != 2){
        cerr << "# [error] strange, this chain links "
             << chain_link_patch[ci].size() << " patches." << endl;
        return __LINE__;
      }
    }
  }

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  matrixst outside_face, outside_face_idx;
  get_outside_face(*fa, outside_face);
  get_outside_face_idx(*fa, outside_face_idx);

  vector<boost::shared_ptr<vertex_connection<UNDIRECT> > > vvc;
  vvc.resize(surface_patches.size());
  for(size_t pi = 0; pi < surface_patches.size(); ++pi){
    map<pair<size_t,size_t>, double> edge_weight;
    unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
          jtf::mesh::edge2cell_adjacent::create(surface_patches[pi]));
    if(!ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }
    for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
      const pair<size_t,size_t> & one_edge = ea->edges_[ei];
      edge_weight[one_edge] = norm(node(colon(), one_edge.first) -
                                   node(colon(), one_edge.second));
    }
    vvc[pi].reset(vertex_connection<UNDIRECT>::create(edge_weight));
  }

  vector<vector<size_t> > modified_boundary(chain_link_patch.size());
  for(size_t ci = 0; ci < chain_link_patch.size(); ++ci){
    const deque<pair<size_t,size_t> > & one_chain = chains[ci];
    if(!is_large_distortion_chain(node, one_chain)) {
      for(size_t ei = 0; ei < one_chain.size(); ++ei){
        modified_boundary[ci].push_back(one_chain[ei].first);
      }
      modified_boundary[ci].push_back(one_chain.back().second);
      continue;
    }

    const size_t & patch_0 = chain_link_patch[ci].front();
    const size_t & patch_1 = chain_link_patch[ci].back();

    vector<size_t> path_0, path_1;

    vvc[patch_0]->get_shortest_path(
          one_chain.front().first, one_chain.back().second, path_0);
    vvc[patch_1]->get_shortest_path(
          one_chain.front().first, one_chain.back().second, path_1);
    double len_0 = calc_path_length(node, path_0);
    double len_1 = calc_path_length(node, path_1);

    if(len_0 < len_1){
      modified_boundary[ci] = path_0;
    }else
      modified_boundary[ci] = path_1;
  }


  update_surface_type_according_to_boundary(surface_type, modified_boundary,
                                            outside_face, outside_face_idx);

  return 0;
}


int cal_face_normal2D(const matrix<size_t> &tri,
                      const matrix<double> &uv,
                      matrix<double> &face_normal)
{
  face_normal = zeros<double>(3, tri.size(2));

  matrix<double> fake_uv = zeros<double>(3,2);
  for(size_t fi = 0; fi < tri.size(2); ++fi){
    fake_uv(colon(0,1),0) = uv(colon(),tri(1,fi)) - uv(colon(),tri(0,fi));
    fake_uv(colon(0,1),1) = uv(colon(),tri(2,fi)) - uv(colon(),tri(0,fi));

    face_normal(colon(),fi) =
        cross(fake_uv(colon(),0), fake_uv(colon(),1));
    double len = norm(face_normal(colon(),fi));
    if(len < 1e-6) len = 1;
    face_normal(colon(),fi) /= len;
  }
  return 0;
}

class node_equal_func : public function_t<double, int32_t>
{
public:
  node_equal_func(const matrix<double> &uv,
                  const size_t & p0,
                  const matrix<double> & des,
                  const double weight)
    :node_num_(uv.size(2)), p0_(p0), des_(des), w_(weight){}
  virtual size_t dim_of_x()const{
    return 2 * node_num_;
  }
  virtual size_t dim_of_f()const{
    return 2;
  }
  virtual int val(const double *x, double *f,
                  hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(2, node_num_, x);
    zjucad::matrix::itr_matrix<double*> f0(2,1,f);
    f0 = (T(colon(), p0_) - des_) * w_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    for(size_t di = 0; di < 2; ++di){
      ptr[di+1] = ptr[di] + 1;
      idx[ptr[di]+0] = 2 * p0_ + di;
      val[ptr[di]+0] = w_;
    }

    return 0;
  }
  virtual size_t jac_nnz() const{
    return 1 * dim_of_f();
  }

private:
  const size_t node_num_;
  const size_t p0_;
  const matrix<double> des_;
  const double w_;
};


hj::function::function_t<double, int32_t> *
build_anti_flip_laplacian_func(
    const matrix<double> &uv,
    const matrix<size_t> &tri,
    const jtf::mesh::edge2cell_adjacent & ea,
    const jtf::mesh::one_ring_point_at_point & p2p,
    const double weight)
{
  matrix<double> face_normal;

  cal_face_normal2D(tri,uv, face_normal);

  matrix<double> average_node = zeros<double>(2,1);
  matrix<double> edge_weight = zeros<double>(ea.edges_.size(), 1);

  boost::shared_ptr<vector<boost::shared_ptr<function_t<double,int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double,int32_t> > >);
  for(jtf::mesh::one_ring_point_at_point::p2p_type::const_iterator cit =
      p2p.p2p_.begin(); cit != p2p.p2p_.end(); ++cit){

    //const boost::unordered_set<size_t> & link_point = cit->second;
      const vector<size_t> & link_point = cit->second;

    average_node = zeros<double>(2,1);
//    for(boost::unordered_set<size_t>::const_iterator lit = link_point.begin();
//        lit != link_point.end(); ++lit){
    for(size_t pi = 0; pi < link_point.size(); ++pi){
        const size_t edge_idx = ea.get_edge_idx(cit->first, link_point[pi]);
      if(edge_idx == -1){
        cerr << "# [error] can not find edge " << cit->first
             << " " << link_point[pi] << endl;
        return 0;
      }
      if(fabs(edge_weight[edge_idx]) < 1e-6){
        const pair<size_t,size_t> & tri_pair = ea.edge2cell_[edge_idx];
        if(ea.is_boundary_edge(tri_pair)) {
          edge_weight[edge_idx] = 1.0;
          continue;
        }
        if(dot(face_normal(colon(),tri_pair.first),
               face_normal(colon(), tri_pair.second)) < 0)
          edge_weight[edge_idx] = -1.0;
        else
          edge_weight[edge_idx] = 1.0;
      }
      average_node += edge_weight[edge_idx] * uv(colon(),link_point[pi]);
    }
    average_node /= link_point.size();

    boost::shared_ptr<node_equal_func> nef(
          new node_equal_func(uv,cit->first, average_node, sqrt(weight)));
    funcs->push_back(nef);
  }

  return new_catenated_function<double, int32_t>(funcs);
}

int mapping_uv_to_triangle(matrix<double> &node,
                           const matrix<size_t> &tri,
                           const matrix<double> &uv,
                           const matrix<double> &orig_node,
                           const matrix<double> &uv_basis)
{
  for(size_t pi = 0; pi < node.size(2); ++pi){
    node(colon(), pi) = orig_node + uv(0,pi) * uv_basis(colon(),0)
                        + uv(1,pi) * uv_basis(colon(),1);
  }

  {
    jtf::mesh::save_obj("after_optimize.obj",tri, node);
  }
  return 0;
}

hj::function::function_t<double, int32_t> *
build_fix_position_func(const matrix<double> & uv,
                        const jtf::mesh::edge2cell_adjacent & ea,
                        const double w)
{
  boost::unordered_set<size_t> boundary_point;
  for(size_t ei = 0; ei < ea.edges_.size(); ++ei){
    const pair<size_t,size_t> & tri_pair = ea.edge2cell_[ei];
    const pair<size_t,size_t> & one_edge = ea.edges_[ei];
    if(ea.is_boundary_edge(tri_pair)){
      boundary_point.insert(one_edge.first);
      boundary_point.insert(one_edge.second);
    }
  }

  boost::shared_ptr<vector<boost::shared_ptr<function_t<double, int32_t> > > >
      funcs(new  vector<boost::shared_ptr<function_t<double, int32_t> > >);
  for(boost::unordered_set<size_t>::const_iterator cit = boundary_point.begin();
      cit != boundary_point.end(); ++cit){
    boost::shared_ptr<node_equal_func> nef(
          new node_equal_func(uv, *cit, uv(colon(),*cit), sqrt(w)));
    funcs->push_back(nef);
  }
  return new_catenated_function<double, int32_t>(funcs);
}


class unsigned_area_func : public hj::function::function_t<double, int32_t>
{
public:
  unsigned_area_func(const matrix<double> &uv,
                     const matrix<size_t> &tri,
                     const jtf::mesh::edge2cell_adjacent &ea,
                     const jtf::mesh::one_ring_face_at_point &orfap)
    :node_num_(uv.size(2)), tri_(tri), ea_(ea), orfap_(orfap){
    z_dir_ = zeros<double>(3,1);
    z_dir_[2] = 1.0;
  }

  virtual size_t dim_of_x()const{
    return 2 * node_num_;
  }
  virtual size_t dim_of_f()const{
    return 1;
  }
  virtual int val(const double *x, double *f,
                  hj::function::func_ctx *ctx = 0) const {
    zjucad::matrix::itr_matrix<const double *> T(2, node_num_, x);
    double area = 0.0;
    zjucad::matrix::matrix<double> e(2,3);
    for(size_t fi = 0; fi < tri_.size(2); ++fi){
      e(colon(),0) = T(colon(), tri_(1,fi)) - T(colon(), tri_(0,fi));
      e(colon(),1) = T(colon(), tri_(2,fi)) - T(colon(), tri_(0,fi));
      e(colon(),2) = T(colon(), tri_(2,fi)) - T(colon(), tri_(1,fi));

      const double a = norm(e(colon(),0));
      const double b = norm(e(colon(),1));
      const double c = norm(e(colon(),2));
      const double p = (a+b+c)/2.0;
      area += sqrt(p * (p-a) * (p-b) * (p-c));
    }

    *f = area;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double*> T(2, node_num_, x);
    matrix<double> tri_valid(tri_.size(2),1);
    matrix<double> fake_normal = zeros<double>(3,1), e1 = zeros<double>(3,1),
        e2 = zeros<double>(3,1);
    for(size_t fi = 0; fi < tri_.size(2); ++fi){
      e1(colon(0,1)) = T(colon(), tri_(1,fi)) - T(colon(), tri_(0,fi));
      e2(colon(0,1)) = T(colon(), tri_(2,fi)) - T(colon(), tri_(0,fi));
      fake_normal = cross(e1, e2);
      if(dot(fake_normal, z_dir_) > 0) tri_valid[fi] = 1.0;
      else
        tri_valid[fi] = -1.0;
    }

    const size_t positive_num = count(tri_valid.begin(), tri_valid.end(), 1.0);
    if(2*positive_num < tri_valid.size()) tri_valid *= -1.0;

    ptr[1] = ptr[0] + node_num_ * 2;
    matrix<double> edge_diff = zeros<double>(2,1);
    for(size_t pi = 0; pi < node_num_; ++pi){
      edge_diff = zeros<double>(2,1);

      jtf::mesh::one_ring_face_at_point::p2f_type::const_iterator cit
          = orfap_.p2f_.find(pi);
      if(cit == orfap_.p2f_.end()) continue;
      const vector<size_t> & one_ring_face = cit->second;
      for(size_t fi = 0; fi < one_ring_face.size(); ++fi){
        if(one_ring_face[fi] == -1) continue;
        const matrix<size_t> & one_face = tri_(colon(),one_ring_face[fi]);
        matrix<size_t>::const_iterator pcit = find(one_face.begin(), one_face.end(), pi);
        const size_t p_idx_in_tri = pcit - one_face.begin();
        const size_t next_idx = (p_idx_in_tri + 1)%one_face.size();
        const size_t nnext_idx = (p_idx_in_tri + 2)%one_face.size();

        edge_diff +=  tri_valid[fi] * (T(colon(),tri_(nnext_idx, one_ring_face[fi])) -
                                       T(colon(), tri_(next_idx, one_ring_face[fi])));
      }
      edge_diff /= -2.0;
      swap(edge_diff[0], edge_diff[1]);

      for(size_t di = 0; di < 2; ++di){
        idx[ptr[0]+ 2 * pi + di] = 2 * pi + di;
        val[ptr[0]+ 2 * pi + di] = edge_diff[di];
      }
    }

    return 0;
  }
  virtual size_t jac_nnz() const{
    return 2 * node_num_ * dim_of_f();
  }

private:
  matrix<double> z_dir_;
  const size_t node_num_;
  const matrix<size_t> & tri_;
  const jtf::mesh::edge2cell_adjacent & ea_;
  const jtf::mesh::one_ring_face_at_point & orfap_;
};

hj::function::function_t<double, int32_t>*
build_unsigned_area_func(const matrix<double> & uv,
                         const matrix<size_t> & tri,
                         const jtf::mesh::edge2cell_adjacent & ea,
                         const jtf::mesh::one_ring_face_at_point & orfap)
{
  return new unsigned_area_func(uv, tri, ea, orfap);
}


hj::function::function_t<double, int32_t> *
build_laplacian_func(const matrix<double> & uv,
                     const matrix<size_t> & tri,
                     const jtf::mesh::edge2cell_adjacent & ea,
                     const jtf::mesh::one_ring_point_at_point & orpap,
                     const jtf::mesh::one_ring_face_at_point & orfap)
{
  boost::shared_ptr<function_t<double,int32_t> > anti_flip_lap_func(
        build_anti_flip_laplacian_func(uv, tri, ea, orpap, lap_w));

  boost::shared_ptr<function_t<double, int32_t> > unsigned_area_func(
        build_unsigned_area_func(uv, tri, ea, orfap));

  boost::shared_ptr<function_t<double,int32_t> > fix_position_func(
        build_fix_position_func(uv, ea, fix_w_3));

  boost::shared_ptr<vector<boost::shared_ptr<function_t<double,int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double,int32_t> > >);

  if(fix_position_func.get()){
    funcs->push_back(fix_position_func);
  }

  //  if(anti_flip_lap_func.get()){
  //    funcs->push_back(anti_flip_lap_func);
  //  }

  if(unsigned_area_func.get()){
    funcs->push_back(unsigned_area_func);
  }

  return new_catenated_function<double,int32_t>(funcs);
}

int optimize_patch_boundary(ptree &pt)
{
  jtf::mesh::meshes trm;
  if(jtf::mesh::load_obj(pt.get<string>("obj.value").c_str(), trm.mesh_, trm.node_))
    return __LINE__;

  if(jtf::mesh::reorder_face(trm.mesh_)) {
    cerr << "# [error] can not reorder face." << endl;
    return __LINE__;
  }

  matrix<double> orig_node, uv_basis;
  matrix<double> uv;
  if(load_from_uv(pt.get<string>("uv.value").c_str(), trm, uv, &orig_node, &uv_basis))
    return __LINE__;

  double lap_w = 1;
  double fix_w_2 = 1;

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(trm.mesh_));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  unique_ptr<jtf::mesh::one_ring_point_at_point> orpap(
        jtf::mesh::one_ring_point_at_point::create(trm.mesh_));
  if(!orpap.get()){
    cerr << "# [error] can not build one_ring_point_at_point." << endl;
    return __LINE__;
  }

  jtf::mesh::one_ring_face_at_point orfap;
  orfap.add_all_faces(trm.mesh_, *ea);
  orfap.sort_int_loop(trm.mesh_, trm.node_);

  for(size_t iter = 0; iter < 5; ++iter){
    shared_ptr<function_t<double,int32_t> > anti_flip_lap_func(
          build_laplacian_func(uv, trm.mesh_,*ea, *orpap, orfap));

    unique_ptr<jtf::function::functionN1_t<double,int32_t> > lsgn(
          jtf::function::least_square_warpper(anti_flip_lap_func));

    jtf::optimize(*lsgn, uv,pt, nullptr,nullptr,nullptr);
  }


  {
    mapping_uv_to_triangle(trm.node_, trm.mesh_, uv, orig_node, uv_basis);
  }
  //tetmesh tm;
  //  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),
  //                               &tm.node, &tm.mesh))
  //    return __LINE__;

  //  boost::unordered_map<size_t,size_t> surface_type;
  //  if(load_surface_restricted_type_static(
  //       pt.get<string>("surface_type.value").c_str(), surface_type))
  //    return __LINE__;

  //  vector<matrixst> surface_patches;
  //  boost::unordered_set<pair<size_t,size_t> > patch_boundary;
  //  convert_surface_type_to_surface_patches(tm.mesh, surface_type, surface_patches,
  //                                          patch_boundary);


  //  optimize_patch_boundary(
  //        tm.mesh, tm.node,surface_patches,patch_boundary,surface_type);

  //  {
  //    unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh));
  //    dump_surface_restricted_type_to_vtk_static(
  //          "optimized_surface_type.vtk", "restriced_type", tm.node, *fa, surface_type);
  //  }
  return 0;
}
