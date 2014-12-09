/* detail control param
/usr/bin/time ../../bin/polycube prog=polycube_trimesh linear_solver/type=direct linear_solver/name=cholmod iter_w=20 output="a.tet" iter=100 tet=../../dat/kitty-4.8k.mesh-split-surface-tet.mesh normal_align_w=5e-2 L1_sqrt_eps=5e-1 epsg=1e-3 adj_normal_w=2e1
 */
#include <boost/property_tree/ptree.hpp>

#include <string>
#include <fstream>
#include <numeric>
#include <hjlib/function/func_aux.h>
#include <hjlib/math/polar.h>
#include <hjlib/sparse/sparse.h>


#include <zjucad/optimizer/optimizer.h>
#include <zjucad/ptree/ptree.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include "../tetmesh/util.h"
#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"
#include "../common/util.h"
#include "../common/IO.h"
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/optimizer/optimizer.h>

#include "polycube_surface_func.h"
#include "smooth_L1.h"
#include "polycube_surface_func.h"

#include "../mesh_func/tri-area.h"
#include "tet_func.h"
#include "../mesh_func/tri-normal.h"
#include "../mesh_func/tri-area-normal.h"

#include "util.h"
#include "polycube_surface_func.h"
#include "quality.h"
#include "../mesh_func/tri-area-normal-func.h"

using namespace std;
using boost::property_tree::ptree;
using namespace hj::function;
using namespace zjucad::matrix;

double boundary_align_w;
double laplacian_w;
extern double L1_sqrt_eps;
extern double arap_w;
extern double adj_normal_w;
double boundary_smooth_w;
double fix_inner_point_w;

function_t<double, int32_t> *trimesh_arap_for_diagnose = 0;

class cot_lap_function : public hj::function::function_t<double,int32_t>
{
public:
  cot_lap_function(const matrix<double> & node,
                   const size_t p0,
                   const vector<size_t> & adj_point,
                   const vector<double> & adj_weight,
                   const double weight = 1.0)
    :node_num_(node.size(2)), p0_(p0) , adj_point_(adj_point),
      adj_weight_(adj_weight), w_(weight){
    total_weight_ = std::accumulate(adj_weight_.begin(), adj_weight_.end(),0.0);
  }
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }
  virtual int val(const double*x, double *f, hj::function::func_ctx *ctx = 0)const{
    const zjucad::matrix::itr_matrix<const double*> T(3, node_num_, x);
    itr_matrix<double*> f0(3,1,f);

    f0 = T(colon(), p0_);

    for(size_t pi = 0; pi < adj_point_.size(); ++pi){
        f0 -= adj_weight_[pi] * T(colon(), adj_point_[pi])/total_weight_;
      }

    f0 *= w_;
    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    for(size_t di = 0; di < 3; ++di){
        ptr[di + 1] = ptr[di] + (1 + adj_point_.size());
        idx[ptr[di] + 0] = 3 * p0_ + di;
        val[ptr[di] + 0] = w_;

        for(size_t pi = 0; pi < adj_point_.size(); ++pi){
            idx[ptr[di] + pi + 1] = 3 * adj_point_[pi] + di;
            val[ptr[di] + pi + 1] = -1.0 * w_ * adj_weight_[pi]/total_weight_;
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return (1 + adj_point_.size()) * dim_of_f();
  }
private:
  const size_t node_num_;
  const size_t p0_;
  const double w_;
  const vector<size_t> adj_point_;
  double total_weight_;
  const vector<double> adj_weight_;
};

double calc_cot_weight(const size_t & p0,
                       const size_t & p1,
                       const matrix<size_t> &faces,
                       const matrix<double> &node,
                       const jtf::mesh::edge2cell_adjacent &ea)
{
  const size_t edge_idx = ea.get_edge_idx(p0, p1);
  if(edge_idx == -1) {
      cerr << "# [error] can not find edge " << p0 << " " << p1 << endl;
      return 0;
    }
  const pair<size_t,size_t> & tri_pair = ea.edge2cell_[edge_idx];
  vector<size_t> tri_pair_vec;
  tri_pair_vec.push_back(tri_pair.first);
  tri_pair_vec.push_back(tri_pair.second);
  matrix<double> e(3,2);
  vector<double> cot_val_vec;
  for(size_t i = 0; i < tri_pair_vec.size(); ++i){
      if(tri_pair_vec[i] == -1) continue;
      const size_t other_idx =
          std::accumulate(faces(colon(),tri_pair.first).begin(),
                          faces(colon(),tri_pair.first).end(),0)
          - p0 - p1;
      e(colon(),0) = node(colon(),p0) - node(colon(),other_idx);
      e(colon(),1) = node(colon(),p1) - node(colon(),other_idx);
      double len_0 = norm(e(colon(),0));
      if(len_0 < 1e-6) len_0 = 1.0;
      double len_1 = norm(e(colon(),1));
      if(len_1 < 1e-6) len_1 = 1.0;
      e(colon(),0) /= len_0;
      e(colon(),1) /= len_1;
      const double cos_val = dot(e(colon(),0), e(colon(),1));
      const double sin_val = sqrt(1 - cos_val * cos_val);
      const double cot_val = cos_val / sin_val;
      cot_val_vec.push_back(cot_val);
    }

  return std::accumulate(cot_val_vec.begin(), cot_val_vec.end(), 0.0) / cot_val_vec.size();
}

hj::function::function_t<double, int32_t> *
build_surface_smooth_func(const matrix<double> & node,
                          const matrix<size_t> &faces,
                          const double lap_w)
{
  //return new smooth_surface_func(node, faces, sqrt(lap_w));

  unique_ptr<jtf::mesh::one_ring_point_at_point> orpap(
        jtf::mesh::one_ring_point_at_point::create(faces));
  if(!orpap.get()){
      cerr << "# [error] can not build one_ring_point_at_point." << endl;
      return 0;
    }

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(faces));
  if(!ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return 0;
    }

  jtf::mesh::one_ring_face_at_point orfap;
  orfap.add_all_faces(faces, *ea);

  matrix<double> face_area = zeros<double>(faces.size(2),1);
  matrix<double> face_node = zeros<double>(3,3);
  for(size_t fi = 0; fi < face_area.size(); ++fi){
      face_node = node(colon(), faces(colon(),fi));
      calc_tri_area_(&face_area[fi], &face_node[0]);
    }

  const double total_area = std::accumulate(face_area.begin(), face_area.end(), 0.0);

  //    vector<size_t> arounding_points;
  vector<double> cot_weight;
  boost::shared_ptr<vector<boost::shared_ptr<function_t<double,int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double,int32_t> > >);

  for(boost::unordered_map<size_t, vector<size_t> >::const_iterator
      cit = orpap->p2p_.begin(); cit != orpap->p2p_.end(); ++cit){

      const vector<size_t> & arounding_points = cit->second;
      //        arounding_points.resize(link_point.size());

      //        copy(link_point.begin(), link_point.end(), arounding_points.begin());
      cot_weight.resize(arounding_points.size());

      for(size_t p = 0; p < arounding_points.size(); ++p){
          cot_weight[p] =  calc_cot_weight(cit->first, arounding_points[p],faces, node, *ea);
        }

      map<size_t,vector<size_t> >::const_iterator p2f_it =
          orfap.p2f_.find(cit->first);
      if(p2f_it == orfap.p2f_.end()) continue;

      const vector<size_t> & adj_faces = p2f_it->second;

      double around_face_area = 0;
      for(size_t fi = 0; fi < adj_faces.size(); ++fi){
          if(adj_faces[fi] == -1) continue;
          around_face_area += face_area[adj_faces[fi]];
        }
      around_face_area /= 3.0;

      boost::shared_ptr<function_t<double, int32_t> > pf(
            new cot_lap_function(node, cit->first, arounding_points,
                                 cot_weight,sqrt(around_face_area/total_area *lap_w)));
      funcs->push_back(pf);
#if 0 // validate jac
      {
        const double err = jac_err(*pf, &node[0]);
        if(err > 1e-5) {
            cerr << "large error in jac: " << err << endl;
          }
      }
#endif
      //boost::shared_ptr<function_t<double, int32_t> >
    }
  return new_catenated_function<double,int32_t>(funcs);
}


hj::function::function_t<double, int32_t> *
build_fix_inner_point(const matrix<double> &node,
                      const matrix<size_t> &faces,
                      const double fix_inner_point_w)
{
  set<size_t> surface_node(faces.begin(), faces.end());
  vector<bool> node_inner_flag(node.size(2), true);
  for(set<size_t>::const_iterator cit = surface_node.begin();
      cit != surface_node.end(); ++cit){
      node_inner_flag[*cit] = false;
    }

  boost::shared_ptr<vector<boost::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double, int32_t> > >);
  size_t p = 0;
  for(size_t pi = 0; pi < node_inner_flag.size(); ++pi){
      if(node_inner_flag[pi]){
          funcs->push_back(boost::shared_ptr<function_t<double, int32_t> >(
                             new fix_node_func(node, pi, sqrt(fix_inner_point_w))));
#if 0 // validate jac
          {
            const double err = jac_err(*(funcs->back()), &node[0]);
            if(err > 1e-5) {
                cerr << "large error in jac: " << err << endl;
              }
          }
#endif
        }
    }
  return new_catenated_function<double,int32_t>(funcs);
}


class edge_dir_smooth_func : public function_t<double,int32_t>
{
public:
  edge_dir_smooth_func(const matrix<double> & node,
                       const size_t &source_point,
                       const size_t pa,
                       const size_t pb,
                       const double weight)
    :node_num_(node.size(2)), point_source_(source_point), point_left_(pa), point_right_(pb), w_(weight){}
  virtual size_t dim_of_x()const{
    return 3 * node_num_;
  }
  virtual size_t dim_of_f() const{
    return 3;
  }
  virtual int val(const value_type *x, value_type *f, func_ctx *ctx) const{
    itr_matrix<const value_type* > x0(3, node_num_,x);
    itr_matrix<value_type*> f0(3,1,f);
    f0 = x0(colon(),point_source_) - (x0(colon(),point_left_) + x0(colon(),point_right_))/2.0;
    f0 *= w_;
    return 0;
  }
  virtual int jac(const value_type *x, value_type *val, int_type *ptr, int_type *idx, func_ctx *ctx) const {
    for(size_t  fi = 0; fi < 3; ++fi){
        ptr[fi+1] = ptr[fi] + 3;

        idx[ptr[fi] + 0] = 3 * point_source_ + fi;
        val[ptr[fi] + 0] = w_;

        idx[ptr[fi] + 1] = 3 * point_left_ + fi;
        val[ptr[fi] + 1] = -1/2.0*w_;

        idx[ptr[fi] + 2] = 3 * point_right_ + fi;
        val[ptr[fi] + 2] = -1/2.0*w_;
      }
    return 0;
  }
private:
  const size_t node_num_;
  const size_t point_source_;
  const size_t point_left_;
  const size_t point_right_;
  const size_t w_;
};


hj::function::function_t<double,int32_t> *
build_adj_edge_smooth_func(const matrix<double> & node,
                           const matrix<size_t> & boundary_edges,
                           const double adj_edge_w)
{
  boost::shared_ptr<vector<boost::shared_ptr<function_t<double,int32_t> > > > funcs(
        new vector<boost::shared_ptr<function_t<double,int32_t> > >());
  map<size_t,vector<size_t> > point2edge_idx;
  matrix<double> dis(boundary_edges.size(2),1);

  for(size_t ei = 0; ei < boundary_edges.size(2); ++ei){
      point2edge_idx[boundary_edges(0,ei)].push_back(ei);
      point2edge_idx[boundary_edges(1,ei)].push_back(ei);
      dis[ei] = norm(node(colon(), boundary_edges(0,ei)) - node(colon(), boundary_edges(1,ei)));
    }

  const double total_edge_len = std::accumulate(dis.begin(), dis.end(), 0.0);
  dis /= total_edge_len;

  vector<size_t> connect_points;
  size_t other_point = -1;
  for(map<size_t,vector<size_t> >::const_iterator cit = point2edge_idx.begin(); cit != point2edge_idx.end(); ++cit){
      const size_t source_point = cit->first;
      const vector<size_t> & connect_edges = cit->second;
      connect_points.resize(connect_edges.size());
      double weight = 0;

      for(size_t ei = 0; ei < connect_edges.size(); ++ei){
          assert(boundary_edges(0, connect_edges[ei]) == source_point ||
                 boundary_edges(1, connect_edges[ei]) == source_point);
          other_point = boundary_edges(0, connect_edges[ei]) + boundary_edges(1, connect_edges[ei])
              - source_point;
          connect_points[ei] = other_point;
          weight += dis[connect_edges[ei]];
        }
      weight *= adj_edge_w/2.0;
      for(size_t pi = 1; pi < connect_points.size(); ++pi){
          funcs->push_back(boost::shared_ptr<function_t<double,int32_t> >(
                             new edge_dir_smooth_func(node,source_point, connect_points[pi-1], connect_points[pi], weight)));
        }

    }
  if(funcs->empty()) return 0;

  return new_catenated_function<double,int32_t>(funcs);
}

hj::function::function_t<double, int32_t> *
build_polycube_function_pure_trimesh(const matrix<double> &node,
                                     const matrix<size_t> &faces,
                                     const matrix<size_t> &boundary_edges,
                                     const matrix<double> * iter_node = 0)
{
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(faces));
  if(!ea.get()){
      cerr  << "# [error] can not build edge2triangle adjacent." << endl;
      return 0;
    }


  //////////////////////////////////////////////////////////////////
  ////////////////  add function term  /////////////////////////////
  //  boost::shared_ptr<function_t<double, int32_t> > surface_smooth_func
  //      (build_surface_smooth_func(node, faces, laplacian_w));

  //////////////////////////////////////////////////////////////////
  ////////////////  add function term  /////////////////////////////

  boost::shared_ptr<function_t<double, int32_t> > surface_arap_func
      (build_surface_arap_func((iter_node?*iter_node:node), faces));
  //(build_surface_meshless_func(node, faces));

  trimesh_arap_for_diagnose = surface_arap_func.get();

  //    boost::shared_ptr<function_t<double, int32_t> > fix_first_node_func
  //            (build_fix_zero_node_func(node, node(colon(),0), 1.0));

  boost::shared_ptr<function_t<double, int32_t> > adj_normal(
        build_adj_normal_sub_func(node, faces, adj_normal_w));

  boost::shared_ptr<function_t<double, int32_t> > adj_edge_smooth(
        build_adj_edge_smooth_func(node, boundary_edges, boundary_smooth_w));


  boost::shared_ptr<vector<boost::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double, int32_t> > >);

  if(surface_arap_func.get() ){
      funcs->push_back(surface_arap_func);
      cerr << "# [info] add surface arap_func, weight " << endl;
    }

  //  if(fix_first_node_func.get() ){
  //    funcs->push_back(fix_first_node_func);
  //    cerr << "# [info] add fix_zeros node, weight " << 1.0 << endl;
  //  }

  if(adj_normal.get() && adj_normal_w > 0){
      funcs->push_back(adj_normal);
      cerr << "# [info] add adj_normal function, weight "
           << adj_normal_w << endl;
    }
  //    if(fix_first_node_func.get() ){
  //        funcs->push_back(fix_first_node_func);
  //        cerr << "# [info] add fix_zeros node, weight " << 1.0 << endl;
  //    }

  if(adj_edge_smooth.get() && boundary_smooth_w > 0){
      funcs->push_back(adj_edge_smooth);
      cerr << "# [info] add boundary_edge_smooth function, weight "
           << boundary_smooth_w << endl;
    }
  //  if(adj_normal.get() && adj_normal_w > 0){
  //    funcs->push_back(adj_normal);
  //    cerr << "# [info] add adj_normal function, weight "
  //         << adj_normal_w << endl;
  //  }

  if(funcs->empty()) return 0;
  return new_catenated_function<double, int32_t>(funcs);
}

hj::function::function_t<double, int32_t> *
build_fix_first_node_func(const matrix<double> &orig_node,
                          const size_t idx,
                          const double weight)
{
  return new fix_node(orig_node, idx, weight);
}

//hj::function::function_t<double, int32_t> *
//build_polycube_function_trimesh(
//        //    const matrix<size_t> & tet,
//        const matrix<double> & orig_node,
//        const matrix<size_t> &faces)
//{
//    unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
//                jtf::mesh::edge2cell_adjacent::create(faces));
//    if(!ea.get()){
//        cerr  << "# [error] can not build edge2triangle adjacent." << endl;
//        return __LINE__;
//    }


//    //////////////////////////////////////////////////////////////////
//    ////////////////  add function term  /////////////////////////////
//    //  boost::shared_ptr<function_t<double, int32_t> > tetmesh_arap(
//    //        build_tetmesh_arap_func(node, tet, arap_w));


//    boost::shared_ptr<function_t<double, int32_t> > surface_smooth_func
//            (build_surface_smooth_func(orig_node, faces, laplacian_w));

//    boost::shared_ptr<function_t<double, int32_t> > fix_first_node_func(
//                build_fix_first_node_func(orig_node, faces[0], 1.0));

//    // fix all inner points
//    boost::shared_ptr<function_t<double, int32_t> > fix_inner_point
//            (build_fix_inner_point(orig_node, faces, fix_inner_point_w));

//    matrix<double> weight;
//    boost::shared_ptr<function_t<double, int32_t> > adj_normal
//            (build_adj_normal_func(orig_node, faces, weight));

//    boost::shared_ptr<vector<boost::shared_ptr<function_t<double, int32_t> > > >
//            funcs(new vector<boost::shared_ptr<function_t<double, int32_t> > >);

//    if(surface_smooth_func.get()){
//        funcs->push_back(surface_smooth_func);
//        cerr << "# [info] add surface smooth_func, weight " << laplacian_w << endl;
//    }

//    if(fix_first_node_func.get()){
//        funcs->push_back(fix_first_node_func);
//        cerr << "# [info] add fix first node func, weight " << 1.0 << endl;
//    }

//    if(fix_inner_point.get() && non_zero(fix_inner_point_w)){
//        funcs->push_back(fix_inner_point);
//        cerr << "# [info] add fix inner node function, weight "
//             << fix_inner_point_w << endl;
//    }

//    if(adj_normal.get() && adj_normal_w > 0){
//        funcs->push_back(adj_normal);
//        cerr << "# [info] add adj_normal function, weight "
//             << adj_normal_w << endl;
//    }

//    return new_catenated_function<double, int32_t>(funcs);
//}

class edge_dir : public jtf::function::functionN1
{
public:
  edge_dir(const zjucad::matrix::matrix<double> &node,
           const zjucad::matrix::matrix<size_t> &edge,
           size_t xyz,
           const double weight)
    :edge_(edge), node_num_(node.size(2)), weight_(weight), xyz_(xyz) {
  }
  virtual size_t dim(void) const {
    return node_num_*3;
  }
  virtual int val(const double *x, double &v) {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> edge = T(zjucad::matrix::colon(), edge_);
    double f[3];

    calc_edge_dir_(f, &edge[0]);
    jtf::math::erase_nan_inf(f[xyz_]);
    v += f[xyz_] * weight_;
    return 0;
  }
  virtual int gra(const double *x, double *g) {
    double sp_g[6];
    int32_t idx[6];
    size_t nnz = 6;
    gra(x, nnz, sp_g, idx);
    for(size_t i = 0; i < 6; ++i)
      g[idx[i]] += sp_g[i];
    return 0;
  }
  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx) {
    if(g == 0) {
        nnz = 6;
        return 0;
      }
    if(nnz != 6) {
        cerr << "strange." << endl;
      }
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> edge =
        T(zjucad::matrix::colon(), edge_), jac(3, 6);
    //calc_tri_area_normal_jac_(&jac[0], &tri[0]);
    calc_edge_dir_jac_(&jac[0], &edge[0]);
    for(size_t i = 0; i < 6; ++i) {
        idx[i] = edge_[i/3]*3+i%3;
        g[i] = jac(xyz_, i)*weight_;
        jtf::math::erase_nan_inf(g[i]);
      }
    return 0;
  }
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha = 1) {
    format = 1;
    if(h == 0 && ptr == 0 && idx == 0) { // query nnz
        nnz = 6*6;
        return 0;
      }
    if(h == 0 && ptr != 0 && idx != 0) { // query pattern
        return 0;
        for(size_t ci = 0; ci < 6; ++ci) {
            const size_t var_ci = edge_[ci/3]*3+ci%3;
            ptr[var_ci+1] = ptr[var_ci]+6;
            for(size_t ri = 0; ri < 6; ++ri) {
                const size_t var_ri = edge_[ri/3]*3+ri%3;
                idx[ptr[var_ci]+ri] = var_ri;
              }
          }
        return 0;
      }
    if(h != 0 && ptr != 0 && idx != 0) { // accumulate
        return 0;
      }
    return __LINE__;
  }
  virtual int hes_block(const double *x, double *h, double alpha = 1) {return 0;}
private:
  const size_t node_num_, xyz_;
  const zjucad::matrix::matrix<size_t> edge_;
  const double weight_;
};

class edge_len_dir : public jtf::function::functionN1_t<double,int32_t>
{
public:
  edge_len_dir(const zjucad::matrix::matrix<double> &node,
               const zjucad::matrix::matrix<size_t> &edge,
               size_t xyz,
               const double weight)
    :edge_(edge), node_num_(node.size(2)), weight_(weight), xyz_(xyz) {
  }
  virtual size_t dim(void) const {
    return node_num_*3;
  }
  virtual int val(const double *x, double &v) {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> edge = T(zjucad::matrix::colon(), edge_);
    double f[3];

    calc_edge_len_dir_(f, &edge[0]);
    jtf::math::erase_nan_inf(f[xyz_]);
    v += f[xyz_] * weight_ ;
    return 0;
  }
  virtual int gra(const double *x, double *g) {
    double sp_g[6];
    int32_t idx[6];
    size_t nnz = 6;
    gra(x, nnz, sp_g, idx);
    for(size_t i = 0; i < 6; ++i)
      g[idx[i]] += sp_g[i];
    return 0;
  }
  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx) {
    if(g == 0) {
        nnz = 6;
        return 0;
      }
    if(nnz != 6) {
        cerr << "strange." << endl;
      }
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> edge =
        T(zjucad::matrix::colon(), edge_), jac(3, 6);
    //calc_tri_area_normal_jac_(&jac[0], &tri[0]);
    calc_edge_len_dir_jac_(&jac[0], &edge[0]);
    for(size_t i = 0; i < 6; ++i) {
        idx[i] = edge_[i/3]*3+i%3;
        g[i] = jac(xyz_, i)*weight_;
        jtf::math::erase_nan_inf(g[i]);
      }
    return 0;
  }
  virtual int hes(const double *x, size_t &nnz, size_t&format, double *h, int32_t *ptr, int32_t *idx, double alpha = 1) {
    format = 1;
    if(h == 0 && ptr == 0 && idx == 0) { // query nnz
        nnz = 0;
        return 0;
      }
    if(h == 0 && ptr != 0 && idx != 0) { // query pattern
        return 0;
        for(size_t ci = 0; ci < 6; ++ci) {
            const size_t var_ci = edge_[ci/3]*3+ci%3;
            ptr[var_ci+1] = ptr[var_ci]+6;
            for(size_t ri = 0; ri < 6; ++ri) {
                const size_t var_ri = edge_[ri/3]*3+ri%3;
                idx[ptr[var_ci]+ri] = var_ri;
              }
          }
        return 0;
      }
    if(h != 0 && ptr != 0 && idx != 0) { // accumulate
        return 0;
      }
    return __LINE__;
  }
  virtual int hes_block(const double *x, double *h, double alpha = 1) {
    return 0;
  }
private:
  const size_t node_num_, xyz_;
  const zjucad::matrix::matrix<size_t> edge_;
  const double weight_;
};

jtf::function::functionN1_t<double,int32_t> *
build_smooth_L1_boundary_trimesh(
    const matrix<double> &node,
    const matrix<size_t> &edges,
    const matrix<double> &length,
    double boundary_align_w)
{
  assert(edges.size(1) == 2);
  if(edges.size(2) == 0) return 0;
  const matrix<double> & edge_len = length;

  const double total_length = std::accumulate(edge_len.begin(), edge_len.end(), 0.0);
  cerr << "# [info] total boundary_length in build boundary L1: " << total_length << endl;

  const size_t en  = edges.size(2);
  vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > sum(en*3);
  const double weight = boundary_align_w / total_length;
  for(size_t ei = 0; ei < en; ++ei) {
      for(size_t i = 0; i < 3; ++i) {
          shared_ptr<jtf::function::functionN1_t<double,int32_t> > normal_func(new edge_len_dir(node, edges(colon(),ei), i, 1));
          sum[ei*3+i].reset(new smooth_L1(normal_func, weight, L1_sqrt_eps*edge_len[ei]));
        }
    }
  unique_ptr<jtf::function::functionN1_t<double,int32_t> > rtn(new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD>(sum));
  return rtn.release();
}

class edge_sum : public jtf::function::functionN1_t<double,int32_t>
{
public:
  edge_sum(size_t node_num,
           const matrix<size_t> &edges, double orig_len)
    :node_num_(node_num), edges_(edges), orig_len_(orig_len) {
  }
  virtual size_t dim(void) const { return node_num_*3; }
  //NOTE: add to v
  virtual int val(const double *x, double &v) {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    matrix<double> edge(3, 2);
    v = -orig_len_;
    for(size_t ei = 0; ei < edges_.size(2); ++ei) {
        edge = T(colon(),edges_(colon(),ei));
        double len;
        calc_edge_len_(&len, &edge[0]);
        v += len;
      }
    return 0;
  }
  //NOTE: add to v
  virtual int gra(const double *x, double *g) {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    matrix<double> edges(3, 2), g6(6);
    for(size_t ei = 0; ei < edges_.size(2); ++ei) {
        edges = T(colon(), edges_(colon(), ei));
        g6(colon()) = 0;
        calc_edge_len_jac_(&g6[0], &edges[0]);
        //calc_tri_area_jac_(&g9[0], &tri[0]);
        for(size_t i = 0; i < 6; ++i) {
            g[i%3+edges_(i/3, ei)*3] += g6[i];
          }
      }
    return 0;
  }
  virtual int gra(const double * x, size_t & nnz, double * val, int32_t * idx){
    cerr << "This gra is time consuming." << endl;
    if(val == 0 && idx == 0){
        nnz = node_num_ * 3;
      }else{
        itr_matrix<double*> val_(node_num_*3,1, val);
        val_ *= 3;
        gra(x,val);
        for(size_t i = 0 ;i < node_num_ * 3; ++i)
          idx[i] = i;
      }
    return 0;
  }
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha = 1) {
    cerr << "should not be here." << endl;
    exit(1);
    return 1;
  }
  virtual int hes_block(const double *x, double *h, double alpha = 1) {return 0;}
protected:
  const zjucad::matrix::matrix<size_t> &edges_;
  const size_t node_num_;
  double orig_len_;
};

// move surface point to the center of one-ring neighborhood by a little
void pre_smooth_trimesh(matrix<double> &node, const matrix<double> &faces)
{
  matrix<size_t> nb_num = zeros<size_t>(node.size(2), 1);
  matrix<double> ct = zeros<double>(3, node.size(2));
  for(size_t fi = 0; fi < faces.size(2); ++fi) {
      nb_num(faces(colon(), fi)) += 1;
      for(int i = 0; i < 3; ++i) {
          ct(colon(), faces(i, fi)) += node(colon(), faces((i+1)%3, fi)) + node(colon(), faces((i+2)%3, fi));
        }
    }
  const double alpha = 0.3;
  for(size_t ni = 0; ni < node.size(2); ++ni) {
      if(nb_num[ni] == 0) continue;
      node(colon(), ni) = (1-alpha)*node(colon(), ni) + alpha*ct(colon(), ni) / nb_num[ni] / 2;
    }
}

int polycube_obj(ptree & pt)
{
  using namespace jtf::mesh;
  meshes trm;
  if(jtf::mesh::load_obj(pt.get<string>("obj.value").c_str(), trm.mesh_, trm.node_))
    return __LINE__;

  remove_extra_node(trm.mesh_, trm.node_);

  cerr << "# read in obj " << trm.node_.size(2) << " " << trm.mesh_.size() << endl;

  pt.put("L1_sqrt_eps.desc", "L1 sqrt eps, default is 1e-1.");
  pt.put("normal_align_w.desc", "surface L1 weight, default is 1e2.");
  pt.put("adj_normal_w.desc", "surface smoothness weight, default is 0.");
  pt.put("iter_w.desc", "iterative weight, default is 1");
  pt.put("fix_zero_w.desc", "fix first node to zero, default is 1");
  pt.put("anti_flip_w.desc", "anti_flip_w, default is 1e-2.");
  pt.put("div_L.desc", "L1_sqrt_eps div L, default is 1.0");
  pt.put("boundary_smooth_w.desc", "smooth boundary points");

  L1_sqrt_eps = pt.get<double>("L1_sqrt_eps.value", 1);
  double normal_align_w = pt.get<double>("normal_align_w.value", 0.05);
  double anti_flip_w = pt.get<double>("anti_flip_w.value",1e-2);
  adj_normal_w = pt.get<double>("adj_normal_w.value", 0.0);
  const double div_L = pt.get<double>("div_L.value",1.0);
  laplacian_w = pt.get<double>("lap_w.value",1.0);
  boundary_align_w = pt.get<double>("boundary_w.value",1.0);
  const size_t iter_w = pt.get<size_t>("iter_w.value",1);
  boundary_smooth_w = pt.get<double>("boundary_smooth_w.value",0.0);

  const double epsg = pt.get<double>("epsg.value", 1e-6);
  const double epsf = pt.get<double>("epsf.value", 1e-3);
  matrix<double> areas = zeros<double>(trm.mesh_.size(2),1);
  for(size_t fi = 0; fi < trm.mesh_.size(2); ++fi){
      areas[fi] = jtf::mesh::cal_face_area(trm.mesh_(colon(),fi),trm.node_);
    }

  const double total_area = std::accumulate(areas.begin(), areas.end(), 0.0);

  cerr << "total_area: " << total_area << endl;
  matrix<double> R = eye<double>(3);
  ostringstream vtk_path_pref;
  vtk_path_pref << normal_align_w << "-" << L1_sqrt_eps << "-"
                << adj_normal_w << "-" << fix_inner_point_w << "-";

  matrix<double> node = trm.node_;

  if(zjucad::has("init.value", pt)){
      matrix<size_t> init_obj;
      matrix<double> init_obj_node;
      if(jtf::mesh::load_obj(pt.get<string>("init.value").c_str(), init_obj,
                             init_obj_node)) {
          cerr << "# [error] can not read init obj." << endl;
          return __LINE__;
        }
      if(init_obj.size() != trm.mesh_.size()){
          cerr << "# [error] mesh incompatible." << endl;
          return __LINE__;
        }
      if(norm(init_obj - trm.mesh_) > 1e-6){
          cerr << "# [error] mesh incompatible." << endl;
          return __LINE__;
        }
      node = init_obj_node;
      cerr << "# [info] use init value." << endl;
    }
  matrix<size_t> boundary_edges;
  matrix<double> boundary_edge_length;

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(trm.mesh_));
  if(!ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }

  get_boundary_edge(*ea, boundary_edges);

  matrix<double> edge_len = zeros<double>(boundary_edges.size(2),1);

  for(size_t ei = 0; ei < boundary_edges.size(2); ++ei){
      edge_len[ei] = norm(trm.node_(colon(), boundary_edges(0,ei))
                          - trm.node_(colon(), boundary_edges(1,ei)));
    }
  const double total_edge_len = std::accumulate(edge_len.begin(), edge_len.end(),0.0);
  cerr << "# [info] total edge length " << total_edge_len << endl;
  {
    if(boundary_edges.size()){
        ofstream ofs_boundary("init_boundary.vtk");
        line2vtk(ofs_boundary, &trm.node_[0], trm.node_.size(2),
            &boundary_edges[0], boundary_edges.size(2));
      }
  }

  const matrix<double> ct = trm.node_*ones<double>(trm.node_.size(2), 1)/trm.node_.size(2);

  for(size_t i = 0; i < iter_w; ++i){
      cerr << "# [info] iteration " << i << endl;
      matrix<double> areas_k = zeros<double>(trm.mesh_.size(2),1);
      for(size_t fi = 0; fi < trm.mesh_.size(2); ++fi){
          areas_k[fi] = jtf::mesh::cal_face_area(trm.mesh_(colon(),fi),node);
        }

      if(1) { // opt global R
          unique_ptr<jtf::function::functionN1_t<double,int32_t> > func(build_polycube_rot_func2(node, trm.mesh_,areas_k));
          jtf::optimize(*func, R, pt, nullptr, nullptr, nullptr);

          cout << R << endl;
          hj::polar3d p;
          p(R);
          cout << R << endl;
          node = temp(R*node);
        }

      shared_ptr<function_t<double, int32_t> > func
          (build_polycube_function_pure_trimesh(trm.node_, trm.mesh_, boundary_edges));

      cerr << "# [info] add polycube surface normal L1 function, wegiht "
           << normal_align_w << endl;
      cerr << "# [info] L1_sqrt_eps " << L1_sqrt_eps  << endl;


      vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > sum;
      vector<pair<jtf::function::functionN1_t<double,int32_t> *, double> > wf;
      sum.push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(jtf::function::least_square_warpper(func)));
      shared_ptr<jtf::function::functionN1_t<double,int32_t> > tri_for_diagnose
          (jtf::function::least_square_warpper(*trimesh_arap_for_diagnose));

      wf.push_back(make_pair(tri_for_diagnose.get(), 1));
      if(normal_align_w > 0) {
          matrix<double> areas_k = zeros<double>(trm.mesh_.size(2),1);
          for(size_t fi = 0; fi < trm.mesh_.size(2); ++fi){
              areas_k[fi] = jtf::mesh::cal_face_area(trm.mesh_(colon(),fi),node);
            }
          sum.push_back(
                shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                  build_smooth_L1_area_normal(node, trm.mesh_, areas_k, normal_align_w)));
          wf.push_back(make_pair(sum[1].get(), normal_align_w));
        }
      if(boundary_align_w > 0 && boundary_edges.size() > 0){
          matrix<double> edge_k = zeros<double>(boundary_edges.size(2),1);
          for(size_t ei = 0; ei < boundary_edges.size(2); ++ei){
              edge_k[ei] = norm(node(colon(), boundary_edges(0,ei)) -
                                node(colon(), boundary_edges(1,ei)));
            }
          shared_ptr<jtf::function::functionN1_t<double,int32_t> > pf(
                build_smooth_L1_boundary_trimesh(node, boundary_edges, edge_k, boundary_align_w));
          if(pf.get()){
              sum.push_back(pf);
              cerr << "# [info] boundary align weight " << boundary_align_w << endl;
            }
        }

      if(anti_flip_w > 0) {// use unflip
          matrix<double> weight;
          shared_ptr<function_t<double, int32_t> > adj_normal
              (build_adj_normal_func(trm.node_, trm.mesh_, weight));
          weight *= anti_flip_w;
          sum.push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                          jtf::function::neg_log_warpper(adj_normal, weight)));
          cerr << "# use anti-flip with weight: " << anti_flip_w << endl;
        }

      shared_ptr<jtf::function::functionN1_t<double,int32_t> > target(new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD>(sum));

      wf.push_back(make_pair(target.get(), (1+normal_align_w)));

      vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > constraint_sum;
      constraint_sum.push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> > (new area_sum(node.size(2), trm.mesh_, total_area)));
      if(boundary_align_w > 0 && boundary_edges.size() > 0)
        constraint_sum.push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> > (new edge_sum(node.size(2), boundary_edges, total_edge_len)));
      shared_ptr<jtf::function::functionN1_t<double,int32_t> > constraint(new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD>(constraint_sum));

      unnormalized_normal_quality_checker cb(node, trm.mesh_, wf, *constraint);

      pt.put("epsg.value", epsg*(1+normal_align_w));

      if(normal_align_w > 0){
        vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > constraint_vec;
        constraint_vec.push_back(constraint);
          jtf::optimize(*target, node, pt, nullptr, &constraint_vec, &cb);
      }else
        jtf::optimize(*target, node, pt, nullptr, nullptr, nullptr);

      const matrix<double> cur_ct = node*ones<double>(trm.node_.size(2), 1)/trm.node_.size(2);
      node += (ct-cur_ct)*ones<double>(1, trm.node_.size(2));

      { // visualize
        ostringstream vtk_surface_path;
        vtk_surface_path << vtk_path_pref.str() << laplacian_w << "-" << i << ".surface.vtk";
        ofstream ofs(vtk_surface_path.str().c_str());
        tri2vtk(ofs, &node[0], node.size(2), &trm.mesh_[0], trm.mesh_.size(2));

        if(boundary_edges.size()){
            ostringstream vtk_bounary_path;
            vtk_bounary_path << vtk_path_pref.str() << laplacian_w << "-" << i << ".boundary.vtk";
            ofstream ofs_b(vtk_bounary_path.str().c_str());
            line2vtk(ofs_b, &node[0], node.size(2), &boundary_edges[0], boundary_edges.size(2));
          }
      }

      if(polycube_L1_area_quality(&node[0], node.size(2), trm.mesh_) < 1e-3) { // the global condition
          cout << "global converge." << endl;
          break;
        }

      normal_align_w *= 1.5;
      boundary_align_w *=1.5;

      L1_sqrt_eps /= sqrt(1.5);
      if(L1_sqrt_eps < 1e-2)
        L1_sqrt_eps = 1e-2;

    }

  if(jtf::mesh::save_obj(pt.get<string>("output.value").c_str(), trm.mesh_ ,node))
    return __LINE__;

  cerr << "success." << endl;
  return 0;
}
