#include <iostream>
#include <set>
#include <fstream>
#include <numeric>

#include <zjucad/linear_solver/linear_solver.h>
#include <hjlib/sparse/sparse.h>
#include <hjlib/sparse/operation.h>

#include <zjucad/ptree/ptree.h>
#include <zjucad/optimizer/optimizer.h>
#include <zjucad/matrix/matrix.h>

#include <hjlib/function/function.h>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/math/math.h>
#include <jtflib/algorithm/gauss_elimination.h>

#include "../tetmesh/util.h"

#include "../numeric/util.h"
#include "../tetmesh/tetmesh.h"
#include "../common/transition_type.h"
#include "../common/vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace hj::function;
using namespace boost::property_tree;

static bool check_quaternion(const matrix<double> & quat)
{
  for(size_t i = 0 ; i < quat.size(2); ++i)
    {
      if(quat(3,i) < 0) {
          cerr << "# [error] quat " << i << " has negative qw." << endl;
          return false;
        }
    }
  return true;
}
static void fix_rot(matrix<double> & rot)
{
  assert(rot.size(1) == 3);
  assert(rot.size(2) == 3);
  vector<pair<double, int> > trace(24);
  for(size_t i = 0; i < 24; ++i){ // 24 transition type
      matrix<double> transition_rot = type_transition2(i);
      transition_rot = temp(transition_rot * rot);
      trace[i].first = transition_rot(0,0) + transition_rot(1,1) + transition_rot(2,2);
      trace[i].second = i;
    }
  sort(trace.begin(), trace.end());
  rot = temp(type_transition2(trace.back().second) * rot);
}

int cal_best_aligned_axes(const matrix<double> & point_normal,
                          matrix<double> &axis)
{
  static const matrix<double> eye_axes = eye<double>(3);
  vector<pair<double, int> > err_idx(6);
  for(size_t i = 0; i < err_idx.size(); ++i){
      err_idx[i].first = dot(point_normal, eye_axes(colon(), i /2)) * (i%2==0?1.0:-1.0);
      err_idx[i].second =i;
    }
  sort(err_idx.begin(), err_idx.end());
  axis = eye_axes(colon(), err_idx.back().second/2)
      * (err_idx.back().second%2==0?1.0:-1.0);
  return 0;
}

int cal_quaternions(const matrix<size_t> & outside_face,
                    const matrix<double> & point_normal,
                    matrix<double> & quaternions)
{
  // init quat as (0,0,0,1)^T for each colon
  quaternions = ones<double>(4, point_normal.size(2));
  quaternions(colon(0,2),colon()) *= 0;

  set<size_t> outside_points(outside_face.begin(), outside_face.end());
  matrix<double> axis = zeros<double>(3,1), rot = eye<double>(3);
  double zyx[3];

  for(set<size_t>::const_iterator cit = outside_points.begin();
      cit != outside_points.end(); ++cit){
      cal_best_aligned_axes(point_normal(colon(),*cit),axis);
      matrix<double> rotate_axis = cross(point_normal(colon(),*cit), axis);
      if(norm(rotate_axis) < 1e-6){ // point_normal is parallel with axis
          if(dot(point_normal(colon(),*cit), axis) > 0){ // the same direction, rot = I
              rot = eye<double>(3);
            }else{
              cerr << "# [error] It's impossible a closest axis will"
                      "be opposite with this direction." << endl;
              return __LINE__;
            }
        }else{
          // assume point_normal, axis are normalized
          const double angle = jtf::math::cos2deg(dot(point_normal(colon(),*cit), axis));
          from_angle_to_rotation_matrix(angle*boost::math::constants::pi<double>()
                                        /180.0, rotate_axis, rot);
          fix_rot(rot);
          convert_rotation_matrix_to_quat(rot, &quaternions(0,*cit));
          //          if(check_quaternion(quaternions(colon(), *cit)) == false)
          //            cerr << "# [error] has not fix" << endl;
        }
    }

  // check_quaternion(quaternions);
  return 0;
}

class lap_function : public hj::function::function_t<double,int32_t>
{
public:
  lap_function(const set<size_t>  & one_ring_p2p,
               const size_t idx,
               const size_t num_x)
    : one_ring_p2p_(one_ring_p2p), idx_(idx), num_x_(num_x){}

  virtual size_t dim_of_x(void) const {
    return num_x_;
  }
  virtual size_t dim_of_f(void) const {
    return 1;
  }

  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const itr_matrix<const double*> x0(num_x_,1, x);
    *f= 0;
    for(set<size_t>::const_iterator scit = one_ring_p2p_.begin();
        scit != one_ring_p2p_.end(); ++scit){
        *f-= x0[*scit];
      }
    *f /= one_ring_p2p_.size();
    *f += x0[idx_];
    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {

    ptr[1] = ptr[0] + (1 + one_ring_p2p_.size());

    idx[ptr[0] + 0] = idx_;
    val[ptr[0] + 0] = 1.0;

    size_t i = 0;
    for(set<size_t>::const_iterator scit = one_ring_p2p_.begin();
        scit != one_ring_p2p_.end(); ++scit,++i){
        idx[ptr[0] + i + 1] = *scit;
        val[ptr[0] + i + 1] = -1.0/one_ring_p2p_.size();
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return one_ring_p2p_.size() + 1 ;
  }
private:
  const set<size_t> & one_ring_p2p_;
  const size_t idx_;
  const size_t num_x_;
};

class fix_function : public hj::function::function_t<double,int32_t>
{
public:
  fix_function(const double & fix_c,
               const size_t & num_x,
               size_t idx,
               const double w = 1.0)
    : fix_c_(fix_c), idx_(idx), num_x_(num_x), w_(w){}

  virtual size_t dim_of_x(void) const {
    return num_x_;
  }
  virtual size_t dim_of_f(void) const {
    return 1;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const itr_matrix<const double*> x0(num_x_,1, x);
    *f = (x0[idx_] - fix_c_) * w_;
    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    ptr[1] = ptr[0] + 1;
    idx[ptr[0] + 0] = idx_;
    val[ptr[0] + 0] = 1.0 * w_;
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 1 ;
  }
private:
  const double fix_c_;
  const size_t idx_;
  const size_t num_x_;
  const double w_;
};

hj::function::function_t<double,int32_t> *
build_lap_function(const map<size_t, set<size_t> > &one_ring_point_of_p)
{
  std::shared_ptr<vector<std::shared_ptr<function_t<double,int32_t> > > >
      func(new vector<std::shared_ptr<function_t<double,int32_t> > >);

  for(map<size_t, set<size_t> >::const_iterator cit = one_ring_point_of_p.begin();
      cit != one_ring_point_of_p.end(); ++cit){
      func->push_back(std::shared_ptr<function_t<double,int32_t> >(
                        new lap_function(cit->second, cit->first, one_ring_point_of_p.size())));
    }
  return new_catenated_function<double,int32_t>(func);
}

hj::function::function_t<double,int32_t> *
build_fix_function(const map<size_t,matrix<double> > & fix_c,
                   const size_t di, const size_t num_x, const double w = 1.0)
{
  std::shared_ptr<vector<std::shared_ptr<function_t<double,int32_t> > > >
      func(new vector<std::shared_ptr<function_t<double,int32_t> > >());
  for(map<size_t, matrix<double> >::const_iterator cit = fix_c.begin();
      cit != fix_c.end(); ++cit){
      func->push_back(std::shared_ptr<function_t<double,int32_t> >(
                        new fix_function(cit->second[di], num_x, cit->first, sqrt(w))));
    }
  return new_catenated_function<double,int32_t>(func);
}

int propagate_quaternious_linear_solver(const matrix<size_t> & tet,
                                        const matrix<double> & node,
                                        const matrix<size_t> &outside_face,
                                        const map<size_t, set<size_t> > & one_ring_point_of_p,
                                        matrix<double> &quaternions)
{
  // all surface quaternious are fixed, and inner are propagated
  set<size_t> surface_points(outside_face.begin(), outside_face.end());

  map<size_t, matrix<double> > fix_c;
  {
    std::shared_ptr<jtf::mesh::edge2cell_adjacent> ea(
          jtf::mesh::edge2cell_adjacent::create(outside_face));
    if(!ea.get()){
        cerr << "# [error] can not build edge2cell_adjacent." << endl;
        return __LINE__;
      }
    matrix<double> face_normal;
    jtf::mesh::cal_face_normal(outside_face, node, face_normal);
    set<size_t> feature_points;
    for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
        const pair<size_t,size_t> & one_edge = ea->edges_[ei];
        const pair<size_t,size_t> & face_pair = ea->edge2cell_[ei];
        assert(face_pair.first != -1 && face_pair.second != -1);
        if(jtf::math::cos2deg(dot(face_normal(colon(), face_pair.first),
                       face_normal(colon(), face_pair.second))) > 45){ // if angle > 45
            feature_points.insert(one_edge.first);
            feature_points.insert(one_edge.second);
          }
      }

    for(set<size_t>::const_iterator cit = feature_points.begin();
        cit != feature_points.end(); ++cit){
        set<size_t>::iterator it =
            find(surface_points.begin(), surface_points.end(), *cit) ;
        if(it != surface_points.end())
          surface_points.erase(it);
      }

    for(set<size_t>::const_iterator cit = surface_points.begin();
        cit != surface_points.end(); ++cit)
      fix_c[*cit] = quaternions(colon(), *cit);
  }

  ptree pt;

  // Lx=0
  // Bx=C ==> Ax=C
  //A^TAx=A^TC

  matrix<double> C;

  hj::sparse::csc<double,int32_t> AT;
  const size_t num_x = quaternions.size(2);
  const size_t num_f = one_ring_point_of_p.size() + fix_c.size();
  C = zeros<double>(num_f, 4);
  size_t num_nnz = fix_c.size();
  for(map<size_t, set<size_t> >::const_iterator cit = one_ring_point_of_p.begin();
      cit != one_ring_point_of_p.end(); ++cit){
      num_nnz += cit->second.size() + 1;
    }
  AT.resize(num_x, num_f, num_nnz);
  size_t i = 0;
  for(map<size_t, set<size_t> >::const_iterator cit = one_ring_point_of_p.begin();
      cit != one_ring_point_of_p.end(); ++cit, ++i){
      const set<size_t> & linked_points = cit->second;
      AT.ptr()[i+1] = AT.ptr()[i] + linked_points.size() + 1;
      AT.idx()[AT.ptr()[i] + 0] = cit->first;
      AT.val()[AT.ptr()[i] + 0] = 1.0;

      size_t j = 1;
      for(set<size_t>::const_iterator scit = linked_points.begin();
          scit != linked_points.end(); ++scit, ++j){
          AT.idx()[AT.ptr()[i] + j] = *scit;
          AT.val()[AT.ptr()[i] + j] = -1.0/linked_points.size();
        }
    }
  for(map<size_t, matrix<double> >::const_iterator cit = fix_c.begin();
      cit != fix_c.end(); ++cit, ++i){
      AT.ptr()[i+1] = AT.ptr()[i] + 1;
      AT.idx()[AT.ptr()[i]] = cit->first;
      AT.val()[AT.ptr()[i]] = 1.0;

      C(cit->first,colon()) = trans(cit->second);
    }

  hj::sparse::csc<double,int32_t> ATA;
  hj::sparse::AAT<hj::sparse::map_by_unsorted_vector > aat(AT,ATA);
  matrix<double> ATC = zeros<double>(AT.rows_, C.size(2));
  hj::sparse::mm(false, AT, C, ATC);

  std::shared_ptr<linear_solver> ls(linear_solver::create(&ATA.val()[0],
                                    &ATA.idx()[0], &ATA.ptr()[0], ATA.nnz(), ATA.size(1), ATA.size(2), pt));
  if(!ls.get()){
      cerr << "# [error] can not build linear solver." << endl;
      return __LINE__;
    }

  matrix<double> x = zeros<double>(quaternions.size(2),1);
  for(size_t i = 0; i < quaternions.size(1); ++i){
      const matrix<double> ci = ATC(colon(),i);
      ls->solve(&ci[0],&x[0],1,pt);
      quaternions(i,colon()) = trans(x);
    }

  pt.put("linear_solver/type.value", "direct");
  pt.put("linear_solver/name.value", "chomod");

  for(size_t i = 0; i < quaternions.size(2); ++i){
      const double len = norm(quaternions(colon(),i));
      if(len > 1e-6)
        quaternions(colon(),i) /= len;
      else
        {
          cerr << "# [error] degenerated quaternions." << endl;
        }
    }
  return 0;
}

int propagate_quaternious(const matrix<size_t> & tet,
                          const matrix<double> & node,
                          const matrix<size_t> &outside_face,
                          const map<size_t, set<size_t> > & one_ring_point_of_p,
                          matrix<double> &quaternions)
{
  // all surface quaternious are fixed, and inner are propagated
  set<size_t> surface_points(outside_face.begin(), outside_face.end());

  map<size_t, matrix<double> > fix_c;
  {
    std::shared_ptr<jtf::mesh::edge2cell_adjacent> ea(
          jtf::mesh::edge2cell_adjacent::create(outside_face));
    if(!ea.get()){
        cerr << "# [error] can not build edge2cell_adjacent." << endl;
        return __LINE__;
      }
    matrix<double> face_normal;
    jtf::mesh::cal_face_normal(outside_face, node, face_normal);
    set<size_t> feature_points;
    for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
        const pair<size_t,size_t> & one_edge = ea->edges_[ei];
        const pair<size_t,size_t> & face_pair = ea->edge2cell_[ei];
        assert(face_pair.first != -1 && face_pair.second != -1);
        if(jtf::math::cos2deg(dot(face_normal(colon(), face_pair.first),
                       face_normal(colon(), face_pair.second))) > 45){ // if angle > 45
            feature_points.insert(one_edge.first);
            feature_points.insert(one_edge.second);
          }
      }

    for(set<size_t>::const_iterator cit = feature_points.begin();
        cit != feature_points.end(); ++cit){
        set<size_t>::iterator it =
            find(surface_points.begin(), surface_points.end(), *cit) ;
        if(it != surface_points.end())
          surface_points.erase(it);
      }

    for(set<size_t>::const_iterator cit = surface_points.begin();
        cit != surface_points.end(); ++cit)
      fix_c[*cit] = quaternions(colon(), *cit);
  }

  ptree pt;
  //    pt.put("package.value", "hj");
  //    pt.put("alg.value", "More");
  //    pt.put("GS_LM_mu.value", "0.05");
  //    pt.put("iter.value", 30);
  pt.put("package.value","alglib");
  pt.put("alg.value", "lbfgs");
  pt.put("lbfgs-len.value", "7");
  pt.put("iter.value", 500);

  std::shared_ptr<hj::function::function_t<double,int32_t> > lap_func(
        build_lap_function(one_ring_point_of_p));

  const double fix_w = 1;
  matrix<double> temp_node;
  for(size_t di = 0; di < 4; ++di){
      std::shared_ptr<hj::function::function_t<double,int32_t> > fix_func(
            build_fix_function(fix_c, di, quaternions.size(2), fix_w));

      std::shared_ptr<vector< std::shared_ptr<const function_t<double,int32_t> > > >
          all_funcs(new vector< std::shared_ptr<const function_t<double,int32_t> > >);

      all_funcs->push_back(lap_func);
      all_funcs->push_back(fix_func);

      std::shared_ptr<const function_t<double,int32_t> >
          func(new_catenated_function<double, int32_t>(all_funcs));

      temp_node = quaternions(di,colon());

      matrix<double> residual = zeros<double>(func->dim_of_f());
      zjucad::optimize(*func, temp_node, residual, pt);

      quaternions(di,colon()) = temp_node;
    }

  for(size_t i = 0; i < quaternions.size(2); ++i){
      const double len = norm(quaternions(colon(),i));
      if(len > 1e-6)
        quaternions(colon(),i) /= len;
      else
        {
          cerr << "# [error] degenerated quaternions." << endl;
        }
    }
  return 0;
}



int quaternions2rot_matrix(const matrix<double> & quat,
                           matrix<matrix<double> > &rot)
{
  rot.resize(quat.size(2),1);
  for(size_t i = 0; i < quat.size(2); ++i){
      rot[i].resize(3,3);
      convert_quat_to_rotation_matrix(&quat(0,i), rot[i]);
    }
  return 0;
}

class vertex_update_func: public hj::function::function_t<double,int32_t>
{
public:
  vertex_update_func(const matrix<double> & orig_node,
                     const size_t & idx,
                     const set<size_t> & connect_points,
                     const matrix<matrix<double> > &rot)
    :orig_node_(orig_node), idx_(idx), connect_points_(connect_points),
      rot_(rot){}
  virtual ~vertex_update_func(){}

  virtual size_t dim_of_x(void) const {
    return orig_node_.size();
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }

  // vi - \sum vj/N = \sum (\nabla fi+\nabla fj)/2(\bar{vi}-\bar{vj})/N
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const itr_matrix<const double*> x0(orig_node_.size(1), orig_node_.size(2), x);
    itr_matrix<double*> f0(3,1,f);
    f0 = x0(colon(), idx_);

    matrix<double> avg_connect_node = zeros<double>(3,1);
    matrix<double> avg_orig_node = zeros<double>(3,1);

    for(set<size_t>::const_iterator cit = connect_points_.begin();
        cit != connect_points_.end(); ++cit){
        avg_connect_node += x0(colon(), *cit);

        avg_orig_node += (rot_[idx_] + rot_[*cit])/2.0
            * (orig_node_(colon(), idx_) - orig_node_(colon(), *cit));
      }
    avg_connect_node /= connect_points_.size();
    avg_orig_node /= connect_points_.size();

    f0 -= avg_connect_node;
    f0 -= avg_orig_node;

    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    for(size_t fi = 0; fi < 3; ++fi){
        ptr[fi+1] = ptr[fi] + connect_points_.size() + 1;
        idx[ptr[fi]+0] = 3 * idx_ + fi;
        val[ptr[fi]+0] = 1.0;

        size_t add_item = 1;
        for(set<size_t>::const_iterator cit = connect_points_.begin();
            cit != connect_points_.end(); ++cit, ++add_item){
            idx[ptr[fi] + add_item] = 3 * (*cit) + fi;
            val[ptr[fi] + add_item] = -1.0 / connect_points_.size();
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 3* (1 + connect_points_.size());
  }

private:
  const matrix<double> & orig_node_;
  const size_t idx_;
  const set<size_t> & connect_points_;
  const matrix<matrix<double> > & rot_;
};

class fix_point: public hj::function::function_t<double,int32_t>
{
public:
  fix_point(const size_t & node_number,
            const size_t & idx,
            const matrix<double> & fix_position)
    :node_number_(node_number), idx_(idx), fix_position_(fix_position){}

  virtual ~fix_point(){}

  virtual size_t dim_of_x(void) const {
    return 3 * node_number_;
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const itr_matrix<const double*> x0(3, node_number_, x);
    itr_matrix<double*> f0(3,1,f);
    f0 = x0(colon(), idx_) - fix_position_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    for(size_t fi = 0; fi < 3; ++fi){
        ptr[fi+1] = ptr[fi] + 1;
        idx[ptr[fi]+0] = 3 * idx_ + fi;
        val[ptr[fi]+0] = 1.0;
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 3;
  }

private:
  const matrix<double> fix_position_;
  const size_t idx_;
  const size_t node_number_;
};

class surface_smooth: public hj::function::function_t<double,int32_t>
{
public:
  surface_smooth(const size_t & node_number,
                 const size_t & idx,
                 const vector<size_t> & connected_points,
                 const double weight)
    :node_number_(node_number), idx_(idx), cp_(connected_points),
      w_(weight){}

  virtual size_t dim_of_x(void) const {
    return 3 * node_number_;
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const itr_matrix<const double*> x0(3, node_number_, x);
    itr_matrix<double*> f0(3,1,f);
    f0 = x0(colon(),idx_);
    for(size_t pi = 0; pi < cp_.size(); ++pi){
        f0 -= x0(colon(), cp_[pi])/cp_.size();
      }
    f0 *= w_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    for(size_t di = 0; di < 3; ++di){
        ptr[di+1] = ptr[di] + cp_.size() + 1;
        idx[ptr[di]] = 3 * idx_ + di;
        val[ptr[di]] = 1.0 * w_;

        for(size_t pi = 0; pi < cp_.size(); ++pi){
            idx[ptr[di] + pi + 1] = 3 * cp_[pi] + di;
            val[ptr[di] + pi + 1] = -1.0/cp_.size() * w_;
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 3*(1+cp_.size());
  }

private:
  const vector<size_t> cp_;
  const size_t idx_;
  const size_t node_number_;
  const double w_;
};

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


hj::function::function_t<double,int32_t> *
build_surface_smooth_func(
    const matrix<double> & node,
    const boost::unordered_map<size_t, boost::unordered_set<size_t> > & p2p,
    const double weight)
{ // smooth surface
  std::shared_ptr<vector<std::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<std::shared_ptr<function_t<double, int32_t> > >);

  vector<size_t> connect_points;
  for(boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
      cit = p2p.begin();  cit != p2p.end(); ++cit){
      connect_points.resize(cit->second.size());

      copy(cit->second.begin(), cit->second.end(), connect_points.begin());

      funcs->push_back(
            std::shared_ptr<function_t<double,int32_t> >(
              new surface_smooth(node.size(2), cit->first, connect_points, sqrt(weight))));
    }
  return new_catenated_function<double,int32_t>(funcs);
}
int update_vertex(const matrix<size_t> &tet,
                  const matrix<double> &orig_node,
                  matrix<double> &node,
                  const matrix<matrix<double> > &rotation_gradient,
                  const boost::unordered_map<size_t, boost::unordered_set<size_t> > &orpap)
{
  map<size_t, set<size_t> > one_ring_point_of_p;
  for(size_t ti = 0; ti < tet.size(2); ++ti)
    for(size_t pi = 0; pi < tet.size(1); ++pi){
        for(size_t pj =1 ; pj < tet.size(1); ++pj)
          one_ring_point_of_p[tet(pi,ti)].insert(tet((pi+pj)%tet.size(1),ti));
      }

  assert(one_ring_point_of_p.begin()->first == 0);

  std::shared_ptr<vector< std::shared_ptr<const function_t<double,int32_t> > > >
      all_funcs(new vector< std::shared_ptr<const function_t<double,int32_t> > >);

  all_funcs->push_back(std::shared_ptr<function_t<double,int32_t> >(
                         new fix_point(orig_node.size(2),0,orig_node(colon(),0))));

  for(map<size_t, set<size_t> >::const_iterator cit = one_ring_point_of_p.begin();
      cit != one_ring_point_of_p.end(); ++cit){
      const set<size_t> & linked_points = cit->second;
      all_funcs->push_back(
            std::shared_ptr<function_t<double,int32_t> >(
              new vertex_update_func(orig_node, cit->first,
                                     linked_points, rotation_gradient)));
    }

//  {
//    std::shared_ptr<function_t<double,int32_t> > smooth_surface_func(
//          build_surface_smooth_func(orig_node, orpap, 1e-3));
//    all_funcs->push_back(smooth_surface_func);
//  }

  std::shared_ptr<const function_t<double,int32_t> >
      func(new_catenated_function<double, int32_t>(all_funcs));
  node = orig_node;
  ptree pt;
  pt.put("package.value","alglib");
  pt.put("alg.value", "lbfgs");
  pt.put("lbfgs-len.value", "7");
  pt.put("iter.value", 500);

  //    pt.put("package.value", "hj");
  //    pt.put("alg.value", "More");
  //    pt.put("GS_LM_mu.value", "0.05");
  //    pt.put("iter.value", 30);
  matrix<double> residual = zeros<double>(func->dim_of_f(),1);
  zjucad::optimize(*func, node, residual, pt);

//  matrix<double> avg_point = zeros<double>(3,1);
//  for(size_t i = 0; i < 1; ++i){
//      for(boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
//          cit = orpap.begin();  cit != orpap.end(); ++cit){
//          const boost::unordered_set<size_t> & connect_points = cit->second;
//          avg_point = zeros<double>(3,1);
//          for(boost::unordered_set<size_t>::const_iterator scit = connect_points.begin();
//              scit != connect_points.end(); ++scit){
//              avg_point += node(colon(),*scit);
//            }
//          avg_point /= connect_points.size();
//          node(colon(),cit->first) = (node(colon(),cit->first) + avg_point)/2.0;
//        }
//    }

  return 0;
}

int deform_tet(jtf::mesh::meshes & tm)
{
  matrix<size_t> outside_face;
  std::shared_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }

  get_outside_face(*fa, outside_face, true);

  matrix<double> point_normal;
  matrix<double> quaternions;
  matrix<matrix<double> > rotation_gradient;
  matrix<double> orig_node = tm.node_;

  map<size_t, set<size_t> > one_ring_point_of_p;
  boost::unordered_map<size_t, boost::unordered_set<size_t> > one_ring_surface_point_of_p;
  for(size_t ti = 0; ti < tm.mesh_.size(2); ++ti)
    for(size_t pi = 0; pi < tm.mesh_.size(1); ++pi){
        for(size_t pj =1 ; pj < tm.mesh_.size(1); ++pj)
          one_ring_point_of_p[tm.mesh_(pi,ti)].insert(tm.mesh_((pi+pj)%tm.mesh_.size(1),ti));
      }

  for(size_t fi = 0; fi < outside_face.size(2); ++fi)
    for(size_t pi = 0; pi < outside_face.size(1); ++pi){
        for(size_t pj =1 ; pj < outside_face.size(1); ++pj)
          one_ring_surface_point_of_p[outside_face(pi,fi)].insert(
                outside_face((pi+pj)%outside_face.size(1),fi));
      }

  for(size_t i = 0; i < 6; ++i){ // based on paper, 5 iterations are used

      jtf::mesh::cal_point_normal(outside_face, tm.node_, point_normal);
      jtf::tetmesh::smooth_one_ring_point_normal(point_normal, one_ring_surface_point_of_p,5);

      cal_quaternions(outside_face, point_normal, quaternions);
      //      propagate_quaternious_linear_solver(tm.mesh_, tm.node_, outside_face,
      //                                          one_ring_point_of_p, quaternions);

      propagate_quaternious(tm.mesh_, tm.node_, outside_face,
                            one_ring_point_of_p, quaternions);
      {
        ostringstream os;
        os << "quat_" << i << ".vtk";
        ofstream ofs(os.str().c_str());
        tet2vtk(ofs, &tm.node_[0], tm.node_.size(2), &tm.mesh_[0], tm.mesh_.size(2));
        matrix<double> x = quaternions(0,colon());
        matrix<double> y = quaternions(1,colon());
        matrix<double> z = quaternions(2,colon());
        matrix<double> w = quaternions(3,colon());

        point_data(ofs, &x[0],x.size(),"x");
        vtk_data(ofs, &y[0],y.size(),"y");
        vtk_data(ofs, &z[0],z.size(),"z");
        vtk_data(ofs, &w[0],w.size(),"w");
      }
      quaternions2rot_matrix(quaternions, rotation_gradient);
      update_vertex(tm.mesh_, orig_node, tm.node_, rotation_gradient,
                    one_ring_surface_point_of_p);

      orig_node = tm.node_;

      ostringstream os;
      os << "output_" << i << ".vtk";
      ofstream ofs(os.str().c_str());
      tet2vtk(ofs, &orig_node[0], orig_node.size(2), &tm.mesh_[0], tm.mesh_.size(2));
    }

  return 0;
}

int sgp_deform(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [error] sgp_deform input_tet output_tet" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  jtf::mesh::meshes output_tm = tm;
  deform_tet(output_tm);

  if(jtf::mesh::tet_mesh_write_to_zjumat(argv[2], &output_tm.node_, &output_tm.mesh_))
    return __LINE__;

  return 0;
}
