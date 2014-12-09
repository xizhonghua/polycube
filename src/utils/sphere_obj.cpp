#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/optimizer/optimizer.h>
#include <zjucad/ptree/ptree.h>
#include <fstream>
#include <iostream>
#include <hjlib/math_func/math_func.h>
#include <hjlib/math_func/func_aux.h>
#include <hjlib/math_func/operation.h>

#include "../tetmesh/tetmesh.h"
#include "../common/util.h"
#include "../vol_param/descriptor/func_terms/arap.h"
#include <hjlib/math/polar.h>

using namespace std;
using namespace zjucad::matrix;
using namespace hj::math_func;
using namespace boost::property_tree;


// dot(p,p)-1
template <typename val_type, typename int_type>
class dis_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  dis_func(const matrix<val_type> &nodes) :nodes_(nodes){}
  virtual ~dis_func(){}
  virtual size_t nx(void)const{
    return nodes_.size();
  }
  virtual size_t nf(void)const{
    return nodes_.size(2);
  }
  virtual int eval(size_t k , const val_type *x, const coo2val_t<val_type, int_type> &cv,
                   func_ctx *ctx = 0) const
  {
    itr_matrix<const val_type*> x0(3, nodes_.size(2) ,x);
    if(k == 0){
        int_type c[1];
        for(size_t pi = 0; pi < nodes_.size(2); ++pi){
            c[0] = pi;
            matrix<val_type> xp = x0(colon(),pi);
            cv[c] += dot(xp,xp) -1.0;
          }
      }
    if(k == 1){
        int_type c[2];
        for(size_t pi = 0; pi < nodes_.size(2); ++pi){
            for(size_t di = 0; di < nodes_.size(1); ++di){
                c[0] = pi;
                c[1] = 3*pi+di;
                cv[c] += 2*x0(di,pi);
              }
          }
      }if(k == 2){
        int_type c[3];
        for(size_t pi = 0; pi < nodes_.size(2); ++pi){
            for(size_t di = 0; di < nodes_.size(1); ++di){
                c[0] = pi;
                c[1] = 3*pi+di;
                c[2] = 3*pi+di;
                cv[c] += 2;
              }
          }
      }
    return 0;
  }
  virtual int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        int_type c[2];
        for(size_t pi = 0; pi < nodes_.size(2); ++pi){
            c[0] = pi;
            for(size_t di = 0; di < nodes_.size(1); ++di){
                c[1] = 3*pi+di;
                l2g.add(cs, c);
              }
          }
      }
    if(k == 2){
        int_type c[3];
        for(size_t pi = 0; pi < nodes_.size(2); ++pi){
            c[0] = pi;
            for(size_t di = 0; di < nodes_.size(1); ++di){
                c[1] = 3*pi+di;
                c[2] = 3*pi+di;
                l2g.add(cs, c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1){
        return nodes_.size();
      }
    if(k == 2){
        return nodes_.size();
      }
  }
private:
  const matrix<val_type> & nodes_;
};


template <typename val_type, typename int_type>
class dis_func2 : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  dis_func2(const matrix<val_type> &nodes,const set<size_t> & points,
            const val_type weight)
    :nodes_(nodes),points_(points),weight_(weight){}
  virtual ~dis_func2(){}
  virtual size_t nx(void)const{
    return nodes_.size();
  }
  virtual size_t nf(void)const{
    return points_.size();
  }
  virtual int eval(size_t k , const val_type *x, const coo2val_t<val_type, int_type> &cv,
                   func_ctx *ctx = 0) const
  {
    itr_matrix<const val_type*> x0(nodes_.size(1), nodes_.size(2) ,x);
    if(k == 0){
        int_type c[1];
        int_type idx = 0;
        for(const auto &pi : points_){
            c[0] = idx++;
            matrix<val_type> xp = x0(colon(),pi);
            cv[c] += weight_*(dot(xp,xp) -1.0);
          }
      }
    if(k == 1){
        int_type c[2];
        int_type idx = 0;
        for(const auto & pi : points_){
            for(size_t di = 0; di < nodes_.size(1); ++di){
                c[0] = idx;
                c[1] = nodes_.size(1)*pi+di;
                cv[c] += weight_*2*x0(di,pi);
              }
            ++idx;
          }
      }if(k == 2){
        int_type c[3];
        int_type idx = 0;
        for(const auto & pi : points_){
            for(size_t di = 0; di < nodes_.size(1); ++di){
                c[0] = idx;
                c[1] = nodes_.size(1)*pi+di;
                c[2] = nodes_.size(1)*pi+di;
                cv[c] += 2*weight_;
              }
            ++idx;
          }
      }
    return 0;
  }
  virtual int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        int_type c[2];
        int_type idx = 0;
        for(const auto & pi : points_){
            for(size_t di = 0; di < nodes_.size(1); ++di){
                c[0] = idx;
                c[1] = nodes_.size(1)*pi+di;
                l2g.add(cs, c);
              }
            ++idx;
          }
      }if(k == 2){
        int_type c[3];
        int_type idx = 0;
        for(const auto & pi : points_){
            for(size_t di = 0; di < nodes_.size(1); ++di){
                c[0] = idx;
                c[1] = nodes_.size(1)*pi+di;
                c[2] = nodes_.size(1)*pi+di;
                l2g.add(cs, c);
              }
            ++idx;
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1){
        return nodes_.size();
      }
    if(k == 2){
        return nodes_.size();
      }
  }
private:
  const matrix<val_type> & nodes_;
  const set<size_t> points_;
  const val_type weight_;
};

// \sum_k p_k - p_i
template <typename val_type, typename int_type>
class smooth_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  smooth_func(const matrix<val_type> &node,
              const jtf::mesh::one_ring_point_at_point &orpap,
              const double weight= 1)
    :node_(node), orpap_(orpap), weight_(weight){}
  virtual ~smooth_func(){}
  virtual size_t nx(void)const{
    return node_.size();
  }
  virtual size_t nf(void) const
  {
    return node_.size();
  }
  virtual int eval(size_t k , const val_type *x, const coo2val_t<val_type, int_type> &cv,
                   func_ctx *ctx = 0) const
  {
    itr_matrix<const val_type*> x0(node_.size(1), node_.size(2), x);
    if(k == 0){
        int_type c[0];
        matrix<double> sum_nodes(node_.size(1),1);
        for(const auto & one_point : orpap_.p2p_){
            const size_t & this_point = one_point.first;
            const vector<size_t> & connected_points = one_point.second;
            sum_nodes *= 0;
            for(size_t vi = 0; vi < connected_points.size(); ++vi){
                sum_nodes += x0(colon(), connected_points[vi]);
              }
            sum_nodes /= connected_points.size();
            sum_nodes -= x0(colon(), this_point);
            sum_nodes *= 1;
            for(size_t di = 0; di < node_.size(1); ++di){
                c[0] = node_.size(1) * this_point + di;
                cv[c] += weight_*sum_nodes[di];
              }
          }
      }
    if(k == 1){
        int_type c[1];
        for(const auto & one_point : orpap_.p2p_){
            const size_t & this_point = one_point.first;
            const vector<size_t> &connected_points = one_point.second;
            for(size_t di = 0; di < node_.size(1); ++di){
                c[0] = node_.size(1)*this_point+di;
                for(size_t i = 0; i < connected_points.size(); ++i){
                    c[1] = node_.size(1)*connected_points[i]+di;
                    cv[c] += 1.0/connected_points.size();
                  }
                c[1] = node_.size(1)*this_point + di;
                cv[c] -= 1.0*weight_;
              }
          }
      }
    return 0;
  }
  virtual int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        int_type c[1];
        for(const auto & one_point : orpap_.p2p_){
            const size_t & this_point = one_point.first;
            const vector<size_t> &connected_points = one_point.second;
            for(size_t di = 0; di < node_.size(1); ++di){
                c[0] = node_.size(1)*this_point+di;
                for(size_t i = 0; i < connected_points.size(); ++i){
                    c[1] = node_.size(1)*connected_points[i]+di;
                    l2g.add(cs,c);
                  }
                c[1] = node_.size(1)*this_point + di;
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }

  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1){
        size_t n = 0;
        for(const auto & one_point : orpap_.p2p_){
            n += one_point.second.size();
            ++n;
          }
        return n*node_.size(1);
      }
  }
private:
  const jtf::mesh::one_ring_point_at_point &orpap_;
  const zjucad::matrix::matrix<val_type> & node_;
  const double weight_;
};

// x
template <typename val_type, typename int_type>
class point_fix_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  point_fix_func(const matrix<val_type> &node, const double weight)
    :node_(node), weight_(weight){}
  virtual ~point_fix_func(){}
  virtual size_t nx(void)const{
    return node_.size();
  }
  virtual size_t nf(void) const
  {
    return node_.size();
  }
  virtual int eval(size_t k , const val_type *x, const coo2val_t<val_type, int_type> &cv,
                   func_ctx *ctx = 0) const
  {
    itr_matrix<const val_type*> x0(3, node_.size(2), x);
    if(k == 0){
        int_type c[1]= {0};
        for(size_t pi = 0; pi < node_.size(2); ++pi){
            for(size_t di = 0; di < node_.size(1); ++di){
                c[0] =  3*pi+di;
                cv[c] += weight_*(x0(di,pi) - node_(di,pi));
              }
          }
      }
    if(k == 1){
        int_type c[2];
        for(size_t pi = 0; pi < node_.size(2); ++pi){
            for(size_t di = 0; di < node_.size(1); ++di){
                c[0] = 3*pi+di;
                c[1] = 3*pi+di;
                cv[c] += weight_;
              }
          }
      }
    return 0;
  }
  virtual int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        int_type c[2];
        for(size_t pi = 0; pi < node_.size(2); ++pi){
            for(size_t di = 0; di < node_.size(1); ++di){
                c[0] = 3*pi+di;
                c[1] = 3*pi+di;
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }

  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1){
        return node_.size();
      }
  }
private:
  const zjucad::matrix::matrix<val_type> & node_;
  const double weight_;
};

static void optimize_model(jtf::mesh::meshes & trim, const double weight)
{
  shared_ptr<jtf::mesh::one_ring_point_at_point> orpap(jtf::mesh::one_ring_point_at_point::create(trim.mesh_));
  if(!orpap.get()){
      throw std::invalid_argument("invalid trimesh.");
    }

  typedef math_func_t<double,int32_t> math_func_T;
  typedef shared_ptr<math_func_T> math_func_ptr;

  shared_ptr<vector<math_func_ptr> >  funcs(new vector<math_func_ptr>);

  math_func_ptr dis_func_ptr(new dis_func<double,int32_t>(trim.node_));
  funcs->push_back(math_func_ptr(new sumsqr<double,int32_t>(dis_func_ptr)));

  math_func_ptr smooth_func_ptr(new smooth_func<double,int32_t>(trim.node_, *orpap));
  funcs->push_back(math_func_ptr(new sumsqr<double,int32_t>(smooth_func_ptr)));

  math_func_ptr point_fix_func_ptr(new point_fix_func<double,int32_t>(trim.node_, weight));
  funcs->push_back(math_func_ptr(new sumsqr<double,int32_t>(point_fix_func_ptr)));

  //  math_func_ptr test_func_ptr(new test_func<double,int32_t>(trim.node_));
  //  funcs->push_back(math_func_ptr(new sumsqr<double,int32_t>(test_func_ptr)));

  math_func_ptr all_funcs_cat( new fcat<double, int32_t, vector<math_func_ptr> >(funcs));
  math_func_ptr all_funcs( new hj::math_func::sum<double,int32_t>(all_funcs_cat));

  boost::property_tree::ptree pt;
  pt.put("package.value","jtf");
  pt.put("alg.value","SQP");
  pt.put("iter.value",10);
  jtf::optimize(*all_funcs, trim.node_, pt, nullptr, nullptr, nullptr);
}

static void optimize_model_with_boundary(jtf::mesh::meshes & trim, const double weight)
{
  shared_ptr<jtf::mesh::one_ring_point_at_point> orpap(jtf::mesh::one_ring_point_at_point::create(trim.mesh_));
  if(!orpap.get()){
      throw std::invalid_argument("invalid trimesh.");
    }

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(trim.mesh_));

  matrix<size_t> edges;
  jtf::mesh::get_boundary_edge(*ea, edges);
  set<size_t> boundary_points(edges.begin(), edges.end());

  typedef math_func_t<double,int32_t> math_func_T;
  typedef shared_ptr<math_func_T> math_func_ptr;

  shared_ptr<vector<math_func_ptr> >  funcs(new vector<math_func_ptr>);

  math_func_ptr dis_func_ptr(new dis_func2<double,int32_t>(trim.node_, boundary_points, sqrt(weight)));
  funcs->push_back(math_func_ptr(new sumsqr<double,int32_t>(dis_func_ptr)));

  math_func_ptr smooth_func_ptr(new smooth_func<double,int32_t>(trim.node_, *orpap));
  funcs->push_back(math_func_ptr(new sumsqr<double,int32_t>(smooth_func_ptr)));

  math_func_ptr all_funcs_cat( new fcat<double, int32_t, vector<math_func_ptr> >(funcs));
  math_func_ptr all_funcs( new hj::math_func::sum<double,int32_t>(all_funcs_cat));

  boost::property_tree::ptree pt;
  pt.put("package.value","jtf");
  pt.put("alg.value","SQP");
  pt.put("iter.value",100);
  jtf::optimize(*all_funcs, trim.node_, pt, nullptr, nullptr, nullptr);
}

int sphere_obj(int argc, char * argv[])
{
  if(argc != 4){
      std::cerr << "# [usage] sphere_obj input_obj output_obj smooth_weight" << std::endl;
      return __LINE__;
    }

  jtf::mesh::meshes trim;
  if(jtf::mesh::load_obj(argv[1], trim.mesh_, trim.node_)){
      cerr << "# [error] can not open input obj." << endl;
      return __LINE__;
    }

  optimize_model(trim, atof(argv[3]));

  if(jtf::mesh::save_obj(argv[2], trim.mesh_, trim.node_)){
      cerr << "# [error] can not save output obj." << endl;
      return __LINE__;
    }

  return 0;
}

int sphere_obj_with_boundary(int argc, char * argv[])
{
  if(argc != 4){
      std::cerr << "# [usage] sphere_obj input_obj output_obj smooth_weight" << std::endl;
      return __LINE__;
    }

  jtf::mesh::meshes trim;
  if(jtf::mesh::load_obj(argv[1], trim.mesh_, trim.node_)){
      cerr << "# [error] can not open input obj." << endl;
      return __LINE__;
    }

  matrix<double> new_node = trim.node_(colon(0,1),colon());
  trim.node_ = new_node;


  optimize_model_with_boundary(trim, atof(argv[3]));
  new_node = trim.node_;

  trim.node_.resize(3, trim.node_.size(2));
  trim.node_(colon(0,1),colon()) = new_node;

  if(jtf::mesh::save_obj(argv[2], trim.mesh_, trim.node_)){
      cerr << "# [error] can not save output obj." << endl;
      return __LINE__;
    }

  return 0;
}

template <typename val_type, typename int_type>
class fix_position_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  fix_position_func(const size_t points_num, const size_t idx,
                    const zjucad::matrix::matrix<double> & target,const double w)
    :num_(points_num), idx_(idx), target_(target), w_(w) {}
  virtual ~fix_position_func(){}
  virtual size_t nx() const{
    return 3*num_;
  }
  virtual size_t nf() const{
    return 3;
  }
  virtual int eval(size_t k, const val_type *x,
                   const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type*> x0(3, num_, &x[0]);
        matrix<val_type> diff = x0(colon(), idx_) - target_;
        for(int_type fi = 0; fi < nf(); ++fi){
            int_type c[1] = {fi};
            cv[c] += w_*diff[fi];
          }
      }
    if(k==1){
        for(int_type fi = 0; fi < nf(); ++fi){
            int_type c[2] = {fi, 3*idx_+fi};
            cv[c] += w_;
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==1){
        for(int_type fi = 0; fi < nf(); ++fi){
            int_type c[2] = {fi, 3*idx_+fi};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k == 0) return -1;
    if(k == 1) return 3;
  }
private:
  const size_t num_;
  const size_t idx_;
  const zjucad::matrix::matrix<double> target_;
  const double w_;
};

shared_ptr<hj::math_func::math_func_t<double, int32_t> >
build_fix_position_func(const matrix<double>  &nodes,
                        const zjucad::matrix::matrix<size_t> &v2s,
                        const zjucad::matrix::matrix<double> &target,
                        double w)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  for(size_t ti = 0; ti < v2s.size(); ++ti){
      if(v2s[ti] == -1) continue;
      all_func->push_back(math_func_ptr(
                            new fix_position_func<double,int32_t>(nodes.size(2), ti, target(colon(),v2s[ti]), sqrt(w))));
    }

  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

template <typename val_type, typename int_type>
class math_func_arap_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  math_func_arap_func(const zjucad::matrix::matrix<size_t> & one_tet,
                      const zjucad::matrix::matrix<double> & node, const double w)
    :one_tet_(one_tet), node_(node), w_(w), grad_op_(4,3) {
    matrix<double> tet_node = node(colon(), one_tet);
    if(calc_tet_def_grad_op(&tet_node[0], &grad_op_[0])){
        std::cerr << "# degenerated tet"  << std::endl;
      }
  }
  virtual ~math_func_arap_func(){}
  virtual size_t nx() const{
    return 3*node_.size(2);
  }
  virtual size_t nf() const{
    return 3*3;
  }
  virtual int eval(size_t k, const val_type *x,
                   const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type*> x0(3, node_.size(2), &x[0]);
        matrix<val_type> f = x0(colon(), one_tet_) * grad_op_;
        matrix<val_type> R = f;
        hj::polar3d p;
        p(R,2);
        f -= R;

        for(int_type fi = 0; fi < nf(); ++fi){
            int_type c[1] = {fi};
            cv[c] += w_*f[fi];
          }
      }
    if(k==1){
        for(int_type fi = 0; fi < nf(); ++fi){
            int_type r = fi%3;
            int_type ci = fi/3;
            for(int_type pi = 0; pi < 4; ++pi){
                int_type c[2] = {fi, one_tet_[pi]*3+r};
                cv[c] += grad_op_(pi,ci) *w_;
              }
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==1){
        for(int_type fi = 0; fi < nf(); ++fi){
            int_type r = fi%3;
            int_type ci = fi/3;
            for(int_type pi = 0; pi < 4; ++pi){
                int_type c[2] = {fi, one_tet_[pi]*3+r};
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k == 0) return -1;
    if(k == 1) return 4*9;
  }
private:
  const zjucad::matrix::matrix<double> & node_;
  const zjucad::matrix::matrix<size_t> one_tet_;
  const double w_;
  zjucad::matrix::matrix<double> grad_op_;
};

shared_ptr<hj::math_func::math_func_t<double, int32_t> >
build_arap_func(const jtf::tet_mesh &tm)
{

  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  double total_vol = 0;
  for(size_t ti = 0; ti < tm.vol_.size(); ++ti)
    total_vol += tm.vol_[ti];

  for(size_t ti = 0; ti < tm.tetmesh_.mesh_.size(2); ++ti){
      const double w = fabs(tm.vol_[ti]/ total_vol);
      all_func->push_back(math_func_ptr(
                            new math_func_arap_func<double,int32_t>(
                              tm.tetmesh_.mesh_(colon(),ti), tm.tetmesh_.node_, sqrt(w))));
    }

  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}


int sphere_tet(int argc, char *argv[])
{
  if(argc != 4){
      std::cerr << "# [usage] sphere_tet input_tet output_tet smooth_weight" << std::endl;
      return __LINE__;
    }

  jtf::tet_mesh tm(argv[1]);

  matrix<size_t> mapping_point;
  matrix<size_t> face = tm.outside_face_;
  matrix<double> node = tm.tetmesh_.node_;

  remove_extra_node(face, node, &mapping_point);

  jtf::mesh::meshes trim;
  trim.mesh_ = face;
  trim.node_ = node;

  optimize_model(trim, atof(argv[3]));

  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);
  func->push_back(math_func_ptr(
                    new hj::math_func::sumsqr<double,int32_t>(
                      build_fix_position_func(tm.tetmesh_.node_, mapping_point, trim.node_, 1))
                    ));
  func->push_back(math_func_ptr(
                    new hj::math_func::sumsqr<double,int32_t>(
                      build_arap_func(tm))));

  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  boost::property_tree::ptree pt;
  pt.put("package.value", "jtf");
  pt.put("alg.value", "SQP");
  pt.put("iter.value",10);

  jtf::optimize(*obj, tm.tetmesh_.node_, pt, nullptr, nullptr, nullptr);

  jtf::mesh::tet_mesh_write_to_zjumat(argv[2], &tm.tetmesh_.node_, &tm.tetmesh_.mesh_);

  return 0;
}
