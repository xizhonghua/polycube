#ifndef MESH_FRAME_FUNC_H
#define MESH_FRAME_FUNC_H

#include <memory>
#include <vector>
#include <hjlib/math_func/math_func.h>
#include <hjlib/math_func/operation.h>
#include <jtflib/mesh/mesh.h>
#include <zjucad/matrix/itr_matrix.h>
#include "../spherical_harmonics/rot_cubic_f_SH.h"
#include "../common/zyz.h"
#include "../common/util.h"
#include "../common/vtk.h"

template <typename val_type, typename int_type>
class inner_smooth_func_sh : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  inner_smooth_func_sh(const size_t tet_num, const size_t face_num,
                       const size_t i, const size_t j,  const double w,
                       const double k0 = 1, const double k1 = 1)
    :i_(i), j_(j), tet_num_(tet_num), face_num_(face_num), w_(w), k0_(k0), k1_(k1){
  }
  virtual ~inner_smooth_func_sh(){}
  virtual size_t nx() const{
    return 9 * tet_num_ + 2*face_num_;
  }
  virtual size_t nf() const{
    return 9;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        zjucad::matrix::itr_matrix<const val_type *> x0(9,tet_num_,x);
        for(int_type i = 0; i < 9; ++i){
            int_type c[1] = {i};
            cv[c] += w_*(k0_*x0(i,i_)-k1_*x0(i,j_));
          }
      }
    if(k == 1){
        for(int_type i = 0; i < 9; ++i){
            int_type c1[2] = {i,9*i_+i};
            cv[c1] += w_*k0_;
            int_type c2[2] = {i,9*j_+i};
            cv[c2] += -w_*k1_;
          }
      }
    return 0;
  }

  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        for(int_type i = 0; i < 9; ++i){
            int_type c1[2] = {i,9*i_+i};
            l2g.add(cs,c1);
            int_type c2[2] = {i,9*j_+i};
            l2g.add(cs,c2);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 2*9;
  }
private:
  const int_type i_;
  const int_type j_;
  const int_type tet_num_;
  const int_type face_num_;
  const double w_;
  const double k0_;
  const double k1_;
};


template <typename val_type, typename int_type>
class inner_smooth_func_sh_L1 : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  inner_smooth_func_sh_L1(const size_t tet_num,  const size_t i, const size_t j,
                          const size_t di , const double w)
    :i_(i), j_(j), tet_num_(tet_num), di_(di), w_(w){
  }
  virtual ~inner_smooth_func_sh_L1(){}
  virtual size_t nx() const{
    return 9 * tet_num_;
  }
  virtual size_t nf() const{
    return 1;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        zjucad::matrix::itr_matrix<const val_type *> x0(9,tet_num_,x);
        val_type a = 0.5;
        val_type diff = sqrt( (x0(di_,i_)-x0(di_,j_),2)+a);
        int_type c[1] = {0};
        cv[c] += w_*diff;
      }
    if(k == 1){
        zjucad::matrix::itr_matrix<const val_type *> x0(9,tet_num_,x);
        val_type diff = x0(di_,i_)-x0(di_,j_);
        val_type a = 0.5;
        int_type c[2] = {0,9*i_+di_};
        cv[c] += w_*diff/sqrt(diff*diff+a);

        int_type c2[2] = {0,9*j_+di_};
        cv[c2] += -w_*diff/sqrt(diff*diff+a);
      }
    if(k == 2){
        zjucad::matrix::itr_matrix<const val_type *> x0(9,tet_num_,x);
        val_type a = 0.5;
        val_type diff0 = x0(di_,i_)-x0(di_,j_);
        val_type diff1 = sqrt(diff0*diff0+a);
        val_type diff2 = a/pow(diff1,3);

        int_type c1[3] = {0,9*i_+di_, 9*i_+di_};
        cv[c1] += w_*diff2;

        int_type c2[3] = {0,9*j_+di_, 9*j_+di_};
        cv[c2] += w_*diff2;

        int_type c3[3] = {0,9*i_+di_, 9*j_+di_};
        cv[c3] += -w_*diff2;


        int_type c4[3] = {0,9*j_+di_, 9*i_+di_};
        cv[c4] += -w_*diff2;
      }
    return 0;
  }

  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        int_type c[2] = {0,9*i_+di_};
        l2g.add(cs,c);

        int_type c2[2] = {0,9*j_+di_};
        l2g.add(cs,c2);
      }
    if(k == 2){
        int_type c1[3] = {0,9*i_+di_, 9*i_+di_};
        l2g.add(cs,c1);

        int_type c2[3] = {0,9*j_+di_, 9*j_+di_};
        l2g.add(cs,c2);

        int_type c3[3] = {0,9*i_+di_, 9*j_+di_};
        l2g.add(cs,c3);

        int_type c4[3] = {0,9*j_+di_, 9*i_+di_};
        l2g.add(cs,c4);
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 2;
    if(k == 2) return 4;
  }
private:
  const int_type i_;
  const int_type j_;
  const int_type tet_num_;
  const int_type di_;
  const double w_;
};

template <typename val_type, typename int_type>
class inner_smooth_func_sh_size : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  inner_smooth_func_sh_size(const size_t tet_num, const size_t i, const size_t j,
                            const zjucad::matrix::matrix<double> & weight_edge)
    :i_(i), j_(j), tet_num_(tet_num), weight_edge_(weight_edge){}
  virtual ~inner_smooth_func_sh_size(){}
  virtual size_t nx() const{
    return 9 * tet_num_ ;
  }
  virtual size_t nf() const{
    return 9*3;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        zjucad::matrix::itr_matrix<const val_type *> x0(9,tet_num_,x);
        zjucad::matrix::matrix<val_type> grad_f = (x0(zjucad::matrix::colon(), i_) - x0(zjucad::matrix::colon(), j_))*zjucad::matrix::trans(weight_edge_);
        for(int_type i = 0; i < 9*3; ++i){
            int_type c[1] = {i};
            cv[c] += grad_f[0];
          }
      }
    if(k == 1){
        zjucad::matrix::itr_matrix<const val_type *> x0(9,tet_num_,x);
        for(int_type i = 0; i < 9; ++i){
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i*3+j, 9*i_+i};
                cv[c] += weight_edge_[j];
              }
          }
      }
    return 0;
  }

  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        for(int_type i = 0; i < 9; ++i){
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i*3+j, 9*i_+i};
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 27*2;
  }
private:
  const int_type i_;
  const int_type j_;
  const int_type tet_num_;
  const zjucad::matrix::matrix<double> weight_edge_;
};

template <typename val_type, typename int_type>
class inner_smooth_func_sh_laplace : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  inner_smooth_func_sh_laplace(const size_t tet_num, const size_t tet_idx,
                               const std::vector<size_t> &adj_tets, const double w )
    : tet_num_(tet_num), tet_idx_(tet_idx), w_(w), adj_tets_(adj_tets){}
  virtual ~inner_smooth_func_sh_laplace(){}
  virtual size_t nx() const{
    return 9 * tet_num_ ;
  }
  virtual size_t nf() const{
    return 9;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        zjucad::matrix::itr_matrix<const val_type *> x0(9,tet_num_,x);
        zjucad::matrix::matrix<val_type> diff = x0(zjucad::matrix::colon(), tet_idx_);
        for(size_t j = 0; j < adj_tets_.size(); ++j){
            diff -= x0(zjucad::matrix::colon(), adj_tets_[j])/adj_tets_.size();
          }
        for(int_type i = 0; i < 9; ++i){
            int_type c[1] = {i};
            cv[c] += w_*diff[i];
          }
      }
    if(k == 1){
        for(int_type i = 0; i < 9; ++i){
            int_type c1[2] = {i,9*tet_idx_+i};
            cv[c1] += w_;
            for(int_type j = 0; j < adj_tets_.size(); ++j){
                int_type c2[2] = {i, int_type(9*adj_tets_[j]+i)};
                cv[c2] += -w_/adj_tets_.size();
              }
          }
      }
    return 0;
  }

  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        for(int_type i = 0; i < 9; ++i){
            int_type c1[2] = {i,9*tet_idx_+i};
            l2g.add(cs,c1);
            for(int_type j = 0; j < adj_tets_.size(); ++j){
                int_type c2[2] = {i, int_type(9*adj_tets_[j]+i)};
                l2g.add(cs,c2);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return (adj_tets_.size()+1)*9;
  }
private:
  const int_type tet_num_;
  const int_type tet_idx_;
  const double w_;
  const std::vector<size_t> adj_tets_;
};

template <typename val_type, typename int_type>
class inner_smooth_func_zyz : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  inner_smooth_func_zyz(const size_t tet_num, const size_t face_num,
                        const size_t i, const size_t j,  const double w,
                        const double k0 = 1, const double k1 = 1)
    :i_(i), j_(j), tet_num_(tet_num), face_num_(face_num), w_(w), k0_(k0), k1_(k1){}
  virtual ~inner_smooth_func_zyz(){}
  virtual size_t nx() const{
    return 3 * tet_num_ + 2*face_num_;
  }
  virtual size_t nf() const{
    return 9;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        zjucad::matrix::itr_matrix<const val_type *> x0(3,tet_num_,x);
        zjucad::matrix::matrix<val_type> sh_i(9,1), sh_j(9,1);
        calc_rot_cubic_f_sh_(&sh_i[0], &x0(0,i_));
        calc_rot_cubic_f_sh_(&sh_j[0], &x0(0,j_));
        for(int_type i = 0; i < 9; ++i){
            int_type c[1] = {i};
            cv[c] += w_*(k0_*sh_i[i]-k1_*sh_j[i]);
          }
      }
    if(k == 1){
        zjucad::matrix::itr_matrix<const val_type *> x0(3,tet_num_,x);
        zjucad::matrix::matrix<val_type> jac_sh_i(9,3), jac_sh_j(9,3);
        calc_jac_rot_cubic_f_sh_(&jac_sh_i[0], &x0(0,i_));
        calc_jac_rot_cubic_f_sh_(&jac_sh_j[0], &x0(0,j_));

        for(int_type i = 0; i < 9; ++i){
            for(int_type di = 0; di< 3; ++di){
                int_type c1[2] = {i, 3*i_+di};
                cv[c1] += w_*k0_*jac_sh_i(i,di);
                int_type c2[2] = {i, 3*j_+di};
                cv[c2] += -w_*k1_*jac_sh_j(i,di);
              }
          }
      }
    return 0;
  }

  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        for(int_type i = 0; i < 9; ++i){
            for(int_type di = 0; di< 3; ++di){
                int_type c1[2] = {i, 3*i_+di};
                l2g.add(cs,c1);
                int_type c2[2] = {i, 3*j_+di};
                l2g.add(cs,c2);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 3*2*9;
  }
private:
  const int_type i_;
  const int_type j_;
  const int_type tet_num_;
  const int_type face_num_;
  const double w_;
  const double k0_;
  const double k1_;
};

template <typename val_type>
std::shared_ptr<const hj::math_func::math_func_t<val_type, int32_t> >
build_inner_smooth_func(const zjucad::matrix::matrix<size_t> &tet_mesh,
                        const zjucad::matrix::matrix<val_type> &node,
                        const jtf::mesh::face2tet_adjacent &fa,
                        const jtf::mesh::one_ring_tet_at_edge &ortae,
                        const zjucad::matrix::matrix<val_type> &vol,
                        const val_type w,
                        const std::string strategy = "face",
                        const std::string opt_type = "init",
                        const zjucad::matrix::matrix<val_type> * size_field = 0)
{
  using namespace std;
  using namespace zjucad::matrix;

  typedef hj::math_func::math_func_t<val_type,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_v = std::accumulate(zjucad::matrix::fabs(vol).begin(),
                                         zjucad::matrix::fabs(vol).end(), 0.0);

  matrix<double> bb(3,2);
  calc_bounding_box(node, &bb[0]);

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  if(strategy == "face"){
      vector<double> cell_weight(tet_mesh.size(2),1);
      matrix<double> center(3,1);
      const double lambda = 1;
      for(size_t ti = 0; ti < tet_mesh.size(2); ++ti){
          center *= 0;
          for(size_t pi = 0; pi < tet_mesh.size(1); ++pi){
              center += node(colon(), tet_mesh(pi,ti));
            }
          center /= tet_mesh.size(1);
          const double radius = 1e-5*sqrt(center[0]*center[0]+center[2]*center[2]);
          if(center[1] < 1 && radius < 1)
            cell_weight[ti] = lambda*radius+(1-lambda)*1;
        }
      {
        ofstream ofs("cell_size.vtk");
        tet2vtk(ofs, &node[0], node.size(2), &tet_mesh[0], tet_mesh.size(2));
        cell_data(ofs, &cell_weight[0], cell_weight.size(), "cell_size");
      }
      for(size_t fi = 0; fi < fa.face2tet_.size(); ++fi){
          const pair<size_t,size_t> & tet_pair = fa.face2tet_[fi];
          if(fa.is_outside_face(tet_pair)) continue;
          double size = 1;
          if(size_field){
              const vector<size_t> & one_face = fa.faces_[fi];
              size = (*size_field)[one_face[0]]+(*size_field)[one_face[1]]+(*size_field)[one_face[2]];
              size /= 3;
              size = 1.0/size;
            }
          {
            size = (cell_weight[tet_pair.first] + cell_weight[tet_pair.second])*0.5;
          }
          if(opt_type == "init")
            all_func->push_back(math_func_ptr(
                                  new inner_smooth_func_sh<val_type,int32_t>(
                                    tet_mesh.size(2), 0,
                                    tet_pair.first, tet_pair.second,
                                    sqrt(size*w*(fabs(vol[tet_pair.first])+fabs(vol[tet_pair.second]))/(4.0*total_v)))));
          else if(opt_type == "zyz")
            all_func->push_back(math_func_ptr(
                                  new inner_smooth_func_zyz<val_type,int32_t>(
                                    tet_mesh.size(2), 0,
                                    tet_pair.first, tet_pair.second,
                                    sqrt(size*w*(fabs(vol[tet_pair.first])+fabs(vol[tet_pair.second]))/(4.0*total_v)))));
          else throw std::invalid_argument("# [error] can not recognize opt type");
        }
    }else if(strategy == "gradient"){
      for(size_t fi = 0; fi < fa.face2tet_.size(); ++fi){
          const pair<size_t,size_t> &  tet_pair = fa.face2tet_[fi];
          if(fa.is_outside_face((tet_pair))) continue;
          const size_t other_v_0 = std::accumulate(tet_mesh(colon(),tet_pair.first).begin(),
                                                   tet_mesh(colon(),tet_pair.first).end(),0)-
              std::accumulate(fa.faces_[fi].begin(),fa.faces_[fi].end(),0);
          const size_t other_v_1 = std::accumulate(tet_mesh(colon(),tet_pair.second).begin(),
                                                   tet_mesh(colon(),tet_pair.second).end(),0)-
              std::accumulate(fa.faces_[fi].begin(),fa.faces_[fi].end(),0);

          const double dis = norm(node(colon(),other_v_0)/4 - node(colon(), other_v_1)/4);
          const double weight = sqrt(w*(fabs(vol[tet_pair.first]) + fabs(vol[tet_pair.second]) )/4)
              * pow(total_v,-1.0/6.0) /dis  ;

          if(opt_type == "init")
            all_func->push_back(math_func_ptr(
                                  new inner_smooth_func_sh<val_type,int32_t>(
                                    tet_mesh.size(2), 0, tet_pair.first,
                                    tet_pair.second, weight)));
          else if(opt_type == "zyz")
            all_func->push_back(math_func_ptr(
                                  new inner_smooth_func_zyz<val_type,int32_t>(
                                    tet_mesh.size(2), 0, tet_pair.first,
                                    tet_pair.second, weight)));
          else throw std::invalid_argument("# [error] can not recognize opt type");
        }
    }else if(strategy == "liuyang"){
      for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator it
          = ortae.e2t_.begin(); it != ortae.e2t_.end(); ++it){
          const vector<size_t> & around_tets = it->second;
          if(!ortae.is_inner_edge(around_tets)) continue;
          for(size_t i = 0; i < around_tets.size()-1; ++i){
              for(size_t j = i + 1; j < around_tets.size()-1; ++j){
                  const double weight = sqrt(w*(fabs(vol[around_tets[i]])+fabs(vol[around_tets[j]]))/(6*total_v));
                  if(opt_type == "init")
                    all_func->push_back(math_func_ptr(
                                          new inner_smooth_func_sh<val_type,int32_t>(
                                            tet_mesh.size(2), 0, around_tets[i],
                                            around_tets[j], weight)));
                  else if(opt_type == "zyz")
                    all_func->push_back(math_func_ptr(
                                          new inner_smooth_func_zyz<val_type,int32_t>(
                                            tet_mesh.size(2), 0, around_tets[i],
                                            around_tets[j], weight)));
                  else throw std::invalid_argument("# [error] can not recognize opt type");
                }
            }
        }
    }else if(strategy == "laplace"){
      matrixst tet_ptr;
      matrixst tet_idx;
      tet2tet_adj_matrix(tet_mesh, tet_ptr, tet_idx, &fa);

      vector<size_t> adj_tets;
      for(size_t ti = 0; ti < tet_ptr.size()-1; ++ti){
          adj_tets.clear();
          for(size_t off = tet_ptr[ti]; off < tet_ptr[ti+1]; ++off){
              adj_tets.push_back(tet_idx[off]);
            }
          if(opt_type == "init")
            all_func->push_back(math_func_ptr(
                                  new inner_smooth_func_sh_laplace<val_type,int32_t>(
                                    tet_mesh.size(2), ti, adj_tets,
                                    sqrt(fabs(vol[ti]/total_v)))));
          else if(opt_type == "zyz"){
              throw std::invalid_argument("not finished");
            }
        }
    }else{
      throw std::invalid_argument("wrong inner smooth strategy,[face/gradient/liuyang/laplace]");
    }
  cerr << "# [info] inner smooth func num " << all_func->size() << endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<val_type,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

template <typename val_type, typename int_type>
class surface_normal_align2_sh: public hj::math_func::math_func_t<val_type, int32_t>
{
public:
  surface_normal_align2_sh(const size_t tet_num, const size_t tet_idx,
                           const zjucad::matrix::matrix<val_type> &rot,
                           const size_t basis_idx, const double v,
                           const val_type w, const val_type ks = 1)
    :tet_num_(tet_num), tet_idx_(tet_idx), w_(w), v_(v), ks_(ks){
    basis_ = rot(basis_idx,zjucad::matrix::colon());
  }
  virtual ~surface_normal_align2_sh(){}
  virtual size_t nx() const{
    return 9*tet_num_;
  }
  virtual size_t nf() const{
    return 1;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==0){
        zjucad::matrix::itr_matrix<const val_type*> x0(9,tet_num_,&x[0]);
        int_type c[1] = {0};
        cv[c] += w_*(ks_*zjucad::matrix::dot(basis_, x0(zjucad::matrix::colon(),tet_idx_))-v_);
      }
    if(k==1){
        for(int_type i = 0; i < 9; ++i){
            int_type c[2] = {0,9*tet_idx_+i};
            cv[c] += w_*ks_*basis_[i];
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==1){
        for(int_type i = 0; i < 9; ++i){
            int_type c[2] = {0,9*tet_idx_+i};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k==0) return -1;
    if(k==1) return 9;
  }
private:
  const int_type tet_num_;
  const int_type tet_idx_;
  zjucad::matrix::matrix<val_type> basis_;
  const val_type w_;
  const val_type v_;
  const val_type ks_;
};

template <typename val_type, typename int_type>
class surface_normal_align2_zyz: public hj::math_func::math_func_t<val_type, int32_t>
{
public:
  surface_normal_align2_zyz(const size_t tet_num, const size_t tet_idx,
                            const zjucad::matrix::matrix<val_type> &rot,
                            const size_t basis_idx, const double v,
                            const val_type w, const val_type ks = 1)
    :tet_num_(tet_num), tet_idx_(tet_idx), w_(w), v_(v), ks_(ks){
    basis_ = rot(basis_idx, zjucad::matrix::colon());
  }
  virtual ~surface_normal_align2_zyz(){}
  virtual size_t nx() const{
    return 3*tet_num_;
  }
  virtual size_t nf() const{
    return 1;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==0){
        zjucad::matrix::itr_matrix<const val_type*> x0(3,tet_num_,&x[0]);
        zjucad::matrix::matrix<val_type> sh(9,1);
        calc_rot_cubic_f_sh_(&sh[0], &x0(0,tet_idx_));

        int_type c[1] = {0};
        cv[c] += w_*(ks_*zjucad::matrix::dot(basis_, sh)-v_);
      }
    if(k==1){
        zjucad::matrix::itr_matrix<const val_type*> x0(3,tet_num_,&x[0]);
        zjucad::matrix::matrix<val_type> jac_sh(9,3);
        calc_jac_rot_cubic_f_sh_(&jac_sh[0], &x0(0,tet_idx_));

        zjucad::matrix::matrix<val_type> jac = basis_ * jac_sh;
        for(int_type i = 0; i < 3; ++i){
            int_type c[2] = {0,3*tet_idx_+i};
            cv[c] += w_*ks_*jac[i];
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==1){
        for(int_type i = 0; i < 3; ++i){
            int_type c[2] = {0,3*tet_idx_+i};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k==0) return -1;
    if(k==1) return 3;
  }
private:
  const int_type tet_num_;
  const int_type tet_idx_;
  zjucad::matrix::matrix<val_type> basis_;
  const val_type w_;
  const val_type v_;
  const val_type ks_;
};

template <typename val_type>
std::shared_ptr<const hj::math_func::math_func_t<val_type, int32_t> >
build_normal_align_func(const zjucad::matrix::matrix<size_t> &tet_mesh,
                        const jtf::mesh::face2tet_adjacent &fa,
                        const zjucad::matrix::matrix<size_t> &tri_idx,
                        const zjucad::matrix::matrix<val_type> &surface_normal,
                        const zjucad::matrix::matrix<val_type> &area,
                        const double w_normal_align,
                        const std::string opt_type = "init")
{
  typedef hj::math_func::math_func_t<val_type,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  using namespace std;
  using namespace zjucad::matrix;

  const val_type total_area = std::accumulate(area.begin(), area.end(), 0.0);
  std::shared_ptr<std::vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  matrix<val_type> rot(9,9),zyz(3,1);
  val_type v = 0;
  for(size_t fi = 0; fi < tri_idx.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[tri_idx[fi]];
      assert(fa.is_outside_face(tet_pair));
      const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);

      rot_n_2_z_by_zyz(&surface_normal(0,fi), &zyz[0]);
      calc_rot_cubic_f_sh_mat_(&rot[0], &zyz[0]);

      for(size_t di = 1; di < 8; ++di){
          if(di == 4) v = sqrt(7.0);
          else v = 0;

          if(opt_type == "init")
            all_func->push_back(math_func_ptr(
                                  new surface_normal_align2_sh<val_type,int32_t>(
                                    tet_mesh.size(2), tet_idx, rot, di, v,
                                    sqrt(w_normal_align * area[fi]/total_area))));
          else if(opt_type == "zyz")
            all_func->push_back(math_func_ptr(
                                  new surface_normal_align2_zyz<val_type,int32_t>(
                                    tet_mesh.size(2), tet_idx, rot, di, v,
                                    sqrt(w_normal_align * area[fi]/total_area))));
          else throw std::invalid_argument("# [error] can not recognize opt type.");
        }
    }
  cerr << "# [info] normal_align_func number " << all_func->size() << endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<val_type,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

template <typename val_type>
std::shared_ptr<const hj::math_func::math_func_t<val_type, int32_t> >
build_full_surface_cons_align(const zjucad::matrix::matrix<size_t> &tet_mesh,
                              const jtf::mesh::face2tet_adjacent &fa,
                              const zjucad::matrix::matrix<size_t> &tri_idx,
                              const zjucad::matrix::matrix<zjucad::matrix::matrix<val_type> > &full_tri_cons,
                              const zjucad::matrix::matrix<val_type> &area,
                              const double w_normal_align,
                              const std::string opt_type = "init")
{
  assert(full_tri_cons.size() == tri_idx.size());
  typedef hj::math_func::math_func_t<val_type,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  using namespace std;
  using namespace zjucad::matrix;

  const val_type total_area = std::accumulate(area.begin(), area.end(), 0.0);
  std::shared_ptr<std::vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  matrix<val_type> rot(9,9),zyz(3,1);
  val_type v = 0;
  for(size_t fi = 0; fi < tri_idx.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[tri_idx[fi]];
      assert(fa.is_outside_face(tet_pair));
      const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);

      for(size_t di = 0; di < full_tri_cons[fi].size(2); ++di){
          rot_n_2_z_by_zyz(&full_tri_cons[fi](0,di), &zyz[0]);
          calc_rot_cubic_f_sh_mat_(&rot[0], &zyz[0]);

          for(size_t di = 1; di < 8; ++di){
              if(di == 4) v = sqrt(7.0);
              else v = 0;

              if(opt_type == "init")
                all_func->push_back(math_func_ptr(
                                      new surface_normal_align2_sh<val_type,int32_t>(
                                        tet_mesh.size(2), tet_idx, rot, di, v,
                                        sqrt(w_normal_align * area[fi]/(total_area*full_tri_cons[fi].size(2))))));
              else if(opt_type == "zyz")
                all_func->push_back(math_func_ptr(
                                      new surface_normal_align2_zyz<val_type,int32_t>(
                                        tet_mesh.size(2), tet_idx, rot, di, v,
                                        sqrt(w_normal_align * area[fi]/(total_area*full_tri_cons[fi].size(2))))));
              else throw std::invalid_argument("# [error] can not recognize opt type.");
            }
        }
    }
  cerr << "# [info] normal_align_func number " << all_func->size() << endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<val_type,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

template <typename val_type>
std::shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_feature_line_align_func(const zjucad::matrix::matrix<size_t> &tet_mesh,
                              const zjucad::matrix::matrix<double> & node,
                              const std::vector<std::vector<size_t> > & feature_lines,
                              const jtf::mesh::one_ring_tet_at_edge &ortae,
                              const double weight_feature,
                              const std::string opt_type = "init",
                              const bool use_all_one_ring_tets = false)
{
  using namespace std;
  using namespace zjucad::matrix;
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  double total_len = 0;
  vector<vector<double> > lens(feature_lines.size());
  for(size_t ei = 0; ei < feature_lines.size(); ++ei){
      const vector<size_t> & one_line = feature_lines[ei];
      for(size_t pi = 1; pi < one_line.size(); ++pi){
          lens[ei].push_back(norm(node(colon(), one_line[pi-1])- node(colon(), one_line[pi])));
          total_len += lens[ei].back();
        }
    }

  matrix<double> rot(9,9), zyz(3,1);

  matrix<double> dir(3,1);
  double v = 0;
  set<size_t> used_tets;
  for(size_t ei = 0; ei < feature_lines.size(); ++ei){
      const vector<size_t> & one_line = feature_lines[ei];
      for(size_t pi = 1; pi < one_line.size(); ++pi){
          dir = node(colon(), one_line[pi]) - node(colon(), one_line[pi-1]);
          dir /= norm(dir);
          auto it = ortae.e2t_.find(make_pair(one_line[pi], one_line[pi-1]));
          if(it == ortae.e2t_.end()){
              it = ortae.e2t_.find(make_pair(one_line[pi-1], one_line[pi]));
            }
          if(it == ortae.e2t_.end()){
              throw std::logic_error("can not find edge in ortae.");
            }

          rot_n_2_z_by_zyz(&dir[0], &zyz[0]);
          calc_rot_cubic_f_sh_mat_(&rot[0], &zyz[0]);

          const vector<size_t> & around_tets = it->second;
          used_tets.clear();
          if(use_all_one_ring_tets){
              used_tets.insert(around_tets.begin(), around_tets.end());
            }else{
              assert(around_tets.size() > 1);
              if(around_tets.front() == -1) used_tets.insert(around_tets[1]);
              if(around_tets.back() == -1) used_tets.insert(around_tets[around_tets.size()-2]);
            }

          for(const auto & one_tet_idx : used_tets){
              if(one_tet_idx == -1) continue;
              for(size_t di = 1; di < 8; ++di){
                  if(di == 4) {v = sqrt(7.0);}
                  else {v = 0;}
                  if(opt_type == "init")
                    all_func->push_back(math_func_ptr(
                                          new surface_normal_align2_sh<double,int32_t>(
                                            tet_mesh.size(2), one_tet_idx, rot, di, v,
                                            sqrt(weight_feature*lens[ei][pi]/total_len))));
                  else if(opt_type == "zyz")
                    all_func->push_back(math_func_ptr(
                                          new surface_normal_align2_zyz<double,int32_t>(
                                            tet_mesh.size(2), one_tet_idx, rot, di, v,
                                            sqrt(weight_feature*lens[ei][pi]/total_len))));
                  else throw std::invalid_argument("# [error] can not recognize opt type.");
                }
            }
        }
    }

  cerr << "# [info] feature_align_func number " << all_func->size() << endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}
#endif // FRAME_FUNC_H
