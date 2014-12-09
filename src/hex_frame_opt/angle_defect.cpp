#include "angle_defect.h"
#include <vector>
#include <map>
#include <fstream>
#include <jtflib/math/math.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include <jtflib/function/func_aux.h>
#include <hjlib/function/function.h>
#include <zjucad/optimizer/optimizer.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <boost/property_tree/ptree.hpp>
#include <jtflib/function/function.h>
#include <jtflib/optimizer/optimizer.h>
#include "../common/solver_ipopt.h"
#include "../common/vtk.h"
#include "../common/util.h"
#include "../numeric/util.h"

using namespace std;
using boost::property_tree::ptree;

class theta_func : public hj::function::function_t<double,int32_t>
{
public:
  theta_func(const size_t num,
             const size_t idx)
    :num_(num), idx_(idx){}
  virtual ~theta_func(){}
  virtual size_t dim_of_f() const
  {
    return 1;
  }
  virtual size_t dim_of_x() const
  {
    return num_;
  }
  virtual int val(const value_type *x, value_type *f, hj::function::func_ctx *ctx) const
  {
    *f = x[idx_];
    return 0;
  }
  virtual int jac(const value_type *x, value_type *val, int_type *ptr, int_type *idx, hj::function::func_ctx *ctx) const
  {
    ptr[1] = ptr[0] + 1;
    idx[ptr[0]] = idx_;
    val[ptr[0]] = 1;
    return 0;
  }
  virtual size_t jac_nnz() const
  {
    return 1;
  }
private:
  const size_t num_;
  const size_t idx_;
};

class defect_func : public hj::function::function_t<double,int32_t>
{
public:
  defect_func(const size_t num,
              const vector<size_t> &edges,
              const double angle_defect,
              const double w)
    :num_(num), edges_(edges), angle_(angle_defect), w_(w){}
  virtual ~defect_func(){}
  virtual size_t dim_of_f() const
  {
    return 1;
  }
  virtual size_t dim_of_x() const
  {
    return num_;
  }
  virtual int val(const value_type *x, value_type *f, hj::function::func_ctx *ctx) const
  {
    *f = -angle_;
    for(size_t i = 0 ; i< edges_.size(); ++i){
        *f += x[edges_[i]];
      }
    *f *= w_;
    return 0;
  }
  virtual int jac(const value_type *x, value_type *val, int_type *ptr, int_type *idx,
                  hj::function::func_ctx *ctx) const
  {
    ptr[1] = ptr[0] + edges_.size();
    for(size_t i = 0; i < edges_.size(); ++i){
        idx[ptr[0] + i] = edges_[i];
        val[ptr[0] + i] = w_;
      }
    return 0;
  }
  virtual size_t jac_nnz() const
  {
    return edges_.size();
  }
private:
  const size_t num_;
  const vector<size_t> edges_;
  const double angle_;
  const double w_;
};


class quadratic_sum : public jtf::function::functionN1_t<double,int32_t>
{
public:
  quadratic_sum(const size_t num,
                const zjucad::matrix::matrix<double> & area):num_(num), area_(area){}

  virtual ~quadratic_sum(){}

  virtual size_t dim() const
  {
    return num_;
  }
  virtual int val(const double *x, double &v)
  {
    zjucad::matrix::itr_matrix<const double *> x0(num_,1,x);
    for(size_t i = 0; i < num_; ++i)
      v += area_[i]*0.5*x[i] * x[i];
    return 0;
  }
  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx)
  {
    if(g == 0 && idx == 0){
        nnz = num_;
        return 0;
      }
    assert(g != 0 && idx != 0);
    for(size_t i = 0; i < num_; ++i){
        g[i] = x[i] * area_[i];
        idx[i] = i;
      }
    return 0;
  }
  virtual int gra(const double *x, double *g)
  {
    zjucad::matrix::itr_matrix<const double*> x0(num_,1,x);
    zjucad::matrix::itr_matrix<double*> g0(num_,1,g);

    for(size_t i = 0; i < num_; ++i){
        g0[i] += x0[i] * area_[i];
      }
    return 0;
  }
  virtual int hes(const double *x, size_t &nnz, size_t &format,
                  double *h, int32_t *ptr, int32_t *idx, double alpha)
  {
    if(h == 0 && ptr == 0 && idx == 0){
        nnz = num_;
        format = 1;
        return 0;
      }
    if(h == 0 && ptr != 0 && idx != 0){
        for(size_t i = 0 ; i < num_; ++i){
            ptr[i+1] = ptr[i] + 1;
            idx[ptr[i] ] = i;
          }
      }
    if(h != 0 && ptr != 0 && idx != 0){
        for(size_t i = 0; i < num_; ++i){
            jtf::function::add_to_csc(h,ptr,idx, i,i,area_[i]);
          }
      }
    return 0;
  }
  virtual int hes_block(const double *x, double *h, double alpha)
  {
    return 0;
  }
private:
  const size_t num_;
  const zjucad::matrix::matrix<double> & area_;
};

class angle_defect_sum_func: public jtf::function::functionN1_t<double,int32_t>
{
public:
  angle_defect_sum_func(const size_t num, const map<size_t,double> & edge_idx_with_order,
                        double angle, double weight = 1, int order = 1)
    :num_(num), edge_idx_with_order_(edge_idx_with_order), angle_(angle),order_(order),weight_(weight){
    assert(order_ == 1 || order_ == 2);
  }

  virtual ~angle_defect_sum_func(){}

  virtual size_t dim()const{
    return num_;
  }

  virtual int val(const double *x, double &v)
  {
    if(order_ == 1){
        double sum_ = sum_angle(x);
        v += weight_*(sum_ - angle_);
      }
    if(order_ == 2){
        double sum_ = sum_angle(x);
        v += 0.5*weight_*pow(sum_ - angle_,2);
      }
    return 0;
  }

  virtual int gra(const double *x, double *g)
  {
    if(order_==1){
        for(const auto & one_edge_idx : edge_idx_with_order_){
            g[one_edge_idx.first] += weight_*one_edge_idx.second;
          }
      }else{
        double sum_ = sum_angle(x);
        for(const auto & one_edge_idx : edge_idx_with_order_){
            g[one_edge_idx.first] += weight_*one_edge_idx.second*(sum_-angle_);
          }
      }
    return 0;
  }

  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx)
  {
    if(g == 0 && idx == 0){
        nnz = edge_idx_with_order_.size();
        return 0;
      }
    assert(g != 0 && idx != 0);
    size_t i = 0;
    if(order_ == 1){
        for(const auto & one_edge_idx : edge_idx_with_order_){
            g[i] = weight_*one_edge_idx.second;
            idx[i] = one_edge_idx.first;
            ++i;
          }
      }else{
        double sum_ = sum_angle(x);
        for(const auto & one_edge_idx : edge_idx_with_order_){
            g[i] = weight_*one_edge_idx.second*(sum_-angle_);
            idx[i] = one_edge_idx.first;
            ++i;
          }
      }
    return 0;
  }

  virtual int hes(const double *x, size_t &nnz, size_t &format,
                  double *h, int32_t *ptr, int32_t *idx, double alpha)
  {
    if(order_==1)
      nnz = 0;
    else{
        if(h == 0 && ptr == 0 && idx == 0){
            nnz = edge_idx_with_order_.size() * edge_idx_with_order_.size();
            format = 1;
            return 0;
          }
        if(h == 0 && ptr != 0 && idx != 0){
            for(const auto & one_edge : edge_idx_with_order_){
                ptr[one_edge.first+1] = ptr[one_edge.first]+edge_idx_with_order_.size();
                size_t i = 0;
                for(const auto & one_edge_i : edge_idx_with_order_){
                    idx[ptr[one_edge.first]+i] = one_edge_i.first;
                    ++i;
                  }
              }
            return 0;
          }
        if(h != 0 && ptr != 0 && idx != 0){
            for(const auto & one_edge_i : edge_idx_with_order_){
                for(const auto & one_edge_j : edge_idx_with_order_)
                  jtf::function::add_to_csc(h, ptr, idx, one_edge_i.first, one_edge_j.first,
                                            weight_*one_edge_i.second*one_edge_j.second);
              }
            return 0;
          }
      }
    return 0;
  }

  double sum_angle(const double *x) const{
    double sum_ = 0;
    for(const auto & one_edge_idx : edge_idx_with_order_)
      sum_ += x[one_edge_idx.first] * one_edge_idx.second;
    return sum_;
  }

  virtual int hes_block(const double *x, double *h, double alpha)
  {
    return 0;
  }
private:
  const size_t num_;
  const map<size_t,double> edge_idx_with_order_;
  const double angle_;
  const int order_;
  const double weight_;
};

void naive_angle_defect::cal_angle_defect(
    const jtf::mesh::one_ring_face_at_point & orfap,
    const zjucad::matrix::matrix<size_t> & face,
    const zjucad::matrix::matrix<double> & node,
    map<size_t,double> & angle_defect_map)
{
  vector<double> angles;
  size_t other_points[2];
  zjucad::matrix::matrix<double> two_edge(3,2);
  for(jtf::mesh::one_ring_face_at_point::p2f_type::const_iterator cit
      = orfap.p2f_.begin(); cit != orfap.p2f_.end(); ++cit){
      angles.clear();
      const vector<size_t> &face_vec = cit->second;
      set<size_t> face_set(face_vec.begin(), face_vec.end());

      if(face_vec.front() != -1){
          for(const auto & face_idx : face_set){
              size_t idx = -1;
              for(size_t j = 0; j< face.size(1); ++j){
                  if(face(j, face_idx) == cit->first) {
                      idx = j;
                      break;
                    }
                }
              if(idx == -1){

                  throw std::logic_error("one_ring_face_at_point is not compatiable with face");
                }
              other_points[0] = face((idx+1)%face.size(1), face_idx);
              other_points[1] = face((idx+2)%face.size(1), face_idx);
              two_edge(zjucad::matrix::colon(),0) = node(zjucad::matrix::colon(), other_points[0]) - node(zjucad::matrix::colon(), cit->first);
              two_edge(zjucad::matrix::colon(),1) = node(zjucad::matrix::colon(), other_points[1]) - node(zjucad::matrix::colon(), cit->first);
              const double len0 = zjucad::matrix::norm(two_edge(zjucad::matrix::colon(),0));
              const double len1 = zjucad::matrix::norm(two_edge(zjucad::matrix::colon(),1));

              angles.push_back(jtf::math::safe_acos(zjucad::matrix::dot(two_edge(zjucad::matrix::colon(),0), two_edge(zjucad::matrix::colon(),1))/(len0*len1)));
            }
          angle_defect_map[cit->first] = 2*My_PI()-std::accumulate(angles.begin(),angles.end(),0.0);
        }else{
          angle_defect_map[cit->first] = 0;
        }
    }
}

void naive_angle_defect::opt(const jtf::mesh::edge2cell_adjacent &ea,
                             const zjucad::matrix::matrix<size_t> &face,
                             const zjucad::matrix::matrix<double> &node,
                             zjucad::matrix::matrix<double> &angle,
                             boost::property_tree::ptree &opt)
{
  angle = zjucad::matrix::zeros(ea.edges_.size(),1);

  std::shared_ptr<vector<shared_ptr<hj::function::function_t<double,int32_t> > > > funcs(
        new vector<shared_ptr<hj::function::function_t<double,int32_t> > >);

  {// add function of \|\theta\|^2
    for(size_t i = 0 ; i < ea.edges_.size(); ++i)
      funcs->push_back(std::shared_ptr<hj::function::function_t<double,int32_t> >(new theta_func(ea.edges_.size(), i)));
  }

  jtf::mesh::one_ring_face_at_point orfap;
  orfap.add_all_faces(face, ea);
  map<size_t,double> angle_defect_map;
  cal_angle_defect(orfap, face, node, angle_defect_map);

  double angle_defect_total = 0;
  for(const auto & one_point : angle_defect_map){
      angle_defect_total += one_point.second;
    }
  cerr << "# before calc, total angle defect " << angle_defect_total << endl;

  std::unique_ptr<jtf::mesh::one_ring_point_at_point> orpap(
        jtf::mesh::one_ring_point_at_point::create(face));
  if(!orpap.get()){
      throw std::invalid_argument("can not build one_ring_point_at_point");
    }

  {// add defect for each point
    vector<size_t> edge_idx;
    for(jtf::mesh::one_ring_point_at_point::p2p_type::const_iterator
        cit = orpap->p2p_.begin(); cit != orpap->p2p_.end(); ++cit){
        edge_idx.clear();
        const vector<size_t> &points= cit->second;
        for(size_t i = 0 ; i < points.size(); ++i){
            edge_idx.push_back(ea.get_edge_idx(cit->first, points[i]));
          }
        funcs->push_back(std::shared_ptr<hj::function::function_t<double,int32_t> >(
                           new defect_func(ea.edges_.size(),edge_idx, angle_defect_map.at(cit->first),10)));
      }
  }

  shared_ptr<hj::function::function_t<double,int32_t> > all_func(
        hj::function::new_catenated_function<double,int32_t>(funcs));

  zjucad::matrix::matrix<double> residual(all_func->dim_of_f(),1);
  boost::property_tree::ptree pt;
  pt.put("package.value", "hj");
  pt.put("iter.value", "20");
  zjucad::optimize(*all_func, angle, residual, pt);
}

int manually_set_angle_defect::cal_angle_defect2(
    const jtf::mesh::one_ring_face_at_point & orfap,
    const zjucad::matrix::matrix<size_t> & face,
    const zjucad::matrix::matrix<double> & node,
    const zjucad::matrix::matrix<double> & normal,
    map<size_t,double> & angle_defect)
{
  using namespace zjucad::matrix;
  cerr << "# [warning] one_ring_face_at_point should be ordered." << endl;

  for(jtf::mesh::one_ring_face_at_point::p2f_type::const_iterator cit =
      orfap.p2f_.begin(); cit != orfap.p2f_.end(); ++cit){
      const vector<size_t> & face_loop = cit->second;

      angle_defect[cit->first] = cal_angle_defect_along_circle(face_loop, node, face,normal);

      //      double test_angle = cal_angle_defect_directly(face_loop, node, face, cit->first);

      //      if(std::fabs(angle_defect[cit->first]-test_angle) > 1e-12){
      //          cerr << "# [error] wrong point angle defect " << cit->first << " "
      //               << angle_defect[cit->first] << " " << test_angle << endl;
      //         // exit(0);
      //        }

      // angle_defect[cit->first] = test_angle;
    }

  return 0;
}

double manually_set_angle_defect::cal_cota_angle(
    const size_t face_idx,
    const zjucad::matrix::matrix<size_t> &face,
    const pair<size_t,size_t> & one_edge,
    const zjucad::matrix::matrix<double> &node)
{
  size_t other_idx = std::accumulate(face(zjucad::matrix::colon(), face_idx).begin(),
                                     face(zjucad::matrix::colon(), face_idx).end(), int(0))
      - one_edge.first - one_edge.second;
  zjucad::matrix::matrix<double> dir0 = node(zjucad::matrix::colon(), one_edge.first) - node(zjucad::matrix::colon(), other_idx);
  zjucad::matrix::matrix<double> dir1 = node(zjucad::matrix::colon(), one_edge.second) - node(zjucad::matrix::colon(), other_idx);
  double cos_angle = zjucad::matrix::dot(dir0,dir1)/(norm(dir0)*norm(dir1));
  double sin_angle = sqrt(1-cos_angle*cos_angle);
  return cos_angle / sin_angle;
}

double manually_set_angle_defect::cal_angle_defect_directly(
    const std::vector<size_t> & face_loop,
    const zjucad::matrix::matrix<double> & node,
    const zjucad::matrix::matrix<size_t> & face,
    const size_t this_point)
{
  using namespace zjucad::matrix;
  if(face_loop.front() == -1 || face_loop.back() == -1) return 0;
  double angle = 0;
  zjucad::matrix::matrix<double> two_dir(3,2);
  for(size_t fi = 0; fi < face_loop.size()-1; ++fi){
      auto it = std::find(face(colon(),face_loop[fi]).begin(),
                          face(colon(),face_loop[fi]).end(),
                          this_point);
      assert(it != face(colon(), face_loop[fi]).end());
      size_t idx = it-face(colon(),face_loop[fi]).begin();
      size_t next = (idx+1)%face.size(1);
      size_t next_next = (next+1)%face.size(1);

      two_dir(colon(), 0) = node(colon(), face(next,face_loop[fi]))
          - node(colon(), face(idx,face_loop[fi]));
      two_dir(colon(),0) /= norm(two_dir(colon(),0));

      two_dir(colon(), 1) = node(colon(), face(next_next,face_loop[fi]))
          - node(colon(), face(idx,face_loop[fi]));
      two_dir(colon(),1) /= norm(two_dir(colon(),1));

      angle += jtf::math::safe_acos(dot(two_dir(colon(),0), two_dir(colon(),1)));
    }
  angle = 2*My_PI()-angle;
  return angle;
}

double manually_set_angle_defect::cal_angle_defect_along_circle(
    const vector<size_t> & face_loop,
    const zjucad::matrix::matrix<double> & node,
    const zjucad::matrix::matrix<size_t> & face,
    const zjucad::matrix::matrix<double> & normal)
{
  if(face_loop.front() != face_loop.back()){
      cerr << "# [error] it's a sharp edge" << endl;
      return 0;
    }
  if(face_loop.front() == -1 || face_loop.back() == -1){
      cerr << "# [error] it's an open loop."<< endl;
      return 0;
    }

  double init_angle = 0;
  double next_angle = 0;

  for(size_t fi = 0; fi < face_loop.size()-1; ++fi){
      next_angle = calc_next_angle_along_circle(
            face_loop[fi], face_loop[fi+1], face, normal, node, next_angle);
    }

  return next_angle - init_angle;
}

double manually_set_angle_defect::calc_next_angle_along_circle(
    const size_t fi, const size_t fj,
    const zjucad::matrix::matrix<size_t> & face,
    const zjucad::matrix::matrix<double> & normal,
    const zjucad::matrix::matrix<double> & node,
    const double init_angle) // init_angle is based on the first edge of fi
{
  using namespace zjucad::matrix;
  size_t common_edge[2];
  jtf::mesh::find_common_edge(face(colon(),fi), face(colon(),fj), common_edge);
  orient_edge(common_edge, face(colon(),fi));

  matrix<double> di = zeros<double>(3,1);
  find_safe_reference(face(colon(),fi), node, common_edge, di);

  matrix<double> d_common = node(colon(), common_edge[1]) - node(colon(), common_edge[0]);
  double angle_dij = calc_angle_from_d0_to_d1(di,d_common,normal(colon(),fi));

  matrix<double> dj = zeros<double>(3,1);
  find_safe_reference(face(colon(),fj), node, common_edge, dj);
  dj *= -1;

  double angle_dji = calc_angle_from_d0_to_d1(dj, d_common, normal(colon(),fj));

  double next_angle = float_mod(init_angle - angle_dij + angle_dji, 2*My_PI());
  if(fabs(next_angle)> My_PI() && next_angle < 0)
    next_angle += 2*My_PI();
  else if(fabs(next_angle) > My_PI() && next_angle > 0)
    next_angle -= 2*My_PI();

  return next_angle;
}

void manually_set_angle_defect::orient_edge(
    size_t * common_edge,
    const zjucad::matrix::matrix<size_t> & one_face)
{
  assert(common_edge);
  auto it = find(one_face.begin(), one_face.end(), common_edge[0]);
  assert(it != one_face.end());
  size_t idx = it-one_face.begin();
  if(one_face[(idx+1)%one_face.size()] != common_edge[1]){
      if(one_face[(idx-1+one_face.size())%one_face.size()] == common_edge[1]){
          swap(common_edge[0], common_edge[1]);
        }else{
          throw std::invalid_argument("this common edge does not belong to face");
        }
    }
}

void manually_set_angle_defect::find_safe_reference(
    const zjucad::matrix::matrix<size_t> & one_face,
    const zjucad::matrix::matrix<double> & node,
    const size_t *one_edge,
    zjucad::matrix::matrix<double> & ri)
{
  using namespace zjucad::matrix;
  assert(one_edge);
  assert(find(one_face.begin(), one_face.end(), one_edge[0]) != one_face.end());
  assert(find(one_face.begin(), one_face.end(), one_edge[1]) != one_face.end());

  size_t other_point_idx = std::accumulate(one_face.begin(), one_face.end(), 0) - one_edge[0] - one_edge[1];
  ri = node(colon(), one_edge[0]) - node(colon(), other_point_idx);
  const double len = norm(ri);
  if(len > 1e-6)
    ri /= len;
}

double manually_set_angle_defect::calc_angle_from_d0_to_d1(
    const zjucad::matrix::matrix<double> & d0,
    const zjucad::matrix::matrix<double> & d1,
    const zjucad::matrix::matrix<double> & normal)
{
  double acos_01 = jtf::math::safe_acos(zjucad::matrix::dot(d0,d1)/(norm(d0)*norm(d1)));
  if(zjucad::matrix::dot(zjucad::matrix::cross(d0,d1),normal) < 0)
    acos_01 *= -1;
  return acos_01;
}

void manually_set_angle_defect::opt(
    const jtf::mesh::edge2cell_adjacent &ea,
    const zjucad::matrix::matrix<size_t> &face,
    const zjucad::matrix::matrix<double> &node,
    zjucad::matrix::matrix<double> &angle,
    boost::property_tree::ptree &pt)
{
  map<size_t,double> angle_defect_map, nsys_angle_defect_map;
  if(load_angle_defect_set_file(pt.get<string>("angle_defect_set.value").c_str(), nsys_angle_defect_map))
    throw std::invalid_argument("wrong angle defect file.");

  zjucad::matrix::matrix<double> face_normal(3, face.size(2));
  jtf::mesh::cal_face_normal(face, node, face_normal);

  jtf::mesh::one_ring_face_at_point orfap;
  orfap.add_all_faces(face, ea);
  orfap.sort_int_loop_with_normal_info(face, node, ea, face_normal);
  //cal_angle_defect(orfap, face, node, nsys_angle_defect_map, angle_defect_map);
  cal_angle_defect2(orfap, face, node, face_normal, angle_defect_map);


  for(map<size_t,double>::const_iterator cit = nsys_angle_defect_map.begin();
      cit != nsys_angle_defect_map.end(); ++cit){
      angle_defect_map.at(cit->first) -= cit->second;
    }

  zjucad::matrix::matrix<double> face_area(face.size(2),1);
  for(size_t fi = 0; fi < face.size(2); ++fi)
    face_area[fi] = jtf::mesh::cal_face_area(&face(0,fi), face.size(1), node);

  zjucad::matrix::matrix<double> edge_weight(ea.edges_.size(),1);
  if(1){
      for(size_t ei = 0; ei < ea.edges_.size(); ++ei){
          const pair<size_t,size_t> & face_pair = ea.edge2cell_[ei];
          assert(ea.is_boundary_edge(face_pair) == false);
          edge_weight[ei] =  (face_area[face_pair.first] + face_area[face_pair.second])/3.0;
        }
      const double total_area = std::accumulate(edge_weight.begin(), edge_weight.end(), 0.0);
      edge_weight /= total_area;
    }else{
      for(size_t ei = 0; ei < ea.edges_.size(); ++ei){
          const pair<size_t,size_t> & face_pair = ea.edge2cell_[ei];
          assert(ea.is_boundary_edge(face_pair) == false);
          edge_weight[ei] = 2/(cal_cota_angle(face_pair.first, face, ea.edges_[ei], node)
                               +cal_cota_angle(face_pair.first, face, ea.edges_[ei], node));
        }
    }

  shared_ptr<jtf::function::functionN1_t<double,int32_t> > obj_func(
        new quadratic_sum(ea.edges_.size(), edge_weight));

  shared_ptr<vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > > eqn_cons(
        new vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > >);

  shared_ptr<jtf::mesh::one_ring_point_at_point> orpap(jtf::mesh::one_ring_point_at_point::create(face));
  if(!orpap.get()){
      throw std::logic_error("# [error] can not build orpap.");
    }

  {
    size_t common_edge[2];
    for(jtf::mesh::one_ring_face_at_point::p2f_type::const_iterator cit = orfap.p2f_.begin();
        cit != orfap.p2f_.end(); ++cit){
        const vector<size_t> & faces_loop = cit->second;
        assert(faces_loop.front() == faces_loop.back());
        map<size_t,double> edge_with_order;
        for(size_t fi = 0; fi < faces_loop.size()-1;++fi){
            jtf::mesh::find_common_edge(face(zjucad::matrix::colon(), faces_loop[fi]),
                                        face(zjucad::matrix::colon(), faces_loop[fi+1]),
                common_edge);
            const size_t edge_idx = ea.get_edge_idx(common_edge);
            const pair<size_t,size_t> & face_pair = ea.edge2cell_[edge_idx];
            if(face_pair.first == faces_loop[fi] && face_pair.second == faces_loop[fi+1]){
                edge_with_order[edge_idx] = 1.0;
              }else if(face_pair.first == faces_loop[fi+1] && face_pair.second == faces_loop[fi]){
                edge_with_order[edge_idx] = -1.0;
              }
          }
        double angle_defect_this_point = 0;
        auto it = angle_defect_map.find(cit->first);
        if(it != angle_defect_map.end()) {
            angle_defect_this_point = -1*it->second;
          }
#define use_ipopt
#ifdef use_ipopt
        eqn_cons->push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                              new angle_defect_sum_func(ea.edges_.size(),edge_with_order, angle_defect_this_point)));
#else
        eqn_cons->push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                              new angle_defect_sum_func(ea.edges_.size(),edge_with_order, angle_defect_this_point, 100,2)));
#endif
      }
  }

  if(angle.size() != ea.edges_.size())
    angle = zjucad::matrix::zeros<double>(ea.edges_.size(),1);
  pt.put("iter.value",10);
  pt.put("package.value","jtf");
  pt.put("alg.value","SQP");
#ifdef use_ipopt
  ipopt_solve(angle, *obj_func,eqn_cons,pt);
#else
  eqn_cons->push_back(obj_func);
  shared_ptr<jtf::function::functionN1_t<double,int32_t> > all_f(
        new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD_CONS>(*eqn_cons));

  jtf::optimize(*all_f, angle, pt,nullptr, nullptr, nullptr);
#endif
}


int manually_set_angle_defect::load_angle_defect_set_file(
    const char *file,
    map<size_t,double> & angle_defect_map)
{
  ifstream ifs(file);
  if(ifs.fail()){
      std::cerr << "# [error] can not open angle defect set file." << std::endl;
      return __LINE__;
    }
  angle_defect_map.clear();
  size_t fix_num = 0;
  ifs >> fix_num;
  size_t p_idx;
  double angle;
  int ki;
  const double Nsys = 4; // 4-cross field
  for(size_t i = 0; i < fix_num; ++i){
      ifs >> p_idx >> ki;
      angle = 2 * My_PI() * ki/ Nsys;
      angle_defect_map[p_idx] = angle;
    }
  return 0;
}


void even_angle_defect::opt(const jtf::mesh::edge2cell_adjacent &ea, const zjucad::matrix::matrix<size_t> &face,
                            const zjucad::matrix::matrix<double> &node,
                            zjucad::matrix::matrix<double> &angle,
                            boost::property_tree::ptree &pt)
{
  map<size_t,double> angle_defect_map;

  zjucad::matrix::matrix<double> face_normal(3, face.size(2));
  jtf::mesh::cal_face_normal(face, node, face_normal);

  jtf::mesh::one_ring_face_at_point orfap;
  orfap.add_all_faces(face, ea);
  orfap.sort_int_loop_with_normal_info(face, node, ea, face_normal);
  cal_angle_defect2(orfap, face, node, face_normal, angle_defect_map);

  zjucad::matrix::matrix<double> face_area(face.size(2),1);
  for(size_t fi = 0; fi < face.size(2); ++fi)
    face_area[fi] = std::fabs(jtf::mesh::cal_face_area(&face(0,fi), face.size(1), node));
  const double total_area = std::accumulate(face_area.begin(), face_area.end(), 0.0);

  zjucad::matrix::matrix<double> edge_weight(ea.edges_.size(),1);
  for(size_t ei = 0; ei < ea.edges_.size(); ++ei){
      const pair<size_t,size_t> & face_pair = ea.edge2cell_[ei];
      assert(ea.is_boundary_edge(face_pair) == false);
      edge_weight[ei] =  (face_area[face_pair.first] + face_area[face_pair.second])/3.0;
    }
  edge_weight /= total_area;

  shared_ptr<jtf::function::functionN1_t<double,int32_t> > obj_func(
        new quadratic_sum(ea.edges_.size(), edge_weight));

  shared_ptr<vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > > eqn_cons(
        new vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > >);

  shared_ptr<jtf::mesh::one_ring_point_at_point> orpap(jtf::mesh::one_ring_point_at_point::create(face));
  if(!orpap.get()){
      throw std::logic_error("# [error] can not build orpap.");
    }

  const set<size_t> face_nodes(face.begin(), face.end());
  const double total_curvature = 2*My_PI() * (static_cast<int>(face_nodes.size()) + static_cast<int>(face.size(2)) - static_cast<int>(ea.edges_.size()));

  cerr << "# [info] total_curvature " << total_curvature << endl;
  const double percent = total_curvature/total_area;

  {
    size_t common_edge[2];
    double total_connection_angle = 0;
    static int vertex_num = 0;
    for(jtf::mesh::one_ring_face_at_point::p2f_type::const_iterator cit = orfap.p2f_.begin();
        cit != orfap.p2f_.end(); ++cit){
        const vector<size_t> & faces_loop = cit->second;
        assert(faces_loop.front() == faces_loop.back());
        map<size_t,double> edge_with_order;
        for(size_t fi = 0; fi < faces_loop.size()-1;++fi){
            jtf::mesh::find_common_edge(face(zjucad::matrix::colon(), faces_loop[fi]),
                                        face(zjucad::matrix::colon(), faces_loop[fi+1]),
                common_edge);
            const size_t edge_idx = ea.get_edge_idx(common_edge);
            const pair<size_t,size_t> & face_pair = ea.edge2cell_[edge_idx];
            if(face_pair.first == faces_loop[fi] && face_pair.second == faces_loop[fi+1]){
                edge_with_order[edge_idx] = 1.0;
              }else if(face_pair.first == faces_loop[fi+1] && face_pair.second == faces_loop[fi]){
                edge_with_order[edge_idx] = -1.0;
              }
          }
        double voroni_cell_area = 0;
        for(size_t fi = 0; fi < faces_loop.size()-1; ++fi){
            if(faces_loop[fi] == -1) continue;
            voroni_cell_area += face_area[faces_loop[fi]]/3.0;
          }

        double angle_defect_this_point = 0;
        auto it = angle_defect_map.find(cit->first);
        if(it != angle_defect_map.end()) {
            ++vertex_num;
            angle_defect_this_point = -1*it->second;
            angle_defect_this_point += percent * voroni_cell_area;
            total_connection_angle += angle_defect_this_point;
          }else {
            cerr << "# [error] strange can not find points in angle_defect_map"<< endl;
          }
        //#define use_ipopt
#ifndef use_ipopt
        eqn_cons->push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                              new angle_defect_sum_func(ea.edges_.size(),edge_with_order, angle_defect_this_point,100,2)));
#else
        eqn_cons->push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                              new angle_defect_sum_func(ea.edges_.size(),edge_with_order, angle_defect_this_point)));

#endif
      }
    cerr << "# [warning] total_connection_angle (should be 0) " << total_connection_angle << endl;
    assert(std::fabs(total_connection_angle) < 1e-8);
  }

  cerr << "# [info] eqn_cons number " << eqn_cons->size() << endl;

  if(angle.size() != ea.edges_.size())
    angle = zjucad::matrix::zeros<double>(ea.edges_.size(),1);
  boost::property_tree::ptree pt_opt;
  pt_opt.put("iter.value",10);
  pt_opt.put("package.value","jtf");
  pt_opt.put("alg.value","SQP");
#ifdef use_ipopt
  ipopt_solve(angle, *obj_func,eqn_cons,pt_opt);
#else
  eqn_cons->push_back(obj_func);
  shared_ptr<jtf::function::functionN1_t<double,int32_t> > all_f(
        new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD_CONS>(*eqn_cons));
  jtf::optimize(*all_f, angle, pt,nullptr, nullptr, nullptr);
#endif
}

bool cmp_pair_large(const std::pair<size_t, double> &a,
                    const std::pair<size_t, double> &b)
{
  return a.second > b.second;
}

void geometry_aware_angle_defect::cal_geodesic_distance(
    vector<vector<pair<size_t,double> > > & distance,
    const zjucad::matrix::matrix<double> &node,
    shared_ptr<const jtf::mesh::one_ring_point_at_point> orpap,
    const double max_dis,
    set<size_t> * feature_points )
{
  using namespace zjucad::matrix;
  if(distance.size() != node.size(2))
    distance.resize(node.size(2));

  vector<bool> feature_point_flag(node.size(2), false);
  if(feature_points){
      for(set<size_t>::const_iterator it = feature_points->begin();
          it != feature_points->end(); ++it){
          feature_point_flag[*it] = true;
        }
    }
#pragma omp parallel
  {
#pragma omp single
    {
      cerr << "# [info] starting calculate geodesic distance " << endl;
      int number = 0;
      int pre_tag = 0;
      cerr << "# [info] ------ 0%";
      for(const auto & one_point : orpap->p2p_) {
          ++number;
          if(int(number*100.0/orpap->p2p_.size()) > pre_tag){
              pre_tag = int(number*100.0/orpap->p2p_.size());
              if(pre_tag < 10) cerr << "\b\b" << pre_tag << "%";
              else if(pre_tag < 101) cerr << "\b\b\b" << pre_tag << "%";
            }
          const vector<size_t> & one_ring_v = one_point.second;
          vector<std::pair<size_t, double> > gathered_points;
          vector<std::pair<size_t, double> > around_points;
          set<size_t> visited_points;
          visited_points.insert(one_point.first);
          for(const auto & one_ring_p : one_ring_v){
              if(feature_point_flag[one_ring_p]) continue;

              around_points.push_back(make_pair(one_ring_p, norm(node(colon(), one_ring_p) - node(colon(), one_point.first))));
              visited_points.insert(one_ring_p);
            }
          std::make_heap(around_points.begin(), around_points.end(), cmp_pair_large);
          while(1){
              if(around_points.size() == 0) break;
              std::pop_heap(around_points.begin(), around_points.end(), cmp_pair_large);
              const pair<size_t,double> min_p = around_points.back();
              around_points.pop_back();
              if(min_p.second > max_dis) break;
              gathered_points.push_back(min_p);
              const vector<size_t> &connected_v = orpap->p2p_.at(min_p.first);
              for(const auto & one_connected_v: connected_v){
                  if(feature_point_flag[one_connected_v]) continue;

                  size_t visited_points_num = visited_points.size();
                  visited_points.insert(one_connected_v);
                  if(visited_points.size() != visited_points_num){
                      const double len = norm(node(colon(), one_connected_v)  - node(colon(), min_p.first));
                      around_points.push_back(make_pair(one_connected_v, len + min_p.second));
                      std::push_heap(around_points.begin(), around_points.end(), cmp_pair_large);
                    }
                }
            }
          distance[one_point.first] = gathered_points;
        }
      cerr << endl;
    }
  }
}

template <typename val_type, int order>
class tylor_exp{
public:
  val_type operator()(val_type x){//should only be called when x in [-1,1]
  }
};

template <typename val_type>
class tylor_exp<val_type, 5>
{
public:
  val_type operator()(val_type x){
    return (5040+x*(5040+x*(2520+x*(840+x*(210+x*(42+x*(7+x)))))))*0.00019841269f;
  }
};

template <typename val_type>
class tylor_exp<val_type, 6>
{
public:
  val_type operator()(val_type x){
    return (40320+x*(40320+x*(20160+x*(6720+x*(1680+x*(336+x*(56+x*(8+x))))))))*2.4801587301e-5;
  }
};

template <typename val_type>
class tylor_exp<val_type, 7>
{
public:
  val_type operator()(val_type x){
    return (362880+x*(362880+x*(181440+x*(60480+x*(15120+x*(3024+x*(504+x*(72+x*(9+x)))))))))*2.75573192e-6;
  }
};

inline double fast_exp(double a)
{
  int a_i = int(a);
  static double exp_table[] = {1,exp(-1),exp(-2),exp(-3),exp(-4)};
  if(fabs(a_i) < 4 && a_i < 0) return exp_table[-a_i]*tylor_exp<double,6>()(a-a_i);
  if(fabs(a_i) < 4 && a_i > 0) return tylor_exp<double,6>()(a-a_i)/exp_table[int(abs(a_i))];
  return exp(a_i)*tylor_exp<double,6>()(a-a_i);
}

void geometry_aware_angle_defect::filter_point_curvature(
    const double max_dis,
    const zjucad::matrix::matrix<size_t> &face,
    const zjucad::matrix::matrix<double> &node,
    const std::map<size_t,double> &k2,
    std::map<size_t, double> &kcorr,
    set<size_t> *feature_points)
{
  using namespace zjucad::matrix;

  shared_ptr<jtf::mesh::one_ring_point_at_point> orpap(jtf::mesh::one_ring_point_at_point::create(face));
  if(!orpap.get()){
      throw std::logic_error("can not build one_ring_point_at_point");
    }

  vector<vector<pair<size_t,double> > > distance(node.size(2));
  vector<double> dual_cell_area(node.size(2),0);
  {
    for(size_t fi = 0; fi < face.size(2); ++fi){
        double area = jtf::mesh::cal_face_area(face(colon(),fi), node);
        for(size_t pi = 0; pi < face.size(1); ++pi){
            dual_cell_area[face(pi,fi)] += area/3.0;
          }
      }
  }

  vector<double> k2_vec(node.size(2),0);
  for(const auto & one_p : k2){
      k2_vec[one_p.first] = one_p.second;
    }
  cal_geodesic_distance(distance, node, orpap, max_dis, feature_points);

  cerr << "# [info] starting calculate kcorr " << endl;

  cerr << "# [info] ------ 0%";
  size_t vj_idx = 0;
  const double max_dis_2 = 1.0/(max_dis*max_dis);

  //#pragma omp parallel for private(vj_idx)
  for(vj_idx = 0; vj_idx < distance.size(); ++vj_idx){
      const vector<pair<size_t,double> > &vi_vec = distance[vj_idx];
      if(vj_idx %(orpap->p2p_.size()/100) == 0){
          int pre_tag = int(vj_idx*100.0/orpap->p2p_.size());
          if(pre_tag < 10) cerr << "\b\b" << pre_tag << "%";
          else if(pre_tag < 101) cerr << "\b\b\b" << pre_tag << "%";
        }
      if(vi_vec.empty()) {
          kcorr[vj_idx] = k2_vec.at(vj_idx); //float_mod(k2.at(vj.first),My_PI()/2.0);
          continue;
        }else{
          double k = 0;
          for(const auto & vi : vi_vec){
              double sum_cik = 0;
              const vector<pair<size_t,double> > &vk_vec = distance.at(vi.first);
              for(const auto & vk : vk_vec){
                  //                  const double cik = dual_cell_area.at(vk.first)*exp(-1*pow(2*vk.second/max_dis,2));
                  //                  sum_cik += cik;
                  const double cik = dual_cell_area.at(vk.first)*fast_exp(-4*vk.second*vk.second*max_dis_2);
                  sum_cik += cik;
                }
              if(fabs(sum_cik) < 1e-8) {
                  cerr << "# [warring] gauss radio may be two small" << endl;
                  k += k2_vec.at(vi.first);
                }else
                k += dual_cell_area.at(vj_idx)*fast_exp(-4*vi.second*vi.second*max_dis_2)
                    * k2_vec.at(vi.first) /sum_cik;
            }
          kcorr[vj_idx] = k;
        }
    }
  cerr << endl;
}

void geometry_aware_angle_defect::opt(const jtf::mesh::edge2cell_adjacent &ea,
                                      const zjucad::matrix::matrix<size_t> &face,
                                      const zjucad::matrix::matrix<double> &node,
                                      zjucad::matrix::matrix<double> &angle,
                                      boost::property_tree::ptree &pt)
{
  map<size_t,double> angle_defect_map;

  zjucad::matrix::matrix<double> face_normal(3, face.size(2));
  jtf::mesh::cal_face_normal(face, node, face_normal);

  jtf::mesh::one_ring_face_at_point orfap;
  orfap.add_all_faces(face, ea);
  orfap.sort_int_loop_with_normal_info(face, node, ea, face_normal);
  cal_angle_defect2(orfap, face, node, face_normal, angle_defect_map);

  zjucad::matrix::matrix<double> face_area(face.size(2),1);
  for(size_t fi = 0; fi < face.size(2); ++fi)
    face_area[fi] = std::fabs(jtf::mesh::cal_face_area(&face(0,fi), face.size(1), node));
  const double total_area = std::accumulate(face_area.begin(), face_area.end(), 0.0);

  zjucad::matrix::matrix<double> edge_weight(ea.edges_.size(),1);
  for(size_t ei = 0; ei < ea.edges_.size(); ++ei){
      const pair<size_t,size_t> & face_pair = ea.edge2cell_[ei];
      assert(ea.is_boundary_edge(face_pair) == false);
      edge_weight[ei] =  (face_area[face_pair.first] + face_area[face_pair.second])/3.0;
    }
  edge_weight /= total_area;

  shared_ptr<jtf::function::functionN1_t<double,int32_t> > obj_func(
        new quadratic_sum(ea.edges_.size(), edge_weight));

  shared_ptr<vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > > eqn_cons(
        new vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > >);

  shared_ptr<jtf::mesh::one_ring_point_at_point> orpap(jtf::mesh::one_ring_point_at_point::create(face));
  if(!orpap.get()){
      throw std::logic_error("# [error] can not build orpap.");
    }

  map<size_t,double> kcorr;
  {
    const double tor = pt.get<double>("weight/gaussian_radius.value");
    zjucad::matrix::matrix<double> bb(3,2);
    calc_bounding_box(node, &bb[0]);
    bb(zjucad::matrix::colon(),0) -= bb(zjucad::matrix::colon(),1);
    const double len = zjucad::matrix::max(zjucad::matrix::fabs(bb(zjucad::matrix::colon(),0)));

    set<size_t> feature_points;
    if(zjucad::has("input/feature_line.value", pt)){
        vector<vector<size_t> > feature_lines;
        if(zjucad::has("input/s2v.value", pt)){
            if(jtf::mesh::load_feature_line(
                 pt.get<string>("input/feature_line.value").c_str(),
                 pt.get<string>("input/s2v.value").c_str(), feature_lines))
              throw std::invalid_argument("# [error] can not load feature line file.");
          }else{
            if(jtf::mesh::load_feature_line(
                 pt.get<string>("input/feature_line.value").c_str(),feature_lines)){
                throw std::invalid_argument("# [error] can not load feature line file.");
              }
          }
        for(const auto & one_line : feature_lines)
          feature_points.insert(one_line.begin(), one_line.end());
        cerr << "# [info] use feature line." << endl;
      }


    filter_point_curvature(2*tor*len, face, node, angle_defect_map, kcorr, &feature_points);
  }

  {
    set<size_t> face_node(face.begin(), face.end());
    const int euler_character = int(face_node.size()) - int(ea.edges_.size()) + int(face.size(2));
    double b = 0;
    for(const auto & one_point : kcorr)
      b += one_point.second;
    cerr << "# [info] total gauss curvature " << euler_character*2<< " PI." << endl;
    cerr << "# [info] total kcorr " << b/My_PI() << " PI." << endl;
    {
      ofstream ofs("point_curvature");
      ofs << kcorr.size() << endl;
      vector<double> kcorr_vec(kcorr.size());
      for(const auto & one_point: kcorr){
          ofs << one_point.second << endl;
          kcorr_vec[one_point.first] = one_point.second;
        }
      ofstream ofs2("face_point_curvature.vtk");
      tri2vtk(ofs2, &node[0], node.size(2), &face[0],face.size(2));
      point_data(ofs2, &kcorr_vec[0], kcorr_vec.size(), "curvature");
    }

  }

  {
    size_t common_edge[2];
    double total_connection_angle = 0;
    static int vertex_num = 0;
    for(jtf::mesh::one_ring_face_at_point::p2f_type::const_iterator cit = orfap.p2f_.begin();
        cit != orfap.p2f_.end(); ++cit){
        const vector<size_t> & faces_loop = cit->second;
        assert(faces_loop.front() == faces_loop.back());
        map<size_t,double> edge_with_order;
        for(size_t fi = 0; fi < faces_loop.size()-1;++fi){
            jtf::mesh::find_common_edge(face(zjucad::matrix::colon(), faces_loop[fi]),
                                        face(zjucad::matrix::colon(), faces_loop[fi+1]),
                common_edge);
            const size_t edge_idx = ea.get_edge_idx(common_edge);
            const pair<size_t,size_t> & face_pair = ea.edge2cell_[edge_idx];
            if(face_pair.first == faces_loop[fi] && face_pair.second == faces_loop[fi+1]){
                edge_with_order[edge_idx] = 1.0;
              }else if(face_pair.first == faces_loop[fi+1] && face_pair.second == faces_loop[fi]){
                edge_with_order[edge_idx] = -1.0;
              }
          }
        double voroni_cell_area = 0;
        for(size_t fi = 0; fi < faces_loop.size()-1; ++fi){
            if(faces_loop[fi] == -1) continue;
            voroni_cell_area += face_area[faces_loop[fi]]/3.0;
          }

        double angle_defect_this_point = 0;
        auto it = angle_defect_map.find(cit->first);
        if(it != angle_defect_map.end()) {
            ++vertex_num;
            angle_defect_this_point = -1*it->second;
            angle_defect_this_point += kcorr.at(cit->first);
            total_connection_angle += angle_defect_this_point;
          }else {
            cerr << "# [error] strange can not find points in angle_defect_map"<< endl;
          }
        //#define use_ipopt
#ifndef use_ipopt
        eqn_cons->push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                              new angle_defect_sum_func(ea.edges_.size(),edge_with_order, angle_defect_this_point,100,2)));
#else
        eqn_cons->push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                              new angle_defect_sum_func(ea.edges_.size(),edge_with_order, angle_defect_this_point)));

#endif
      }
  }

  cerr << "# [info] eqn_cons number " << eqn_cons->size() << endl;

  if(angle.size() != ea.edges_.size())
    angle = zjucad::matrix::zeros<double>(ea.edges_.size(),1);
  boost::property_tree::ptree pt_opt;
  pt_opt.put("iter.value",10);
  pt_opt.put("package.value","jtf");
  pt_opt.put("alg.value","SQP");
#ifdef use_ipopt
  ipopt_solve(angle, *obj_func,eqn_cons,pt_opt);
#else
  eqn_cons->push_back(obj_func);
  shared_ptr<jtf::function::functionN1_t<double,int32_t> > all_f(
        new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD_CONS>(*eqn_cons));
  jtf::optimize(*all_f, angle, pt,nullptr, nullptr, nullptr);
#endif
}
