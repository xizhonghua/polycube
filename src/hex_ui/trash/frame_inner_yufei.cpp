#include <zjucad/ptree/ptree.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>

#include <jtflib/mesh/util.h>
#include <jtflib/function/function.h>
#include <jtflib/function/func_aux.h>
#include <jtflib/optimizer/optimizer.h>
#include <string>
#include <stack>
#include "../common/vtk.h"
#include "../tetmesh/tetmesh.h"
#include "../common/zyz.h"
#include "../hex_frame_opt/optimize_frame_field_L1.h"
#include "../spherical_harmonics/rot_cubic_f_SH.h"
#include "../hex_frame_opt/yufei_permu.h"
#include "../vol_param/solver/solver_ipopt.h"

using namespace std;
using namespace jtf::mesh;
using boost::property_tree::ptree;
using namespace zjucad::matrix;

enum frame_exp{SH,ZYZ,ROT};

template <typename val_type, frame_exp FE>
struct frame_expression_traits;

template <typename val_type>
struct frame_expression_traits<val_type, SH>
{
  typedef zjucad::matrix::matrix<val_type> frame_type;
};

template <typename val_type>
struct frame_expression_traits<val_type, ZYZ>
{
  typedef zjucad::matrix::matrix<val_type> frame_type;
};

template <typename val_type>
struct frame_expression_traits<val_type, ROT>
{
  typedef zjucad::matrix::matrix<val_type> frame_type;
};

template <typename T1, typename T2>
void frame_convert(T1 & frame_exp1, T2 & frame_exp2);

template <typename val_type>
void frame_convert(frame_expression_traits<val_type, ZYZ> & frame_exp1,
                   frame_expression_traits<val_type, ROT> & frame_exp2)
{
  if(frame_exp1.size(1) != 3 || frame_exp2.size(1) != 9
     || frame_exp1.size(2) != frame_exp2.size(2))
    throw std::invalid_argument("wrong input frame");

  size_t fi = 0;
#pragma omp parallel for private(fi)
  for(fi = 0; fi < frame_exp1.size(2); ++fi){
      zyz_angle_2_rotation_matrix1(&frame_exp1(0,fi), &frame_exp2(0,fi));
    }
}

template <typename val_type>
void frame_convert(frame_expression_traits<val_type, ROT> & frame_exp1,
                   frame_expression_traits<val_type, ZYZ> & frame_exp2)
{
  if(frame_exp2.size(1) != 3 || frame_exp1.size(1) != 9
     || frame_exp1.size(2) != frame_exp2.size(2))
    throw std::invalid_argument("wrong input frame");

  size_t fi;
#pragma omp parallel for private(fi)
  for(fi = 0; fi < frame_exp1.size(2); ++fi){
      rotation_matrix_2_zyz_angle(&frame_exp1(0,fi), &frame_exp2(0,fi));
    }
}

template <typename val_type>
void frame_convert(frame_expression_traits<val_type, SH> & frame_exp1,
                   frame_expression_traits<val_type, ZYZ> & frame_exp2)
{
  if(frame_exp2.size(1) != 3 || frame_exp1.size(1) != 9
     || frame_exp1.size(2) != frame_exp2.size(2))
    throw std::invalid_argument("wrong input frame");

  size_t fi = 0;
#pragma omp parallel for private(fi)
  for(fi = 0; fi < frame_exp1.size(2); ++fi){
      SH2zyz_lmder1(&frame_exp2(0,fi), &frame_exp1(0,fi), 200);
    }
}

template <typename val_type>
void frame_convert(frame_expression_traits<val_type, ZYZ> & frame_exp1,
                   frame_expression_traits<val_type, SH> & frame_exp2)
{
  if(frame_exp1.size(1) != 3 || frame_exp2.size(1) != 9
     || frame_exp1.size(2) != frame_exp2.size(2))
    throw std::invalid_argument("wrong input frame.");

  size_t fi = 0;
#pragma omp parallel for private(fi)
  for(fi = 0; fi < frame_exp1.size(2); ++fi){
      calc_rot_cubic_f_sh_(&frame_exp2(0,fi), &frame_exp1(0,fi));
    }
}

class permu_diff: public jtf::function::functionN1_t<double,int32_t>
{
public:
  permu_diff(size_t number, const size_t t0,const size_t t1,const double w)
    :t0_(t0), t1_(t1), w_(w), num_(number){
    if(t0_>t1_) swap(t0_,t1_);
  }
  virtual ~permu_diff(){}
  virtual size_t dim() const {return 3*num_;}
  virtual int val(const val_type *x, val_type &v) const
  {
    itr_matrix<const val_type*> zyz_x(3, num_, &x[0]);
    matrix<double> two_orth_frame(3,2);
    two_orth_frame(colon(),0) = zyz_x(colon(),t0_);
    two_orth_frame(colon(),1) = zyz_x(colon(),t1_);

    val_type d = 0;
    calc_yufei_permu_(&d, &two_orth_frame[0]);
    v += d * w_;
    return 0;
  }

  virtual int gra(const val_type *x, size_t &nnz, val_type *g, int_type *idx)
  {
    if(g == 0 && idx == 0){
        nnz = 6;
        return 0;
      }

    itr_matrix<const val_type*> zyz_x(3, num_, &x[0]);
    matrix<double> two_orth_frame(3,2), gra_x(6,1);
    two_orth_frame(colon(),0) = zyz_x(colon(),t0_);
    two_orth_frame(colon(),1) = zyz_x(colon(),t1_);
    calc_yufei_permu_jac_(&gra_x[0], &two_orth_frame[0]);
    gra_x *= w_;

    for(size_t i = 0; i < 6; ++i){
        g[i] = gra_x[i];
        idx[i] = (i<3?3*t0_+i%3:3*t1_+i%3);
      }
    return 0;
  }

  virtual int gra(const val_type *x, val_type *g)
  {
    static matrix<val_type> grax(6,1);
    static matrix<int_type> idx(6,1);
    size_t nnz = 6;
    gra(x, nnz, &grax[0], &idx[0]);
    for(size_t i = 0; i < 6; ++i){
        g[idx[i]] += grax[i];
      }
    return 0;
  }

  virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                  val_type *h, int_type *ptr, int_type *idx, val_type alpha)
  {
    if(h == 0 && ptr == 0 && idx == 0){
        nnz = 36;
        format = 1;
        return 0;
      }
    if(h == 0 && ptr != 0 && idx != 0){
        for(size_t ri = 0; ri < 6; ++ri){
            if(ri < 3) ptr[3*t0_+ri+1] = ptr[3*t0_+ri] + 6;
            else{
                ptr[3*t1_] = ptr[3*t0_+3];
                ptr[3*t1_+ri%3+1] = ptr[3*t1_+ri%3] + 6;
              }
            for(size_t ci = 0; ci < 6; ++ci){
                idx[ptr[3*(ri<3?t0_:t1_)+ri%3]+ci] = 3*(ci<3?t0_:t1_)+ci%3;
              }
          }
        return 0;
      }

    itr_matrix<const val_type*> zyz_x(3, num_, &x[0]);
    matrix<double> two_orth_frame(3,2), gra_x(6,1);
    two_orth_frame(colon(),0) = zyz_x(colon(),t0_);
    two_orth_frame(colon(),1) = zyz_x(colon(),t1_);

    matrix<double> hes(6,6);
    calc_yufei_permu_hes_(&hes[0], &two_orth_frame[0]);
    hes *= w_;
    for(size_t ri = 0; ri < 6; ++ri){
        const size_t ridx = 3*(ri<3?t0_:t1_)+ri%3;
        for(size_t ci = 0; ci < 6; ++ci){
            const size_t cidx = 3*(ci<3?t0_:t1_)+ci%3;
            jtf::function::add_to_csc(h,ptr,idx, ridx, cidx, hes(ri,ci));
          }
      }
    return 0;
  }
  virtual int hes_block(const val_type *x, val_type *h, val_type alpha)
  {
    return __LINE__;
  }
private:
  const size_t num_;
  size_t t0_;
  size_t t1_;
  const double w_;
};

class vec_fix : public jtf::function::functionN1_t<double,int32_t>
{
public:
  vec_fix(const size_t num,
          const size_t idx,
          const matrix<double> & target,
          const double w)
    :num_(num), idx_(idx), target_(target), w_(w){}
  virtual ~vec_fix(){}
  virtual size_t dim() const
  {
    return 3*num_;
  }
  virtual int val(const val_type *x, val_type &v) const
  {
    itr_matrix<const val_type*> x0(3, num_,x);
    const double vec_norm = norm(target_-x0(colon(),idx_));
    v += 0.5*w_*(vec_norm * vec_norm);
    return 0;
  }
  virtual int gra(const val_type *x, val_type *g)
  {
    itr_matrix<const val_type*> x0(3, num_,x);
    for(size_t i = 0 ; i < 3; ++i){
        g[3*idx_+i] += w_*(x0(i,idx_)-target_[i]);
      }
    return 0;
  }
  virtual int gra(const val_type *x, size_t &nnz, val_type * g, int_type *idx)
  {
    if(g == 0 && idx == 0){
        nnz = target_.size();
        return 0;
      }
    itr_matrix<const val_type*> x0(3, num_,x);
    for(size_t i = 0 ; i< target_.size(); ++i){
        g[i] = w_*(x0(i,idx_) - target_[i]);
        idx[i] = 3*idx_+i;
      }
    return 0;
  }
  virtual int hes(const val_type *x, size_t &nnz, size_t &format,
                  val_type *h, int_type *ptr, int_type *idx, val_type alpha)
  {
    if(h == 0 && ptr == 0 && idx == 0){
        nnz = target_.size();
        format = 1;
        return 0;
      }
    if(h == 0 && ptr != 0 && idx != 0){
        for(size_t i = 0 ;i < target_.size(); ++i){
            ptr[3*idx_+i+1] = ptr[3*idx_+i]+1;
            idx[ptr[3*idx_+i]] = 3*idx_+i;
          }
        return 0;
      }
    for(size_t i = 0 ; i < 3; ++i){
        jtf::function::add_to_csc(h, ptr, idx, 3*idx_+i, 3*idx_+i, w_);
      }

    return 0;
  }

  virtual int hes_block(const val_type *x, val_type *h, val_type alpha)
  {
    return __LINE__;
  }
private:
  const size_t num_;
  const size_t idx_;
   matrix<double> target_;
  const double w_;
};

class frame_optimizer
{
public:
  frame_optimizer(const jtf::tet_mesh &tm):tm_(tm){}
  void opt(boost::property_tree::ptree &pt,
           zjucad::matrix::matrix<double> & frame);
  void init_frame(boost::property_tree::ptree &pt,
                  zjucad::matrix::matrix<double> & frame);
  void opt_frame(boost::property_tree::ptree &pt,
                 zjucad::matrix::matrix<double> & frame);
private:

  ///
  /// \brief load_cross_field load surface cross field.
  /// \param ptree pt
  /// \param frame
  /// \param visited_tets
  /// \return
  ///
  int load_cross_field(boost::property_tree::ptree &pt,
                       zjucad::matrix::matrix<double> & frame,
                       std::vector<bool> & visited_tets);

  int propagate_field0(zjucad::matrix::matrix<double> & frame,
                       std::vector<bool> & visited_tets);

  int propagate_field1(zjucad::matrix::matrix<double> & frame,
                       std::vector<bool> & visited_tets);

  ///
  /// \brief init_frame_yufei init frame field by fix surface frames with cross field
  ///                         which is used by yufei.
  /// \param pt
  /// \param frame
  ///
  void init_frame_yufei(boost::property_tree::ptree &pt,
                        zjucad::matrix::matrix<double> &frame)
  {
    // yufei use a surface cross frame field as input, and fix surface tet
    // frame with w as outside normal, and u,v as cross frame field
    frame.resize(9,tm_.tetmesh_.mesh_.size(2));
    vector<bool> visited_tets(frame.size(2),false);
    load_cross_field(pt, frame, visited_tets);
    propagate_field1(frame, visited_tets);
  }

  ///
  /// \brief opt_frame_yufei use yufei's strategy to opt frame
  /// \param pt
  /// \param frame  input init frame and output opt frame
  ///
  void opt_frame_yufei(boost::property_tree::ptree &pt,
                       zjucad::matrix::matrix<double> & frame);

  const jtf::tet_mesh &tm_;
  zjucad::matrix::matrix<double> frame_;
};


int frame_optimizer::load_cross_field(boost::property_tree::ptree & pt,
                                      zjucad::matrix::matrix<double> & frame,
                                      std::vector<bool> & visited_tets)
{
  const string surface_cross_field = pt.get<string>("input/surface_cross_field.value");
  const string obj_file = pt.get<string>("input/obj.value");
  const string s2v_file = pt.get<string>("input/s2v.value");

  ifstream ifs_field(surface_cross_field);
  if(ifs_field.fail()){
      return __LINE__;
    }

  jtf::mesh::meshes obj_m;
  if(jtf::mesh::load_obj(obj_file.c_str(), obj_m.mesh_, obj_m.node_)){
      std::cerr << "# [error] can not open obj file." << std::endl;
      return __LINE__;
    }

  zjucad::matrix::matrix<int32_t> s2v;
  if(jtf::mesh::read_matrix(s2v_file.c_str(), s2v))
    return __LINE__;

  matrix<double> surface_frame;
  {
    size_t num;
    ifs_field >> num;
    surface_frame.resize(6,num);
    for(size_t i = 0; i < num; ++i){
        for(size_t j = 0 ; j < 6; ++j){
            ifs_field >> surface_frame(j,i);
          }
      }
  }

  // after load surface frame, the next step is to generate 3-cross field for tets
  for(size_t fi = 0 ; fi < obj_m.mesh_.size(2); ++fi){
      const size_t face_idx = tm_.fa_->get_face_idx(s2v[obj_m.mesh_(0,fi)],
          s2v[obj_m.mesh_(1,fi)], s2v[obj_m.mesh_(2,fi)]);
      assert(face_idx != -1);
      const pair<size_t,size_t> & tet_pair = tm_.fa_->face2tet_[face_idx];
      assert(tm_.fa_->is_outside_face(tet_pair));
      const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);
      itr_matrix<double*> frame_itr(3,3,&frame[9*tet_idx]);
      frame_itr(colon(),0) = surface_frame(colon(0,2),fi);
      frame_itr(colon(),1) = surface_frame(colon(3,5),fi);
      frame_itr(colon(),2) = cross(frame_itr(colon(),0),frame_itr(colon(),1));
      visited_tets[tet_idx] = true;
    }

  return 0;
}

int frame_optimizer::propagate_field0(zjucad::matrix::matrix<double> & frame,
                                      std::vector<bool> & visited_tets)
{
  vector<set<pair<double,size_t> > > tet2tet_with_dis(tm_.tetmesh_.mesh_.size(2));
  {
    matrix<double> tet_center = zeros<double>(3, tm_.tetmesh_.mesh_.size(2));
    for(size_t ti = 0; ti < tm_.tetmesh_.mesh_.size(2); ++ti){
        for(size_t pi = 0; pi < tm_.tetmesh_.mesh_.size(1); ++pi){
            tet_center(colon(),ti) += tm_.tetmesh_.node_(colon(),tm_.tetmesh_.mesh_(pi,ti));
          }
      }
    tet_center/=4.0;
    for(size_t fi = 0; fi < tm_.fa_->face2tet_.size(); ++fi){
        const pair<size_t,size_t> & tet_pair = tm_.fa_->face2tet_[fi];
        if(tm_.fa_->is_outside_face(tet_pair)) continue;
        const double distance = norm(tet_center(colon(),tet_pair.first) - tet_center(colon(), tet_pair.second));
        tet2tet_with_dis[tet_pair.first].insert(make_pair(distance, tet_pair.second));
        tet2tet_with_dis[tet_pair.second].insert(make_pair(distance, tet_pair.first));
      }
  }

  std::vector<bool>::iterator vit = find(visited_tets.begin(), visited_tets.end(), false);
  while(vit != visited_tets.end()){
      for(size_t ti = 0; ti < visited_tets.size(); ++ti){
          if(visited_tets[ti] == false) continue;
          const set<pair<double,size_t> > & adj_tets = tet2tet_with_dis[ti];
          for(const auto &one_tet : adj_tets){
              if(visited_tets[one_tet.second]) continue;
              frame(colon(),one_tet.second) = frame(colon(), ti);
              visited_tets[one_tet.second] = true;
            }
        }
      vit = find(visited_tets.begin(), visited_tets.end(), false);
    }

  return 0;
}
int frame_optimizer::propagate_field1(zjucad::matrix::matrix<double> & frame,
                                      std::vector<bool> & visited_tets)
{
  matrix<size_t> nearest_tet(visited_tets.size(),1);
  set<size_t> surface_tets;
  for(size_t ti = 0; ti < visited_tets.size(); ++ti){
      if(visited_tets[ti] == true) surface_tets.insert(ti);
    }

  matrix<double> tet_center = zeros<double>(3, tm_.tetmesh_.mesh_.size(2));
  for(size_t ti = 0; ti < tm_.tetmesh_.mesh_.size(2); ++ti){
      for(size_t pi = 0; pi < tm_.tetmesh_.mesh_.size(1); ++pi){
          tet_center(colon(),ti) += tm_.tetmesh_.node_(colon(),tm_.tetmesh_.mesh_(pi,ti));
        }
    }
  tet_center/=4.0;

  size_t ti = 0;
#pragma omp parallel for private(ti)
  for(ti = 0;ti < visited_tets.size(); ++ti){
      if(find(surface_tets.begin(), surface_tets.end(), ti) == surface_tets.end()){
          double dis = std::numeric_limits<double>::infinity();
          size_t target_tet = -1;
          for(const auto & one_surf_tet : surface_tets){
              double dis_temp = norm(tet_center(colon(), ti) - tet_center(colon(), one_surf_tet));
              if(dis > dis_temp){
                  dis = dis_temp;
                  target_tet = one_surf_tet;
                }
            }
          frame(colon(),ti) = frame(colon(), target_tet);
        }
    }

  return 0;
}

void frame_optimizer::opt(boost::property_tree::ptree &pt,
                          zjucad::matrix::matrix<double> & frame)
{
  init_frame(pt,frame);
  opt_frame(pt,frame);
}

void frame_optimizer::init_frame(boost::property_tree::ptree &pt,
                                 zjucad::matrix::matrix<double> &frame)
{
  pt.put("init_strategy.desc", "yufei");
  const string strategy = pt.get<string>("input/init_strategy.value");
  if(strategy == "yufei"){
      init_frame_yufei(pt,frame);
    }else if(strategy == "io"){
      matrix<double> zyz = zeros<double>(3,tm_.tetmesh_.mesh_.size(2));
      if(read_matrix(pt.get<string>("input/init_zyz.value").c_str(), zyz))
        throw std::invalid_argument("invalid init zyz");
      if(zyz.size(2) != tm_.tetmesh_.mesh_.size(2)){
          throw std::invalid_argument("wrong zyz size.");
        }
      cerr << "# [info] use init_zyz." << endl;
      frame.resize(9,zyz.size(2));
      for(size_t i = 0; i < tm_.tetmesh_.mesh_.size(2); ++i){
          zyz_angle_2_rotation_matrix1(&zyz(0,i), &frame(0,i));
        }
    }else{
      throw std::invalid_argument("only support yufei init_strategy");
    }
}

void frame_optimizer::opt_frame(boost::property_tree::ptree &pt,
                                zjucad::matrix::matrix<double> &frame)
{
  pt.put("opt_strategy.desc", "yufei");
  const string strategy = pt.get<string>("input/opt_strategy.value");
  if(strategy == "yufei"){
      opt_frame_yufei(pt,frame);
    }else{
      throw std::invalid_argument("only support yufei opt_strategy");
    }
}

void frame_optimizer::opt_frame_yufei(boost::property_tree::ptree &pt,
                                      zjucad::matrix::matrix<double> & frame)
{
  zjucad::matrix::matrix<double> zyz(3,frame.size(2));
  size_t ti;
#pragma omp parallel for private(ti)
  for(ti = 0; ti < frame.size(2); ++ti){
      rotation_matrix_2_zyz_angle(&frame(0,ti), &zyz(0,ti), 0);
    }

  shared_ptr<vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > > func_vec_ptr(
        new vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > >);
  // inner smooth
  for(const auto &one_edge : tm_.ortae_.e2t_){
      const vector<size_t> & adj_tets = one_edge.second;
      for(size_t i = 0; i < adj_tets.size(); ++i){
          if(adj_tets[i] == -1) continue; // meet outside tets
          for(size_t j = i + 1; j < adj_tets.size(); ++j){
              if(adj_tets[j] == -1 || adj_tets[j] == adj_tets[i]) continue;
              double w = 1.0;
              func_vec_ptr->push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                                        new permu_diff(frame.size(2), adj_tets[i], adj_tets[j],w)));
            }
        }
    }

  // surface fix
  for(size_t fi = 0; fi < tm_.outside_face_idx_.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = tm_.fa_->face2tet_[tm_.outside_face_idx_[fi]];
      assert(tm_.fa_->is_outside_face(tet_pair));
      const size_t surface_tet_idx = tet_pair.first==-1?tet_pair.second:tet_pair.first;
      const double w = 20.0;
      func_vec_ptr->push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                                new vec_fix(frame.size(2), surface_tet_idx,
                                            zyz(colon(), surface_tet_idx), w)));
    }

  shared_ptr<jtf::function::functionN1_t<double,int32_t> > all_func(
        new jtf::function::sum_function<double,int32_t, jtf::function::SMART_STD>(*func_vec_ptr));

  jtf::optimize(*all_func, zyz, pt, nullptr, nullptr, nullptr);
}


int frame_yufei(ptree &pt)
{
  pt.put("input/tet.desc", "input tetmesh.");
  pt.put("input/init_zyz.desc", "input init zyz.");
  pt.put("output/zyz.desc" , "output zyz");

  jtf::tet_mesh tm(pt.get<string>("input/tet.value").c_str());

  matrix<double> zyz = zeros<double>(3,tm.tetmesh_.mesh_.size(2));
  matrix<double> frame_init(9, zyz.size(2));

  frame_optimizer fo(tm);

  fo.opt(pt,frame_init);

  for(size_t i = 0; i < tm.tetmesh_.mesh_.size(2); ++i){
      rotation_matrix_2_zyz_angle(&frame_init(0,i), &zyz(0,i),0);
    }

  if(write_matrix(pt.get<string>("output/zyz.value").c_str(), zyz))
    return __LINE__;
  return 0;
}

