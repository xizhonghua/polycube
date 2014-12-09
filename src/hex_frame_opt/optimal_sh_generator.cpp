#include "optimal_sh_generator.h"
#include "angle_defect.h"
#include "../numeric/util.h"
#include "../common/vtk.h"
//#include "../numeric/func2opt.h"
#include "../spherical_harmonics/rot_cubic_f_SH.h"
#include <hjlib/sparse/operation.h>
#include <hjlib/sparse/sparse.h>
#include <jtflib/math/math.h>
#include <jtflib/function/func_aux.h>
#include <hjlib/math_func/math_func.h>
#include <hjlib/math_func/operation.h>
#include <hjlib/math_func/func_aux.h>
#include <hjlib/math_func/coo.h>
#include <zjucad/matrix/io.h>
#include <jtflib/optimizer/optimizer.h>
#include "../common/solver_ipopt.h"
#include "frame_function.h"
#include "../common/zyz.h"
#include "relation.h"
#include "../common/util.h"
#include "connection_generator.h"
#include "../vol_param/descriptor/func_terms/arap.h"
#include "../common/transition.h"
#include "../common/visualize_tool.h"
#include "smooth_L1.h"
#include "../numeric/csc_filler.h"
#include "../mesh_func/frame_func.h"

using namespace std;
using namespace zjucad::matrix;

template <typename val_type, typename int_type>
void generate_eye_csc(hj::sparse::csc<val_type, int_type> & A,
                      const size_t dim)
{
  A.resize(dim,dim,dim);
  A.val() = zjucad::matrix::ones<val_type>(dim,1);
  A.ptr()(colon(1,dim),0) = zjucad::matrix::colon(1,dim);
  A.idx() = zjucad::matrix::colon(0,dim-1);
}

template <typename val_type>
val_type get_opposite_angle(const pair<size_t,size_t> & one_edge,
                            const zjucad::matrix::matrix<size_t> & face,
                            const zjucad::matrix::matrix<val_type> & node)
{
  const size_t other_point = std::accumulate(face.begin(), face.end(), static_cast<size_t>(0))
      - one_edge.first - one_edge.second;
  matrix<double> dir1 = node(colon(), one_edge.first) - node(colon(), other_point);
  matrix<double> dir2 = node(colon(), one_edge.second) - node(colon(), other_point);
  const double len1 = norm(dir1);
  const double len2 = norm(dir2);
  if(len1 > 1e-12){ dir1/=len1;}
  if(len2 > 1e-12){ dir2/=len2;}
  return jtf::math::safe_acos(dot(dir1,dir2));
}

int optimal_sh_generator::zyz2frame(const zjucad::matrix::matrix<double> & zyz,
                                    zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & frame)
{
  if(zyz.size(1) != 3){
      cerr << "# [error] convert zyz to frame fail: wrong zyz." << endl;
      return __LINE__;
    }
  if(frame.size() != zyz.size(2)) frame.resize(zyz.size(2),1);

  {
    size_t ti = 0;
#pragma omp parallel for private(ti)
    for(ti = 0;ti < zyz.size(2); ++ti){
        if(frame[ti].size() == 0) frame[ti].resize(3,3);
        zyz_angle_2_rotation_matrix1(&zyz(0,ti), &frame[ti][0]);
      }
  }
  return 0;
}

void optimal_sh_generator::project2zyz(const zjucad::matrix::matrix<double> &sh,
                                       zjucad::matrix::matrix<double> &zyz)const
{
  itr_matrix<const double*> sh_m(9, sh.size()/9, &sh[0]);
  if(zyz.size(1) != 3 || zyz.size(2) != sh_m.size(2))
    zyz.resize(3,sh_m.size(2));
  zyz *= 0;
  size_t i = 0;
#pragma omp parallel for private(i)
  for(i = 0; i < sh_m.size(2); ++i){
      SH2zyz(&zyz(0,i), &sh_m(0,i),100);
    }
}

void optimal_sh_generator::init_connection(boost::property_tree::ptree &pt)
{
  rotation_angle_ = zeros<double>(tm_.ea_outside_->edges_.size(),1);

  pt.put("angle_defect_strategy.desc","manually/even/null");
  const string strategy = pt.get<string>("angle_defect_strategy.value");
  if(strategy == "manually"){
      cerr << "# [info] use manually angle defect." << endl;
      std::unique_ptr<angle_defect> ad(new manually_set_angle_defect);
      ad->opt(*tm_.ea_outside_,tm_.outside_face_,tm_.tetmesh_.node_, rotation_angle_, pt);
    }else if(strategy == "even"){
      cerr << "# [info] use even angle_defect " << endl;
      std::unique_ptr<angle_defect> ad(new even_angle_defect);
      ad->opt(*tm_.ea_outside_,tm_.outside_face_,tm_.tetmesh_.node_, rotation_angle_, pt);
    }else if(strategy == "geometry"){
      cerr << "# [info] use geometry aware angle_defect " << endl;
      std::shared_ptr<angle_defect> ad(new geometry_aware_angle_defect);
      ad->opt(*tm_.ea_outside_,tm_.outside_face_,tm_.tetmesh_.node_, rotation_angle_, pt);
      pt.put("input/need_smooth_normal.desc","smooth surface normal according to gaussian radius, [Y/N]");
      if(zjucad::has("input/need_smooth_normal.value",pt)){
          string need_smooth_normal = pt.get<string>("input/need_smooth_normal.value","y");
          if(need_smooth_normal == "YES" || need_smooth_normal == "yes" ||
             need_smooth_normal == "Y" || need_smooth_normal == "y"){
              cerr << "# [info] smooth surface normal according to geometry aware filter radius." << endl;
              std::shared_ptr<geometry_aware_angle_defect> ad_g =dynamic_pointer_cast<geometry_aware_angle_defect>(ad);
              smooth_surface_normal(tm_.outside_face_normal_, tm_.outside_face_, tm_.tetmesh_.node_, ad_g, pt, reference_dir_);
            }
        }
    }else if(strategy == "null"){
      //levi-civitia connection, other connection can be built
      cerr << "# [info] use null angle_defect." << endl;
    }

  rij_ = zeros<double>(3,tm_.ea_outside_->edges_.size());
  for(size_t ei = 0; ei < tm_.ea_outside_->edges_.size(); ++ei){
      const std::pair<size_t,size_t> & face_pair = tm_.ea_outside_->edge2cell_[ei];
      if(tm_.ea_outside_->is_boundary_edge(face_pair)) continue;
      double angle = get_rij_angle(face_pair,tm_.outside_face_, tm_.tetmesh_.node_, tm_.outside_face_normal_);
      //double angle2 = get_rij_angle2(face_pair,tm_.outside_face_, tm_.tetmesh_.node_, tm_.outside_face_normal_);

      rij_(0,ei) = -1*(angle - rotation_angle_[ei]); // -1 means this angle should be minused during measuring smoothness
    }
}

#define debug
#ifdef debug
static void cal_ref_align_rot(const zjucad::matrix::matrix<size_t> & one_face,
                              const zjucad::matrix::matrix<double> & node,
                              const zjucad::matrix::matrix<double> &rz2n,
                              zjucad::matrix::matrix<double> &rot_ref)
{
  zjucad::matrix::matrix<double> x = rz2n(colon(),0);
  zjucad::matrix::matrix<double> dir = node(colon(), one_face[1])
      - node(colon(), one_face[0]);
  const double len_x = zjucad::matrix::norm(x);
  const double len_dir = zjucad::matrix::norm(dir);
  assert(len_x > 1e-18 && len_dir > 1e-18);
  x/=len_x;
  dir /= len_dir;
  double angle = jtf::math::safe_acos(dot(x,dir));

  const matrix<double> normal = cross(dir, node(colon(), one_face[2]) - node(colon(), one_face[1]));
  if(dot(normal, cross(x,dir))<0)
    angle *= -1;
  from_angle_to_rotation_matrix(angle, normal,rot_ref);
}
#endif

void optimal_sh_generator::smooth_surface_normal(
    zjucad::matrix::matrix<double> & face_normal,
    const zjucad::matrix::matrix<size_t> & face,
    const zjucad::matrix::matrix<double> & node,
    std::shared_ptr<geometry_aware_angle_defect> ga,
    boost::property_tree::ptree &pt,
    zjucad::matrix::matrix<double> & reference_dir)
{
  double radius = pt.get<double>("weight/gaussian_radius.value");
  if(fabs(radius) < 1e-8) return;

  matrix<double> bb(3,2);
  calc_bounding_box(tm_.tetmesh_.node_, &bb[0]);
  bb(colon(),0) -= bb(colon(),1);
  const double max_len = zjucad::matrix::max(zjucad::matrix::fabs(bb(colon(),0)));
  radius *= 2*max_len;

  matrix<double> point_normal(3, node.size(2));
  jtf::mesh::cal_point_normal(face, node, point_normal);

  vector<vector<std::pair<size_t,double> > > distance;
  std::shared_ptr<const jtf::mesh::one_ring_point_at_point> orpap(
        jtf::mesh::one_ring_point_at_point::create(face));
  if(!orpap.get())
    throw std::logic_error("can not build one_ring_point_at_point.");

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

  ga->cal_geodesic_distance(distance, node, orpap, radius, &feature_points);

  matrix<double> this_point_normal = zeros<double>(3,1);
  for(size_t pi = 0; pi < distance.size(); ++pi){
      const vector<std::pair<size_t,double> > & connect_points = distance[pi];
      if(connect_points.empty()) continue;
      this_point_normal *= 0;
      for(size_t pj = 0; pj < connect_points.size(); ++pj){
          this_point_normal += point_normal(colon(), connect_points[pj].first)/connect_points[pj].second;
        }
      this_point_normal /= norm(this_point_normal);
      point_normal(colon(),pi) = this_point_normal;
    }

  // convert face normal to point normal
  zjucad::matrix::matrix<double> old_face_normal = face_normal;
  for(size_t fi = 0; fi < face.size(2); ++fi){
      face_normal(colon(),fi) =
          (point_normal(colon(), face(0,fi)) +
           point_normal(colon(), face(1,fi)) +
           point_normal(colon(), face(2,fi)))/3.0;
    }
  assert(reference_dir.size(1) == 3 && reference_dir.size(2) == old_face_normal.size(2));

  zjucad::matrix::matrix<double> rot(3,3);
  for(size_t fi = 0; fi < old_face_normal.size(2); ++fi){
      zjucad::matrix::matrix<double> cross_edge = cross(old_face_normal(colon(),fi), face_normal(colon(),fi));
      const double len = norm(cross_edge);
      if(len < 1e-8) continue;
      cross_edge /= len;
      const double angle = jtf::math::safe_asin(len/(norm(old_face_normal(colon(),fi))*norm(face_normal(colon(),fi))));
      from_angle_to_rotation_matrix(angle, cross_edge, rot);
      reference_dir(colon(),fi) = rot*reference_dir(colon(),fi);
    }
}

void optimal_sh_generator::load_edge_dense_field(map<pair<size_t,size_t>, double> & edge_dense,
                                                 boost::property_tree::ptree &pt)const
{
  ifstream ifs(pt.get<string>("input/edge_dense_field.value").c_str());
  if(ifs.fail())
    throw std::invalid_argument("can not open edge dense field file.");

  size_t edge_num = 0;
  ifs >> edge_num;
  if(edge_num != tm_.ortae_.e2t_.size())
    throw std::invalid_argument("edge number is not compatible with tet mesh edges");

  pair<size_t,size_t> one_edge;
  double d;
  for(size_t ei = 0; ei < edge_num; ++ei){
      ifs >> one_edge.first >> one_edge.second >> d;
      if(one_edge.first > one_edge.second) swap(one_edge.first, one_edge.second);
      edge_dense[one_edge] = d;
    }
}

void optimal_sh_generator::load_surface_charts(
    const char * chart_file,
    std::map<std::pair<size_t,bool>, std::set<size_t> > & charts)
{
  ifstream ifs(chart_file);
  if(ifs.fail()){
      throw std::invalid_argument("can not open chart file.");
    }
  size_t chart_num, face_num, face_idx_temp;
  int use_real_normal;
  ifs >> chart_num;
  set<size_t> face_idx;
  for(size_t ci = 0; ci < chart_num; ++ci){
      ifs >> face_num >> use_real_normal; // use_real_normal : 1 true; 0 false
      face_idx.clear();
      for(size_t fi = 0; fi < face_num; ++fi){
          ifs >> face_idx_temp;
          face_idx.insert(face_idx_temp);
        }
      if(use_real_normal){
          charts[make_pair(face_num, true)] = face_idx;
        }else{
          charts[make_pair(face_num, false)] = face_idx;
        }
    }
}

void optimal_sh_generator::load_guiding_field(matrix<double> &gv,
                                              boost::property_tree::ptree &pt)
{
  ifstream ifs(pt.get<string>("input/guiding_field.value").c_str());
  if(ifs.fail()){
      throw std::invalid_argument("can not open guiding field file.");
    }
  size_t face_num ;
  ifs >> face_num;
  if(face_num != gv.size(2)){
      gv.resize(3,face_num);
      cerr << "# [warning] input guding field number is not same as premalloc memory." << endl;
    }
  double trash;
  for(size_t fi = 0; fi < face_num; ++fi){
      ifs >> gv(0,fi) >> gv(1,fi) >> gv(2,fi);
      ifs >> trash >> trash >> trash;
    }
}

void optimal_sh_generator::init(boost::property_tree::ptree &pt)
{
  edge_area_weight_ = zjucad::matrix::matrix<double>(tm_.ea_outside_->edges_.size(),1);
  for(size_t ei = 0; ei < tm_.ea_outside_->edge2cell_.size(); ++ei){
      const pair<size_t,size_t> & face_pair = tm_.ea_outside_->edge2cell_[ei];
      if(tm_.ea_outside_->is_boundary_edge(face_pair)) continue;
      edge_area_weight_[ei] = (std::fabs(tm_.outside_face_area_[face_pair.first])
          +std::fabs(tm_.outside_face_area_[face_pair.second]))/3.0;
    }

  local2global_rot_ = zjucad::matrix::matrix<double>(3,tm_.outside_face_.size(2));

  reference_dir_ = zjucad::matrix::matrix<double>(3, tm_.outside_face_.size(2));

  for(size_t fi = 0; fi < tm_.outside_face_.size(2); ++fi){
      reference_dir_(colon(),fi) = tm_.tetmesh_.node_(colon(), tm_.outside_face_(1,fi))
          - tm_.tetmesh_.node_(colon(), tm_.outside_face_(0,fi));
      reference_dir_(colon(),fi) /= norm(reference_dir_(colon(),fi));
    }
  init_connection(pt);

  matrix<double> rot(3,3);
  for(size_t fi = 0; fi < tm_.outside_face_.size(2); ++fi){
      rot(colon(),0) = reference_dir_(colon(),fi);
      rot(colon(),2) = tm_.outside_face_normal_(colon(),fi);
      rot(colon(),1) = cross(rot(colon(),2), rot(colon(),0));
      rot(colon(),1) /= norm(rot(colon(),1));
      rotation_matrix_2_zyz_angle(&rot[0], &local2global_rot_(0,fi),nullptr);
    }

  if(zjucad::has("input/guiding_field.value",pt)){
      matrix<double> gf(3, tm_.outside_face_.size(2));
      load_guiding_field(gf,pt);
      surface_target_idx_on_boundary_.resize(tm_.outside_face_idx_.size(),1);
      surface_target_dir_angle_.resize(tm_.outside_face_idx_.size(),1);
      for(size_t fi = 0; fi < tm_.outside_face_.size(2); ++fi){
          surface_target_idx_on_boundary_[fi] = fi;
          gf(colon(),fi) /= norm(gf(colon(),fi));
          double angle = jtf::math::safe_acos(dot(gf(colon(),fi), reference_dir_(colon(),fi)));
          if(dot(tm_.outside_face_normal_(colon(),fi), cross(reference_dir_(colon(),fi), gf(colon(),fi))) < 0)
            angle *= -1;
          surface_target_dir_angle_[fi] = 4.*angle;
        }
      surface_target_idx_ = tm_.outside_face_idx_;
    }
  if(zjucad::has("input/edge_dense_field.value",pt)){
      load_edge_dense_field(edge_dense_field_, pt);
    }
  if(zjucad::has("input/feature_line.value",pt) &&
     zjucad::has("input/s2v.value",pt)){
      if(jtf::mesh::load_feature_line(pt.get<string>("input/feature_line.value").c_str(),
                                      pt.get<string>("input/s2v.value").c_str(), feature_lines_))
        throw std::invalid_argument("can not load feature line file.");
    }
  if(zjucad::has("input/surface_chart.value",pt)){
      map<pair<size_t,bool>, set<size_t> > charts; // chart: <idx,use_real_normal>, <face_idx_in_this_chart>
      load_surface_charts(pt.get<string>("input/surface_chart.value").c_str(),charts);
      matrix<double> normal(3,1);
      for(const auto & one_chart : charts){
          const pair<size_t,bool> & chart_idx = one_chart.first;
          if(chart_idx.second == true) continue; // use real normal
          normal *= 0;
          const set<size_t>&one_chart_face_idx = one_chart.second;\
          for(const auto & one_face : one_chart_face_idx){
              normal += tm_.outside_face_normal_(colon(), one_face)*fabs(tm_.outside_face_area_[one_face]);
            }
          normal /= norm(normal);
          for(const auto & one_face : one_chart_face_idx){
              tm_.outside_face_normal_(colon(), one_face) = normal;
            }
        }
      cerr << "# [info] use surface chart." << endl;

      {
        matrix<size_t> face_tag(tm_.outside_face_.size(2),1);
        for(const auto & one_chart : charts){
            const set<size_t> & one_chart_set = one_chart.second;
            for(const auto & one_face : one_chart_set){
                face_tag[one_face] = one_chart.first.first;
              }
          }
        ofstream ofs("chart.vtk");
        tri2vtk(ofs, &tm_.tetmesh_.node_[0], tm_.tetmesh_.node_.size(2), &tm_.outside_face_[0], tm_.outside_face_.size(2));
        cell_data(ofs, &face_tag[0], face_tag.size(), "chart");
      }
    }
  if(zjucad::has("input/normal.value",pt)){
      zjucad::matrix::matrix<double> normal;
      if(jtf::mesh::read_matrix(pt.get<string>("input/normal.value").c_str(), normal))
        throw std::invalid_argument("can not open normal file");
      if(normal.size(2) != tm_.outside_face_normal_.size(2))
        throw std::invalid_argument("wrong normal size.");
      tm_.outside_face_normal_ = normal;
      cerr << "# [info] use given normal." << endl;
    }
}


template <typename val_type, typename int_type>
class inner_smooth_func_with_connection_sh : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  inner_smooth_func_with_connection_sh(const size_t tet_num,
                                       const size_t i, const size_t j,  const double w,
                                       const zjucad::matrix::matrix<double> & connection)
    :i_(i), j_(j), tet_num_(tet_num), w_(w), connection_(connection){
  }
  virtual ~inner_smooth_func_with_connection_sh(){}
  virtual size_t nx() const{
    return 9 * tet_num_;
  }
  virtual size_t nf() const{
    return 9;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type *> x0(9,tet_num_,x);
        matrix<val_type> diff = x0(colon(), i_) - connection_ * x0(colon(),j_);
        for(int_type i = 0; i < 9; ++i){
            int_type c[1] = {i};
            cv[c] += w_*diff[i];
          }
      }
    if(k == 1){
        for(int_type i = 0; i < 9; ++i){
            int_type c1[2] = {i,9*i_+i};
            cv[c1] += w_;
            for(int_type j = 0; j < 9; ++j){
                int_type c2[2] = {i,9*j_+j};
                cv[c2] += -w_*connection_(i,j);
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
            int_type c1[2] = {i,9*i_+i};
            l2g.add(cs,c1);
            for(int_type j = 0; j < 9; ++j){
                int_type c2[2] = {i,9*j_+j};
                l2g.add(cs,c2);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 10*9;
  }
private:
  const int_type i_;
  const int_type j_;
  const int_type tet_num_;
  const double w_;
  const zjucad::matrix::matrix<double> connection_;
};


template <typename val_type, typename int_type>
class inner_smooth_func_with_connection_zyz : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  inner_smooth_func_with_connection_zyz(const size_t tet_num,
                                        const size_t i, const size_t j,  const double w,
                                        const zjucad::matrix::matrix<double> & connection)
    :i_(i), j_(j), tet_num_(tet_num), w_(w), connection_(connection){
  }
  virtual ~inner_smooth_func_with_connection_zyz(){}
  virtual size_t nx() const{
    return 3 * tet_num_;
  }
  virtual size_t nf() const{
    return 9;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type *> x0(3,tet_num_,x);
        matrix<val_type> sh0(9,1), sh1(9,1);
        calc_rot_cubic_f_sh_(&sh0[0], &x0(0,i_));
        calc_rot_cubic_f_sh_(&sh1[0], &x0(0,j_));

        matrix<val_type> diff = sh0 - connection_ * sh1;
        for(int_type i = 0; i < 9; ++i){
            int_type c[1] = {i};
            cv[c] += w_*diff[i];
          }
      }
    if(k == 1){
        itr_matrix<const val_type*> x0(3,tet_num_, x);
        matrix<val_type> jac_sh0(9,3), jac_sh1(9,3);
        calc_jac_rot_cubic_f_sh_(&jac_sh0[0], &x0(0,i_));
        calc_jac_rot_cubic_f_sh_(&jac_sh1[0], &x0(0,j_));

        matrix<val_type> jac1 = connection_*jac_sh1;

        for(int_type i = 0; i < 9; ++i){
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i, 3*i_+j};
                cv[c] += w_*jac_sh0(i,j);
                int_type c2[2] = {i, 3*j_+j};
                cv[c2] += -w_*jac1(i,j);
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
                int_type c[2] = {i, 3*i_+j};
                l2g.add(cs,c);
                int_type c2[2] = {i, 3*j_+j};
                l2g.add(cs,c2);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 6*9;
  }
private:
  const int_type i_;
  const int_type j_;
  const int_type tet_num_;
  const double w_;
  const zjucad::matrix::matrix<double> connection_;
};



shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_inner_smooth_func_with_connection3(const zjucad::matrix::matrix<size_t> &tet_mesh,
                                         const zjucad::matrix::matrix<double> &node,
                                         const jtf::mesh::face2tet_adjacent &fa,
                                         const jtf::mesh::one_ring_tet_at_edge &ortae,
                                         const zjucad::matrix::matrix<double> &vol,
                                         const double w,
                                         const map<pair<size_t,size_t>,zjucad::matrix::matrix<double> > & connection_matrix,
                                         const string strategy = "face",
                                         const string opt_type = "init")
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  zjucad::matrix::matrix<double> connection(3,3), connection_zyz(3,1), connection_sh(9,9);

  const double total_v = std::accumulate(zjucad::matrix::fabs(vol).begin(),
                                         zjucad::matrix::fabs(vol).end(), 0.0);

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  for(size_t fi = 0; fi < fa.face2tet_.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[fi];
      if(fa.is_outside_face(tet_pair)) continue;

      auto it  = connection_matrix.find(tet_pair);
      if(it != connection_matrix.end()) connection = it->second;
      else connection = eye<double>(3);

      rotation_matrix_2_zyz_angle(&connection[0], &connection_zyz[0],0);
      calc_rot_cubic_f_sh_mat_(&connection_sh[0], &connection_zyz[0]);
      if(strategy == "face"){
          if(opt_type == "init")
            all_func->push_back(math_func_ptr(
                                  new inner_smooth_func_with_connection_sh<double,int32_t>(
                                    tet_mesh.size(2), tet_pair.first, tet_pair.second,
                                    sqrt(w*(fabs(vol[tet_pair.first])+fabs(vol[tet_pair.second]))/(4.0*total_v))
                                ,connection_sh
                                )));
          else if(opt_type == "zyz")
            all_func->push_back(math_func_ptr(
                                  new inner_smooth_func_with_connection_zyz<double,int32_t>(
                                    tet_mesh.size(2), tet_pair.first, tet_pair.second,
                                    sqrt(w*(fabs(vol[tet_pair.first])+fabs(vol[tet_pair.second]))/(4.0*total_v))
                                ,connection_sh
                                )));
          else throw std::invalid_argument("# [error] can not recognize opt type");
        }
    }

  cerr << "# [info] inner smooth func with_connection " << all_func->size() << endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}



shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_inner_smooth_func_L1(const zjucad::matrix::matrix<size_t> &tet_mesh,
                           const zjucad::matrix::matrix<double> &node,
                           const jtf::mesh::face2tet_adjacent &fa,
                           const jtf::mesh::one_ring_tet_at_edge &ortae,
                           const zjucad::matrix::matrix<double> &vol,
                           const double w,
                           const string strategy = "face",
                           const string opt_type = "init",
                           const matrix<double> * size_field = 0)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_v = std::accumulate(zjucad::matrix::fabs(vol).begin(),
                                         zjucad::matrix::fabs(vol).end(), 0.0);

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  if(strategy == "face"){
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
          if(opt_type == "init")
            for(size_t di = 0; di < 9; ++di)
              all_func->push_back(math_func_ptr(
                                    new inner_smooth_func_sh_L1<double,int32_t>(
                                      tet_mesh.size(2),
                                      tet_pair.first, tet_pair.second, di,
                                      sqrt(size*w*(fabs(vol[tet_pair.first])+fabs(vol[tet_pair.second]))/(4.0*total_v)))));
          else if(opt_type == "zyz")
            all_func->push_back(math_func_ptr(
                                  new inner_smooth_func_zyz<double,int32_t>(
                                    tet_mesh.size(2), 0,
                                    tet_pair.first, tet_pair.second,
                                    sqrt(size*w*(fabs(vol[tet_pair.first])+fabs(vol[tet_pair.second]))/(4.0*total_v)))));
          else throw std::invalid_argument("# [error] can not recognize opt type");
        }
    }else{
      throw std::invalid_argument("wrong inner smooth strategy,[face/gradient/liuyang/laplace]");
    }
  cerr << "# [info] inner smooth func_L1 num " << all_func->size() << endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}


template <typename val_type, typename int_type>
class fix_length_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  fix_length_func(const size_t tet_num, const size_t idx, const double target, const double w)
    :tet_num_(tet_num), idx_(idx), target_(target), w_(w){}
  virtual ~fix_length_func(){}
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
        itr_matrix<const val_type *> x0(9,tet_num_,x);
        matrix<val_type> xd = x0(colon(),idx_);
        val_type d = dot(xd, xd) - target_;
        int_type c[1]={0};
        cv[c] += w_*d;
      }
    if(k == 1){
        itr_matrix<const val_type *> x0(9,tet_num_,x);
        for(int_type i = 0; i < 9; ++i){
            int_type c[2] = {0,9*idx_+i};
            cv[c] += 2*w_*x0(i,idx_);
          }
      }
    if(k == 2){
        for(int_type i = 0; i < 9; ++i){
            int_type c[3] = {0,9*idx_+i, 9*idx_+i};
            cv[c] += 2*w_;
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
            int_type c[2] = {0,9*idx_+i};
            l2g.add(cs,c);
          }
      }
    if(k == 2){
        for(int_type i = 0; i < 9; ++i){
            int_type c[3] = {0,9*idx_+i, 9*idx_+i};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 9;
    if(k == 2) return 9;
  }
private:
  const size_t tet_num_;
  const int_type idx_;
  const double target_;
  const double w_;
};

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_length_fix_func(const zjucad::matrix::matrix<size_t> &tet_mesh,
                      const zjucad::matrix::matrix<double> &node,
                      const zjucad::matrix::matrix<double> &vol,
                      const double target,
                      const double w)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  const double total_vol = std::accumulate(fabs(vol).begin(), fabs(vol).end(), 0.0);

  for(size_t ti = 0; ti < tet_mesh.size(2); ++ti){
      const double weight = sqrt(w * fabs(vol[ti])/total_vol);
      all_func->push_back(math_func_ptr(new fix_length_func<double,int32_t>(tet_mesh.size(2), ti, target, weight)));
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

template <typename val_type, typename int_type>
class surface_smooth_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  surface_smooth_func(const size_t tet_num, const size_t tri_num,
                      const size_t fi, const size_t fj,
                      const zjucad::matrix::matrix<double> rot, const double w)
    :tet_num_(tet_num), tri_num_(tri_num), fi_(fi), fj_(fj), rot_(rot), w_(w){}
  virtual ~surface_smooth_func(){}
  virtual size_t nx() const
  {
    return 9*tet_num_+2*tri_num_;
  }
  virtual size_t nf() const
  {
    return 2;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type*> x0(2,tri_num_, &x[9*tet_num_]);
        matrix<val_type> diff = rot_ * x0(colon(),fi_) - x0(colon(), fj_);
        for(int_type i = 0; i < 2; ++i){
            int_type c[1] = {i};
            cv[c] += w_*diff[i];
          }
      }
    if(k == 1){
        for(int_type i = 0; i < 2; ++i){
            for(int_type j = 0; j < 2; ++j){
                int_type c[2] = {i,9*tet_num_+2*fi_+j};
                cv[c] += w_*rot_(i,j);
              }
            int_type c[2] = {i,9*tet_num_+2*fj_+i};
            cv[c] += -w_;
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        for(int_type i = 0; i < 2; ++i){
            for(int_type j = 0; j < 2; ++j){
                int_type c[2] = {i,9*tet_num_+2*fi_+j};
                l2g.add(cs,c);
              }
            int_type c[2] = {i,9*tet_num_+2*fj_+i};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 3 * 2;
  }
private:
  const int_type tet_num_;
  const int_type tri_num_;
  const int_type fi_;
  const int_type fj_;
  const zjucad::matrix::matrix<val_type> rot_;
  const double w_;
};

template <typename val_type, typename int_type>
class surface_dir_fix_sh : public hj::math_func::math_func_t<val_type, int_type>

{
public:
  surface_dir_fix_sh(const size_t tet_num, const size_t tet_idx,
                     const matrix<double> & rl2g, const double target_angle,
                     const double w, const double ks = 1)
    :tet_num_(tet_num) , tet_idx_(tet_idx), rl2g_(rl2g), w_(w), ks_(ks){
    zjucad::matrix::matrix<val_type> rc_ = zeros<val_type>(2,9);
    rc_(0,8) = 1.0;
    rc_(1,0) = 1.0;
    rc_rl2g_T__ = ks_*rc_*trans(rl2g);
    target_.resize(2,1);
    target_[0] = sqrt(5.0) * cos(target_angle);
    target_[1] = sqrt(5.0) * sin(target_angle);
  }
  virtual ~surface_dir_fix_sh(){}
  virtual size_t nx() const
  {
    return 9*tet_num_;
  }
  virtual size_t nf() const
  {
    return 2;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type*> x0(9,tet_num_, x);
        matrix<val_type> diff = rc_rl2g_T__ * x0(colon(),tet_idx_) - target_;
        for(int_type i = 0; i < 2; ++i){
            int_type c[1] = {i};
            cv[c] += w_*diff[i];
          }
      }
    if(k == 1){
        for(int_type i = 0; i < 2; ++i){
            for(int_type j = 0; j < 9; ++j){
                int_type c[2] = {i,9*tet_idx_+j};
                cv[c] += w_*rc_rl2g_T__(i,j);
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
        for(int_type i = 0; i < 2; ++i){
            for(int_type j = 0; j < 9; ++j){
                int_type c[2] = {i,9*tet_idx_+j};
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return  18;
  }
private:
  zjucad::matrix::matrix<double> rc_rl2g_T__;
  zjucad::matrix::matrix<double> target_;
  const size_t tet_num_;
  const int_type tet_idx_;
  const matrix<double> & rl2g_;
  const double w_;
  const double ks_;
};


template <typename val_type, typename int_type>
class surface_dir_fix_zyz : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  surface_dir_fix_zyz(const size_t tet_num, const size_t tet_idx,
                      const matrix<double> & rl2g, const double target_angle,
                      const double w, const double ks = 1)
    :tet_num_(tet_num), tet_idx_(tet_idx), w_(w), ks_(ks){
    zjucad::matrix::matrix<val_type> rc_ = zeros<val_type>(2,9);
    rc_(0,8) = 1.0;
    rc_(1,0) = 1.0;
    rc_rl2g_T_ = ks_*rc_*trans(rl2g);
    target_.resize(2,1);
    target_[0] = sqrt(5.0) * cos(target_angle);
    target_[1] = sqrt(5.0) * sin(target_angle);
  }
  virtual ~surface_dir_fix_zyz(){}
  virtual size_t nx() const
  {
    return 3*tet_num_;
  }
  virtual size_t nf() const
  {
    return 2;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type*> x0(3,tet_num_, x);
        matrix<val_type> sh0(9,1);
        calc_rot_cubic_f_sh_(&sh0[0], &x0(0,tet_idx_));

        matrix<val_type> diff = rc_rl2g_T_* sh0 - target_;
        for(int_type i = 0; i < 2; ++i){
            int_type c[1] = {i};
            cv[c] += w_*diff[i];
          }
      }
    if(k == 1){
        itr_matrix<const val_type*> x0(3,tet_num_, x);
        matrix<val_type> jac_sh0(9,3);
        calc_jac_rot_cubic_f_sh_(&jac_sh0[0], &x0(0,tet_idx_));

        matrix<val_type> jac0 = rc_rl2g_T_ * jac_sh0;

        for(int_type i = 0; i < 2; ++i){
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i,3*tet_idx_+j};
                cv[c] += w_*jac0(i,j);
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
        for(int_type i = 0; i < 2; ++i){
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i,3*tet_idx_+j};
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 6;
  }
private:
  const size_t tet_num_;
  const int_type tet_idx_;
  zjucad::matrix::matrix<val_type> rc_rl2g_T_;
  const double w_;
  zjucad::matrix::matrix<val_type>  target_;
  const double ks_;
};


template <typename val_type, typename int_type>
class surface_smooth_func2_sh : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  surface_smooth_func2_sh(const size_t tet_num, const size_t t0, const size_t t1,
                          const zjucad::matrix::matrix<double> &rot,
                          const zjucad::matrix::matrix<double> &rl2g0,
                          const zjucad::matrix::matrix<double> &rl2g1,
                          const double w, const double k0 = 1, const double k1 = 1)
    :tet_num_(tet_num), t0_(t0), t1_(t1), w_(w), k0_(k0), k1_(k1){
    zjucad::matrix::matrix<val_type> rc_ = zeros<val_type>(2,9);
    rc_(0,8) = 1.0;
    rc_(1,0) = 1.0;
    rij_rc_rl2g0_T_ = rot * rc_*trans(rl2g0);
    rc_rl2g1_T_ = rc_*trans(rl2g1);
  }
  virtual ~surface_smooth_func2_sh(){}
  virtual size_t nx() const
  {
    return 9*tet_num_;
  }
  virtual size_t nf() const
  {
    return 2;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type*> x0(9,tet_num_, x);
        matrix<val_type> diff = k0_ * rij_rc_rl2g0_T_ * x0(colon(),t0_) - k1_*rc_rl2g1_T_*x0(colon(), t1_);
        for(int_type i = 0; i < 2; ++i){
            int_type c[1] = {i};
            cv[c] += w_*diff[i];
          }
      }
    if(k == 1){
        for(int_type i = 0; i < 2; ++i){
            for(int_type j = 0; j < 9; ++j){
                int_type c[2] = {i,9*t0_+j};
                cv[c] += w_*k0_*rij_rc_rl2g0_T_(i,j);
              }
            for(int_type j = 0; j < 9; ++j){
                int_type c[2] = {i,9*t1_+j};
                cv[c] += -w_*k1_*rc_rl2g1_T_(i,j);
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
        for(int_type i = 0; i < 2; ++i){
            for(int_type j = 0; j < 9; ++j){
                int_type c[2] = {i,9*t0_+j};
                l2g.add(cs,c);
              }
            for(int_type j = 0; j < 9; ++j){
                int_type c[2] = {i,9*t1_+j};
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 2*18;
  }
private:
  const size_t tet_num_;
  const int_type t0_;
  const int_type t1_;
  zjucad::matrix::matrix<val_type> rij_rc_rl2g0_T_;
  zjucad::matrix::matrix<val_type> rc_rl2g1_T_;
  const double w_;
  const double k0_, k1_;
};

template <typename val_type, typename int_type>
class surface_smooth_func2_zyz : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  surface_smooth_func2_zyz(const size_t tet_num, const size_t t0, const size_t t1,
                           const zjucad::matrix::matrix<double> &rot,
                           const zjucad::matrix::matrix<double> &rl2g0,
                           const zjucad::matrix::matrix<double> &rl2g1,
                           const double w, const double k0 = 1, const double k1 = 1)
    :tet_num_(tet_num), t0_(t0), t1_(t1), w_(w), k0_(k0), k1_(k1){
    zjucad::matrix::matrix<val_type> rc_ = zeros<val_type>(2,9);
    rc_(0,8) = 1.0;
    rc_(1,0) = 1.0;
    rij_rc_rl2g0_T_ = rot * rc_*trans(rl2g0);
    rc_rl2g1_T_ = rc_*trans(rl2g1);
  }
  virtual ~surface_smooth_func2_zyz(){}
  virtual size_t nx() const
  {
    return 3*tet_num_;
  }
  virtual size_t nf() const
  {
    return 2;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type*> x0(3,tet_num_, x);
        matrix<val_type> sh0(9,1), sh1(9,1);
        calc_rot_cubic_f_sh_(&sh0[0], &x0(0,t0_));
        calc_rot_cubic_f_sh_(&sh1[0], &x0(0,t1_));

        matrix<val_type> diff = k0_ * rij_rc_rl2g0_T_ * sh0 - k1_*rc_rl2g1_T_*sh1;
        for(int_type i = 0; i < 2; ++i){
            int_type c[1] = {i};
            cv[c] += w_*diff[i];
          }
      }
    if(k == 1){
        itr_matrix<const val_type*> x0(3,tet_num_, x);
        matrix<val_type> jac_sh0(9,3), jac_sh1(9,3);
        calc_jac_rot_cubic_f_sh_(&jac_sh0[0], &x0(0,t0_));
        calc_jac_rot_cubic_f_sh_(&jac_sh1[0], &x0(0,t1_));

        matrix<val_type> jac0 = rij_rc_rl2g0_T_ * jac_sh0;
        matrix<val_type> jac1 = rc_rl2g1_T_*jac_sh1;

        for(int_type i = 0; i < 2; ++i){
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i,3*t0_+j};
                cv[c] += w_*k0_*jac0(i,j);
              }
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i,3*t1_+j};
                cv[c] += -w_*k1_*jac1(i,j);
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
        for(int_type i = 0; i < 2; ++i){
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i,3*t0_+j};
                l2g.add(cs,c);
              }
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i,3*t1_+j};
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 2*6;
  }
private:
  const size_t tet_num_;
  const int_type t0_;
  const int_type t1_;
  zjucad::matrix::matrix<val_type> rij_rc_rl2g0_T_;
  zjucad::matrix::matrix<val_type> rc_rl2g1_T_;
  const double w_;
  const double k0_, k1_;
};

template <typename val_type, typename int_type>
class surface_smooth_func3 : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  surface_smooth_func3(const size_t tet_num, const size_t t0, const size_t t1,
                       const zjucad::matrix::matrix<double> &rot, const double w)
    :tet_num_(tet_num), t0_(t0), t1_(t1), w_(w), rot_(rot){}
  virtual ~surface_smooth_func3(){}
  virtual size_t nx() const
  {
    return 9*tet_num_;
  }
  virtual size_t nf() const
  {
    return 9;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type*> x0(9,tet_num_, x);
        matrix<val_type> diff = rot_ * x0(colon(), t0_) - x0(colon(), t1_);
        for(int_type i = 0; i < 9; ++i){
            int_type c[1] = {i};
            cv[c] += w_*diff[i];
          }
      }
    if(k == 1){
        for(int_type i = 0 ; i < 9; ++i){
            for(int_type j = 0; j < 9; ++j){
                int_type c[2] = {i,9*t0_+j};
                cv[c] += w_*rot_(i,j);
              }
            int_type c[2] = {i,9*t1_+i};
            cv[c] += -w_;
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        for(int_type i = 0 ; i < 9; ++i){
            for(int_type j = 0; j < 9; ++j){
                int_type c[2] = {i,9*t0_+j};
                l2g.add(cs,c);
              }
            int_type c[2] = {i,9*t1_+i};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 9*10;
  }
private:
  const size_t tet_num_;
  const int_type t0_;
  const int_type t1_;
  const zjucad::matrix::matrix<val_type> rot_;
  const double w_;
};

template <typename val_type, typename int_type>
class surface_smooth_func2_pow2 : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  surface_smooth_func2_pow2(const size_t tet_num, const size_t t0, const size_t t1,
                            const zjucad::matrix::matrix<double> &rot,
                            const zjucad::matrix::matrix<double> &rl2g0,
                            const zjucad::matrix::matrix<double> &rl2g1,
                            const double w)
    :tet_num_(tet_num), t0_(t0), t1_(t1), w_(w){
    zjucad::matrix::matrix<val_type> rc_ = zeros<val_type>(2,9);
    rc_(0,8) = 1.0;
    rc_(1,0) = 1.0;
    rij_rc_rl2g0_T_ = rot * rc_*trans(rl2g0);
    rc_rl2g1_T_ = rc_*trans(rl2g1);
    ATA_ = trans(rij_rc_rl2g0_T_) * rij_rc_rl2g0_T_;
    BTB_ = trans(rc_rl2g1_T_) * rc_rl2g1_T_;
    ATB_ = trans(rij_rc_rl2g0_T_) * rc_rl2g1_T_;
    BTA_ = trans(ATB_);
  }
  virtual ~surface_smooth_func2_pow2(){}
  virtual size_t nx() const
  {
    return 9*tet_num_;
  }
  virtual size_t nf() const
  {
    return 1;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type*> x0(9,tet_num_, x);
        matrix<val_type> diff = rij_rc_rl2g0_T_ * x0(colon(),t0_) - rc_rl2g1_T_*x0(colon(), t1_);
        int_type c[1] = {0};
        cv[c] += w_*dot(diff,diff);
      }
    if(k == 1){
        itr_matrix<const val_type*> x0(9,tet_num_, x);
        matrix<val_type> graA(9,1), graB(9,1);
        graA = 2*w_*(ATA_*x0(colon(),t0_)-ATB_*x0(colon(),t1_));
        graB = 2*w_*(BTB_*x0(colon(),t1_)-BTA_*x0(colon(),t0_));

        for(int_type j = 0; j < 9; ++j){
            int_type ca[2] = {0,9*t0_+j};
            cv[ca] += graA[j];

            int_type cb[2] = {0,9*t1_+j};
            cv[cb] += graB[j];
          }
      }
    if(k == 2){
        for(int_type i = 0; i < 9; ++i)
          for(int_type j = 0; j < 9; ++j){
              int_type caa[3] = {0,9*t0_+i,9*t0_+j};
              cv[caa] += 2*w_*ATA_(i,j);
              int_type cab[3] = {0,9*t0_+i,9*t1_+j};
              cv[cab] += -2*w_*ATB_(i,j);
              int_type cba[3] = {0,9*t1_+i,9*t0_+j};
              cv[cba] += -2*w_*BTA_(i,j);
              int_type cbb[3] = {0,9*t1_+i,9*t1_+j};
              cv[cbb] += 2*w_*BTB_(i,j);
            }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        for(int_type j = 0; j < 9; ++j){
            int_type ca[2] = {0,9*t0_+j};
            l2g.add(cs,ca);

            int_type cb[2] = {0,9*t1_+j};
            l2g.add(cs,cb);
          }
      }
    if(k == 2){
        for(int_type i = 0; i < 9; ++i)
          for(int_type j = 0; j < 9; ++j){
              int_type caa[3] = {0,9*t0_+i,9*t0_+j};
              l2g.add(cs,caa);

              int_type cab[3] = {0,9*t0_+i,9*t1_+j};
              l2g.add(cs,cab);

              int_type cba[3] = {0,9*t1_+i,9*t0_+j};
              l2g.add(cs,cba);

              int_type cbb[3] = {0,9*t1_+i,9*t1_+j};
              l2g.add(cs,cbb);
            }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 18;
    if(k == 2) return 18*18;
  }
private:
  const size_t tet_num_;
  const int_type t0_;
  const int_type t1_;
  zjucad::matrix::matrix<val_type> rij_rc_rl2g0_T_;
  zjucad::matrix::matrix<val_type> rc_rl2g1_T_;
  zjucad::matrix::matrix<val_type> ATA_, BTB_, ATB_, BTA_;
  const double w_;
};

template <typename val_type, typename int_type>
class surface2inner_func2: public hj::math_func::math_func_t<val_type, int_type>
{
public:
  surface2inner_func2(const size_t tet_num, const size_t tri_num,
                      const size_t tet_idx, const size_t tri_idx,
                      const zjucad::matrix::matrix<val_type> local2global_rot,
                      const val_type w)
    :tet_num_(tet_num), tri_num_(tri_num),tet_idx_(tet_idx), tri_idx_(tri_idx),w_(w),
      local2global_rot_(local2global_rot){
    assert(tri_idx_ < tri_num && tet_idx_ < tet_num);
    assert(local2global_rot.size(1) == 9 && local2global_rot.size(2) == 9);
  }
  virtual ~surface2inner_func2(){}
  virtual size_t nx() const{
    return 9*tet_num_ + 2*tri_num_;
  }
  virtual size_t nf() const{
    return 1;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0)const
  {
    itr_matrix<const val_type*> sh_x(9,tet_num_, &x[0]);
    itr_matrix<const val_type*> f_x(2, tri_num_, &x[9*tet_num_]);

    if(k == 0){
        matrix<val_type> sh_f(11,1);
        itr_matrix<val_type*> sh_f_sh(9,1,&sh_f[0]);
        itr_matrix<val_type*> sh_f_f(2,1, &sh_f[9]);
        sh_f_sh = sh_x(colon(),tet_idx_);
        sh_f_f = f_x(colon(),tri_idx_);
        val_type v;
        relate_diff_energy_(&v, &sh_f[0], &local2global_rot_[0]);
        jtf::math::erase_nan_inf(v);

        int_type c[1] = {0};
        cv[c] += w_*v;

      }
    if(k == 1){
        matrix<val_type> sh_f(11,1);
        itr_matrix<val_type*> sh_f_sh(9,1,&sh_f[0]);
        itr_matrix<val_type*> sh_f_f(2,1, &sh_f[9]);
        sh_f_sh = sh_x(colon(),tet_idx_);
        sh_f_f = f_x(colon(),tri_idx_);
        matrix<val_type> jac(11,1);
        relate_diff_energy_jac_(&jac[0], &sh_f[0], &local2global_rot_[0]);
        for(size_t i = 0; i < jac.size(); ++i) jtf::math::erase_nan_inf(jac[i]);

        for(size_t j = 0; j < 9; ++j){
            int_type c[2] = {0,9*tet_idx_+j};
            cv[c] += w_*jac[j];
          }
        for(size_t j = 0; j < 2; ++j){
            int_type c[2] = {0, 9*tet_num_ + tri_idx_*2+j};
            cv[c] += w_*jac[j+9];
          }
      }
    if(k==2){
        matrix<val_type> sh_f(11,1);
        itr_matrix<val_type*> sh_f_sh(9,1,&sh_f[0]);
        itr_matrix<val_type*> sh_f_f(2,1, &sh_f[9]);
        sh_f_sh = sh_x(colon(),tet_idx_);
        sh_f_f = f_x(colon(),tri_idx_);
        matrix<val_type> hes(11,11);
        relate_diff_energy_hes_(&hes[0], &sh_f[0], &local2global_rot_[0]);
        for(size_t i = 0; i < hes.size(); ++i) jtf::math::erase_nan_inf(hes[i]);

        for(size_t i = 0; i < 11; ++i){
            const size_t jac_idx = (i < 9?9*tet_idx_+i:9*tet_num_+tri_idx_*2+(i-9));
            for(size_t j = 0; j < 11; ++j){
                const size_t hes_idx = (j < 9?9*tet_idx_+j:9*tet_num_+tri_idx_*2+(j-9));
                int_type c[3] = {0,jac_idx, hes_idx};
                cv[c] += w_*hes(i,j);
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
        for(size_t j = 0; j < 9; ++j){
            int_type c[2] = {0,9*tet_idx_+j};
            l2g.add(cs,c);
          }
        for(size_t j = 0; j < 2; ++j){
            int_type c[2] = {0, 9*tet_num_ + tri_idx_*2+j};
            l2g.add(cs,c);
          }
      }
    if(k==2){
        for(size_t i = 0; i < 11; ++i){
            const size_t jac_idx = (i < 9?9*tet_idx_+i:9*tet_num_+tri_idx_*2+(i-9));
            for(size_t j = 0; j < 11; ++j){
                const size_t hes_idx = (j < 9?9*tet_idx_+j:9*tet_num_+tri_idx_*2+(j-9));
                int_type c[3] = {0,jac_idx, hes_idx};
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k==0) return -1;
    if(k==1) return 11;
    if(k==2) return 121;
  }
private:
  const size_t tet_num_;
  const size_t tri_num_;
  const size_t tet_idx_;
  const size_t tri_idx_;
  const double w_;
  const zjucad::matrix::matrix<val_type> local2global_rot_;
};


template <typename val_type, typename int_type>
class surface_fix_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  surface_fix_func(const size_t tet_num, const size_t tri_num,
                   const size_t tri_idx, const zjucad::matrix::matrix<val_type> & cos_sin,
                   const val_type w):tet_num_(tet_num), tri_num_(tri_num), tri_idx_(tri_idx), cos_sin_(cos_sin), w_(w){}
  virtual ~surface_fix_func(){}
  virtual size_t nx() const
  {
    return 9*tet_num_ + 2 * tri_num_;
  }
  virtual size_t nf() const{
    return 2;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==0){
        itr_matrix<const val_type*> x0(2, tri_num_, &x[9*tet_num_]);
        matrix<val_type> diff = x0(colon(), tri_idx_) - cos_sin_;
        for(int_type fi = 0; fi < nf(); ++fi){
            int_type c[1] = {fi};
            cv[c] += w_*diff[fi];
          }
      }
    if(k==1){
        for(int_type fi = 0; fi < nf(); ++fi){
            int_type c[2] = {fi,9*tet_num_+2*tri_idx_+fi};
            cv[c] += w_;
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        for(int_type fi = 0; fi < nf(); ++fi){
            int_type c[2] = {fi,9*tet_num_+2*tri_idx_+fi};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k == 0) return -1;
    if(k == 1) return 2;
  }

private:
  const int_type tet_num_;
  const int_type tri_num_;
  const int_type tri_idx_;
  const zjucad::matrix::matrix<double> cos_sin_;
  const val_type w_;
};

template <typename val_type, typename int_type>
class surface_fix_func2 : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  surface_fix_func2(const size_t tet_num, const size_t tet_idx,
                    const zjucad::matrix::matrix<val_type> &cos_sin,
                    const zjucad::matrix::matrix<val_type> rc_rot,
                    const val_type w)
    :tet_num_(tet_num), tet_idx_(tet_idx), cos_sin_(cos_sin), w_(w), rc_rot_(rc_rot){}
  virtual ~surface_fix_func2(){}
  virtual size_t nx() const
  {
    return 9*tet_num_;
  }
  virtual size_t nf() const{
    return 2;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==0){
        itr_matrix<const val_type*> x0(9, tet_num_, &x[0]);
        matrix<val_type> diff = rc_rot_*x0(colon(), tet_idx_) - cos_sin_;
        for(int_type fi = 0; fi < nf(); ++fi){
            int_type c[1] = {fi};
            cv[c] += w_*diff[fi];
          }
      }
    if(k==1){
        for(int_type fi = 0; fi < nf(); ++fi){
            for(int_type di = 0; di < 9; ++di){
                int_type c[2] = {fi,9*tet_idx_+di};
                cv[c] += w_;
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
            for(int_type di = 0; di < 9; ++di){
                int_type c[2] = {fi,9*tet_idx_+di};
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k == 0) return -1;
    if(k == 1) return 2*9;
  }

private:
  const size_t tet_num_;
  const int_type tet_idx_;
  const zjucad::matrix::matrix<val_type> cos_sin_;
  const zjucad::matrix::matrix<val_type> rc_rot_;
  const val_type w_;
};

template <typename val_type, typename int_type>
class fix_zyz_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  fix_zyz_func(const size_t tet_num, const size_t tet_idx,
               const zjucad::matrix::matrix<double> & target_zyz,
               const val_type w)
    :tet_num_(tet_num), tet_idx_(tet_idx), target_zyz_(target_zyz), w_(w){}
  virtual ~fix_zyz_func(){}
  virtual size_t nx() const
  {
    return 3*tet_num_;
  }
  virtual size_t nf() const{
    return 3;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==0){
        itr_matrix<const val_type*> x0(3, tet_num_, &x[0]);
        matrix<val_type> diff = x0(colon(), tet_idx_) - target_zyz_;
        for(int_type fi = 0; fi < nf(); ++fi){
            int_type c[1] = {fi};
            cv[c] += w_*diff[fi];
          }
      }
    if(k==1){
        for(int_type fi = 0; fi < nf(); ++fi){
            int_type c[2] = {fi,3*tet_idx_+fi};
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
            int_type c[2] = {fi,3*tet_idx_+fi};
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
  const size_t tet_num_;
  const int_type tet_idx_;
  const zjucad::matrix::matrix<val_type> target_zyz_;
  const val_type w_;
};

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_surface_fix_func(const zjucad::matrix::matrix<size_t> & tet_mesh,
                       const zjucad::matrix::matrix<size_t> & outside_face_idx,
                       const jtf::mesh::face2tet_adjacent &fa,
                       const zjucad::matrix::matrix<double> &surface_area,
                       const zjucad::matrix::matrix<double> &init_zyz,
                       double w_surface,
                       const string opt_type)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(surface_area.begin(), surface_area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  for(size_t fi = 0 ; fi < outside_face_idx.size(); ++fi){
      const pair<size_t,size_t> &tet_pair = fa.face2tet_[outside_face_idx[fi]];
      const size_t tet_idx = tet_pair.first == -1?tet_pair.second:tet_pair.first;
      const double w = sqrt(fabs(surface_area[fi]/total_area));
      all_func->push_back(math_func_ptr(
                            new fix_zyz_func<double,int32_t>(
                              tet_mesh.size(2), tet_idx,
                              init_zyz(colon(), tet_idx),
                              w)));
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_surface_fix_func2(const zjucad::matrix::matrix<size_t> &tet_mesh,
                        const zjucad::matrix::matrix<size_t> &tri_mesh_idx,
                        const zjucad::matrix::matrix<double> &area,
                        const jtf::mesh::face2tet_adjacent &fa,
                        const zjucad::matrix::matrix<double> &local2global,
                        const zjucad::matrix::matrix<double> &surface_cos_sin,
                        const double w)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  assert(surface_cos_sin.size() == 2 * tri_mesh_idx.size());
  itr_matrix<const double*> surface_cos_sin_m(2, surface_cos_sin.size()/2, &surface_cos_sin[0]);

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  matrix<double> rc = zeros<double>(2,9);
  rc(0,8) = 1.0; rc(1,0) = 1.0;
  matrix<double> rot(9,9);
  matrix<double> rc_rot(2,9);
  for(size_t fi = 0 ; fi < tri_mesh_idx.size(2); ++fi){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[tri_mesh_idx[fi]];
      const size_t tet_idx = tet_pair.first==-1?tet_pair.second:tet_pair.first;
      calc_rot_cubic_f_sh_mat_(&rot[0], &local2global(0,fi));

      rot = temp(trans(rot));
      rc_rot = rc * rot;
      all_func->push_back(math_func_ptr(
                            new surface_fix_func2<double,int32_t>(
                              tet_mesh.size(2), tet_idx, surface_cos_sin_m(colon(),fi),
                              rc_rot,  sqrt(w*(area[fi])/total_area))));
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}


shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_surface_smooth_func(const zjucad::matrix::matrix<size_t> &tet_mesh,
                          const zjucad::matrix::matrix<size_t> &tri_mesh,
                          const jtf::mesh::edge2cell_adjacent &ea,
                          const zjucad::matrix::matrix<double> &area,
                          const zjucad::matrix::matrix<double> &rij,
                          const double w)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  matrix<double> rot(2,2);
  const double Nsys = 4;
  for(size_t ei = 0; ei < ea.edge2cell_.size(); ++ei){
      const pair<size_t,size_t> & face_pair = ea.edge2cell_[ei];
      if(ea.is_boundary_edge(face_pair)) continue;
      get_2d_rotation(Nsys*rij(0,ei), rot);
      all_func->push_back(math_func_ptr(
                            new surface_smooth_func<double,int32_t>(
                              tet_mesh.size(2), tri_mesh.size(2),
                              face_pair.first, face_pair.second, rot,
                              sqrt(w*(area[face_pair.first]+area[face_pair.second])/(3.0*total_area)))));
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_surface_smooth_func3(const zjucad::matrix::matrix<size_t> &tet_mesh,
                           const zjucad::matrix::matrix<size_t> &tri_mesh,
                           const jtf::mesh::edge2cell_adjacent &ea,
                           const zjucad::matrix::matrix<double> &area,
                           const zjucad::matrix::matrix<double> &rij,
                           const double w)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  matrix<double> rot(2,2);
  const size_t N_sys = 4;
  for(size_t ei = 0; ei < ea.edge2cell_.size(); ++ei){
      const pair<size_t,size_t> & face_pair = ea.edge2cell_[ei];
      if(ea.is_boundary_edge(face_pair)) continue;
      get_2d_rotation(N_sys*rij(0,ei), rot);
      all_func->push_back(math_func_ptr(
                            new surface_smooth_func<double,int32_t>(
                              0, tri_mesh.size(2),
                              face_pair.first, face_pair.second, rot,
                              sqrt(w*(area[face_pair.first]+area[face_pair.second])/(3.0*total_area)))));
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_surface_smooth_func2(const zjucad::matrix::matrix<size_t> &tet_mesh,
                           const zjucad::matrix::matrix<size_t> &tri_mesh_idx,
                           const jtf::mesh::edge2cell_adjacent &ea,
                           const jtf::mesh::face2tet_adjacent &fa,
                           const zjucad::matrix::matrix<double> &area,
                           const zjucad::matrix::matrix<double> &rij,
                           const zjucad::matrix::matrix<double> &local2global,
                           const double w)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  matrix<double> rot(2,2);
  matrix<double> rl2g0(9,9), rl2g1(9,9);
  const double Nsys = 4;
  for(size_t ei = 0; ei < ea.edge2cell_.size(); ++ei){
      const pair<size_t,size_t> & face_pair = ea.edge2cell_[ei];
      if(ea.is_boundary_edge(face_pair)) continue;
      const pair<size_t,size_t> &tet_pair0 = fa.face2tet_[tri_mesh_idx[face_pair.first]];
      const pair<size_t,size_t> &tet_pair1 = fa.face2tet_[tri_mesh_idx[face_pair.second]];
      const size_t t0 = tet_pair0.first==-1?tet_pair0.second:tet_pair0.first;
      const size_t t1 = tet_pair1.first==-1?tet_pair1.second:tet_pair1.first;

      get_2d_rotation(Nsys*rij(0,ei), rot);
      calc_rot_cubic_f_sh_mat_(&rl2g0[0], &local2global(0, face_pair.first));
      calc_rot_cubic_f_sh_mat_(&rl2g1[0], &local2global(0, face_pair.second));

      all_func->push_back(math_func_ptr(
                            new surface_smooth_func2_sh<double,int32_t>(
                              tet_mesh.size(2), t0,t1, rot, rl2g0, rl2g1,
                              sqrt(w*(area[face_pair.first]+area[face_pair.second])/(3.0*total_area)))));
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_surface_smooth_func4(const zjucad::matrix::matrix<size_t> &tet_mesh,
                           const zjucad::matrix::matrix<size_t> &tri_mesh_idx,
                           const jtf::mesh::edge2cell_adjacent &ea,
                           const jtf::mesh::face2tet_adjacent &fa,
                           const zjucad::matrix::matrix<double> &area,
                           const zjucad::matrix::matrix<double> &rij,
                           const zjucad::matrix::matrix<double> &local2global,
                           const zjucad::matrix::matrix<double> &ks,
                           const double w,
                           const string opt_type = "init")
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  matrix<double> rot(2,2);
  matrix<double> rl2g0(9,9), rl2g1(9,9);
  const double Nsys = 4;
  for(size_t ei = 0; ei < ea.edge2cell_.size(); ++ei){
      const pair<size_t,size_t> & face_pair = ea.edge2cell_[ei];
      if(ea.is_boundary_edge(face_pair)) continue;
      const pair<size_t,size_t> &tet_pair0 = fa.face2tet_[tri_mesh_idx[face_pair.first]];
      const pair<size_t,size_t> &tet_pair1 = fa.face2tet_[tri_mesh_idx[face_pair.second]];
      const size_t t0 = tet_pair0.first==-1?tet_pair0.second:tet_pair0.first;
      const size_t t1 = tet_pair1.first==-1?tet_pair1.second:tet_pair1.first;

      get_2d_rotation(Nsys*rij(0,ei), rot);
      calc_rot_cubic_f_sh_mat_(&rl2g0[0], &local2global(0, face_pair.first));
      calc_rot_cubic_f_sh_mat_(&rl2g1[0], &local2global(0, face_pair.second));

      if(opt_type == "init")
        all_func->push_back(math_func_ptr(
                              new surface_smooth_func2_sh<double,int32_t>(
                                tet_mesh.size(2), t0,t1, rot, rl2g0, rl2g1,
                                sqrt(w*(fabs(area[face_pair.first])+fabs(area[face_pair.second]))/(3.0*total_area)))));
      else if(opt_type == "zyz")
        all_func->push_back(math_func_ptr(
                              new surface_smooth_func2_zyz<double,int32_t>(
                                tet_mesh.size(2), t0,t1, rot, rl2g0, rl2g1,
                                sqrt(w*(fabs(area[face_pair.first])+fabs(area[face_pair.second]))/(3.0*total_area)))));
      else throw std::invalid_argument("# [error] can not recognize opt type.");
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}


shared_ptr<hj::math_func::math_func_t<double, int32_t> >
build_surface_dir_align_func(const zjucad::matrix::matrix<size_t> &tet_mesh,
                             const jtf::mesh::face2tet_adjacent &fa,
                             const zjucad::matrix::matrix<double> &area,
                             const zjucad::matrix::matrix<double> &local2global,
                             const zjucad::matrix::matrix<double> &surface_dir_angle,
                             const zjucad::matrix::matrix<size_t> & surface_dir_idx,
                             const zjucad::matrix::matrix<size_t> & surface_dix_idx_on_boundary,
                             const double w,
                             const string opt_type = "init")
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  double total_area = 0;
  for(size_t fi = 0; fi < surface_dix_idx_on_boundary.size(); ++fi)
    total_area += fabs(area[surface_dix_idx_on_boundary[fi]]);

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  matrix<double> rl2g(9,9);

  for(size_t fi = 0; fi < surface_dir_idx.size(); ++fi){
      const pair<size_t,size_t> &tet_pair = fa.face2tet_[surface_dir_idx[fi]];
      const size_t tet_idx = tet_pair.first==-1?tet_pair.second:tet_pair.first;
      calc_rot_cubic_f_sh_mat_(&rl2g[0], &local2global(0, surface_dix_idx_on_boundary[fi]));
      if(opt_type == "init")
        all_func->push_back(math_func_ptr(
                              new surface_dir_fix_sh<double,int32_t>(
                                tet_mesh.size(2), tet_idx, rl2g, surface_dir_angle[fi],
                                sqrt(w*fabs(area[surface_dix_idx_on_boundary[fi]])/total_area))));
      else if(opt_type == "zyz"){
          all_func->push_back(math_func_ptr(
                                new surface_dir_fix_zyz<double,int32_t>(
                                  tet_mesh.size(2), tet_idx, rl2g, surface_dir_angle[fi],
                                  sqrt(w*fabs(area[surface_dix_idx_on_boundary[fi]])/total_area))));
        }
      else throw std::invalid_argument("# [error] can not recognize opt type.");
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_surface_smooth_func5(const zjucad::matrix::matrix<size_t> &tet_mesh,
                           const zjucad::matrix::matrix<size_t> &tri_mesh_idx,
                           const jtf::mesh::edge2cell_adjacent &ea,
                           const jtf::mesh::face2tet_adjacent &fa,
                           const zjucad::matrix::matrix<double> &area,
                           const zjucad::matrix::matrix<double> &rij,
                           const zjucad::matrix::matrix<double> &local2global,
                           const zjucad::matrix::matrix<double> &ks,
                           const double w)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  matrix<double> rot = zeros<double>(3,1);
  matrix<double> rl2g0(9,9), rl2g1(9,9), rot_sh(9,9),R(9,9);
  const double Nsys = 4;
  for(size_t ei = 0; ei < ea.edge2cell_.size(); ++ei){
      const pair<size_t,size_t> & face_pair = ea.edge2cell_[ei];
      if(ea.is_boundary_edge(face_pair)) continue;
      const pair<size_t,size_t> &tet_pair0 = fa.face2tet_[tri_mesh_idx[face_pair.first]];
      const pair<size_t,size_t> &tet_pair1 = fa.face2tet_[tri_mesh_idx[face_pair.second]];
      const size_t t0 = tet_pair0.first==-1?tet_pair0.second:tet_pair0.first;
      const size_t t1 = tet_pair1.first==-1?tet_pair1.second:tet_pair1.first;

      rot[0] = Nsys*rij(0,ei);
      calc_rot_cubic_f_sh_mat_(&rot_sh[0], &rot[0]);
      calc_rot_cubic_f_sh_mat_(&rl2g0[0], &local2global(0, face_pair.first));
      calc_rot_cubic_f_sh_mat_(&rl2g1[0], &local2global(0, face_pair.second));
      R = rl2g1 * rot_sh * trans(rl2g0);
      all_func->push_back(math_func_ptr(
                            new surface_smooth_func3<double,int32_t>(
                              tet_mesh.size(2), t0,t1, R,
                              sqrt(w*(area[face_pair.first]+area[face_pair.second])/(3.0*total_area)))));
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

template <typename val_type, typename int_type>
class surface2inner_func: public hj::math_func::math_func_t<val_type, int_type>
{
public:
  surface2inner_func(const size_t tet_num, const size_t tri_num,
                     const size_t tet_idx, const size_t tri_idx,
                     const zjucad::matrix::matrix<val_type> local2global_rot,
                     const val_type w)
    :tet_num_(tet_num), tri_num_(tri_num),tet_idx_(tet_idx), tri_idx_(tri_idx),w_(w){
    assert(tri_idx_ < tri_num && tet_idx_ < tet_num);
    assert(local2global_rot.size(1) == 9 && local2global_rot.size(2) == 9);
    zjucad::matrix::matrix<val_type> e4_  = zeros<double>(9,1);
    e4_[4] = 1.0;
    zjucad::matrix::matrix<val_type> rcT_ = zeros<double>(9,2);
    rcT_(0,1) = 1.0;
    rcT_(8,0) = 1.0;
    Rlocal2global_RcT_ = local2global_rot * rcT_ * sqrt(5.0)*w_;
    Rlocal2global_e4_ = local2global_rot * e4_ * sqrt(7.0)*w_;
  }
  virtual size_t nx() const{
    return 9*tet_num_ + 2*tri_num_;
  }
  virtual size_t nf() const{
    return 9;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0)const
  {
    itr_matrix<const val_type*> sh_x(9,tet_num_, &x[0]);
    itr_matrix<const val_type*> f_x(2, tri_num_, &x[9*tet_num_]);

    if(k == 0){
        //matrix<val_type> diff = local2global_rot_ *(sqrt(5.0)*rcT_*f_x(colon(),tri_idx_)+sqrt(7.0)*e4_) - sh_x(colon(),tet_idx_);
        matrix<val_type> diff = Rlocal2global_RcT_ * f_x(colon(),tri_idx_) + Rlocal2global_e4_ - w_*sh_x(colon(),tet_idx_);
        for(size_t i = 0; i < 9; ++i){
            int_type c[1] = {i};
            cv[c] += diff[i];
          }
      }
    if(k == 1){
        for(size_t i = 0; i < 9; ++i){
            for(size_t j = 0; j < 2; ++j){
                int_type c[2] = {i, 9*tet_num_ + tri_idx_*2+j};
                cv[c] += Rlocal2global_RcT_(i,j);
              }
            int_type c[2] = {i,9*tet_idx_+i};
            cv[c] += -w_;
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        for(size_t i = 0; i < 9 ; ++i){
            for(size_t j = 0; j < 2; ++j){
                int_type c[2] = {i, 9*tet_num_ + tri_idx_*2+j};
                l2g.add(cs,c);
              }
            int_type c[2] = {i,9*tet_idx_+i};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k==0) return -1;
    if(k==1) return 3*9;
  }
private:
  const size_t tet_num_;
  const size_t tri_num_;
  const size_t tet_idx_;
  const size_t tri_idx_;
  const double w_;
  matrix<val_type> Rlocal2global_RcT_;
  matrix<val_type> Rlocal2global_e4_;
};

template <typename val_type, typename int_type>
class surface_cross_norm_sh: public hj::math_func::math_func_t<val_type, int32_t>
{
public:
  surface_cross_norm_sh(const size_t tet_num, const size_t tet_idx,
                        const zjucad::matrix::matrix<val_type> &rot,
                        const val_type w, const val_type ks = 1)
    :tet_num_(tet_num), tet_idx_(tet_idx), w_(w), ks_(ks){
    basis0_ = rot(0,colon());
    basis8_ = rot(8,colon());
    A_ = trans(basis0_)*basis0_  + trans(basis8_)*basis8_ ;
    A_ *= ks_;
  }
  virtual ~surface_cross_norm_sh(){}
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
        itr_matrix<const val_type*> x0(9,tet_num_,&x[0]);
        int_type c[1] = {0};
        cv[c] += w_*(pow(ks_*dot(basis0_, x0(colon(),tet_idx_)),2.0)+
                     pow(ks_*dot(basis8_, x0(colon(),tet_idx_)),2.0)-5.0);
      }
    if(k==1){
        itr_matrix<const val_type*> x0(9,tet_num_,&x[0]);
        matrix<val_type> jac = A_*x0(colon(), tet_idx_);
        jac *= 2.0*w_;
        for(int_type i = 0; i < 9; ++i){
            int_type c[2] = {0,9*tet_idx_+i};
            cv[c] += jac[i];
          }
      }
    if(k == 2){
        for(int_type i = 0; i < 9; ++i)
          for(int_type j = 0 ; j < 9; ++j){
              int_type c[3] = {0, 9*tet_idx_+i, 9*tet_idx_+j};
              cv[c] += 2*w_*A_(i,j);
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
    if(k == 2){
        for(int_type i = 0; i < 9; ++i)
          for(int_type j = 0 ; j < 9; ++j){
              int_type c[3] = {0, 9*tet_idx_+i, 9*tet_idx_+j};
              l2g.add(cs,c);
            }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k==0) return -1;
    if(k==1) return 9;
    if(k==2) return 81;
  }
private:
  const int_type tet_num_;
  const int_type tet_idx_;
  zjucad::matrix::matrix<val_type> basis0_, basis8_;
  zjucad::matrix::matrix<val_type> A_;
  const val_type w_;
  const val_type ks_;
};

template <typename val_type, typename int_type>
class surface_normal_align2_pow2: public hj::math_func::math_func_t<val_type, int32_t>
{
public:
  surface_normal_align2_pow2(const size_t tet_num, const size_t tet_idx,
                             const zjucad::matrix::matrix<val_type> &rot,
                             const size_t basis_idx, const double v,
                             const val_type w)
    :tet_num_(tet_num), tet_idx_(tet_idx), w_(w), v_(v){
    basis_ = rot(basis_idx,colon());
    basisT_basis_ = trans(basis_) * basis_;
  }
  virtual ~surface_normal_align2_pow2(){}
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
        itr_matrix<const val_type*> x0(9,tet_num_,&x[0]);
        int_type c[1] = {0};
        cv[c] += w_*pow(dot(basis_, x0(colon(),tet_idx_))-v_,2.0);
      }
    if(k==1){
        itr_matrix<const val_type*> x0(9,tet_num_,&x[0]);
        zjucad::matrix::matrix<val_type> bbx = basisT_basis_*x0(colon(), tet_idx_);
        for(int_type i = 0; i < 9; ++i){
            int_type c[2] = {0,9*tet_idx_+i};
            cv[c] += 2*w_*(bbx[i]-v_*basis_[i]);
          }
      }
    if(k==2){
        const val_type bb = dot(basis_,basis_);
        for(int_type i = 0; i < 9; ++i){
            for(int_type j = 0; j < 9; ++j){
                int_type c[3] = {0, 9 *tet_idx_+i, 9*tet_idx_+j} ;
                cv[c] += 2*w_*basisT_basis_(i,j);
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
        for(int_type i = 0; i < 9; ++i){
            int_type c[2] = {0,9*tet_idx_+i};
            l2g.add(cs,c);
          }
      }
    if(k==2){
        for(int_type i = 0; i < 9; ++i){
            for(int_type j = 0; j < 9; ++j){
                int_type c[3] = {0, 9 *tet_idx_+i, 9*tet_idx_+j} ;
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k==0) return -1;
    if(k==1) return 9;
    if(k==2) return 81;
  }
private:
  const size_t tet_num_;
  const int_type tet_idx_;
  zjucad::matrix::matrix<val_type> basis_;
  zjucad::matrix::matrix<val_type> basisT_basis_;
  const val_type w_;
  const val_type v_;
};

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_normal_align_func2(const zjucad::matrix::matrix<size_t> &tet_mesh,
                         const jtf::mesh::face2tet_adjacent &fa,
                         const zjucad::matrix::matrix<size_t> &tri_idx,
                         const zjucad::matrix::matrix<double> &local2global,
                         const zjucad::matrix::matrix<double> &area,
                         const double w_normal_align)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  matrix<double> rot(9,9);
  const double v = sqrt(7.0); const size_t basis = 4;
  for(size_t fi = 0; fi < tri_idx.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[tri_idx[fi]];
      assert(fa.is_outside_face(tet_pair));
      const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);

      calc_rot_cubic_f_sh_mat_(&rot[0], &local2global(0,fi));
      rot = temp(trans(rot));

      all_func->push_back(math_func_ptr(
                            new surface_normal_align2_sh<double,int32_t>(
                              tet_mesh.size(2), tet_idx, rot, basis, v,
                              sqrt(w_normal_align * area[fi]/total_area))));
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_normal_align_func4(const zjucad::matrix::matrix<size_t> &tet_mesh,
                         const jtf::mesh::face2tet_adjacent &fa,
                         const zjucad::matrix::matrix<size_t> &tri_idx,
                         const zjucad::matrix::matrix<double> &local2global,
                         const zjucad::matrix::matrix<double> &area,
                         const zjucad::matrix::matrix<double> &ks,
                         const double w_normal_align,
                         const bool only_basis_4 = true,
                         const string opt_type = "init")
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  matrix<double> rot(9,9);
  double v = 0;
  for(size_t fi = 0; fi < tri_idx.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[tri_idx[fi]];
      assert(fa.is_outside_face(tet_pair));
      const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);

      calc_rot_cubic_f_sh_mat_(&rot[0], &local2global(0,fi));
      rot = temp(trans(rot));

      for(size_t di = 1; di < 8; ++di){
          if(di == 4) {
              v = sqrt(7.0);
            } else {
              if(only_basis_4) continue;
              v = 0;
            }
          if(opt_type == "init")
            all_func->push_back(math_func_ptr(
                                  new surface_normal_align2_sh<double,int32_t>(
                                    tet_mesh.size(2), tet_idx, rot, di, v,
                                    sqrt(w_normal_align * area[fi]/total_area))));
          else if(opt_type == "zyz")
            all_func->push_back(math_func_ptr(
                                  new surface_normal_align2_zyz<double,int32_t>(
                                    tet_mesh.size(2), tet_idx, rot, di, v,
                                    sqrt(w_normal_align * area[fi]/total_area))));
          else throw std::invalid_argument("# [error] can not recognize opt type.");
        }
    }
  cerr << "# [info] normal_align_func number " << all_func->size() << endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_normal_align_func6(const zjucad::matrix::matrix<size_t> &tet_mesh,
                         const jtf::mesh::face2tet_adjacent &fa,
                         const zjucad::matrix::matrix<size_t> &tri_idx,
                         const zjucad::matrix::matrix<double> &surface_normal,
                         const zjucad::matrix::matrix<double> &area,
                         const zjucad::matrix::matrix<double> &ks,
                         const double w_normal_align,
                         const bool only_basis_4 = true,
                         const string opt_type = "init")
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  matrix<double> rot(9,9),zyz(3,1);
  double v = 0;
  for(size_t fi = 0; fi < tri_idx.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[tri_idx[fi]];
      assert(fa.is_outside_face(tet_pair));
      const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);

      rot_n_2_z_by_zyz(&surface_normal(0,fi), &zyz[0]);
      calc_rot_cubic_f_sh_mat_(&rot[0], &zyz[0]);

      for(size_t di = 1; di < 8; ++di){
          if(di == 4) {
              v = sqrt(7.0);
            } else {
              if(only_basis_4) continue;
              v = 0;
            }
          if(opt_type == "init")
            all_func->push_back(math_func_ptr(
                                  new surface_normal_align2_sh<double,int32_t>(
                                    tet_mesh.size(2), tet_idx, rot, di, v,
                                    sqrt(w_normal_align * area[fi]/total_area))));
          else if(opt_type == "zyz")
            all_func->push_back(math_func_ptr(
                                  new surface_normal_align2_zyz<double,int32_t>(
                                    tet_mesh.size(2), tet_idx, rot, di, v,
                                    sqrt(w_normal_align * area[fi]/total_area))));
          else throw std::invalid_argument("# [error] can not recognize opt type.");
        }
    }
  cerr << "# [info] normal_align_func number " << all_func->size() << endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}




shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_normal_align_func4_with_deformation_gradient(
    const zjucad::matrix::matrix<size_t> &tet_mesh,
    const jtf::mesh::face2tet_adjacent &fa,
    const zjucad::matrix::matrix<size_t> &tri_idx,
    const zjucad::matrix::matrix<double> &normal,
    const zjucad::matrix::matrix<zjucad::matrix::matrix<double> > &deformation_gradient,
    const zjucad::matrix::matrix<double> &area,
    const double w_normal_align,
    const bool only_basis_4 = true,
    const string opt_type = "init")
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  matrix<double> rot(9,9);
  double v = 0;
  matrix<double> deformed_normal(3,1), zyz(3,1);
  matrix<double> grad;
  for(size_t fi = 0; fi < tri_idx.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[tri_idx[fi]];
      assert(fa.is_outside_face(tet_pair));
      const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);
      grad = deformation_gradient[tet_idx];
      inv(grad);
      deformed_normal = trans(grad) * normal(colon(),fi);
      deformed_normal/=norm(deformed_normal);
      rot_n_2_z_by_zyz(&deformed_normal[0], &zyz[0]);
      calc_rot_cubic_f_sh_mat_(&rot[0], &zyz[0]);

      for(size_t di = 1; di < 8; ++di){
          if(di == 4) {v = sqrt(7.0);
            } else {
              if(only_basis_4) continue;
              v = 0;
            }
          if(opt_type == "init")
            all_func->push_back(math_func_ptr(
                                  new surface_normal_align2_sh<double,int32_t>(
                                    tet_mesh.size(2), tet_idx, rot, di, v,
                                    sqrt(w_normal_align * area[fi]/total_area))));
          else if(opt_type == "zyz")
            all_func->push_back(math_func_ptr(
                                  new surface_normal_align2_zyz<double,int32_t>(
                                    tet_mesh.size(2), tet_idx, rot, di, v,
                                    sqrt(w_normal_align * area[fi]/total_area))));
          else throw std::invalid_argument("# [error] can not recognize opt type.");
        }
    }
  cerr << "# [info] normal_align_func number " << all_func->size() << endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_normal_align_func5(const zjucad::matrix::matrix<size_t> &tet_mesh,
                         const jtf::mesh::face2tet_adjacent &fa,
                         const zjucad::matrix::matrix<size_t> &tri_idx,
                         const zjucad::matrix::matrix<double> &local2global,
                         const zjucad::matrix::matrix<double> &area,
                         const zjucad::matrix::matrix<double> &ks,
                         const double w_normal_align,
                         const bool only_basis_4 = true,
                         const string opt_type = "init")
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  matrix<double> rot(9,9);
  double v = 0;
  for(size_t fi = 0; fi < tri_idx.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[tri_idx[fi]];
      assert(fa.is_outside_face(tet_pair));
      const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);

      calc_rot_cubic_f_sh_mat_(&rot[0], &local2global(0,fi));
      rot = temp(trans(rot));

      for(size_t di = 1; di < 8; ++di){
          if(di == 4) {
              v = sqrt(7.0);
            } else {
              if(only_basis_4) continue;
              v = 0;
            }
          if(opt_type == "init")
            all_func->push_back(math_func_ptr(
                                  new surface_normal_align2_sh<double,int32_t>(
                                    tet_mesh.size(2), tet_idx, rot, di, v,
                                    sqrt(w_normal_align * area[fi]/total_area))));
          else if(opt_type == "zyz")
            all_func->push_back(math_func_ptr(
                                  new surface_normal_align2_zyz<double,int32_t>(
                                    tet_mesh.size(2), tet_idx, rot, di, v,
                                    sqrt(w_normal_align * area[fi]/total_area))));
          else throw std::invalid_argument("# [error] can not recognize opt type.");
        }
      if(opt_type == "init")
        all_func->push_back(math_func_ptr(
                              new surface_cross_norm_sh<double,int32_t>(
                                tet_mesh.size(2), tet_idx, rot,
                                sqrt(w_normal_align * area[fi]/total_area))));
      //      else if(opt_type == "zyz")
      //        all_func->push_back(math_func_ptr(
      //                              new surface_cross_norm_zyz<double,int32_t>(
      //                                tet_mesh.size(2), tet_idx, rot,
      //                                sqrt(w_normal_align * area[fi]/total_area))));
      //      else throw std::invalid_argument("# [error] can not recognize opt type.");
    }
  cerr << "# [info] normal_align_func number " << all_func->size() << endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

template <typename val_type, typename int_type>
class surface_norm_func : public hj::math_func::math_func_t<val_type,int_type>
{
public:
  surface_norm_func(const size_t tri_num, const size_t fi, const double w)
    :tri_num_(tri_num), fi_(fi), w_(w){}
  virtual ~surface_norm_func(){}
  virtual size_t nx() const{
    return 2*tri_num_;
  }
  virtual size_t nf() const{
    return 1;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==0){
        itr_matrix<const val_type*> x0(2, tri_num_,x);
        matrix<val_type> one_x = x0(colon(),fi_);
        int_type c[1] = {0};
        cv[c] += w_*dot(one_x,one_x)-w_;
      }
    if(k==1){
        itr_matrix<const val_type*> x0(2, tri_num_,x);
        for(int_type i = 0; i < 2; ++i){
            int_type c[2] = {0,2*fi_+i};
            cv[c] += 2*w_*x0(i,fi_);
          }
      }
    if(k==2){
        for(int_type i = 0; i < 2; ++i){
            int_type c[3] = {0,2*fi_+i,2*fi_+i};
            cv[c] += 2*w_;
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k==1){
        for(int_type i = 0; i < 2; ++i){
            int_type c[2] = {0,2*fi_+i};
            l2g.add(cs,c);
          }
      }
    if(k==2){
        for(int_type i = 0; i < 2; ++i){
            int_type c[3] = {0,2*fi_+i,2*fi_+i};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k == 0) return -1;
    if(k == 1) return 2;
    if(k == 2) return 2;
  }

private:
  const size_t tri_num_;
  const int_type fi_;
  const double w_;
};

shared_ptr<const hj::math_func::math_func_t<double,int32_t> >
build_surface_constraint_func3(
    const zjucad::matrix::matrix<size_t> &tri_mesh,
    const zjucad::matrix::matrix<double> &area,
    const double w)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  for(size_t fi = 0; fi < tri_mesh.size(2); ++fi){
      all_func->push_back(math_func_ptr(
                            new surface_norm_func<double,int32_t>(
                              tri_mesh.size(2), fi, sqrt(w*area[fi]/total_area))));
    }

  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

template<typename val_type, typename int_type>
class sh_align_surface_field_func: public hj::math_func::math_func_t<val_type, int_type>
{
public:
  sh_align_surface_field_func(const size_t tet_num, const size_t tet_idx,
                              const zjucad::matrix::matrix<val_type> &surface_field,
                              const zjucad::matrix::matrix<val_type> &local2global,
                              const val_type w)
    :tet_num_(tet_num), tet_idx_(tet_idx), surface_field_(surface_field),
      global2local_(trans(local2global)),w_(w){
    ref_ = zeros<val_type>(9,1);
    const double len = norm(surface_field);
    ref_[4] = sqrt(7.0/5.0)*len;
    ref_[0] = surface_field[1];
    ref_[8] = surface_field[0];
  }
  virtual ~sh_align_surface_field_func(){}
  virtual size_t nx() const{
    return 9*tet_num_;
  }
  virtual size_t nf() const{
    return 9;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const{
    if(k == 0){
        itr_matrix<const val_type*> x0(9,tet_num_, x);
        matrix<val_type> diff = global2local_ * x0(colon(), tet_idx_) - ref_;
        for(size_t i = 0; i < 9 ; ++i){
            int_type c[1] = {i};
            cv[c] += w_ * diff[i];
          }
      }
    if(k==1){
        for(size_t i =  0; i < 9; ++i)
          for(size_t j = 0; j < 9; ++j){
              int_type c[2] = {i, 9 * tet_idx_ + j};
              cv[c] += w_*global2local_(i,j);
            }
      }
    return 0;
  }
  virtual int patt(size_t k , hj::math_func::coo_set<int_type> & cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx * ctx = 0)const {
    if(k == 1){
        for(size_t i = 0; i < 9; ++i)
          for(size_t j = 0; j < 9; ++j){
              int_type c[2] = {i, 9 * tet_idx_+j};
              l2g.add(cs,c);
            }
      }
    return 0;
  }
  virtual size_t nnz(size_t k)const{
    if(k == 0) return -1;
    if(k == 1) return 9*9;
  }
private:
  const size_t tet_num_;
  const size_t tet_idx_;
  const zjucad::matrix::matrix<val_type> surface_field_;
  const zjucad::matrix::matrix<val_type> global2local_;
  zjucad::matrix::matrix<val_type> ref_;
  const val_type w_;
};

template<typename val_type, typename int_type>
class sh_align_surface_field_func2 : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  sh_align_surface_field_func2(const size_t tet_num, const size_t tet_idx,
                               const zjucad::matrix::matrix<val_type> & surface_field,
                               const zjucad::matrix::matrix<val_type> & local2global,
                               const val_type w, const val_type s = 1)
    :tet_num_(tet_num), tet_idx_(tet_idx), surface_field_(surface_field), w_(w){
    ref_ = zeros<val_type>(9,1);
    const double len = norm(surface_field);
    ref_[4] = sqrt(7.0/5.0)*len;
    ref_[0] = surface_field[1];
    ref_[8] = surface_field[0];
    ref_ = temp(local2global * ref_);
    s_v_ = {1,1,1,1,1,1,1,1,1};
  }
  virtual ~sh_align_surface_field_func2(){}
  virtual size_t nx()const{
    return 9 * tet_num_;
  }
  virtual size_t nf()const{
    return 9;
  }
  virtual int eval(size_t k , const val_type *x,
                   const hj::math_func::coo2val_t<val_type, int_type> & cv,
                   hj::math_func::func_ctx * ctx = 0 )const{
    if(k==0){
        itr_matrix<const val_type*> x0(9, tet_num_, &x[0]);
        matrix<val_type> diff = x0(colon(), tet_idx_) - ref_;
        for(int_type i = 0; i < 9; ++i){
            int_type c[1] = {i};
            cv[c] += w_*(diff[i]) * s_v_[i];
          }
      }
    if(k==1){
        for(int_type i = 0; i < 9; ++i){
            int_type c[2] = {i,9*tet_idx_+i};
            cv[c] += w_*s_v_[i];
          }
      }
    return 0;
  }

  virtual int patt(size_t k, hj::math_func::coo_set<int_type> & cs,
                   const hj::math_func::coo_l2g & l2g,
                   hj::math_func::func_ctx * ctx = 0)const{
    if(k == 1){
        for(int_type i = 0; i < 9; ++i){
            int_type c[2] = {i,9*tet_idx_+i};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k)const{
    if(k == 0) return -1;
    if(k == 1) return 9;
  }
private:
  const int_type tet_num_;
  const int_type tet_idx_;
  const zjucad::matrix::matrix<val_type> surface_field_;
  zjucad::matrix::matrix<val_type> ref_;
  const val_type w_;
  vector<val_type> s_v_;
};

shared_ptr<const hj::math_func::math_func_t<double,int32_t> >
build_surface_constraint_func4(
    const zjucad::matrix::matrix<size_t> &tet_mesh,
    const zjucad::matrix::matrix<size_t> &tri_mesh_idx,
    const zjucad::matrix::matrix<double> & local2global,
    const zjucad::matrix::matrix<double> &area,
    const jtf::mesh::face2tet_adjacent & fa,
    const zjucad::matrix::matrix<double> &surface_field,
    const double w, const double s)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  matrix<double> rot(9,9);

  itr_matrix<const double*> surface_field_m(2, tri_mesh_idx.size(), &surface_field[0]);
  for(size_t fi = 0; fi < tri_mesh_idx.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[tri_mesh_idx[fi]];
      const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);
      calc_rot_cubic_f_sh_mat_(&rot[0], &local2global(0,fi));
      const double weight = sqrt(w*area[fi]/total_area);

      all_func->push_back(math_func_ptr(
                            new sh_align_surface_field_func2<double,int32_t>(
                              tet_mesh.size(2), tet_idx, surface_field_m(colon(),fi), rot,
                              weight, s)));
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

template <typename val_type, typename int_type>
class fix_surface_field_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  fix_surface_field_func(const size_t tet_num, const size_t variable_idx,
                         const val_type target, const val_type weight)
    :tet_num_(tet_num), variable_idx_(variable_idx), target_(target), weight_(weight){}
  virtual ~fix_surface_field_func(){}
  virtual size_t nx() const{
    return 9*tet_num_;
  }
  virtual size_t nf()const{
    return 1;
  }
  virtual int eval(size_t k , const val_type *x,
                   const hj::math_func::coo2val_t<val_type, int_type> & cv,
                   hj::math_func::func_ctx * ctx = 0 ) const
  {
    if(k == 0)  {
        int_type c[1] = {0};
        cv[c] += weight_*(x[variable_idx_]-target_);
      }
    if(k == 1){
        int_type c[2] = {0, variable_idx_};
        cv[c] += weight_;
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> & cs,
                   const hj::math_func::coo_l2g & l2g,
                   hj::math_func::func_ctx * ctx = 0) const
  {
    if(k == 1){
        int_type c[2] = {0, variable_idx_};
        l2g.add(cs,c);
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k==0) return -1;
    if(k==1) return 1;
  }
private:
  const size_t tet_num_;
  const size_t variable_idx_;
  const val_type target_;
  const val_type weight_;
};

shared_ptr<const hj::math_func::math_func_t<double,int32_t> >
build_surface_constraint_func5(
    const zjucad::matrix::matrix<size_t> &tet_mesh,
    const zjucad::matrix::matrix<size_t> &tri_mesh_idx,
    const zjucad::matrix::matrix<double> & local2global,
    const zjucad::matrix::matrix<double> &area,
    const jtf::mesh::face2tet_adjacent & fa,
    const zjucad::matrix::matrix<double> &surface_field,
    const double w, const double s, const string weight_scheme = "w1")
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  matrix<double> rot(9,9);

  itr_matrix<const double*> surface_field_m(2, tri_mesh_idx.size(), &surface_field[0]);
  vector<double> weight_s = {1+s,0,0,0,1-s,0,0,0,1+s};
  if(weight_scheme == "w2"){
      weight_s[1] = 1;
      weight_s[2] = 1;
      weight_s[3] = 1;
      weight_s[5] = 1;
      weight_s[6] = 1;
      weight_s[7] = 1;
    }

  zjucad::matrix::matrix<double> ref = zeros<double>(9,1);
  matrix<double> one_surface_field;
  for(size_t fi = 0; fi < tri_mesh_idx.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[tri_mesh_idx[fi]];
      const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);
      calc_rot_cubic_f_sh_mat_(&rot[0], &local2global(0,fi));
      rot = temp(trans(rot));
      //const double weight = sqrt(w*area[fi]/total_area);
      const double weight = w*area[fi]/total_area;
      one_surface_field = surface_field_m(colon(),fi);
      const double len = norm(one_surface_field);
      ref *= 0;
      ref[4] = sqrt(7.0/5.0)*len;
      ref[0] = surface_field_m(1,fi);
      ref[8] = surface_field_m(0,fi);

      for(size_t j = 0; j< 9; ++j){
          if(fabs(weight_s[j]) > 1e-8){
              all_func->push_back(math_func_ptr(
                                    new surface_normal_align2_pow2<double,int32_t>(
                                      tet_mesh.size(2), tet_idx,
                                      rot, j, ref[j],weight*weight_s[j])));
            }
        }
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

void optimal_sh_generator::opt(zjucad::matrix::matrix<double> &sh,
                               boost::property_tree::ptree &pt)
{
  const int opt_num = pt.get<int>("input/opt_num.value",7);
  switch (opt_num)
    {
    case 7:   opt7(sh,pt);  break;
    case 13:  opt13(sh,pt); break;
    default:  cerr << "# [error] wrong opt_num." << endl;  break;
    }
}

static void convert_local_rot_2_global_rot(zjucad::matrix::matrix<double> & zyz,
                                           const zjucad::matrix::matrix<size_t> & trimesh,
                                           const zjucad::matrix::matrix<double> & node,
                                           const zjucad::matrix::matrix<double> & face_normal)
{
  assert(zyz.size(2) == trimesh.size(2));
  zjucad::matrix::matrix<double> rz2n(3,1);
  matrix<double> rot(3,3), rz2n_mat(3,3), rot_ref = eye<double>(3);
  matrix<double> error(3,1);
  for(size_t fi = 0; fi < face_normal.size(2); ++fi){
      rot_n_2_z_by_zyz(&face_normal(0,fi), &rz2n[0]);
      swap(rz2n[0],rz2n[2]);
      rz2n *= -1;
      assert(fabs(zyz(1,fi)) < 1e-6);
      from_angle_to_rotation_matrix(zyz(0,fi)+zyz(2,fi), face_normal(colon(),fi), rot);
      zyz_angle_2_rotation_matrix1(&rz2n[0], &rz2n_mat[0]);
      cal_ref_align_rot(trimesh(colon(),fi),node,rz2n_mat, rot_ref);
      rot = temp(rot*rot_ref*rz2n_mat);
      rotation_matrix_2_zyz_angle(&rot[0], &zyz(0,fi),&error[0]);
    }
}

static void convert_local_rot_2_global_rot(zjucad::matrix::matrix<double> & zyz,
                                           jtf::mesh::tri_mesh &tri_mesh)
{
  convert_local_rot_2_global_rot(zyz, tri_mesh.trimesh_.mesh_, tri_mesh.trimesh_.node_,
                                 tri_mesh.face_normal_);
}

void optimal_sh_generator::opt6(zjucad::matrix::matrix<double> &sh,
                                boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  zjucad::matrix::matrix<double> tradition_cos_sin;

  {
    string input_cos_sin_file;
    if(zjucad::has("input/cos_sin_mat.value",pt)){
        input_cos_sin_file = pt.get<string>("input/cos_sin_mat.value");
      }
    ifstream ifs(input_cos_sin_file);
    if(ifs.fail()){
        cerr << "# [warning] can not read surface field " << input_cos_sin_file << endl;
        jtf::mesh::meshes triangle_mesh;
        triangle_mesh.mesh_ = tm_.outside_face_;
        triangle_mesh.node_ = tm_.tetmesh_.node_;
        jtf::mesh::tri_mesh tri_mesh(triangle_mesh);
        optimal_surface_sh_generator ossg(tri_mesh);
        ossg.opt_cos_sin(tradition_cos_sin,pt);
        jtf::mesh::write_matrix("cos_sin.mat", tradition_cos_sin);
        cerr << "# [info] generate surface field cos_sin.mat " << endl;
      }else{
        jtf::mesh::read_matrix(input_cos_sin_file.c_str(), tradition_cos_sin);
        cerr << "# [info] read surface field " << input_cos_sin_file << endl;
      }
  }

  init(pt);

  pt.put("weight/need_normalized.desc","need normalized or not [yes/no]");
  pt.put("weight/s.desc", "to control the surface field alignment,[-1,0,1]");
  const string need_normalized = pt.get<string>("weight/need_normalized.value");
  const double s=pt.get<double>("weight/s.value");
  itr_matrix<double*> cos_sin_m(2, tradition_cos_sin.size()/2, &tradition_cos_sin[0]);
  if(need_normalized == "YES" | need_normalized == "yes" |
     need_normalized == "Y" | need_normalized == "y"){
      for(size_t fi = 0; fi < cos_sin_m.size(2); ++fi){
          cos_sin_m(colon(),fi) /= norm(cos_sin_m(colon(),fi));
        }
      cerr << "# [info] normalized surface field." << endl;
    }else
    cerr << "# [info] do not normalized surface field" << endl;

  pt.put("weight/relationship.desc","weight to constraint surface and inner relationship");
  const double w_relationship = pt.get<double>("weight/relationship.value");

  pt.put("input/smooth_scheme.desc", "inner field smooth scheme [face/gradient]");
  const string smooth_scheme = pt.get<string>("input/smooth_scheme.value");
  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);
  {
    func->push_back(math_func_ptr(
                      new hj::math_func::sumsqr<double,int32_t>(
                        build_inner_smooth_func(tm_.tetmesh_.mesh_, tm_.tetmesh_.node_,
                                                *tm_.fa_, tm_.ortae_, tm_.vol_, 1.0, smooth_scheme))));
    cerr << "# [info] add inner smooth constraints. smooth scheme " << smooth_scheme << endl;
  }

  shared_ptr<vector<math_func_ptr> > eqns(new vector<math_func_ptr>);

  pt.put("input/weight_scheme.desc","[w1/w2]");
  const string weight_scheme = pt.get<string>("input/weight_scheme.value");
  if(w_relationship > 1e-8){
      eqns->push_back(math_func_ptr(
                        new hj::math_func::sum<double,int32_t>(
                          build_surface_constraint_func5(
                            tm_.tetmesh_.mesh_, tm_.outside_face_idx_,
                            local2global_rot_, tm_.outside_face_area_, *tm_.fa_,
                            tradition_cos_sin, w_relationship,s,weight_scheme))));
      //                                    new hj::math_func::sumsqr<double,int32_t>(
      //                                      build_surface_constraint_func4(
      //                                        tm_.tetmesh_.mesh_, tm_.outside_face_idx_,
      //                                        local2global_rot_, tm_.outside_face_area_, *tm_.fa_,
      //                                        tradition_cos_sin, w_relationship,s))));
      cerr << "# [info] add surface constraint function." << endl;
      cerr << "# [info] weight scheme " << weight_scheme << endl;
    }
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  for(size_t i = 0; i < func->size(); ++i) all_func->push_back((*func)[i]);
  for(size_t i = 0; i < eqns->size(); ++i) all_func->push_back((*eqns)[i]);

  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(all_func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));
  //shared_ptr<func2opt> jtf_func(new func2opt(obj));

  sh =zeros<double>(9,tm_.tetmesh_.mesh_.size(2));
  {
    sh(4,colon()) = sqrt(7.0);
    sh(8,colon()) = sqrt(5.0);
  }
  jtf::optimize(*obj, sh, pt, nullptr, nullptr, nullptr);

  for(size_t i = 0; i < sh.size(2); ++i){
      sh(colon(),i) /= norm(sh(colon(),i));
    }
  sh *= sqrt(12.0);
}

void optimal_sh_generator::opt5(zjucad::matrix::matrix<double> &sh,
                                boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  init(pt);

  pt.put("weight/surface_smooth.desc", "surface smoothing weight.");
  pt.put("weight/normal_align.desc","weight to add surface normal alignment");
  pt.put("input/smooth_scheme.desc", "input smooth scheme strategy [face/gradient]");

  const double w_surface = pt.get<double>("weight/surface_smooth.value");
  const double w_normal_align = pt.get<double>("weight/normal_align.value");
  const string smooth_scheme = pt.get<string>("input/smooth_scheme.value");

  matrix<double> ks = ones<double>(tm_.outside_face_.size(2),1);

  sh = zeros<double>(9,tm_.tetmesh_.mesh_.size(2));
  sh(4,colon()) += sqrt(7.0);
  sh(8,colon()) += sqrt(5.0);
  if(zjucad::has("input/init.value" ,pt)){
      matrix<double> zyz;
      jtf::mesh::read_matrix(pt.get<string>("input/init.value").c_str(), zyz);
      if(zyz.size(2) != tm_.tetmesh_.mesh_.size(2)){
          cerr << "# [error] wrong zyz file." << endl;
          return ;
        }
      size_t i = 0;
#pragma omp parallel for private(i)
      for(i = 0; i < tm_.tetmesh_.mesh_.size(2); ++i){
          calc_rot_cubic_f_sh_(&sh(0,i),&zyz(0,i));
        }
      cerr << "# [info] load init zyz." << endl;
    }

  pt.put("weight/iteration.desc", "iterations for two step optimization");
  const size_t iteration_num = pt.get<size_t>("weight/iteration.value",5);
  ofstream ofs("ks");
  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);
  for(size_t it = 0; it < iteration_num; ++it){
      cerr << "# iter " << it << endl;
      func->clear();
      {
        func->push_back(math_func_ptr(
                          new hj::math_func::sumsqr<double,int32_t>(
                            build_inner_smooth_func<double>(tm_.tetmesh_.mesh_, tm_.tetmesh_.node_,
                                                    *tm_.fa_, tm_.ortae_, tm_.vol_, 1, smooth_scheme))));
        cerr << "# [info] add inner smooth function." << endl;
      }
      if(w_surface > 1e-8){
          if(zjucad::has("input/cos_sin.value",pt)){
              matrix<double> tradition;
              jtf::mesh::read_matrix(pt.get<string>("input/cos_sin.value").c_str(),tradition);
              func->push_back(math_func_ptr(
                                new hj::math_func::sumsqr<double,int32_t>(
                                  build_surface_fix_func2(tm_.tetmesh_.mesh_, tm_.outside_face_idx_,tm_.outside_face_area_,
                                                          *tm_.fa_, local2global_rot_, tradition, w_surface))));
              cerr << "# [info] add surface fix function: " << w_surface << endl;
            }else{
              func->push_back(math_func_ptr(
                                new hj::math_func::sumsqr<double,int32_t>(
                                  build_surface_smooth_func4(tm_.tetmesh_.mesh_, tm_.outside_face_idx_, *tm_.ea_outside_,
                                                             *tm_.fa_, tm_.outside_face_area_,rij_,
                                                             local2global_rot_,ks, w_surface))));
              cerr << "# [info] add surface smooth function: " << w_surface << endl;
            }
        }
      if(w_normal_align > 1e-8){
          func->push_back(math_func_ptr(
                            new hj::math_func::sumsqr<double,int32_t>(
                              build_normal_align_func4(tm_.tetmesh_.mesh_, *tm_.fa_, tm_.outside_face_idx_, local2global_rot_,
                                                       tm_.outside_face_area_,  ks, w_normal_align, true))));
          cerr << "# [info] add normal align function: " << w_normal_align << endl;
        }

      math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
      math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

      jtf::optimize(*obj, sh, pt, nullptr, nullptr, nullptr);

      cerr << "# average ks " << std::accumulate(ks.begin(), ks.end(),0.0)/ks.size() << endl;
      ofs << std::accumulate(ks.begin(), ks.end(),0.0)/ks.size() << endl;
      for(size_t ti = 0; ti < tm_.outside_face_idx_.size(); ++ti){
          const pair<size_t,size_t> & tet_pair = tm_.fa_->face2tet_[tm_.outside_face_idx_[ti]];
          const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);
          const double len = norm(sh(colon(),tet_idx));
          if(len > 1e-18)
            ks[ti] = sqrt(12.0)/len;
        }
    }

  cerr << sh(colon(),colon(0,5)) << endl;
  for(size_t i = 0; i < sh.size(2); ++i){
      sh(colon(),i) /= norm(sh(colon(),i));
    }
  sh *= sqrt(12.0);
}

inline double face_len(const size_t face_idx0, const size_t face_idx1,
                       const zjucad::matrix::matrix<double> & node,
                       const vector<vector<size_t> > &faces)
{
  itr_matrix<const size_t*> face0(3,1, &faces[face_idx0][0]);
  itr_matrix<const size_t*> face1(3,1, &faces[face_idx1][0]);
  matrix<double> diff = node(colon(), face0) - node(colon(), face1);
  return norm(diff);
}

void optimal_sh_generator::compute_dis2surface(zjucad::matrix::matrix<double> &sf,
                                               boost::property_tree::ptree &pt)
{
  cerr << "# [info] using size_field type : dis2surfaces." << endl;
  if(sf.size() != tm_.tetmesh_.node_.size(2))
    sf.resize(tm_.tetmesh_.node_.size(2),1);

  sf *= 0;
  matrix<double> fsf = -1*ones<double>(tm_.fa_->faces_.size(),1);
  deque<size_t> que;
  for(size_t fi = 0; fi < tm_.outside_face_idx_.size(); ++fi){
      que.push_back(tm_.outside_face_idx_[fi]);
      fsf[tm_.outside_face_idx_[fi]] = 0;
    }
  vector<bool> visited_face_flag(tm_.fa_->faces_.size(), false);
  while(!que.empty()){
      const size_t face_idx = que.front();
      que.pop_front();
      visited_face_flag[face_idx] = true;
      const pair<size_t,size_t> & tet_pair = tm_.fa_->face2tet_[face_idx];
      vector<size_t> tet_pair_candidate;

      if(tet_pair.first != -1) tet_pair_candidate.push_back(tet_pair.first);
      if(tet_pair.second != -1) tet_pair_candidate.push_back(tet_pair.second);

      for(const auto & tet_idx : tet_pair_candidate){
          for(size_t i = 0; i < tm_.tetmesh_.mesh_.size(1); ++i){
              const size_t other_face_idx =
                  tm_.fa_->get_face_idx(tm_.tetmesh_.mesh_(i, tet_idx),
                                        tm_.tetmesh_.mesh_((i+1)%4, tet_idx),
                                        tm_.tetmesh_.mesh_((i+2)%4, tet_idx));
              if(visited_face_flag[other_face_idx]) continue;
              fsf[other_face_idx] = fsf[face_idx] +
                  face_len(other_face_idx, face_idx, tm_.tetmesh_.node_, tm_.fa_->faces_);
              que.push_back(other_face_idx);
            }
        }
    }

  { // convert face size to point;
    vector<vector<double> > point_size(tm_.tetmesh_.node_.size(2));
    for(size_t fi = 0; fi < tm_.fa_->faces_.size(); ++fi){
        const vector<size_t> & one_face = tm_.fa_->faces_[fi];
        for(size_t i = 0 ; i < one_face.size(); ++i){
            point_size[one_face[i]].push_back(fsf[fi]);
          }
      }
    for(size_t fi = 0; fi < point_size.size(); ++fi){
        sf[fi] = std::accumulate(point_size[fi].begin(), point_size[fi].end(), 0.0)/point_size[fi].size();
      }
  }
}

void optimal_sh_generator::compute_dis2x(zjucad::matrix::matrix<double> &sf,
                                         boost::property_tree::ptree &pt)
{
  cerr << "# [info] using size_field type : dis2x." << endl;
  if(sf.size() != tm_.tetmesh_.node_.size(2))
    sf.resize(tm_.tetmesh_.node_.size(2),1);

  matrix<double> bb(3,2);
  calc_bounding_box(tm_.tetmesh_.node_, &bb[0]);
  const double min_x = bb(0,0), max_x = bb(0,1);

  const double bound = 2*My_PI();
  for(size_t pi = 0 ; pi < tm_.tetmesh_.node_.size(2); ++pi){
      double dis_x = fabs(tm_.tetmesh_.node_(0,pi));
      sf[pi] = 10*((sin((dis_x-min_x)*bound/(max_x-min_x)-My_PI()/2.0)+2));//+(sin((dis_z-min_z)*bound/(max_z-min_z)-My_PI()/2.0)+2));
    }

  const double max = zjucad::matrix::max(sf);
  const double min = zjucad::matrix::min(sf);

  const double bound_max = 100;
  const double bound_min = 1;
  for(size_t pi = 0; pi < sf.size(); ++pi){
      sf[pi] = (sf[pi]-min)/(max-min)*(bound_max-bound_min)+bound_min;
    }

  {
    matrix<double> tet_size = zeros<double>(tm_.tetmesh_.mesh_.size(2),1);
    for(size_t ti = 0; ti < tm_.tetmesh_.mesh_.size(2); ++ti){
        for(size_t pi = 0; pi < tm_.tetmesh_.mesh_.size(1); ++pi){
            tet_size[ti] += sf[tm_.tetmesh_.mesh_(pi,ti)];
          }
      }
    tet_size /= 4;
    ofstream ofs("tet_size.vtk");
    tet2vtk(ofs, &tm_.tetmesh_.node_[0], tm_.tetmesh_.node_.size(2),
        &tm_.tetmesh_.mesh_[0], tm_.tetmesh_.mesh_.size(2));
    cell_data(ofs, &tet_size[0], tet_size.size(), "tet_size");
  }
}

void optimal_sh_generator::compute_dis2param_edge(zjucad::matrix::matrix<double> &sf,
                                                  boost::property_tree::ptree &pt)
{
  cerr << "# [info] using size_field type : dis2param_edge." << endl;

  jtf::tet_mesh tm_param(pt.get<string>("input/param_tet.value").c_str());

  if(sf.size() != tm_.tetmesh_.node_.size(2))
    sf.resize(tm_.tetmesh_.node_.size(2),1);

  vector<vector<double> > edge_ratio(sf.size());
  for(const auto & one_edge : tm_param.ortae_.e2t_){
      const pair<size_t,size_t> & edge_point = one_edge.first;
      const double ratio = norm(tm_param.tetmesh_.node_(colon(),edge_point.first)-tm_param.tetmesh_.node_(colon(),edge_point.second))
          /norm(tm_.tetmesh_.node_(colon(),edge_point.first)-tm_.tetmesh_.node_(colon(),edge_point.second));
      edge_ratio[edge_point.first].push_back(ratio);
      edge_ratio[edge_point.second].push_back(ratio);
    }
  for(size_t pi = 0; pi < edge_ratio.size(); ++pi){
      sf[pi] = std::accumulate(edge_ratio[pi].begin(), edge_ratio[pi].end(),0.0)/edge_ratio[pi].size();
    }
  {
    matrix<double> tet_size = zeros<double>(tm_.tetmesh_.mesh_.size(2),1);
    for(size_t ti = 0; ti < tm_.tetmesh_.mesh_.size(2); ++ti){
        for(size_t pi = 0; pi < tm_.tetmesh_.mesh_.size(1); ++pi){
            tet_size[ti] += sf[tm_.tetmesh_.mesh_(pi,ti)];
          }
      }
    tet_size /= 4;
    ofstream ofs("tet_size.vtk");
    tet2vtk(ofs, &tm_.tetmesh_.node_[0], tm_.tetmesh_.node_.size(2),
        &tm_.tetmesh_.mesh_[0], tm_.tetmesh_.mesh_.size(2));
    cell_data(ofs, &tet_size[0], tet_size.size(), "tet_size");
  }
}
void optimal_sh_generator::compute_dis2param_vol(zjucad::matrix::matrix<double> &sf,
                                                 boost::property_tree::ptree &pt)
{
  cerr << "# [info] using size_field type : dis2param_vol." << endl;
  jtf::tet_mesh tm_param(pt.get<string>("input/param_tet.value").c_str());

  if(sf.size() != tm_.tetmesh_.node_.size(2))
    sf.resize(tm_.tetmesh_.node_.size(2),1);

  vector<vector<double> > vol_ratio(tm_.tetmesh_.node_.size(2));
  for(size_t ti = 0; ti < tm_.tetmesh_.mesh_.size(2); ++ti){
      for(size_t pi = 0; pi < tm_.tetmesh_.mesh_.size(1); ++pi)
        vol_ratio[tm_.tetmesh_.mesh_(pi,ti)].push_back(pow(fabs(tm_param.vol_[ti]/tm_.vol_[ti]),1.0/3));
    }
  for(size_t pi = 0; pi < sf.size(); ++pi){
      sf[pi] = std::accumulate(vol_ratio[pi].begin(), vol_ratio[pi].end(), 0.0)/vol_ratio[pi].size();
    }
  {
    matrix<double> tet_size = zeros<double>(tm_.tetmesh_.mesh_.size(2),1);
    for(size_t ti = 0; ti < tm_.tetmesh_.mesh_.size(2); ++ti){
        for(size_t pi = 0; pi < tm_.tetmesh_.mesh_.size(1); ++pi){
            tet_size[ti] += sf[tm_.tetmesh_.mesh_(pi,ti)];
          }
      }
    tet_size /= 4;
    ofstream ofs("tet_size.vtk");
    tet2vtk(ofs, &tm_.tetmesh_.node_[0], tm_.tetmesh_.node_.size(2),
        &tm_.tetmesh_.mesh_[0], tm_.tetmesh_.mesh_.size(2));
    cell_data(ofs, &tet_size[0], tet_size.size(), "tet_size");
  }
}

void optimal_sh_generator::compute_disrand(zjucad::matrix::matrix<double> &sf,
                                           boost::property_tree::ptree &pt)
{
  cerr << "# [info] using size_field type : disrand." << endl;
  if(sf.size() != tm_.tetmesh_.node_.size(2))
    sf.resize(tm_.tetmesh_.node_.size(2),1);

  sf = fabs(zjucad::matrix::rand<double>(tm_.tetmesh_.node_.size(2),1)) + 1;

  const double max = zjucad::matrix::max(sf);
  const double min = zjucad::matrix::min(sf);

  const double bound_max = 100;
  const double bound_min = 1;
  for(size_t pi = 0; pi < sf.size(); ++pi){
      sf[pi] = (sf[pi]-min)/(max-min)*(bound_max-bound_min)+bound_min;
    }

  {
    matrix<double> tet_size = zeros<double>(tm_.tetmesh_.mesh_.size(2),1);
    for(size_t ti = 0; ti < tm_.tetmesh_.mesh_.size(2); ++ti){
        for(size_t pi = 0; pi < tm_.tetmesh_.mesh_.size(1); ++pi){
            tet_size[ti] += sf[tm_.tetmesh_.mesh_(pi,ti)];
          }
      }
    tet_size /= 4;
    ofstream ofs("tet_size.vtk");
    tet2vtk(ofs, &tm_.tetmesh_.node_[0], tm_.tetmesh_.node_.size(2),
        &tm_.tetmesh_.mesh_[0], tm_.tetmesh_.mesh_.size(2));
    cell_data(ofs, &tet_size[0], tet_size.size(), "tet_size");
  }
}
void optimal_sh_generator::compute_dis2orig(zjucad::matrix::matrix<double> &sf,
                                            boost::property_tree::ptree &pt)
{
  cerr << "# [info] using size_field type : dis2orig." << endl;
  if(sf.size() != tm_.tetmesh_.node_.size(2))
    sf.resize(tm_.tetmesh_.node_.size(2),1);

  for(size_t pi = 0 ; pi < tm_.tetmesh_.node_.size(2); ++pi){
      sf[pi] = norm(tm_.tetmesh_.node_(colon(),pi));
    }

  const double max = zjucad::matrix::max(sf);
  const double min = zjucad::matrix::min(sf);

  const double bound_max = 100;
  const double bound_min = 1;
  for(size_t pi = 0; pi < sf.size(); ++pi){
      sf[pi] = (sf[pi]-min)/(max-min)*(bound_max-bound_min)+bound_min;
    }

  {
    matrix<double> tet_size = zeros<double>(tm_.tetmesh_.mesh_.size(2),1);
    for(size_t ti = 0; ti < tm_.tetmesh_.mesh_.size(2); ++ti){
        for(size_t pi = 0; pi < tm_.tetmesh_.mesh_.size(1); ++pi){
            tet_size[ti] += sf[tm_.tetmesh_.mesh_(pi,ti)];
          }
      }
    tet_size /= 4;
    ofstream ofs("tet_size.vtk");
    tet2vtk(ofs, &tm_.tetmesh_.node_[0], tm_.tetmesh_.node_.size(2),
        &tm_.tetmesh_.mesh_[0], tm_.tetmesh_.mesh_.size(2));
    cell_data(ofs, &tet_size[0], tet_size.size(), "tet_size");
  }
}

void optimal_sh_generator::compute_size_field(zjucad::matrix::matrix<double> &sf,
                                              boost::property_tree::ptree &pt)
{
  pt.put("input/size_field_type.desc", "[dis2surface/dis2x]");
  const string size_field_strategy = pt.get<string>("input/size_field_type.value");
  if(size_field_strategy == "dis2surface"){
      compute_dis2surface(sf,pt);
    }else if(size_field_strategy == "dis2x"){
      compute_dis2x(sf,pt);
    }else if(size_field_strategy == "dis2orig"){
      compute_dis2orig(sf,pt);
    }else if(size_field_strategy == "disrand"){
      compute_disrand(sf,pt);
    }else if(size_field_strategy == "dis2param_vol"){
      compute_dis2param_vol(sf,pt);
    }else if(size_field_strategy == "dis2param_edge"){
      compute_dis2param_edge(sf,pt);
    }else throw std::invalid_argument("# [error] can not regonize sf type.");

  //adjust_size_field(sf);
}


void optimal_sh_generator::compute_inner_connection(
    zjucad::matrix::matrix<double> &sf,
    std::map<std::pair<size_t,size_t>, zjucad::matrix::matrix<double> > &inner_face_connection)
{
  // input size field in points, output connection in faces

  fxz::connection_generator cg(tm_, sf, inner_face_connection);
  cg.run();
  cerr << "# [info] inner_face_connection number " << inner_face_connection.size() << endl;
  cerr << "# [info] inner edge number " << tm_.ortae_.e2t_.size() << endl;
}

void check_param_tet_pair_coord(const zjucad::matrix::matrix<size_t> &five_points,
                                const zjucad::matrix::matrix<double> &node,
                                const map<pair<size_t,size_t>,double> &edge_ratio,
                                zjucad::matrix::matrix<double> &new_node)
{
  for(const auto & one_edge : edge_ratio){
      const pair<size_t,size_t> & edge_idx = one_edge.first;
      const double orig_len = norm(node(colon(), five_points[edge_idx.first]) - node(colon(), five_points[edge_idx.second]));
      const double new_len = norm(new_node(colon(), edge_idx.first) - new_node(colon(), edge_idx.second));
      if(fabs(new_len - orig_len * one_edge.second) > 1e-7)
        {
          cerr << "edge <" << edge_idx.first << "," << edge_idx.second << ">\t" <<
                  "orig_len \t" << orig_len << " scale \t" << one_edge.second << " new_len \t" << new_len << endl;
        }
    }
}

void  recover_param_tet_pair_coord(const zjucad::matrix::matrix<size_t> &five_points,
                                   const zjucad::matrix::matrix<double> &node,
                                   const map<pair<size_t,size_t>,double> &edge_ratio,
                                   zjucad::matrix::matrix<double> &new_node)
{
  new_node = zeros<double>(3,5); // nail down p0
  auto it01 = edge_ratio.find(make_pair(0,1));
  assert(it01 != edge_ratio.end());
  const double len01 = norm(node(colon(), five_points[0]) - node(colon(), five_points[1]));
  new_node(0,1) = len01*it01->second; // nail down p1

  auto it12 = edge_ratio.find(make_pair(1,2));
  auto it02 = edge_ratio.find(make_pair(0,2));
  const double new_len12 = norm(node(colon(), five_points[1]) - node(colon(), five_points[2])) * it12->second;
  const double new_len02 = norm(node(colon(), five_points[0]) - node(colon(), five_points[2])) * it02->second;
  new_node(0,2) = 0.5*(new_len02*new_len02-new_len12*new_len12)/new_node(0,1) + 0.5*new_node(0,1);
  new_node(1,2) = sqrt(new_len02*new_len02-new_node(0,2)*new_node(0,2)); // nail down p2

  auto it13 = edge_ratio.find(make_pair(1,3));
  auto it03 = edge_ratio.find(make_pair(0,3));
  auto it23 = edge_ratio.find(make_pair(2,3));
  const double new_len13 = norm(node(colon(), five_points[1]) - node(colon(), five_points[3])) * it13->second;
  const double new_len03 = norm(node(colon(), five_points[0]) - node(colon(), five_points[3])) * it03->second;
  const double new_len23 = norm(node(colon(), five_points[2]) - node(colon(), five_points[3])) * it23->second;
  new_node(0,3) = 0.5*(new_len03*new_len03-new_len13*new_len13)/new_node(0,1)+0.5*new_node(0,1);
  new_node(1,3) = 0.5*(new_len03*new_len03-new_node(0,3)*new_node(0,3)-new_len23*new_len23
                       +(new_node(0,3)-new_node(0,2))*(new_node(0,3)-new_node(0,2))+new_node(1,2)*new_node(1,2))/new_node(1,2);

  if(new_len03*new_len03 - new_node(0,3)*new_node(0,3) - new_node(1,3)*new_node(1,3) < 1e-6){
      throw std::invalid_argument("given size field does not satisfy triangle condition.");
    }
  new_node(2,3) = sqrt(new_len03*new_len03 - new_node(0,3)*new_node(0,3) - new_node(1,3)*new_node(1,3)); // nail down p3

  const double new_len31 = norm(new_node(colon(),3)-new_node(colon(),1));
  auto it34 = edge_ratio.find(make_pair(3,4));
  auto it24 = edge_ratio.find(make_pair(2,4));
  auto it14 = edge_ratio.find(make_pair(1,4));
  const double new_len34 = norm(node(colon(), five_points[3])-node(colon(), five_points[4]))*it34->second;
  const double new_len24 = norm(node(colon(), five_points[2])-node(colon(), five_points[4]))*it24->second;
  const double new_len14 = norm(node(colon(), five_points[1])-node(colon(), five_points[4]))*it14->second;

  const double vol = fabs(tet_vol(new_len12, new_len23, new_len31, new_len34, new_len14, new_len24));

  const double area123 = jtf::mesh::cal_face_area(new_len12, new_len23, new_len31);
  const double height = 3*vol/area123;

  if(new_len14 < height || new_len24 < height || new_len34 < height){
      throw std::invalid_argument("given size field does not satisfy triangle condition.");
    }

  const double new_len1c = sqrt(new_len14*new_len14-height*height);
  const double new_len2c = sqrt(new_len24*new_len24-height*height);
  const double new_len3c = sqrt(new_len34*new_len34-height*height);

  matrix<double> lambda(2,1),b(2,1);
  matrix<double> cc(2,2);

  cc(0,0) = new_len12*new_len12;
  cc(1,0) = dot(new_node(colon(),3)-new_node(colon(),1),new_node(colon(),2)-new_node(colon(),1));
  cc(0,1) = cc(1,0);
  cc(1,1) = new_len13*new_len13;
  b(0,0) = 0.5*(-new_len2c*new_len2c+new_len1c*new_len1c+new_len12*new_len12);
  b(1,0) = 0.5*(-new_len3c*new_len3c+new_len1c*new_len1c+new_len13*new_len13);

  inv(cc);
  lambda = cc*b;
  matrix<double> face123_center = lambda[0]*(new_node(colon(),2)-new_node(colon(),1))+
      lambda[1]*(new_node(colon(),3)-new_node(colon(),1))+new_node(colon(),1);

  matrix<double> face_normal = cross(new_node(colon(),2)-new_node(colon(),1),
                                     new_node(colon(),3)-new_node(colon(),2));
  face_normal /= norm(face_normal);
  face_normal *= height;
  new_node(colon(),4) = face123_center+face_normal; //nail down p4
}

void optimal_sh_generator::compute_inner_connection_jtf3(
    boost::property_tree::ptree &pt,
    zjucad::matrix::matrix<double> & point_size,
    map<std::pair<size_t,size_t>, matrix<double> > & connection_matrix)
{
  assert(point_size.size() == tm_.tetmesh_.node_.size(2));

  matrix<double> R(3,3);
  matrix<double> new_node(3,5);

  matrix<size_t> common_face(3,1);
  matrix<size_t> five_points(5,1);
  matrix<double> tet_node0(3,4), tet_node1(3,4), tet_node0_def(3,4), tet_node1_def(3,4);
  matrix<double> grad_op0(4,3), grad_op1(4,3);
  matrix<double> grad0(3,3), grad1(3,3);
  map<pair<size_t,size_t>,double> edge_ratio;

  jtf::tet_mesh tm(pt.get<string>("input/param_tet.value").c_str());

  for(size_t fi = 0; fi < tm_.fa_->face2tet_.size(); ++fi){
      edge_ratio.clear();
      const pair<size_t,size_t> & tet_pair = tm_.fa_->face2tet_[fi];
      if(tm_.fa_->is_outside_face(tet_pair)) continue;

      jtf::mesh::find_common_face(tm_.tetmesh_.mesh_(colon(), tet_pair.first),
                                  tm_.tetmesh_.mesh_(colon(), tet_pair.second),
                                  &common_face[0]);

      const size_t other_point = std::accumulate(tm_.tetmesh_.mesh_(colon(), tet_pair.second).begin(),
                                                 tm_.tetmesh_.mesh_(colon(), tet_pair.second).end(),0)
          - std::accumulate(common_face.begin(), common_face.end(), 0);

      const size_t this_ext_point = std::accumulate(tm_.tetmesh_.mesh_(colon(), tet_pair.first).begin(),
                                                    tm_.tetmesh_.mesh_(colon(), tet_pair.first).end(),0)
          - std::accumulate(common_face.begin(), common_face.end(), 0);

      five_points(colon(0,3),0) = tm_.tetmesh_.mesh_(colon(), tet_pair.first);

      auto it = find(five_points.begin(), five_points.end(), this_ext_point);
      if(it != five_points.begin()){ // reorder
          const size_t it_idx = it-five_points.begin();
          swap(five_points[0], five_points[it_idx]);
          size_t j = 1;
          for(; j < 4; ++j) if(j != it_idx) break;
          swap(five_points[j], five_points[6-j-it_idx]);
        }
      five_points(4,0) = other_point;

      for(size_t i = 0; i < 4; ++i){
          for(size_t j = i+1; j < 5; ++j){
              if(i == 0 && j == 4) continue;
              edge_ratio[make_pair(i,j)] = (point_size[five_points[i]]+point_size[five_points[j]])/2.0;
              //edge_ratio[make_pair(i,j)] = norm(tm.tetmesh_.node_(colon(), five_points[i]) - tm.tetmesh_.node_(colon(), five_points[j]))/norm(tm_.tetmesh_.node_(colon(), five_points[i]) - tm_.tetmesh_.node_(colon(), five_points[j]));
            }
        }

      recover_param_tet_pair_coord(five_points, tm_.tetmesh_.node_, edge_ratio, new_node);

      check_param_tet_pair_coord(five_points, tm_.tetmesh_.node_, edge_ratio, new_node);

      tet_node0 = tm_.tetmesh_.node_(colon(), five_points(colon(0,3),0));
      tet_node1 = tm_.tetmesh_.node_(colon(), five_points(colon(1,4),0));

      tet_node0_def = new_node(colon(), colon(0,3));
      tet_node1_def = new_node(colon(), colon(1,4));

      calc_tet_def_grad_op(&tet_node0[0], &grad_op0[0]);
      calc_tet_def_grad_op(&tet_node1[0], &grad_op1[0]);

      grad0 = tet_node0_def*grad_op0;
      grad1 = tet_node1_def*grad_op1;

      inv(grad0);
      R = grad0*grad1;
      hj::polar3d p;
      p(R,2);

      connection_matrix[tet_pair] = R;
      connection_matrix[make_pair(tet_pair.second, tet_pair.first)] = trans(R);
    }
}

void optimal_sh_generator::compute_inner_connection_jtf2(
    boost::property_tree::ptree &pt,
    map<std::pair<size_t,size_t>,matrix<double> > &connection_matrix)
{
  cerr << "# [info] using deformation gradient connection." << endl;
  jtf::tet_mesh tm_param(pt.get<string>("input/param_tet.value").c_str());
  if(tm_param.tetmesh_.mesh_.size(2) != tm_.tetmesh_.mesh_.size(2)){
      throw std::invalid_argument("wrong param_tet, tet number is not same.");
    }
  vector<zjucad::matrix::matrix<double> > grad(tm_.tetmesh_.mesh_.size(2));
  vector<zjucad::matrix::matrix<double> > inv_grad(tm_.tetmesh_.mesh_.size(2));
  matrix<double> grad_op(4,3);
  matrix<double> tet_node(3,4), param_tet_node(3,4);
  matrix<double> inv_a(3,3);
  matrix<double> R(3,3);
  matrix<double> S = eye<double>(3);
  for(size_t ti = 0; ti < tm_.tetmesh_.mesh_.size(2); ++ti){
      tet_node = tm_.tetmesh_.node_(colon(), tm_.tetmesh_.mesh_(colon(),ti));
      param_tet_node = tm_param.tetmesh_.node_(colon(), tm_.tetmesh_.mesh_(colon(),ti));
      calc_tet_def_grad_op(&tet_node[0], &grad_op[0]);
      grad[ti] = param_tet_node * grad_op;

      inv_a = grad[ti];
      inv(inv_a);
      inv_grad[ti] = inv_a;
    }

  for(size_t fi = 0; fi < tm_.fa_->face2tet_.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = tm_.fa_->face2tet_[fi];
      if(tm_.fa_->is_outside_face(tet_pair)) continue;
      R = inv_grad[tet_pair.first] * grad[tet_pair.second];

      hj::polar3d p;
      p(R,2);

      connection_matrix[tet_pair] = R;
      connection_matrix[make_pair(tet_pair.second, tet_pair.first)] = trans(R);
    }
}
void optimal_sh_generator::compute_inner_connection_jtf(
    zjucad::matrix::matrix<double> &sf,
    std::map<std::pair<size_t, size_t>, double> &inner_face_connection)
{
  matrix<double> sf_on_tet = zeros<double>(tm_.tetmesh_.mesh_.size(2),1);
  for(size_t fi = 0; fi < tm_.fa_->faces_.size(); ++fi){
      const pair<size_t,size_t>  &tet_pair = tm_.fa_->face2tet_[fi];
      if(tet_pair.first != -1)
        sf_on_tet[tet_pair.first] += sf[fi];
      if(tet_pair.second != -1)
        sf_on_tet[tet_pair.second] += sf[fi];
    }
  sf_on_tet /= 4;

  matrix<double> areaed_normal(3,1), sum_normal(3,1);
  matrix<size_t> common_face(3,1);
  double area = 0;
  for(auto & one_edge : tm_.ortae_.e2t_){
      const pair<size_t,size_t> & edge_idx = one_edge.first;
      const vector<size_t> & around_tets = one_edge.second;
      if(!tm_.ortae_.is_inner_edge(around_tets)) continue;
      areaed_normal *= 0;
      sum_normal *= 0;

      double total_vol = 0;
      for(size_t i = 0; i < around_tets.size()-1; ++i){
          jtf::mesh::find_common_face(tm_.tetmesh_.mesh_(colon(), around_tets[i]),
                                      tm_.tetmesh_.mesh_(colon(), around_tets[i+1]),
              &common_face[0]);
          const size_t other_p = std::accumulate(common_face.begin(), common_face.end(),0)
              - edge_idx.first - edge_idx.second;
          common_face[0] = edge_idx.first;
          common_face[1] = edge_idx.second;
          common_face[2] = other_p;
          jtf::mesh::cal_face_normal(common_face, tm_.tetmesh_.node_, areaed_normal, false);
          sum_normal += areaed_normal *(sf_on_tet[around_tets[i+1]] - sf_on_tet[around_tets[i]]);
          total_vol += tm_.vol_[around_tets[i]];
        }

      area = 3*total_vol/norm(tm_.tetmesh_.node_(colon(), edge_idx.second) - tm_.tetmesh_.node_(colon(), edge_idx.first));
      double angle = norm(sum_normal)/(area*2*My_PI());
      angle = float_mod(angle, 2*My_PI());
      inner_face_connection[edge_idx]= angle;
      inner_face_connection[make_pair(edge_idx.second, edge_idx.first)]= -1*angle;
    }
}


void optimal_sh_generator::opt7(zjucad::matrix::matrix<double> &field,
                                boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  init(pt);

  pt.put("weight/surface_smooth.desc", "surface smoothing weight.");
  pt.put("weight/gf.desc", "guiding field align weight.");
  pt.put("weight/normal_align.desc","weight to add surface normal alignment");
  pt.put("weight/feature_line.desc", "weight to add feature line alignment");
  pt.put("input/smooth_scheme.desc", "input smooth scheme strategy [face/gradient]");
  pt.put("input/opt_type.desc", "init or zyz");


  const double w_surface = pt.get<double>("weight/surface_smooth.value");
  const double w_normal_align = pt.get<double>("weight/normal_align.value");
  const double w_gf = pt.get<double>("weight/gf.value",0);
  const double w_feature = pt.get<double>("weight/feature_line.value",0);

  const string smooth_scheme = pt.get<string>("input/smooth_scheme.value");
  const string opt_type = pt.get<string>("input/opt_type.value");

  matrix<double> ks = ones<double>(tm_.outside_face_.size(2),1);

  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);
  {
    func->push_back(math_func_ptr(
                      new hj::math_func::sumsqr<double,int32_t>(
                        build_inner_smooth_func<double>(
                          tm_.tetmesh_.mesh_, tm_.tetmesh_.node_,
                          *tm_.fa_, tm_.ortae_, tm_.vol_, 1, smooth_scheme,
                          opt_type))));
    cerr << "# [info] add inner smooth function. " << smooth_scheme << endl;
  }
  if(w_surface > 1e-8){
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_surface_smooth_func4(
                            tm_.tetmesh_.mesh_, tm_.outside_face_idx_, *tm_.ea_outside_,
                            *tm_.fa_, tm_.outside_face_area_,rij_,
                            local2global_rot_,ks, w_surface, opt_type))));
      cerr << "# [info] add surface smooth function: " << w_surface << endl;
    }
  if(w_normal_align > 1e-8){
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          //                          build_normal_align_func4(
                          //                            tm_.tetmesh_.mesh_, *tm_.fa_, tm_.outside_face_idx_,
                          //                            local2global_rot_, tm_.outside_face_area_,  ks,
                          //                            w_normal_align, false, opt_type))));
                          build_normal_align_func6(
                            tm_.tetmesh_.mesh_, *tm_.fa_, tm_.outside_face_idx_,
                            tm_.outside_face_normal_, tm_.outside_face_area_,  ks,
                            w_normal_align, false, opt_type))));
      cerr << "# [info] add normal align function: " << w_normal_align << endl;
    }
  if(surface_target_dir_angle_.size() > 0 && w_gf > 1e-8){
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_surface_dir_align_func(
                            tm_.tetmesh_.mesh_, *tm_.fa_, tm_.outside_face_area_,
                            local2global_rot_, surface_target_dir_angle_,
                            surface_target_idx_,  surface_target_idx_on_boundary_,
                            w_gf, opt_type))));
      cerr << "# [info] add surface direction align function: " << w_gf << endl;
    }
  if(feature_lines_.size() > 0 && w_feature > 1e-8){
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_feature_line_align_func<double>(
                            tm_.tetmesh_.mesh_, tm_.tetmesh_.node_, feature_lines_, tm_.ortae_,
                            w_feature, opt_type, true))));
      cerr << "# [info] add feature line align function: " << w_feature << endl;
    }

  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  if(opt_type == "init"){
      field = zeros<double>(9,tm_.tetmesh_.mesh_.size(2));
      field(4,colon()) += sqrt(7.0);
      field(8,colon()) += sqrt(5.0);
    }

  jtf::optimize(*obj, field, pt, nullptr, nullptr, nullptr);

  if(opt_type == "init"){
      for(size_t i = 0; i < field.size(2); ++i){
          field(colon(),i) /= norm(field(colon(),i));
        }
      field *= sqrt(12.0);
    }
}

void optimal_sh_generator::opt13(zjucad::matrix::matrix<double> &field,
                                 boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  init(pt);

  pt.put("weight/surface_smooth.desc", "surface smoothing weight.");
  pt.put("weight/gf.desc", "guiding field align weight.");
  pt.put("weight/normal_align.desc","weight to add surface normal alignment");
  pt.put("input/smooth_scheme.desc", "input smooth scheme strategy [face/gradient/laplace/liuyang]");
  pt.put("input/opt_type.desc", "init or zyz");

  const double w_surface = pt.get<double>("weight/surface_smooth.value");
  const double w_normal_align = pt.get<double>("weight/normal_align.value");
  const double w_gf = pt.get<double>("weight/gf.value",0);
  const string smooth_scheme = pt.get<string>("input/smooth_scheme.value");
  const string opt_type = pt.get<string>("input/opt_type.value");

  matrix<double> ks = ones<double>(tm_.outside_face_.size(2),1);

  if(zjucad::has("input/size_field_type.value",pt)){
      matrix<double> sf; // defined on each face
      compute_size_field(sf,pt);

      compute_inner_connection(sf, connection_matrix_);
      {
        ofstream ofs("tet_point_size.vtk");
        tet2vtk(ofs, &tm_.tetmesh_.node_[0], tm_.tetmesh_.node_.size(2), &tm_.tetmesh_.mesh_[0], tm_.tetmesh_.mesh_.size(2));
        point_data(ofs, &sf[0], sf.size(), "point_size");
      }
    }else{
      compute_inner_connection_jtf2(pt,connection_matrix_);
      extract_singularity_from_connection(connection_matrix_);
    }

  {
    jtf::tet_mesh tm(pt.get<string>("input/param_tet.value").c_str());
    matrix<matrix<double> > deformation_gradient;
    cal_deformation_gradient(tm_.tetmesh_.mesh_, tm_.tetmesh_.node_, tm.tetmesh_.node_, deformation_gradient);

    matrix<double> grad;
    matrix<double>  shear;
    vector<double> shear_energy;
    for(size_t fi = 0; fi < tm_.outside_face_idx_.size(); ++fi){
        const pair<size_t,size_t> & tet_pair = tm_.fa_->face2tet_[tm_.outside_face_idx_[fi]];
        assert(tm_.fa_->is_outside_face(tet_pair));
        const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);
        grad = deformation_gradient[tet_idx];
        hj::polar3d p;
        p(grad,2);
        shear = trans(grad) * deformation_gradient[tet_idx];
        for(size_t i = 0; i < shear.size(1); ++i)
          shear(i,i) = 0;
        shear_energy.push_back(norm(shear));
      }
    ofstream ofs("boundary_shear.vtk");
    tri2vtk(ofs, &tm_.tetmesh_.node_[0], tm_.tetmesh_.node_.size(2), &tm_.outside_face_[0], tm_.outside_face_.size(2));
    cell_data(ofs, &shear_energy[0], shear_energy.size(), "shear");
  }
  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);

  {
    func->push_back(math_func_ptr(
                      new hj::math_func::sumsqr<double,int32_t>(
                        build_inner_smooth_func_with_connection3(
                          tm_.tetmesh_.mesh_, tm_.tetmesh_.node_, *tm_.fa_,
                          tm_.ortae_, tm_.vol_, 1, connection_matrix_,
                          smooth_scheme, opt_type))));
    cerr << "# [info] add inner smooth function. " << smooth_scheme << endl;
  }
  if(w_surface > 1e-8){
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_surface_smooth_func4(
                            tm_.tetmesh_.mesh_, tm_.outside_face_idx_, *tm_.ea_outside_,
                            *tm_.fa_, tm_.outside_face_area_,rij_,
                            local2global_rot_,ks, w_surface, opt_type))));
      cerr << "# [info] add surface smooth function: " << w_surface << endl;
    }
  if(w_normal_align > 1e-8){
      jtf::tet_mesh tm(pt.get<string>("input/param_tet.value").c_str());
      matrix<matrix<double> > deformation_gradient;
      cal_deformation_gradient(tm_.tetmesh_.mesh_, tm_.tetmesh_.node_, tm.tetmesh_.node_, deformation_gradient);

      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_normal_align_func4(
                            tm_.tetmesh_.mesh_, *tm_.fa_, tm_.outside_face_idx_,
                            local2global_rot_, tm_.outside_face_area_,  ks,
                            w_normal_align, false, opt_type))
                        //                          build_normal_align_func4_with_deformation_gradient(
                        //                            tm_.tetmesh_.mesh_, *tm_.fa_, tm_.outside_face_idx_,
                        //                            tm_.outside_face_normal_, deformation_gradient,
                        //                            tm_.outside_face_area_, w_normal_align, false, opt_type)
                        ));
      cerr << "# [info] add normal align function: " << w_normal_align << endl;
    }
  if(surface_target_dir_angle_.size() > 0 && w_gf){
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_surface_dir_align_func(
                            tm_.tetmesh_.mesh_, *tm_.fa_, tm_.outside_face_area_,
                            local2global_rot_, surface_target_dir_angle_,
                            surface_target_idx_,  surface_target_idx_on_boundary_,
                            w_gf, opt_type))));
      cerr << "# [info] add surface direction align function: " << w_gf << endl;
    }

  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  if(opt_type == "init"){
      field = zeros<double>(9,tm_.tetmesh_.mesh_.size(2));
      field(4,colon()) += sqrt(7.0);
      field(8,colon()) += sqrt(5.0);
    }

  jtf::optimize(*obj, field, pt, nullptr, nullptr, nullptr);

  if(opt_type == "init"){
      for(size_t i = 0; i < field.size(2); ++i){
          field(colon(),i) /= norm(field(colon(),i));
        }
      field *= sqrt(12.0);
    }

  //  {
  //    jtf::tet_mesh tm(pt.get<string>("input/param_tet.value").c_str());
  //    matrix<matrix<double> > deformation_gradient;
  //    cal_deformation_gradient(tm_.tetmesh_.mesh_, tm_.tetmesh_.node_, tm.tetmesh_.node_, deformation_gradient);

  //    zjucad::matrix::matrix<double> zyz;
  //    project2zyz(field, zyz);
  //    matrix<matrix<double> > rot(zyz.size(2),1);
  //    matrix<double> rot_m(3,3);
  //    for(size_t ti = 0; ti < zyz.size(2); ++ti){
  //        rot[ti].resize(3,3);
  //        zyz_angle_2_rotation_matrix1(&zyz(0,ti), &rot[ti][0]);
  //        rot_m = deformation_gradient[ti];
  //        inv(rot_m);
  //        hj::polar3d p;
  //        p(rot_m,2);
  //        rot[ti] = temp(rot_m * rot[ti]);
  //        rotation_matrix_2_zyz_angle(&rot[ti][0], &zyz(0,ti),0);
  //        calc_rot_cubic_f_sh_(&field(0,ti), &zyz(0,ti));
  //      }
  //  }
}

void optimal_sh_generator::extract_singularity_from_connection(
    const std::map<std::pair<size_t,size_t>, matrix<double> > &connection_matrix)
{
  matrix<double> rot = eye<double>(3),quat(4,1), sum_rot = eye<double>(3);
  matrix<double> dir(3,1);
  vector<size_t> lines;
  vector<double> line_type,line_type2,line3;
  vector<matrix<double> > Frame,Frame2;
  for(const auto & one_edge : tm_.ortae_.e2t_){
      const pair<size_t,size_t> & edge_points = one_edge.first;
      const vector<size_t> & around_tets = one_edge.second;
      if(!tm_.ortae_.is_inner_edge(around_tets)) continue;
      Frame.clear();
      Frame2.clear();
      sum_rot = eye<double>(3);
      for(size_t i = 0; i < around_tets.size()-1; ++i){
          auto it = connection_matrix.find(make_pair(around_tets[i+1], around_tets[i]));
          sum_rot = temp(it->second*sum_rot);
        }

      {
        double theta;
        convert_rotation_matrix_to_quat(sum_rot,&quat[0]);
        convert_quat_to_axis_angle(&quat[0], &dir[0], theta);

        matrix<double> dir_edge = tm_.tetmesh_.node_(colon(), edge_points.second)
            - tm_.tetmesh_.node_(colon(), edge_points.first);
        dir_edge /= norm(dir_edge);
        const double angle = float_mod(jtf::math::safe_acos(dot(dir_edge, dir)),2*My_PI());
        {
          lines.push_back(edge_points.first);
          lines.push_back(edge_points.second);
          line_type.push_back(dot(dir_edge,dir)*theta);
          //          line_type2.push_back(angle);
          //          line3.push_back(theta);
        }
      }
      //            {
      //              if(norm(sum_rot-eye<double>(3))>1e-6){
      //                  lines.push_back(edge_points.first);
      //                  lines.push_back(edge_points.second);
      //                  line_type.push_back(norm(sum_rot-eye<double>(3)));
      //                }
      //            }
    }

  std::ofstream ofs("test_singularity.vtk");
  line2vtk(ofs, &tm_.tetmesh_.node_[0], tm_.tetmesh_.node_.size(2), &lines[0], lines.size()/2);
  cell_data(ofs, &line_type[0], line_type.size(), "dot");
  //  vtk_data(ofs, &line_type2[0], line_type2.size(), "dot_angle");
  //  vtk_data(ofs, &line3[0], line3.size(), "rot_angle");
}

void optimal_sh_generator::opt8(zjucad::matrix::matrix<double> &field,
                                boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  init(pt);

  pt.put("weight/surface_smooth.desc", "surface smoothing weight.");
  pt.put("weight/normal_align.desc","weight to add surface normal alignment");
  pt.put("input/smooth_scheme.desc", "input smooth scheme strategy [face/gradient]");
  pt.put("input/opt_type.desc", "init or zyz");

  const double w_surface = pt.get<double>("weight/surface_smooth.value");
  const double w_normal_align = pt.get<double>("weight/normal_align.value");
  const string smooth_scheme = pt.get<string>("input/smooth_scheme.value");
  const string opt_type = pt.get<string>("input/opt_type.value");

  matrix<double> ks = ones<double>(tm_.outside_face_.size(2),1);

  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);
  {
    func->push_back(math_func_ptr(
                      new hj::math_func::sumsqr<double,int32_t>(
                        build_inner_smooth_func<double>(tm_.tetmesh_.mesh_, tm_.tetmesh_.node_,
                                                *tm_.fa_, tm_.ortae_, tm_.vol_, 1, smooth_scheme, opt_type))));
    cerr << "# [info] add inner smooth function. " << smooth_scheme << endl;
  }
  if(w_surface > 1e-8){
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_surface_smooth_func4(tm_.tetmesh_.mesh_, tm_.outside_face_idx_, *tm_.ea_outside_,
                                                     *tm_.fa_, tm_.outside_face_area_,rij_,
                                                     local2global_rot_,ks, w_surface, opt_type))));
      cerr << "# [info] add surface smooth function: " << w_surface << endl;
    }
  if(w_normal_align > 1e-8){
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_normal_align_func5(tm_.tetmesh_.mesh_, *tm_.fa_, tm_.outside_face_idx_, local2global_rot_,
                                                   tm_.outside_face_area_,  ks, w_normal_align, false, opt_type)))
                      );
      cerr << "# [info] add normal align function: " << w_normal_align << endl;
    }

  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  if(opt_type == "init" && !zjucad::has("input/init_sh.value",pt)){
      field = zeros<double>(9,tm_.tetmesh_.mesh_.size(2));
      field(4,colon()) += sqrt(7.0);
      field(8,colon()) += sqrt(5.0);
    }else if(opt_type == "init" && zjucad::has("input/init_sh.value", pt)){
      matrix<double> init_zyz_sh = zeros<double>(3, tm_.tetmesh_.mesh_.size(2));
      if(jtf::mesh::read_matrix(pt.get<string>("input/init_sh.value").c_str(), init_zyz_sh)){
          cerr << "# [error] can not load init_zyz file." << endl;
          return ;
        }
      if(init_zyz_sh.size(1) != 3 || init_zyz_sh.size(2) != tm_.tetmesh_.mesh_.size(2)){
          cerr << "# [error] init zyz is not compatible with mesh." << endl;
          return ;
        }
      field = zeros<double>(9,tm_.tetmesh_.mesh_.size(2));
      for(size_t ti = 0; ti < init_zyz_sh.size(2) ; ++ti){
          calc_rot_cubic_f_sh_(&field(0,ti), &init_zyz_sh(0,ti));
        }
    }

  jtf::optimize(*obj, field, pt, nullptr, nullptr, nullptr);

  if(opt_type == "init"){
      for(size_t i = 0; i < field.size(2); ++i){
          field(colon(),i) /= norm(field(colon(),i));
        }
      field *= sqrt(12.0);
    }

  {
    matrix<double> rot(9,9);
    for(size_t fi = 0; fi < tm_.outside_face_idx_.size(); ++fi){
        const pair<size_t,size_t>  & tet_pair = tm_.fa_->face2tet_[tm_.outside_face_idx_[fi]];
        assert(tm_.fa_->is_outside_face(tet_pair));
        const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);
        calc_rot_cubic_f_sh_mat_(&rot[0], &local2global_rot_(0,fi));
        matrix<double> d = (trans(rot)*field(colon(), tet_idx));
        if(fabs(d[4]-sqrt(7.0)) > 1e-8)
          cerr << "# [error] wrong sh on face " << fi << d << endl;
      }
  }
}

template <typename val_type, typename int_type>
class field_fix_func_sh : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  field_fix_func_sh(const size_t tet_num, const size_t idx,
                    const matrix<val_type> &target, const double w)
    :tet_num_(tet_num), idx_(idx), target_(target), w_(w){}
  virtual ~field_fix_func_sh(){}
  virtual size_t nx() const{
    return 9 * tet_num_;
  }
  virtual size_t nf() const{
    return 9;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type *> x0(9,tet_num_,x);
        for(int_type i = 0; i < 9; ++i){
            int_type c[1] = {i};
            cv[c] += w_*(x0(i,idx_)-target_[i]);
          }
      }
    if(k == 1){
        for(int_type i = 0; i < 9; ++i){
            int_type c1[2] = {i,9*idx_+i};
            cv[c1] += w_;
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
            int_type c1[2] = {i,9*idx_+i};
            l2g.add(cs,c1);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 9;
  }
private:
  const size_t tet_num_;
  const int_type idx_;
  const double w_;
  const matrix<val_type> target_;
};

double weight_func( double x)
{
  const double pow_2 = x*x;
  const double pow_4 = pow_2*pow_2;
  return -1e5*(pow_2*pow_4/6.0 - pow_4* x/5.0 +0.0875*pow_4-125.0/7500.0*pow_2*x+3/2500*pow_2)+10;
}

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_fix_inner_field_func(const vector<matrix<double> > &sh,
                           const vector<double> & radius,
                           const matrix<size_t> & mesh,
                           const matrix<double> & node,
                           const matrix<double> & vol,
                           const matrix<double> & size_field,
                           const jtf::mesh::face2tet_adjacent &fa)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  assert(radius.back() == *max_element(radius.begin(), radius.end()));
  const double max_size = *max_element(size_field.begin(), size_field.end());
  const double min_size = *min_element(size_field.begin(), size_field.end());

  matrix<double> size_field_resize = size_field;
  const double max_r = radius.back();
  const double min_r = radius.front();
  const double per = (max_r-min_r)/(max_size-min_size);
  for(size_t i = 0; i < size_field_resize.size(); ++i){
      size_field_resize[i] = per * (size_field_resize[i]-min_size)+min_r;
    }

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  const double total_vol = std::accumulate(vol.begin(), vol.end(), 0.0);

  double delta_step = 1;
  if(radius.size() > 2) delta_step = fabs(radius[1] - radius[0]);
  for(size_t ti = 0; ti < mesh.size(2); ++ti){
      const double w = vol[ti]/total_vol;
      for(size_t i = 0; i < sh.size(); ++i){

          all_func->push_back(math_func_ptr(
                                new field_fix_func_sh<double,int32_t>(mesh.size(2), ti, sh[i](colon(),ti),
                                sqrt(w*exp(-1*pow(2./delta_step*(size_field_resize[ti]-radius[i]),2.0))))
                              //sqrt(w/(pow(size_field_resize[ti]-radius[i],2.0)+1e-5)))
                              ));
        }
    }

  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}


void optimal_sh_generator::opt11(zjucad::matrix::matrix<double> &field,
                                 boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  pt.put("input/opt_type.desc", "init or zyz");
  const string opt_type = pt.get<string>("input/opt_type.value");
  if(opt_type != "zyz")
    throw std::invalid_argument("opt11 only support zyz opt.");

  pt.put("input/smooth_scheme.desc", "input smooth scheme strategy [face/gradient]");
  pt.put("input/opt_type.desc", "init or zyz");

  const double w_surface = pt.get<double>("weight/surface_fix.value",10);
  const string smooth_scheme = pt.get<string>("input/smooth_scheme.value");

  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);
  {
    func->push_back(math_func_ptr(
                      new hj::math_func::sumsqr<double,int32_t>(
                        build_inner_smooth_func<double>(tm_.tetmesh_.mesh_, tm_.tetmesh_.node_,
                                                *tm_.fa_, tm_.ortae_, tm_.vol_, 1, smooth_scheme, opt_type))));
    cerr << "# [info] add inner smooth function. " << smooth_scheme << endl;
  }
  if(w_surface > 1e-8){
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_surface_fix_func(tm_.tetmesh_.mesh_, tm_.outside_face_idx_,
                                                 *tm_.fa_, tm_.outside_face_area_,field, w_surface, opt_type))));
      cerr << "# [info] add surface fix function: " << w_surface << endl;
    }

  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  jtf::optimize(*obj, field, pt, nullptr, nullptr, nullptr);

}

void optimal_sh_generator::opt10(zjucad::matrix::matrix<double> &field,
                                 boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  double r = 0;

  const size_t seq_num = 3;
  vector<matrix<double> > sh(seq_num);
  vector<double> radius(seq_num);
  pt.put("input/opt_type.desc", "init or zyz");
  const string opt_type = pt.get<string>("input/opt_type.value");

  for(size_t i = 0; i < seq_num; ++i){
      pt.put("weight/gaussian_radius.value", r);
      opt7(sh[i], pt);
      radius[i] = r;
      r += 0.05;
    }

  matrix<double> face_size_field;
  compute_size_field(face_size_field, pt);

  matrix<double> tet_size = zeros<double>(tm_.tetmesh_.mesh_.size(2),1);
  for(size_t fi = 0; fi < tm_.fa_->faces_.size(); ++fi){
      const pair<size_t,size_t>  &tet_pair = tm_.fa_->face2tet_[fi];
      if(tet_pair.first != -1)
        tet_size[tet_pair.first] += face_size_field[fi];
      if(tet_pair.second != -1)
        tet_size[tet_pair.second] += face_size_field[fi];
    }

  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);
  func->push_back(math_func_ptr(
                    new hj::math_func::sumsqr<double,int32_t>(
                      build_fix_inner_field_func(sh, radius, tm_.tetmesh_.mesh_,
                                                 tm_.tetmesh_.node_,
                                                 tm_.vol_,tet_size, *(tm_.fa_)))));

  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  field = sh[0];
  jtf::optimize(*obj, field, pt, nullptr, nullptr, nullptr);

  if(opt_type == "init"){
      for(size_t i = 0; i < field.size(2); ++i){
          field(colon(),i) /= norm(field(colon(),i));
        }
      field *= sqrt(12.0);
    }
}


void optimal_sh_generator::opt12(zjucad::matrix::matrix<double> &field,
                                 boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  init(pt);

  pt.put("weight/surface_smooth.desc", "surface smoothing weight.");
  pt.put("weight/normal_align.desc","weight to add surface normal alignment");
  pt.put("input/smooth_scheme.desc", "input smooth scheme strategy [face/gradient]");
  pt.put("input/opt_type.desc", "init or zyz");

  const double w_surface = pt.get<double>("weight/surface_smooth.value");
  const double w_normal_align = pt.get<double>("weight/normal_align.value");
  const string smooth_scheme = pt.get<string>("input/smooth_scheme.value");
  const string opt_type = pt.get<string>("input/opt_type.value");

  matrix<double> ks = ones<double>(tm_.outside_face_.size(2),1);

  //matrix<double> size_field = ones<double>(tm_.fa_->faces_.size(),1);
  // compute_size_field(size_field, pt);

  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);
  {
    func->push_back(math_func_ptr(
                      new hj::math_func::sumsqr<double,int32_t>(
                        build_inner_smooth_func<double>(tm_.tetmesh_.mesh_, tm_.tetmesh_.node_,
                                                *tm_.fa_, tm_.ortae_, tm_.vol_, 1, smooth_scheme, opt_type))));
    cerr << "# [info] add inner smooth function. " << smooth_scheme << endl;
  }
  if(w_surface > 1e-8){
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_surface_smooth_func4(tm_.tetmesh_.mesh_, tm_.outside_face_idx_, *tm_.ea_outside_,
                                                     *tm_.fa_, tm_.outside_face_area_,rij_,
                                                     local2global_rot_,ks, w_surface, opt_type))));
      cerr << "# [info] add surface smooth function: " << w_surface << endl;
    }
  if(w_normal_align > 1e-8){
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_normal_align_func4(tm_.tetmesh_.mesh_, *tm_.fa_, tm_.outside_face_idx_, local2global_rot_,
                                                   tm_.outside_face_area_,  ks, w_normal_align, false, opt_type)))
                      );
      cerr << "# [info] add normal align function: " << w_normal_align << endl;
    }

  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  if(opt_type == "init"){
      field = zeros<double>(9,tm_.tetmesh_.mesh_.size(2));
      field(4,colon()) += sqrt(7.0);
      field(8,colon()) += sqrt(5.0);
    }

  jtf::optimize(*obj, field, pt, nullptr, nullptr, nullptr);

  if(opt_type == "init"){
      const double target = 12.;
      double w = 10;

      for(size_t i = 0; i < 2; ++i){
          cerr << "# [info] fix |sh|^2 = 12, weight " << w << endl;
          math_func_ptr func_length_func(
                new hj::math_func::sumsqr<double,int32_t>(
                  build_length_fix_func(tm_.tetmesh_.mesh_, tm_.tetmesh_.node_, tm_.vol_, target, w)));

          func->push_back(func_length_func);
          math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
          math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

          jtf::optimize(*obj, field, pt, nullptr, nullptr, nullptr);
          func->pop_back();
          w *= 2;
        }
    }

  if(opt_type == "init"){
      for(size_t i = 0; i < field.size(2); ++i){
          field(colon(),i) /= norm(field(colon(),i));
        }
      field *= sqrt(12.0);
    }
}

void optimal_sh_generator::opt9(zjucad::matrix::matrix<double> &sh,
                                boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  init(pt);

  hj::sparse::csc<double, int32_t> A,AT;
  get_jacobian_matrix2(A,pt);
  hj::sparse::trans(A,AT);

  hj::sparse::csc<double,int32_t> ATA;
  hj::sparse::AAT<hj::sparse::map_by_sorted_vector >(AT,ATA);

  matrix<double> eigen_vector;
  solve_smallest_eigenvalue_problem(ATA, eigen_vector);

  itr_matrix<const double*> x0(9,tm_.tetmesh_.mesh_.size(2), &eigen_vector[0]);

  sh = x0;

  for(size_t i = 0; i < sh.size(2); ++i){
      sh(colon(),i) /= norm(sh(colon(),i));
    }
  sh *= sqrt(12.0);

  {
    matrix<double> rot(9,9);
    for(size_t fi = 0; fi < tm_.outside_face_idx_.size(); ++fi){
        const pair<size_t,size_t>  & tet_pair = tm_.fa_->face2tet_[tm_.outside_face_idx_[fi]];
        assert(tm_.fa_->is_outside_face(tet_pair));
        const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);
        calc_rot_cubic_f_sh_mat_(&rot[0], &local2global_rot_(0,fi));
        matrix<double> d = (trans(rot)*sh(colon(), tet_idx));
        if(fabs(d[4]-sqrt(7.0)) > 1e-8)
          cerr << "# [error] wrong sh on face " << fi << d << endl;
      }
  }
}

void optimal_sh_generator::opt3(zjucad::matrix::matrix<double> &sh,
                                boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  init(pt);

  pt.put("weight/surface_smooth.desc", "surface smoothing weight.");
  pt.put("weight/normal_align.desc","weight to add surface normal alignment");
  const double w_surface = pt.get<double>("weight/surface_smooth.value");
  const double w_normal_align = pt.get<double>("weight/normal_align.value");
  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);
  {
    func->push_back(math_func_ptr(
                      new hj::math_func::sumsqr<double,int32_t>(
                        build_inner_smooth_func<double>(tm_.tetmesh_.mesh_, tm_.tetmesh_.node_,
                                                *tm_.fa_, tm_.ortae_, tm_.vol_, 1))));
    cerr << "# [info] add inner smooth function." << endl;
  }
  if(w_surface > 1e-8){
      if(zjucad::has("input/cos_sin.value",pt)){
          matrix<double> tradition;
          jtf::mesh::read_matrix(pt.get<string>("input/cos_sin.value").c_str(),tradition);
          func->push_back(math_func_ptr(
                            new hj::math_func::sumsqr<double,int32_t>(
                              build_surface_fix_func2(tm_.tetmesh_.mesh_, tm_.outside_face_idx_,tm_.outside_face_area_,
                                                      *tm_.fa_, local2global_rot_, tradition, w_surface))));
          cerr << "# [info] add surface fix function: " << w_surface << endl;
        }else{
          func->push_back(math_func_ptr(
                            new hj::math_func::sumsqr<double,int32_t>(
                              build_surface_smooth_func2(tm_.tetmesh_.mesh_, tm_.outside_face_idx_, *tm_.ea_outside_,
                                                         *tm_.fa_, tm_.outside_face_area_,rij_,
                                                         local2global_rot_, w_surface))));
          cerr << "# [info] add surface smooth function: " << w_surface << endl;
        }
    }
  if(w_normal_align > 1e-8){
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_normal_align_func2(tm_.tetmesh_.mesh_, *tm_.fa_, tm_.outside_face_idx_, local2global_rot_,
                                                   tm_.outside_face_area_, w_normal_align))));
      cerr << "# [info] add normal align function: " << w_normal_align << endl;
    }

  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  sh = zeros<double>(9,tm_.tetmesh_.mesh_.size(2));
  sh(4,colon()) += sqrt(7.0);
  sh(8,colon()) += sqrt(5.0);
  if(zjucad::has("input/init.value" ,pt)){
      matrix<double> zyz;
      jtf::mesh::read_matrix(pt.get<string>("input/init.value").c_str(), zyz);
      if(zyz.size(2) != tm_.tetmesh_.mesh_.size(2)){
          cerr << "# [error] wrong zyz file." << endl;
          return;
        }
      size_t i = 0;
#pragma omp parallel for private(i)
      for(i = 0; i < tm_.tetmesh_.mesh_.size(2); ++i){
          calc_rot_cubic_f_sh_(&sh(0,i),&zyz(0,i));
        }
      cerr << "# [info] load init zyz." << endl;
    }
  jtf::optimize(*obj, sh, pt, nullptr, nullptr, nullptr);


  //ipopt_solve(sh, *jtf_func, nullptr, pt);
  cerr << sh(colon(),colon(0,5)) << endl;
  for(size_t i = 0; i < sh.size(2); ++i){
      sh(colon(),i) /= norm(sh(colon(),i));
    }
  sh *= sqrt(12.0);
}

void optimal_sh_generator::get_jacobian_matrix(hj::sparse::csc<double,int32_t> &A)
{
  const double total_v = std::accumulate(tm_.vol_.begin(), tm_.vol_.end(),0.0);
  const double total_area = std::accumulate(edge_area_weight_.begin(), edge_area_weight_.end(), 0.0);

  const size_t inner_face_num = tm_.fa_->faces_.size() - tm_.outside_face_.size(2);
  const size_t edge_num = tm_.ea_outside_->edges_.size();

  csc_filler<double> cf(9*inner_face_num+2*edge_num+6*tm_.outside_face_.size(2), tm_.tetmesh_.mesh_.size(2)*9);

  { // for inner smoothness
    size_t face_i = 0;
    for(size_t fi = 0; fi < tm_.fa_->face2tet_.size(); ++fi){
        const pair<size_t,size_t> & tet_pair = tm_.fa_->face2tet_[fi];
        if(tm_.fa_->is_outside_face(tet_pair)) continue;
        const double w = sqrt((tm_.vol_[tet_pair.first]+tm_.vol_[tet_pair.second])/(4*total_v));
        for(size_t j = 0; j < 9; ++j){
            cf.add(9*face_i+j,9*tet_pair.first+j, w);
            cf.add(9*face_i+j,9*tet_pair.second+j,-w);
          }
        ++face_i;
      }
    assert(face_i == inner_face_num);
  }

  {// for surface smoothness
    matrix<double> RL2G0(9,9), RL2G1(9,9);
    matrix<double> R0(2,9), R1(2,9);
    matrix<double> Rc = zeros<double>(2,9);
    Rc(0,8) = 1.0; Rc(1,0) = 1.0;
    matrix<double> Rij(2,2);
    const double Nsys = 4;
    for(size_t ei = 0; ei < tm_.ea_outside_->edge2cell_.size(); ++ei){
        const pair<size_t,size_t> & face_pair = tm_.ea_outside_->edge2cell_[ei];
        assert(tm_.ea_outside_->is_boundary_edge(face_pair) != true);
        const pair<size_t,size_t> & tet_pair0 = tm_.fa_->face2tet_.at(tm_.outside_face_idx_[face_pair.first]);
        const pair<size_t,size_t> & tet_pair1 = tm_.fa_->face2tet_.at(tm_.outside_face_idx_[face_pair.second]);
        const size_t t0 = (tet_pair0.first==-1?tet_pair0.second:tet_pair0.first);
        const size_t t1 = (tet_pair1.first==-1?tet_pair1.second:tet_pair1.first);

        calc_rot_cubic_f_sh_mat_(&RL2G0[0], &local2global_rot_(0,face_pair.first));
        calc_rot_cubic_f_sh_mat_(&RL2G1[0], &local2global_rot_(0,face_pair.second));
        get_2d_rotation(Nsys*rij_(0,ei), Rij);
        R0 = Rij*Rc*trans(RL2G0);
        R1 = Rc*trans(RL2G1);

        const double w = sqrt(edge_area_weight_[ei]/(total_area));
        for(size_t i = 0; i < 2; ++i){
            for(size_t j = 0;  j < 9; ++j){
                cf.add(9*inner_face_num+2*ei+i, 9*t0+j, w*R0(i,j));
                cf.add(9*inner_face_num+2*ei+i, 9*t1+j, -1.0*w*R1(i,j));
              }
          }
      }
  }

  {
    // for normal alignment
    matrix<double> Rpick = zeros<double>(6,9);
    Rpick(0,1) = 1.0;
    Rpick(1,2) = 1.0;
    Rpick(2,3) = 1.0;
    Rpick(3,5) = 1.0;
    Rpick(4,6) = 1.0;
    Rpick(5,7) = 1.0;

    matrix<double> Rl2g(9,9);
    matrix<double> R(6,9);
    for(size_t fi = 0; fi < tm_.outside_face_idx_.size(); ++fi){
        const pair<size_t,size_t>  & tet_pair = tm_.fa_->face2tet_[tm_.outside_face_idx_[fi]];
        assert(tm_.fa_->is_outside_face(tet_pair));
        const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);
        calc_rot_cubic_f_sh_mat_(&Rl2g[0], &local2global_rot_(0,fi));
        R = Rpick * trans(Rl2g);
        const double w = sqrt(1*tm_.outside_face_area_[fi]/total_area);
        for(size_t j = 0; j < R.size(1); ++j){
            for(size_t i = 0; i < R.size(2); ++i){
                cf.add(9*inner_face_num+2*edge_num+6*fi+j, 9*tet_idx+i, w*R(j,i));
              }
          }
      }
  }
  cf.out(A);
}

void optimal_sh_generator::get_jacobian_matrix2(hj::sparse::csc<double,int32_t> &A,
                                                boost::property_tree::ptree &pt)
{
  const double total_v = std::accumulate(tm_.vol_.begin(), tm_.vol_.end(),0.0);
  const double total_area = std::accumulate(edge_area_weight_.begin(), edge_area_weight_.end(), 0.0);

  const size_t inner_face_num = tm_.fa_->faces_.size() - tm_.outside_face_.size(2);

  const size_t edge_num = tm_.ea_outside_->edges_.size();
  csc_filler<double> cf(9*inner_face_num+2*edge_num+6*tm_.outside_face_.size(2), tm_.tetmesh_.mesh_.size(2)*9);

  { // for inner smoothness
    cerr << "# [info] add inner smoothness. " << endl;
    size_t face_i = 0;
    for(size_t fi = 0; fi < tm_.fa_->face2tet_.size(); ++fi){
        const pair<size_t,size_t> & tet_pair = tm_.fa_->face2tet_[fi];
        if(tm_.fa_->is_outside_face(tet_pair)) continue;
        const double w = sqrt((tm_.vol_[tet_pair.first]+tm_.vol_[tet_pair.second])/(4*total_v));
        for(size_t j = 0; j < 9; ++j){
            cf.add(9*face_i+j,9*tet_pair.first+j, w);
            cf.add(9*face_i+j,9*tet_pair.second+j,-w);
          }
        ++face_i;
      }
    assert(face_i == inner_face_num);
  }

  pt.put("weight/surface_smooth.desc", "surface smooth weighting");
  const double w_s = pt.get<double>("weight/surface_smooth.value");
  if(w_s > 1e-8) {// for surface smoothness
      cerr << "# [info] add surface smoothness constraint " << w_s << endl;
      matrix<double> RL2G0(9,9), RL2G1(9,9);
      matrix<double> Rij = zeros<double>(2,2);
      matrix<double> Rc = zeros<double>(2,9);
      Rc(0,8) = 1.0; Rc(1,0) = 1.0;
      matrix<double> Ri(2,9) , Rj(2,9);
      const double Nsys = 4;
      for(size_t ei = 0; ei < tm_.ea_outside_->edge2cell_.size(); ++ei){
          const pair<size_t,size_t> & face_pair = tm_.ea_outside_->edge2cell_[ei];
          assert(tm_.ea_outside_->is_boundary_edge(face_pair) != true);
          const pair<size_t,size_t> & tet_pair0 = tm_.fa_->face2tet_.at(tm_.outside_face_idx_[face_pair.first]);
          const pair<size_t,size_t> & tet_pair1 = tm_.fa_->face2tet_.at(tm_.outside_face_idx_[face_pair.second]);
          const size_t t0 = (tet_pair0.first==-1?tet_pair0.second:tet_pair0.first);
          const size_t t1 = (tet_pair1.first==-1?tet_pair1.second:tet_pair1.first);

          calc_rot_cubic_f_sh_mat_(&RL2G0[0], &local2global_rot_(0,face_pair.first));
          calc_rot_cubic_f_sh_mat_(&RL2G1[0], &local2global_rot_(0,face_pair.second));
          get_2d_rotation(Nsys*rij_(0,ei), Rij);

          Ri = Rij * (Rc*trans(RL2G0));
          Rj = Rc*trans(RL2G1);

          const double w = sqrt(w_s*edge_area_weight_[ei]/total_area);
          for(size_t i = 0; i < 2; ++i){
              for(size_t j = 0;  j < 9; ++j){
                  cf.add(9*inner_face_num+2*ei+i, 9*t0+j, w*Ri(i,j));
                  cf.add(9*inner_face_num+2*ei+i, 9*t1+j, -1.0*w*Rj(i,j));
                }
            }
        }
    }

  pt.put("weight/normal_align.desc", "normal align weighting");
  const double w_n = pt.get<double>("weight/normal_align.value");
  vector<size_t> zeros_idx = {1,2,3,5,6,7};
  if(w_n > 1e-8) { // normal alignment
      cerr << "# [info] add normal alignment constraint " << w_n << endl;
      matrix<double> Rl2g(9,9);
      matrix<double> R(9,1);
      for(size_t fi = 0; fi < tm_.outside_face_idx_.size(); ++fi){
          const pair<size_t,size_t>  & tet_pair = tm_.fa_->face2tet_[tm_.outside_face_idx_[fi]];
          assert(tm_.fa_->is_outside_face(tet_pair));
          const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);
          calc_rot_cubic_f_sh_mat_(&Rl2g[0], &local2global_rot_(0,fi));

          const double w = sqrt(w_n*tm_.outside_face_area_[fi]/total_area);
          for(size_t idx = 0; idx < zeros_idx.size(); ++idx){
              R = Rl2g(colon(),zeros_idx[idx]);
              for(size_t i = 0; i < 9; ++i){
                  cf.add(9*inner_face_num+2*edge_num+6*fi+idx, 9*tet_idx+i, w*R[i]);
                }
            }
        }
    }
  cf.out(A);
}

double optimal_surface_sh_generator::get_edge_adj_tri_area(
    const pair<size_t,size_t> & one_edge,
    const zjucad::matrix::matrix<double> & node,
    const zjucad::matrix::matrix<size_t> & face)
{
  using namespace zjucad::matrix;
  /*    if(jtf::mesh::is_acute_triange(face,node)){
        matrix<double> center;
        const double r = jtf::mesh::get_triangle_circumcenter(node(colon(), face),center);
        double len = norm(node(colon(),one_edge.first)-node(colon(), one_edge.second));
        return jtf::mesh::cal_face_area(r,r,len);
    }else*/{
    return jtf::mesh::cal_face_area(face, node)/3.0;
  }
}

void optimal_surface_sh_generator::init(boost::property_tree::ptree &pt)
{
  edge_area_weight_ = zjucad::matrix::zeros<double>(tm_.ea_->edges_.size(),1);
  for(size_t ei = 0; ei < tm_.ea_->edge2cell_.size(); ++ei){
      const pair<size_t,size_t> & face_pair = tm_.ea_->edge2cell_[ei];
      if(face_pair.first != -1){
          edge_area_weight_[ei] += fabs(get_edge_adj_tri_area(tm_.ea_->edges_[ei], tm_.trimesh_.node_, tm_.trimesh_.mesh_(colon(), face_pair.first)));
        }
      if(face_pair.second != -1)
        edge_area_weight_[ei] += fabs(get_edge_adj_tri_area(tm_.ea_->edges_[ei], tm_.trimesh_.node_, tm_.trimesh_.mesh_(colon(), face_pair.second)));
    }
  init_connection(pt);
}

template <typename val_type, typename int_type>
class smooth_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  smooth_func(const hj::sparse::csc<val_type, int_type> &AT):AT_(AT){}
  virtual ~smooth_func(){}
  virtual size_t nx() const
  {
    return AT_.size(1);
  }
  virtual size_t nf() const
  {
    return AT_.size(2);
  }
  virtual int eval(size_t k , const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    itr_matrix<const val_type*> x0(AT_.size(2),1,x);
    if(k == 0){
        int_type c[1];
        matrix<val_type> Ax = zeros<val_type>(AT_.size(2),1);
        hj::sparse::mv(true, AT_, x0, Ax);
        for(size_t i = 0; i < AT_.size(1); ++i){
            c[0] = i;
            cv[c] += Ax[i];
          }
      }
    if(k == 1){
        int_type c[2];
        for(size_t i = 0; i < AT_.size(2); ++i){
            c[0] = i;
            for(size_t off = AT_.ptr()[i]; off != AT_.ptr()[i+1]; ++off){
                c[1] = AT_.idx()[off];
                cv[c] += AT_.val()[off];
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
        int_type c[2];
        for(size_t i = 0; i < AT_.size(2); ++i){
            c[0] = i;
            for(size_t off = AT_.ptr()[i]; off != AT_.ptr()[i+1]; ++off){
                c[1] = AT_.idx()[off];
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
        return AT_.nnz();
      }
  }
private:
  const hj::sparse::csc<val_type, int_type> &AT_;
};

template <typename val_type, typename int_type>
class fit_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  fit_func(const zjucad::matrix::matrix<val_type> &target,
           const val_type w)
    :target_(target), w_(w){}
  virtual ~fit_func(){}
  virtual size_t nx() const
  {
    return target_.size();
  }
  virtual size_t nf() const
  {
    return target_.size();
  }
  virtual int eval(size_t k , const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    itr_matrix<const val_type*> x0(target_.size(),1,x);
    if(k == 0){
        int_type c[1];
        for(size_t i = 0; i < target_.size(); ++i){
            c[0] = i;
            cv[c] += w_*(x0[i] - target_[i]);
          }
      }
    if(k == 1){
        int_type c[2];
        for(size_t i = 0; i < target_.size(); ++i){
            c[0] = i; c[1] = i;
            cv[c] += w_;
          }
      }
    return 0;
  }

  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        int_type c[2];
        for(size_t i = 0; i < target_.size(); ++i){
            c[0] = i; c[1] = i;
            l2g.add(cs,c);
          }
      }
    return 0;
  }

  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1){
        return target_.size();
      }
  }
private:
  const zjucad::matrix::matrix<val_type> &target_;
  const double w_;
};


//! @brief \|Af\|^2
//!
shared_ptr<const hj::math_func::math_func_t<double,int32_t> >
build_smooth_energy(const matrix<double> & sh,
                    const hj::sparse::csc<double,int32_t> &AT)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef shared_ptr<const math_func_type> math_func_ptr;

  math_func_ptr f(new smooth_func<double,int32_t>(AT));
  return math_func_ptr(new hj::math_func::sumsqr<double,int32_t>(f));
}

shared_ptr<const hj::math_func::math_func_t<double,int32_t> >
build_fitting_energy(const matrix<double> &target, const double w)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef shared_ptr<const math_func_type> math_func_ptr;

  math_func_ptr f(new fit_func<double,int32_t>(target,w));

  return math_func_ptr(new hj::math_func::sumsqr<double,int32_t>(f));
}

template <typename val_type, typename int_type>
class smooth_func_zyz: public hj::math_func::math_func_t<val_type, int_type>
{
public:
  smooth_func_zyz(const size_t fi, const size_t fj, const size_t num,
                  const zjucad::matrix::matrix<double> &rot,
                  const double w) :fi_(fi), fj_(fj), num_(num), rot_(rot), w_(w){}
  virtual ~smooth_func_zyz(){}
  virtual size_t nx() const
  {
    return 3*num_;
  }
  virtual size_t nf() const
  {
    return 9;
  }
  virtual int eval(size_t k , const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    itr_matrix<const val_type*> x0(3,num_,x);
    if(k == 0){
        matrix<double> sh_i(9,1), sh_j(9,1);

        calc_rot_cubic_f_sh_(&sh_i[0], &x0(0,fi_));
        calc_rot_cubic_f_sh_(&sh_j[0], &x0(0,fj_));
        sh_i = temp(rot_ * sh_i) - sh_j;
        for(int_type i = 0; i < 9; ++i){
            int_type c[1]= {i};
            val_type d = w_*sh_i[i];
            jtf::math::erase_nan_inf(d);
            cv[c] += d;
          }
      }
    if(k == 1){
        matrix<double> jac_i(9,3), jac_j(9,3);

        calc_jac_rot_cubic_f_sh_(&jac_i[0], &x0(0,fi_));
        jac_i = temp(rot_ * jac_i);
        calc_jac_rot_cubic_f_sh_(&jac_j[0], &x0(0,fj_));
        for(int_type i = 0; i < 9; ++i){
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i,3*fi_ + j};
                val_type d = w_*jac_i(i,j);
                jtf::math::erase_nan_inf(d);
                cv[c] += d;
              }
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] ={i, 3*fj_ + j};
                val_type d = -1*w_*jac_j(i,j);
                jtf::math::erase_nan_inf(d);
                cv[c] += d;
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
                int_type c[2] = {i,3*fi_ + j};
                l2g.add(cs,c);
              }
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i,3*fj_ + j};
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
        return 2*3*9;
      }
  }
private:
  const int_type fi_;
  const int_type fj_;
  const size_t num_;
  const zjucad::matrix::matrix<double> rot_;
  const double w_;
};

shared_ptr<const hj::math_func::math_func_t<double,int32_t> >
build_smooth_zyz_energy(const zjucad::matrix::matrix<double> &zyz,
                        const jtf::mesh::edge2cell_adjacent &ea,
                        const zjucad::matrix::matrix<double> &rij,
                        const zjucad::matrix::matrix<double> &edge_area_weight)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef shared_ptr<const math_func_type> math_func_ptr;

  const double total_area = std::accumulate(edge_area_weight.begin(), edge_area_weight.end(),0.0);
  matrix<double> rot(9,9);
  shared_ptr<vector<math_func_ptr> > funcs(new vector<math_func_ptr>);
  for(size_t ei = 0; ei < ea.edges_.size(); ++ei){
      const pair<size_t,size_t> &face_pair = ea.edge2cell_[ei];
      if(ea.is_boundary_edge(face_pair)) continue;
      calc_rot_cubic_f_sh_mat_(&rot[0], &rij(0,ei));
      funcs->push_back(math_func_ptr(
                         new smooth_func_zyz<double,int32_t>(face_pair.first, face_pair.second,
                                                             zyz.size(2), rot, sqrt(edge_area_weight[ei]))));
    }
  return math_func_ptr(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(funcs));
}

int load_fv(const char *filename,
            vector<zjucad::matrix::matrix<double> >  &fv)
{
  ifstream ifs(filename);
  if(ifs.fail()) {
      cerr << "# [error] can not load face vector file." << endl;
      return __LINE__;
    }
  size_t face_num;
  ifs >> face_num;
  fv.resize(face_num);
  for(size_t fi = 0; fi < face_num; ++fi){
      fv[fi].resize(3,2);
      for(size_t i = 0; i < 6; ++i)
        ifs >> fv[fi][i];
    }
  return 0;
}

void optimal_surface_sh_generator::opt_cos_sin(
    zjucad::matrix::matrix<double> &tradition,
    boost::property_tree::ptree &pt)
{
  init(pt);

  hj::sparse::csc<double,int32_t > A,AT,ATA;
  get_tradition_jacobian_matrix(A);

  //  {
  //    if(zjucad::has("input/init.value",pt)){
  //        vector<matrix<double> > fv;
  //        if(load_fv(pt.get<string>("input/init.value").c_str(), fv)){
  //            throw std::invalid_argument("# [error] can not read init surface field.");
  //          }
  //        matrix<double> cos_sin(A.size(2),1);
  //        assert(fv.size()*2 == cos_sin.size());
  //        vector<pair<double,size_t> > err_idx(4);
  //        for(size_t fi = 0 ; fi < fv.size(); ++fi){
  //            matrix<double> ref = tm_.trimesh_.node_(colon(), tm_.trimesh_.mesh_(1,fi))
  //                - tm_.trimesh_.node_(colon(), tm_.trimesh_.mesh_(0,fi));
  //            ref /= norm(ref);
  //            for(size_t di = 0; di < 4; ++di){
  //                err_idx[di] = make_pair(dot(ref, (di%2==0?1.0:-1.0)*fv[fi](colon(),di/2)),di);
  //              }
  //            sort(err_idx.begin(), err_idx.end());
  //            double angle = jtf::math::safe_acos(err_idx.back().first);
  //            cos_sin[2*fi+0] = cos(4*angle);
  //            cos_sin[2*fi+1] = sin(4*angle);
  //          }
  //        matrix<double> zero_m = zeros<double>(A.size(1),1);
  //        hj::sparse::mv(false, A, cos_sin, zero_m);
  //        cerr << "# [info] init_face_vector_smoothness " << norm(zero_m) << endl;
  //      }
  //  }
  hj::sparse::trans(A,AT);

  hj::sparse::AAT<hj::sparse::map_by_sorted_vector >(AT,ATA);

  solve_smallest_eigenvalue_problem(ATA,tradition);
}

void optimal_surface_sh_generator::opt(zjucad::matrix::matrix<double> &sh,
                                       boost::property_tree::ptree &pt)
{
  zjucad::matrix::matrix<double> tradition_cos_sin;
  opt_cos_sin(tradition_cos_sin,pt);

  sh.resize(9, tm_.trimesh_.mesh_.size(2));
  sh *= 0;
  sh(4,colon()) = ones<double>(1,sh.size(2))*sqrt(7.0);

  itr_matrix<double*> sh_tradition_m(2, tradition_cos_sin.size()/2, &tradition_cos_sin[0]);
  for(size_t i = 0; i < sh_tradition_m.size(2); ++i){
      sh_tradition_m(colon(),i) /= norm(sh_tradition_m(colon(),i));
      sh(0,i) = sh_tradition_m(1,i)*sqrt(5.0);
      sh(8,i) = sh_tradition_m(0,i)*sqrt(5.0);
    }
  {
    jtf::mesh::write_matrix("cos_sin.mat", tradition_cos_sin);
  }
  return ;
}



int optimal_sh_generator::reorder_edge_according_face_vertex(
    size_t * common_edge,
    const zjucad::matrix::matrix<size_t> &one_face)
{
  assert(common_edge);
  auto it = std::find(one_face.begin(), one_face.end(), common_edge[0]);
  assert(it != one_face.end());
  const size_t i = it-one_face.begin();
  if(one_face[(i+1)%one_face.size(1)] != common_edge[1]){
      swap(common_edge[0],common_edge[1]);
      return 1;
    }
  return 0;
}

void optimal_surface_sh_generator::init_connection(boost::property_tree::ptree &pt)
{
  //levi-civitia connection, other connection can be built
  rotation_angle_ = zeros<double>(tm_.ea_->edges_.size(),1);

  pt.put("angle_defect_strategy.desc","manually/even/null");
  const string strategy = pt.get<string>("angle_defect_strategy.value");
  if(strategy == "manually"){
      cerr << "# [info] use manually angle defect." << endl;
      std::unique_ptr<angle_defect> ad(new manually_set_angle_defect);
      ad->opt(*tm_.ea_,tm_.trimesh_.mesh_,tm_.trimesh_.node_, rotation_angle_, pt);
    }else if(strategy == "even"){
      cerr << "# [info] use even angle_defect " << endl;
      std::unique_ptr<angle_defect> ad(new even_angle_defect);
      ad->opt(*tm_.ea_,tm_.trimesh_.mesh_,tm_.trimesh_.node_, rotation_angle_, pt);
    }else if(strategy == "geometry"){
      cerr << "# [info] use geometry aware angle_defect " << endl;
      std::unique_ptr<angle_defect> ad(new geometry_aware_angle_defect);
      ad->opt(*tm_.ea_,tm_.trimesh_.mesh_, tm_.trimesh_.node_, rotation_angle_, pt);
    }else if(strategy == "null"){
      cerr << "# [info] use null angle_defect." << endl;
    }

  rij_ = zeros<double>(3,tm_.ea_->edges_.size());
  for(size_t ei = 0; ei < tm_.ea_->edge2cell_.size(); ++ei){
      const std::pair<size_t,size_t> & face_pair = tm_.ea_->edge2cell_[ei];
      if(tm_.ea_->is_boundary_edge(face_pair)) continue;
      double angle = get_rij_angle(face_pair, tm_.trimesh_.mesh_, tm_.trimesh_.node_,tm_.face_normal_);
      rij_(0,ei) = -1*(angle - rotation_angle_[ei]); // -1 means this angle should be minused during measuring smoothness
    }

  //    {
  //        generate_cross_field_according_to_angle_defet(*tm_.ea_, tm_.trimesh_.mesh_, tm_.trimesh_.node_, rotation_angle_, tm_.face_normal_, "test_fv");
  //    }
}

double optimal_sh_generator::get_rij_angle(const std::pair<size_t,size_t> & face_pair,
                                           const zjucad::matrix::matrix<size_t> & faces,
                                           const zjucad::matrix::matrix<double> & node,
                                           const zjucad::matrix::matrix<double> & face_normal)
{
  assert(face_pair.first != -1 && face_pair.second != -1);
  matrix<double> refi =
      node(colon(), faces(1,face_pair.first)) -
      node(colon(), faces(0,face_pair.first));
  refi /= norm(refi);
  matrix<double> refj =
      node(colon(), faces(1,face_pair.second)) -
      node(colon(), faces(0,face_pair.second));
  refj /= norm(refj);


  const matrix<double> ni = face_normal(colon(), face_pair.first);
  const matrix<double> nj = face_normal(colon(), face_pair.second);

  matrix<double> rot(3,3);
  {
    double dihedral_angle = jtf::math::safe_acos(dot(ni,nj));
    // adjust angle sign
    size_t common_edge[2];
    jtf::mesh::find_common_edge(faces(colon(), face_pair.first),
                                faces(colon(), face_pair.second), common_edge);

    reorder_edge_according_face_vertex(common_edge, faces(colon(),face_pair.first));

    matrix<double> dir = node(colon(), common_edge[1]) - node(colon(), common_edge[0]);

    const zjucad::matrix::matrix<double> cross_ni_nj = cross(ni,nj);
    if(norm(cross_ni_nj) > 1e-8) {// not parallel
        if(dot(dir, cross_ni_nj) < 0) dihedral_angle *= -1;
      }
    from_angle_to_rotation_matrix(dihedral_angle, dir, rot);
  }

  matrix<double> refi_rot = rot * refi;

  double angle = jtf::math::safe_acos(dot(refi_rot, refj));

  if(dot(nj, cross(refi_rot,refj)) < 0) angle *= -1;

  return angle;
}

double optimal_sh_generator::get_rij_angle2(const std::pair<size_t,size_t> & face_pair,
                                            const zjucad::matrix::matrix<size_t> & faces,
                                            const zjucad::matrix::matrix<double> & node,
                                            const zjucad::matrix::matrix<double> & face_normal)
{
  assert(face_pair.first != -1 && face_pair.second != -1);

  const matrix<double> ni = face_normal(colon(), face_pair.first);
  const matrix<double> nj = face_normal(colon(), face_pair.second);

  matrix<double> rot(3,3);
  {
    zjucad::matrix::matrix<double> cross_dir = cross(ni,nj);
    cross_dir /= norm(cross_dir);

    matrix<double> ni_old,nj_old;
    jtf::mesh::cal_face_normal(faces(colon(), face_pair.first), node, ni_old);
    jtf::mesh::cal_face_normal(faces(colon(), face_pair.second), node, nj_old);

    double dihedral_angle = jtf::math::safe_acos(dot(ni,nj));
    // adjust angle sign
    size_t common_edge[2];
    jtf::mesh::find_common_edge(faces(colon(), face_pair.first),
                                faces(colon(), face_pair.second), common_edge);

    reorder_edge_according_face_vertex(common_edge, faces(colon(),face_pair.first));

    matrix<double> dir = node(colon(), common_edge[1]) - node(colon(), common_edge[0]);

    const zjucad::matrix::matrix<double> cross_ni_nj = cross(ni_old,nj_old);
    if(norm(cross_ni_nj) > 1e-8) {// not parallel
        if(dot(dir, cross_ni_nj) < 0) dihedral_angle *= -1;
      }
    from_angle_to_rotation_matrix(dihedral_angle, cross_dir, rot);
  }

  matrix<double> refi_rot = rot * reference_dir_(colon(),face_pair.first);

  double angle = jtf::math::safe_acos(dot(refi_rot, reference_dir_(colon(),face_pair.second)));

  if(dot(nj, cross(refi_rot,reference_dir_(colon(), face_pair.second))) < 0) angle *= -1;

  return angle;
}

void optimal_surface_sh_generator::get_tradition_jacobian_matrix(
    hj::sparse::csc<double, int32_t> &A) const

{
  size_t edges_num = 0;

  for(size_t ei = 0; ei < tm_.ea_->edge2cell_.size(); ++ei){
      if(tm_.ea_->is_boundary_edge(tm_.ea_->edge2cell_[ei]))  continue;
      ++edges_num;
    }
  csc_filler<double> cf(2*edges_num, tm_.trimesh_.mesh_.size(2)*2);

  cerr << "# edges num " << edges_num << endl;

  edges_num = 0;

  const size_t N_sys = 4;

  matrix<double> rot(2,2);
  for(size_t ei = 0; ei < tm_.ea_->edge2cell_.size(); ++ei){
      const std::pair<size_t,size_t> & face_pair = tm_.ea_->edge2cell_[ei];
      if(tm_.ea_->is_boundary_edge(face_pair)) continue;
      const size_t &fi_idx = face_pair.first;
      const size_t &fj_idx = face_pair.second;
      const double w = sqrt(edge_area_weight_[ei]);

      get_2d_rotation(N_sys*rij_(0,ei),rot);

      for(size_t i = 0; i< 2; ++i){
          for(size_t j = 0; j < 2; ++j){
              cf.add(2*edges_num+i,2*fi_idx+j,w*rot(i,j));
            }
          cf.add(2*edges_num+i,2*fj_idx+i,-1*w);
        }
      ++edges_num;
    }
  cf.out(A);
}

void optimal_surface_sh_generator::get_jacobian_matrix(hj::sparse::csc<double,int32_t>&A) const
{
  matrix<double> sh_rot(9,9);
  size_t edges_num = 0;

  for(size_t ei = 0; ei < tm_.ea_->edge2cell_.size(); ++ei){
      if(tm_.ea_->is_boundary_edge(tm_.ea_->edge2cell_[ei]))  continue;
      ++edges_num;
    }
  csc_filler<double> cf(9*edges_num, tm_.trimesh_.mesh_.size(2)*9);

  cerr << "# edges num " << edges_num << endl;

  edges_num = 0;

  for(size_t ei = 0; ei < tm_.ea_->edge2cell_.size(); ++ei){
      const std::pair<size_t,size_t> & face_pair = tm_.ea_->edge2cell_[ei];
      if(tm_.ea_->is_boundary_edge(face_pair)) continue;
      const size_t &fi_idx = face_pair.first;
      const size_t &fj_idx = face_pair.second;
      const double w = sqrt(edge_area_weight_[ei]);
      calc_rot_cubic_f_sh_mat_(&sh_rot[0], &rij_(0,ei));
      // fill R_ij, and -I
      for(size_t i = 0; i < sh_rot.size(1); ++i){
          for(size_t j = 0; j < sh_rot.size(2); ++j){
              cf.add(9*edges_num + i, 9*fi_idx+j, w*sh_rot(i,j));
            }
          cf.add(9*edges_num+i,9*fj_idx+i, -1.0*w);
        }
      ++edges_num;
    }
  cf.out(A);
}



template <typename val_type, typename int_type>
class end_fit_func : public hj::math_func::math_func_t<double,int32_t>
{
public:
  end_fit_func(const size_t i, const size_t j, const double w): i_(i), j_(j),w_(w){
    assert(w_ > 0);
  }
  virtual ~end_fit_func(){}
  virtual size_t nx() const
  {
    return 9;
  }
  virtual size_t nf() const
  {
    return 1;
  }
  virtual int eval(size_t k , const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        int_type c[1] = {0};
        cv[c] += w_*(x[i_]*x[i_]+x[j_]*x[j_]-5.0);
      }
    if(k == 1){
        int_type c1[2]={0,i_};
        cv[c1] += w_*2*x[i_];

        int_type c2[2]={0,j_};
        cv[c2] += w_*2*x[j_];
      }
    if(k == 2){
        int_type c1[3]={0,i_,i_};
        cv[c1] += 2*w_;
        int_type c2[3] ={0,j_,j_};
        cv[c2] += 2*w_;
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        int_type c1[2]={0,i_};
        l2g.add(cs, c1);

        int_type c2[2]={0,j_};
        l2g.add(cs, c2);
      }
    if(k == 2){
        int_type c1[3]={0,i_,i_};
        l2g.add(cs, c1);

        int_type c2[3] ={0,j_,j_};
        l2g.add(cs, c2);
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 2;
    if(k == 2) return 2;
  }
private:
  const int_type i_, j_;
  const double w_;
};

template <typename val_type, typename int_type>
class fix_func_sh : public hj::math_func::math_func_t<double,int32_t>
{
public:
  fix_func_sh(const size_t i, val_type w, val_type c):i_(i),w_(w),c_(c){}
  virtual ~fix_func_sh(){}
  virtual size_t nx() const{
    return 9;
  }
  virtual size_t nf()const{
    return 1;
  }
  virtual int eval(size_t k , const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        int_type c[1] = {0};
        cv[c] += c_*(x[i_] - w_);
      }
    if(k == 1){
        int_type c[2] = {0,i_};
        cv[c] += c_;
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        int_type c[2] = {0,i_};
        l2g.add(cs,c);
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 1;
  }
private:
  const int_type i_;
  const double w_;
  const double c_;
};

template <typename val_type, typename int_type>
class fit_init_sh_func : public hj::math_func::math_func_t<double,int32_t>
{
public:
  fit_init_sh_func(const zjucad::matrix::matrix<double> &sh):sh_(sh){}
  virtual ~fit_init_sh_func(){}
  virtual size_t nx() const
  {
    return 9;
  }
  virtual size_t nf() const
  {
    return 9;
  }
  virtual int eval(size_t k , const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        for(int_type i = 0; i < 9; ++i){
            int_type c[1] = {i};
            cv[c] += x[i]-sh_[i];
          }
      }
    if(k == 1){
        for(int_type i = 0; i < 9; ++i){
            int_type c[2] = {i,i};
            cv[c] += 1.0;
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
            int_type c[2] = {i,i};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 9;
  }
private:
  const matrix<double> sh_;
};

void optimal_surface_sh_generator::projection_z_axis_rotaiton_each_cell(
    zjucad::matrix::matrix<double> &sh)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef shared_ptr<math_func_type> math_func_ptr;

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  all_func->push_back(math_func_ptr(new end_fit_func<double,int32_t>(0,8,sqrt(10))));
  for(size_t i = 1; i < 8; ++i){
      if(i != 4)
        all_func->push_back(math_func_ptr(new fix_func_sh<double,int32_t>(i, 0.0, sqrt(50.0))));
      else
        all_func->push_back(math_func_ptr(new fix_func_sh<double,int32_t>(i, sqrt(7.0),sqrt(20))));
    }
  all_func->push_back(math_func_ptr(new fit_init_sh_func<double,int32_t>(sh)));

  math_func_ptr all_func_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  math_func_ptr all_f(new hj::math_func::sumsqr<double,int32_t>(all_func_cat));

  boost::property_tree::ptree pt;
  pt.put("package.value", "jtf");
  pt.put("alg.value","SQP");
  pt.put("iter.value",10);

  jtf::optimize(*all_f, sh, pt, nullptr, nullptr, nullptr);
}

void optimal_surface_sh_generator::projection_z_axis_rotaiton(zjucad::matrix::matrix<double> &sh)
{
  itr_matrix<double*> sh_m(9,sh.size()/9,&sh[0]);

  matrix<double> sh_temp;
  for(size_t i = 0; i < sh_m.size(2); ++i){

      cerr << "# i/total " << i << "/" << sh_m.size(2)<< endl;
      sh_m(colon(),i) /= norm(sh_m(colon(),i));
      sh_temp = sh_m(colon(),i)*sqrt(12);

      projection_z_axis_rotaiton_each_cell(sh_temp);
      sh_m(colon(),i) = sh_temp;
    }
}

void optimal_surface_sh_generator::projection(zjucad::matrix::matrix<double> &sh) const
{
  zjucad::matrix::itr_matrix<double*> sh_m(9,sh.size()/9,&sh[0]);

  zjucad::matrix::matrix<double> f_temp(sh_m.size(1),1);
  double e4_raw[9]= {0,0,0,0,1,0,0,0,0};
  itr_matrix<const double*> e4(9,1,&e4_raw[0]);
  for(size_t i = 0; i < sh_m.size(2); ++i){
      f_temp = sh_m(colon(), i) - dot(sh_m(colon(),i),e4)*e4;
      sh_m(colon(),i) = f_temp + sqrt(7.0/5.0)*zjucad::matrix::norm(f_temp) *e4;
      sh_m(colon(),i) /= norm(sh_m(colon(),i));
      if(fabs(dot(sh_m(colon(),i), e4) - sqrt(7.0/12.0)) > 1e-8){
          cerr << "# [info] sh " << i << " error " << endl;
        }
    }
  sh *= sqrt(12.0);
}

int check_tetmesh(const jtf::tet_mesh &tm)
{
  jtf::tet_mesh a(tm.tetmesh_);

  //  zjucad::matrix::matrix<size_t> outside_face_;
  //  zjucad::matrix::matrix<size_t> outside_face_idx_;
  //  zjucad::matrix::matrix<double> outside_face_normal_;
  //  zjucad::matrix::matrix<double> outside_face_area_;
  //  std::shared_ptr<jtf::mesh::face2tet_adjacent> fa_;
  //  std::shared_ptr<jtf::mesh::edge2cell_adjacent> ea_outside_;
  //  jtf::mesh::one_ring_tet_at_edge ortae_;

  cerr << "# [info] vol diff " << norm(a.vol_ - tm.vol_) << endl;
  cerr << "# [info] outside_face diff " << norm(a.outside_face_ - tm.outside_face_) << endl;
  cerr << "# [info] outside_face_idx diff " << norm(a.outside_face_idx_ - tm.outside_face_idx_) << endl;
  cerr << "# [info] outside_face_normal diff " << norm(a.outside_face_normal_ - tm.outside_face_normal_) << endl;

  for(const auto & e : tm.vol_){
      if(jtf::math::is_nan(e)) {
          cerr << "# vol nan. " << endl;
        }
      if(e < 0 || fabs(e) < 1e-6)
        cerr << "# negative vol" << endl;
    }

  for(const auto & e : tm.outside_face_normal_)
    if(jtf::math::is_nan(e)) {
        cerr << "# normal nan. " << endl;
      }

  return 0;
}

void optimal_sh_generator_seq::deformation_project_sh(
    const jtf::tet_mesh &tm_current,
    const jtf::tet_mesh & tm_target,
    matrix<double> &sh)
{
  // sh in this tm should be mapped to target.
  matrix<matrix<double> > deform_gradient;
  cal_deformation_gradient(tm_current.tetmesh_.mesh_, tm_current.tetmesh_.node_,
                           tm_target.tetmesh_.node_, deform_gradient);

  matrix<double> zyz;
  project2zyz(sh, zyz);
  matrix<double> R;
  matrix<double>  rot_m(3,3);
  for(size_t ti = 0; ti < rot_m.size(); ++ti){
      zyz_angle_2_rotation_matrix1(&zyz(0,ti), &rot_m[0]);
      R = deform_gradient[ti];
      hj::polar3d p;
      p(R,2);
      rot_m = temp(R*rot_m);
      rotation_matrix_2_zyz_angle(&rot_m[0], &zyz(0,ti),0);
      calc_rot_cubic_f_sh_(&sh(0,ti), &zyz(0,ti));
    }
}

void optimal_sh_generator_seq::opt0(zjucad::matrix::matrix<double> &field,
                                    boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  vector<matrix<double> > sh(tm_seq_.size());
  pt.put("input/opt_type.desc", "init or zyz");
  const string opt_type = pt.get<string>("input/opt_type.value");

  for(size_t i = 0; i < tm_seq_.size(); ++i){
      stringstream ss;
      ss << "sh_" << time_step_[i]  << ".sh";
      if(jtf::mesh::read_matrix(ss.str().c_str(), sh[i]) == 0) {
          cerr << "# [info] load " << ss.str() << endl;
          continue;
        }
      optimal_sh_generator osg(tm_seq_[i]);
      osg.opt(sh[i], pt);
      deformation_project_sh(tm_seq_[i], tm_seq_[0], sh[i]);
      osg.write_sh(ss.str().c_str(), sh[i]);
    }

  matrix<double> point_size_field, face_size_field;
  compute_size_field(point_size_field, pt);
  {
    face_size_field = zeros<double>(tm_.fa_->faces_.size(),1);
    for(size_t fi = 0; fi < tm_.fa_->faces_.size(); ++fi){
        const vector<size_t> &one_face = tm_.fa_->faces_[fi];
        for(size_t pi = 0; pi < one_face.size(); ++pi){
            face_size_field[fi] += point_size_field[one_face[pi]];
          }
        face_size_field[fi] /= 3;
      }
  }

  matrix<double> tet_size = zeros<double>(tm_.tetmesh_.mesh_.size(2),1);
  for(size_t fi = 0; fi < tm_.fa_->faces_.size(); ++fi){
      const pair<size_t,size_t>  &tet_pair = tm_.fa_->face2tet_[fi];
      if(tet_pair.first != -1)
        tet_size[tet_pair.first] += face_size_field[fi];
      if(tet_pair.second != -1)
        tet_size[tet_pair.second] += face_size_field[fi];
    }

  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);
  func->push_back(math_func_ptr(
                    new hj::math_func::sumsqr<double,int32_t>(
                      build_fix_inner_field_func(sh, time_step_, tm_.tetmesh_.mesh_,
                                                 tm_.tetmesh_.node_,
                                                 tm_.vol_,tet_size, *(tm_.fa_)))));

  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  field = sh[0];
  jtf::optimize(*obj, field, pt, nullptr, nullptr, nullptr);

  if(opt_type == "init"){
      for(size_t i = 0; i < field.size(2); ++i){
          field(colon(),i) /= norm(field(colon(),i));
        }
      field *= sqrt(12.0);
    }
}
