#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <map>
#include <omp.h>

#include <boost/property_tree/ptree.hpp>
#include <minpack.h>

#include <jtflib/mesh/io.h>
#include <jtflib/math/math.h>
#include <jtflib/function/function.h>
#include <jtflib/function/func_aux.h>
#include <jtflib/optimizer/optimizer.h>
#include <hjlib/sparse/sparse.h>
#include <hjlib/function/function.h>
#include <zjucad/optimizer/optimizer.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/ptree/ptree.h>

#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/lapack.h>

#include "../hex_frame_opt/function_term.h"
#include "../common/zyz.h"
#include "../common/IO.h"
#include "../common/util.h"
#include "../common/vtk.h"
#include "../common/visualize_tool.h"
#include "../common/transition.h"
#include "../common/transition_type.h"
#include "../common/face_orient.h"
#include "../common/solver_ipopt.h"
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
#include "../tetmesh/util.h"
#include "../hex_param/hex_param.h"

#include "../hex_frame_opt/frame_function.h"
#include "../hex_frame_opt/util.h"
#include "../hex_frame_opt/hex_frame_opt.h"
#include "./common.h"
#include "../hex_param/io.h"
#include "../hex_param/cut_tet.h"
#include "../hex_param/hex_param.h"
#include "../hex_param/common.h"
#include "../hex_param/find_singularities.h"
#include "../hex_param/global_alignment.h"
#include "../hex_param/topology_analysis.h"
#include "../hex_param/topology_operation.h"
#include "../hex_param/wyz_format.h"
#include "../hex_param/remove_surface_wedge.h"

#include <jtflib/mesh/mesh.h>
#include "../hexmesh/io.h"
#include "../hex_process/hex_process.h"

#include "../tetmesh_refine/tetmesh_refine.h"
#include "../spherical_harmonics/rot_cubic_f_SH.h"

#include <jtflib/mesh/trimesh.h>
#include "../hex_frame_opt/optimal_sh_generator.h"

#include "../crane_fairing/Application.h"
#include "../crane_fairing/Application2.h"
#include "../crane_fairing/Mesh.h"

#include <hjlib/math_func/math_func.h>
#include <jtflib/optimizer/optimizer.h>
#include <hjlib/math_func/operation.h>

#include "../vol_param/descriptor/func_terms/arap.h"
#include <hjlib/math/polar.h>

#include "../hex_frame_opt/optimize_frame_field_with_type.h"
#include "../hex_frame_opt/map_tets.h"

using namespace std;
using namespace zjucad::matrix;
using namespace hj::sparse;
using namespace hj::function;

using boost::property_tree::ptree;

int init_zyz_inner_with_surface_tri(ptree &pt)
{
  jtf::mesh::meshes tri;

  if(jtf::mesh::load_obj(pt.get<string>("obj.value").c_str(), tri.mesh_, tri.node_))
    return __LINE__;

  shared_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(tri.mesh_));
  if(!ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }
  matrixd zyz = zeros<double>(3, tri.mesh_.size(2));
  matrixd sh(9, tri.mesh_.size(2));

  pt.put("align_w.desc","weight of align the surface normal or surface frame.");
  const double align = pt.get<double>("align_w.value");

  pt.put("init_zyz.desc","tet zyz frame file");
  if(zjucad::has("init_zyz.value",pt)) {
      ifstream ifs(pt.get<string>("init_zyz.value").c_str(), ifstream::binary);
      if(jtf::mesh::read_matrix(ifs, zyz, 3)) {
          cerr << "# not a zyz file." << endl;
          return __LINE__;
        }
      if(zyz.size(2) != tri.mesh_.size(2)) {
          cerr << "incompatible inner zyz file." << endl;
          return __LINE__;
        }
    }

  shared_ptr<vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > > funcs(
        new vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > >);

  matrix<double> normal(3, tri.mesh_.size(2));
  jtf::mesh::cal_face_normal(tri.mesh_, tri.node_, normal);

  {
    // surface smooth
    pt.put("surface_smooth.desc","weight for surface smoothing");
    const double w = pt.get<double>("surface_smooth.value");
    if(w > 1e-6){
        std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > surface_smooth_func
            =  build_frame_surface_smooth_function_jtf_tri(
              tri.mesh_.size(2),tri.mesh_, tri.node_, normal, *ea, "sh", w, pt);
        if(surface_smooth_func.get()){
            funcs->push_back(surface_smooth_func);
            cerr << "# [usage] add surface smooth function. " << endl;
          }
      }
  }

  pt.put("strategy.desc", "surface_normal/surface_frame, [normal/frame]");
  const string strategy = pt.get<string>("strategy.value");

  const size_t normal_smooth_iter = 0;
  matrix<double> face_normal = normal;
  if(normal_smooth_iter > 0){
      matrixd point_normal;
      jtf::mesh::cal_point_normal(tri.mesh_, tri.node_, point_normal);
      std::shared_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(tri.mesh_));
      jtf::tetmesh::one_ring_point_at_point orpap;
      orpap.build(ea->edges_);
      jtf::tetmesh::smooth_one_ring_point_normal(point_normal, orpap.p2p_, normal_smooth_iter);
      for(size_t fi = 0; fi < tri.mesh_.size(2); ++fi){
          face_normal(colon(),fi) *= 0;
          for(size_t pi = 0; pi < tri.mesh_.size(1); ++pi){
              face_normal(colon(),fi) += point_normal(colon(),tri.mesh_(pi,fi));
            }
          const double len = norm(face_normal(colon(),fi));
          if(len > 1e-6)
            face_normal(colon(),fi) /= len;
        }
      cerr << "# [info] finish point normal smooth." << endl;
    }

  if(strategy == "normal"){
      std::shared_ptr<function_t<double,int32_t> > normal_align_func(
            build_normal_align_func_tri(
              tri.mesh_, tri.node_, normal, "sh" ,align));
      funcs->push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                         jtf::function::least_square_warpper(normal_align_func)));
      cerr << "# [info] add normal align function." << endl;
    }

  // initial value
  for(size_t i = 0; i < zyz.size(2); ++i)
    calc_rot_cubic_f_sh_(&sh(0, i), &zyz(0, i));

  cerr << sh(colon(), colon(0,8)) << endl;
  cerr << zyz(colon(), colon(0,8)) << endl;
  cerr << "# begin to compute" << endl;

  clock_t beg = clock();

  shared_ptr<jtf::function::functionN1_t<double,int32_t> > all_f(
        new jtf::function::sum_function<double,int32_t, jtf::function::SMART_STD>(*funcs));
  jtf::optimize(*all_f, sh, pt, nullptr, nullptr, nullptr);
  //ipopt_solve(sh, *all_f, nullptr,pt);


  cout << "time for lsq opt: " << (clock()-beg)/double(CLOCKS_PER_SEC) << endl;

  cerr << "# SH 2 zyz." << endl;

  zyz = zeros<double>(3, sh.size(2));
  cerr << "# zyz size " << zyz.size(1) << "," << zyz.size(2) << endl;

  size_t i;
  beg = clock();

#pragma omp parallel for private(i)
  for(i = 0; i < zyz.size(2); ++i) {
      SH2zyz_lmder1(&zyz(0, i), &sh(0, i), 100);  // return 2 for healthy
    }
  cout << "time for SH2zyz: " << (clock()-beg)/double(CLOCKS_PER_SEC) << endl;

  if(zjucad::has("output.value",pt)){
      ofstream ofs(pt.get<string>("output.value").c_str(), ofstream::binary);
      jtf::mesh::write_matrix(ofs, zyz);
      cout << "zyz: " << zyz(colon(), colon(0, 6));
    }
  return 0;
}

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

int frame_inner_with_surface(ptree &pt)
{
  jtf::tet_mesh tm(pt.get<string>("input/tet.value").c_str());

  optimal_sh_generator osg(tm);
  zjucad::matrix::matrix<double> zyz = zeros<double>(3,tm.tetmesh_.mesh_.size(2));
  if(zjucad::has("input/init_zyz.value",pt)){
      if(jtf::mesh::read_matrix(pt.get<string>("input/init_zyz.value").c_str(), zyz)){
          cerr << "# [error] can not load init_zyz file." << endl;
          return __LINE__;
        }
      if(zyz.size(1) != 3 || zyz.size(2) != tm.tetmesh_.mesh_.size(2)){
          cerr << "# [error] init zyz is not compatible with mesh." << endl;
          return __LINE__;
        }

      cerr << "# [info] load init zyz." << endl;
    }

  osg.opt(zyz,pt);

  matrix<double> zyz_bkp = zyz;
  hj_frame_alignemt(tm.tetmesh_.mesh_, *tm.fa_, zyz_bkp, zyz);

  ofstream ofs(pt.get<string>("output.value").c_str(), ofstream::binary);
  jtf::mesh::write_matrix(ofs, zyz);

  return 0;
}

int frame_inner_with_type(ptree &pt)
{
  jtf::tet_mesh tm(pt.get<string>("input/tet.value").c_str());
  boost::unordered_map<pair<size_t,size_t>,size_t> inner_type;
  if(load_inner_face_jump_type(pt.get<string>("input/inner_type.value").c_str(), inner_type)){
      return __LINE__;
    }

  matrix<double> zyz;
  if(jtf::mesh::read_matrix(pt.get<string>("input/zyz.value").c_str(),zyz)){
      cerr << "# [error] can not open zyz." << endl;
      return __LINE__;
    }else{
      if(zyz.size(2) != tm.tetmesh_.mesh_.size(2)){
          cerr << "# [error] not compatible zyz." << endl;
          return __LINE__;
        }
    }
  optimize_frame_field_with_type(tm, inner_type, zyz);
  jtf::mesh::write_matrix(pt.get<string>("output/zyz.value").c_str(), zyz);
  return 0;
}

int init_zyz_inner_with_surface2(ptree &pt)
{
  jtf::tet_mesh tm(pt.get<string>("input/tet.value").c_str());

  optimal_sh_generator osg(tm);
  zjucad::matrix::matrix<double> sh;
  osg.opt(sh,pt);

  osg.write_zyz(pt.get<string>("output.value").c_str(),sh);

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
            int_type c[2] = {fi, static_cast<int_type>(3*idx_+fi)};
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
            int_type c[2] = {fi, static_cast<int_type>(3*idx_+fi)};
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
                int_type c[2] = {fi, static_cast<int_type>(one_tet_[pi]*3+r)};
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
                int_type c[2] = {fi, static_cast<int_type>(one_tet_[pi]*3+r)};
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


void fairing_surface(jtf::tet_mesh &tm, const double step,
                     boost::property_tree::ptree &pt)
{
  using namespace zjucad::matrix;

  {
    stringstream tet_output;
    tet_output << "deformed_tet_" << step << ".tet";
    if(tm.load(tet_output.str().c_str()) == 0){
        cerr << "# [info] load tet " << tet_output.str() << endl;
        return;
      }
  }
  matrix<size_t> v2s;
  matrix<size_t> face = tm.outside_face_;
  matrix<double> node = tm.tetmesh_.node_;

  remove_extra_node(face,node, &v2s);

  double r1 = calc_bounding_sphere_size(node);
  DDG::Mesh ddg_mesh;
  ddg_mesh.read(&face[0], face.size(2), &node[0], node.size(2));
  DDG::Application a;
  a.run(step, ddg_mesh);

  ddg_mesh.write(&face[0], &node[0]);
  double r2 = calc_bounding_sphere_size(node);
  node *= r1/r2;

  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);
  func->push_back(math_func_ptr(
                    new hj::math_func::sumsqr<double,int32_t>(
                      build_fix_position_func(tm.tetmesh_.node_, v2s, node, 1))
                    ));
  func->push_back(math_func_ptr(
                    new hj::math_func::sumsqr<double,int32_t>(
                      build_arap_func(tm))));

  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  jtf::optimize(*obj, tm.tetmesh_.node_, pt, nullptr, nullptr, nullptr);

  tm.load(tm.tetmesh_);
  {
    stringstream tet_output;
    tet_output << "deformed_tet_" << step;
    jtf::mesh::tet_mesh_write_to_zjumat((tet_output.str()+".tet").c_str(), &tm.tetmesh_.node_, &tm.tetmesh_.mesh_);
  }
}

int init_zyz_inner_with_surface3(ptree &pt)
{
  jtf::tet_mesh tm(pt.get<string>("input/tet.value").c_str());

  pt.put("weight/seq_num.desc", "how many seq shape do you need, default is 3");
  const size_t seq_num = pt.get<size_t>("weight/seq_num.value",3);
  vector<jtf::tet_mesh> tm_seq(seq_num);

  cerr << "# [info] weight/seq_num = " << seq_num << endl;
  double time_step = 0;
  pt.put("weight/delta_time.desc","time step for deformation sequence");
  const double delta_time = pt.get<double>("weight/delta_time.value",0.05);
  cerr << "# [info] wegiht/delta_time  = " << delta_time << endl;
  vector<double> time_step_seq(seq_num);

  for(size_t i = 0; i  < seq_num; ++i){
      tm_seq[i] = tm;
      fairing_surface(tm_seq[i], time_step, pt);
      time_step_seq[i] = time_step;
      time_step += delta_time;
    }

  optimal_sh_generator_seq osg(tm_seq, time_step_seq);
  zjucad::matrix::matrix<double> sh;
  osg.opt(sh,pt);

  osg.write_zyz(pt.get<string>("output.value").c_str(), sh);

  return 0;
}


int load_seq_tet_files(const char * file, vector<jtf::tet_mesh> & tm_seq)
{
  ifstream ifs(file);
  if(ifs.fail()){
      cerr << "# [error] can not load_seq_tet_file." << endl;
      return __LINE__;
    }
  string path;
  while(!ifs.eof()){
      ifs >> path;
      if(ifs.eof()) return 0;
      jtf::tet_mesh a(path.c_str());
      tm_seq.push_back(a);
    }
  return 0;
}

int init_zyz_inner_with_surface4(ptree &pt)
{
  vector<jtf::tet_mesh> tm_seq;

  pt.put("input/seq_tet_path_file.desc", "a file store the seq tet path");
  const string input_seq_tet_file_dir = pt.get<string>("input/seq_tet_path_file.value");

  if(load_seq_tet_files(input_seq_tet_file_dir.c_str(), tm_seq)){
      cerr << "# [error] load seq tet fail." << endl;
      return __LINE__;
    }

  vector<double> time_step_seq(tm_seq.size());

  for(size_t ti = 0; ti < time_step_seq.size(); ++ti)
    time_step_seq[ti] = ti;

  optimal_sh_generator_seq osg(tm_seq, time_step_seq);
  zjucad::matrix::matrix<double> sh;
  osg.opt(sh,pt);

  osg.write_zyz(pt.get<string>("output.value").c_str(), sh);

  return 0;
}

int init_zyz_inner_with_surface_tri2(ptree &pt)
{
  jtf::mesh::tri_mesh tri_mesh(pt.get<string>("input/obj.value").c_str());

  optimal_surface_sh_generator ossg(tri_mesh);
  zjucad::matrix::matrix<double> sh;
  ossg.opt(sh,pt);

  zjucad::matrix::matrix<double> zyz = zeros<double>(3,sh.size()/9);

  ossg.project2zyz(sh,zyz);

  convert_local_rot_2_global_rot(zyz, tri_mesh);

  ofstream ofs(pt.get<string>("output.value").c_str(), ofstream::binary);
  jtf::mesh::write_matrix(ofs, zyz);
  return 0;
}

int fairing_obj(ptree &opt)
{
  DDG::Mesh mesh;
  mesh.read(opt.get<string>("input/obj.value"));
  DDG::Application2 a;
  a.run(0.01, mesh);
  mesh.write(opt.get<string>("output/obj.value"));
  return 0;
}

int fairing_tet(ptree &pt)
{
  pt.put("input/tet.desc", "input tet mesh");
  pt.put("weight/time_step.desc", "delta time step");
  pt.put("output/tet.desc", "output tetmesh");

  jtf::tet_mesh tm(pt.get<string>("input/tet.value").c_str());
  const double time_step = pt.get<double>("weight/time_step.value");
  fairing_surface(tm, time_step, pt);

  jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("output/tet.value").c_str(),
                                      &tm.tetmesh_.node_, &tm.tetmesh_.mesh_);
  return 0;
}

int adaptive_frame_field_generation_with_polycube(ptree &pt)
{
  jtf::tet_mesh tm(pt.get<string>("input/tet.value").c_str());
  jtf::tet_mesh tm_polycube(pt.get<string>("input/polycube_tet.value").c_str());

  optimal_sh_generator_polycube osgp(tm,tm_polycube);
  zjucad::matrix::matrix<double> field;
  osgp.opt(field,pt);

  osgp.write_zyz(pt.get<string>("output/zyz.value").c_str(), field, true);
  return 0;
}

static int load_guiding_field(const char * file, matrix<double> & fv0, matrix<double> &fv1)
{
  ifstream ifs(file);
  if(ifs.fail()){
      cerr << "# [error] can not open guiding field." << endl;
      return __LINE__;
    }
  size_t face_num;
  ifs >> face_num;
  if(face_num != fv0.size(2)){
      fv0.resize(3, face_num);
      fv1.resize(3, face_num);
    }
  for(size_t fi = 0; fi < face_num; ++fi){
      ifs >> fv0(0,fi) >> fv0(1,fi) >> fv0(2,fi);
      ifs >> fv1(0,fi) >> fv1(1,fi) >> fv1(2,fi);
    }
  return 0;
}

int inteploate_frame_field(ptree &pt)
{
  jtf::tet_mesh tm(pt.get<string>("input/tet.value").c_str());
  matrix<double> fv0,fv1;
  load_guiding_field(pt.get<string>("input/guiding_field.value").c_str(), fv0, fv1);

  matrix<matrix<double> > frame(tm.tetmesh_.mesh_.size(2),1);
  for(size_t fi = 0; fi < tm.outside_face_idx_.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = tm.fa_->face2tet_[tm.outside_face_idx_[fi]];
      const size_t tet_idx = tet_pair.first==-1?tet_pair.second:tet_pair.first;
      frame[tet_idx].resize(3,3);
      frame[tet_idx](colon(),0) = fv0(colon(),fi);
      frame[tet_idx](colon(),2) = tm.outside_face_normal_(colon(),fi);
      frame[tet_idx](colon(),1) =
          cross(frame[tet_idx](colon(),2), frame[tet_idx](colon(),0));
    }


  deque<size_t> que;
  for(size_t fi = 0; fi < tm.outside_face_idx_.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = tm.fa_->face2tet_[tm.outside_face_idx_[fi]];
      const size_t tet_idx = tet_pair.first==-1?tet_pair.second:tet_pair.first;
      que.push_back(tet_idx);
    }

  while(!que.empty()){
      const size_t tet_idx = que.front();
      que.pop_front();

      for(size_t i = 0; i < tm.tetmesh_.mesh_.size(1); ++i){
          const size_t other_face_idx =
              tm.fa_->get_face_idx(tm.tetmesh_.mesh_(i, tet_idx),
                                   tm.tetmesh_.mesh_((i+1)%4, tet_idx),
                                   tm.tetmesh_.mesh_((i+2)%4, tet_idx));
          const pair<size_t,size_t> & tet_pair = tm.fa_->face2tet_[other_face_idx];
          const size_t other_tet_idx = tet_pair.first+tet_pair.second-tet_idx;
          if(other_tet_idx == -1) continue;
          if(frame[other_tet_idx].size() != 0) continue;
          frame[other_tet_idx] = frame[tet_idx];
          que.push_back(other_tet_idx);
        }
    }

  matrix<double> zyz(3,frame.size());
  for(size_t ti = 0 ; ti < frame.size(); ++ti){
      rotation_matrix_2_zyz_angle(&frame[ti][0], &zyz(0,ti), 0);
    }

  matrix<double> zyz_bkp = zyz;
  hj_frame_alignemt(tm.tetmesh_.mesh_,*tm.fa_,zyz_bkp,zyz);
  jtf::mesh::write_matrix(pt.get<string>("output/zyz.value").c_str(), zyz);

  return 0;
}

int map_tets(ptree &pt)
{
  pt.put("input/tet0.desc", "map tet0 to tet1, input tet0 name");
  pt.put("input/tet1.desc", "map tet0 to tet1, input tet1 name");

  jtf::tet_mesh tm0(pt.get<string>("input/tet0.value").c_str());
  jtf::tet_mesh tm1(pt.get<string>("input/tet1.value").c_str());

  if(fabs(norm(tm0.tetmesh_.mesh_-tm1.tetmesh_.mesh_)) > 1e-6){
      throw std::invalid_argument("tet mesh is not compatible.");
    }
  map_tets(tm0, tm1, pt);
  jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("output/tet.value").c_str(),
                                      &tm0.tetmesh_.node_, &tm0.tetmesh_.mesh_);
  return 0;
}

int map_tris(ptree &pt)
{
  using namespace zjucad::matrix;

  pt.put("input/tri0.desc", "map tri0 to tri1, input tri0 name");
  pt.put("input/tri1.desc", "map tri0 to tri1, input tri1 name");

  matrix<size_t> tri0, tri1;
  matrix<double> tri0_node, tri1_node;


  if(jtf::mesh::load_obj(pt.get<string>("input/tri0.value").c_str(), tri0, tri0_node))
    return __LINE__;

  if(jtf::mesh::load_obj(pt.get<string>("input/tri1.value").c_str(), tri1, tri1_node))
    return __LINE__;

  if(fabs(norm(tri0-tri1)) > 1e-6){
      throw std::invalid_argument("tet mesh is not compatible.");
    }
  map_tris(tri0, tri0_node, tri1_node, pt);

  return 0;
}
