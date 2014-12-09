#include <zjucad/ptree/ptree.h>
#include <zjucad/matrix/io.h>

#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/function/function.h>
#include <jtflib/optimizer/optimizer.h>

#include <hjlib/function/function.h>

#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"
#include "../mesh_func/tri-area-normal-func.h"
#include "polycube_surface_func.h"
#include "quality.h"

#include <utility>
#include <map>

extern double arap_w;
extern double adj_normal_w;

using namespace zjucad::matrix;
using namespace std;

using boost::property_tree::ptree;

class MSG_LOG
{
public:
  template <typename OS>
  void out(OS &os){os << std::endl;}

  template <typename OS, typename T, typename... Param>
  void out(OS &os, T &&val, Param && ... param){
    os << val << " ";
    out(os, param...);
  }
};

static inline string suffix(const string& path)
{
  size_t sufix_p = path.find_last_of('.');
  if(sufix_p == std::string::npos) return "";
  return path.substr(sufix_p+1);
}

class polycube_generator{
public:

  typedef hj::function::function_t<double, int32_t> hj_func;
  typedef shared_ptr<hj_func> hj_func_ptr;
  typedef jtf::function::functionN1_t<double,int32_t> jtf_func;
  typedef shared_ptr<jtf_func> jtf_func_ptr;

  polycube_generator(ptree &pt):pt_(pt),diagnose_(nullptr){}

public:
  int generate_tet();
  int generate_tri();
  int output(const string & path) const;

protected:
  int regist_weight() const;
  int global_rotate_mesh(const matrix<size_t> & surface, const matrix<double> & areas_k);
  int load_input_model(std::string path);
  int set_init_node(std::string path);
  int translate_node();
  void draw_temp_model(string path);

  hj_func_ptr set_polycube_func_tet(const matrix<double> & original_node,
                                    const matrix<double> & deform_node,
                                    const matrix<double> & zero_pos,
                                    const matrix<size_t> & surface,
                                    const matrix<size_t> & tet_mesh,
                                    const double adj_normal_w) ;

  hj_func_ptr set_polycube_func_tri(const matrix<double> & original_node,
                                    const matrix<double> & deform_node,
                                    const matrix<double> & zero_pos,
                                    const matrix<size_t> & surface)const{}
private:
  ptree &pt_;
  jtf::mesh::meshes m_;
  matrix<double> original_node_;
  hj_func * diagnose_;
  MSG_LOG log_;
};

void polycube_generator::draw_temp_model(string path)
{ // visualize
  ofstream ofs(path.c_str());
  tet2vtk(ofs, &m_.node_[0], m_.node_.size(2), &m_.mesh_[0], m_.mesh_.size(2));
}

polycube_generator::hj_func_ptr polycube_generator::set_polycube_func_tet(
    const matrix<double> & original_node,
    const matrix<double> & deform_node,
    const matrix<double> & zero_pos,
    const matrix<size_t> & tet_mesh,
    const matrix<size_t> & surface,
    const double adj_normal_w)
{
  assert(surface.size(1) == 3 && tet_mesh.size(1) == 4);
  return hj_func_ptr(build_polycube_function(original_node, deform_node, zero_pos, tet_mesh, surface, diagnose_, nullptr, nullptr, adj_normal_w, 0,0));
}

int polycube_generator::translate_node()
{
  matrix<double> cur_ct = m_.node_*ones<double>(m_.node_.size(2), 1)/(1.0*m_.node_.size(2));
  matrix<double> ct = original_node_*ones<double>(original_node_.size(2),1)/(1.0*original_node_.size(2));
  m_.node_ += (ct-cur_ct)*ones<double>(1, m_.node_.size(2));
  return 0;
}

int polycube_generator::output(const string & path)const
{
  const string sufix = suffix(path);
  if(sufix == "obj" || sufix == "OBJ"){
      if(jtf::mesh::save_obj(pt_.get<string>("output.value").c_str(), m_.mesh_ ,m_.node_))
        return __LINE__;
    }else if(sufix == "tet" || sufix == "TET"){
      if(jtf::mesh::tet_mesh_write_to_zjumat(pt_.get<string>("output.value").c_str(), &m_.node_, &m_.mesh_))
        return __LINE__;
      if(jtf::mesh::tet_mesh_write_to_zjumat("after_collapse_output.tet", &original_node_, &m_.mesh_))
        return __LINE__;
    }

  return 0;
}
int polycube_generator::regist_weight()const
{
  pt_.put("L1_sqrt_eps.desc", "L1 sqrt eps, default is 1e-1.");
  pt_.put("normal_align_w.desc", "surface L1 weight, default is 1e2.");
  pt_.put("adj_normal_w.desc", "surface smoothness weight, default is 0.");
  pt_.put("anti_flip_w.desc", "surface anti-flip weight, default is 1e-2.");
  pt_.put("iter_w.desc", "iterative weight, default is 1");
  pt_.put("fix_zero_w.desc", "fix first node to zero, default is 1");
  pt_.put("div_L_w.desc", "div eps by L, default is sqrt2");
  pt_.put("area_preserve_w.desc", "polycube area preserving, default is 1.0");
  pt_.put("alpha_times.desc", "alpha times x, x default is 2");

  return 0;
}

int polycube_generator::load_input_model(string path)
{
  const string sufix = suffix(path);
  if(sufix == "obj" || sufix == "OBJ"){
      if(jtf::mesh::load_obj(path.c_str(), m_.mesh_, m_.node_))
        return __LINE__;
    }else if(sufix == "tet" || sufix == "TET"){
      if(jtf::mesh::tet_mesh_read_from_zjumat(path.c_str(), &m_.node_, &m_.mesh_))
        return __LINE__;
    }else{
      log_.out(std::cerr, "# [error] I don't know the type.");
      return __LINE__;
    }
  return 0;
}

int polycube_generator::global_rotate_mesh(const matrix<size_t> & surface,const matrix<double> & areas_k)
{
  static matrix<double> R = eye<double>(3);
  shared_ptr<jtf::function::functionN1_t<double,int32_t> > func(build_polycube_rot_func2(m_.node_, surface,areas_k));
  jtf::optimize(*func, R, pt_, nullptr, nullptr, nullptr);
  cout << R << endl;
  hj::polar3d p;
  p(R);
  cout << R << endl;
  m_.node_ = temp(R*m_.node_);
  return 0;
}

int polycube_generator::generate_tet()
{
  if(load_input_model(pt_.get<string>("tet.value"))) return __LINE__;
  if(zjucad::has("init.value",pt_)) {
      if(set_init_node(pt_.get<string>("init.value")))
        return __LINE__;
    }

  //////////////////////////////////////////////////////////
  const double alpha_times = pt_.get<double>("alpha_times.value",2);
  L1_sqrt_eps = pt_.get<double>("L1_sqrt_eps.value", 0.5);
  arap_w = pt_.get<double>("arap_w.value",1);
  double normal_align_w = pt_.get<double>("normal_align_w.value", 5e-2);
  double anti_flip_w = pt_.get<double>("anti_flip_w.value", 1e-2);
  adj_normal_w = pt_.get<double>("adj_normal_w.value", 0.0);

  const double div_L = pt_.get<double>("div_L.value",sqrt(2.0));
  const double epsg = pt_.get<double>("epsg.value", 1e-6);
  //////////////////////////////////////////////////////////

  jtf::tet_mesh tm(m_);
  original_node_ = m_.node_;
  const double total_area = std::accumulate(tm.outside_face_area_.begin(), tm.outside_face_area_.end(), 0.0);

  log_.out(std::cerr, "# [info] total_area: ", total_area);

  const size_t iter = pt_.get<size_t>("iter_w.value",1);

  for(size_t i = 0; i < iter; ++i){
      matrix<double> areas_k = zeros<double>(tm.outside_face_.size(2),1);
      for(size_t fi = 0; fi < tm.outside_face_.size(2); ++fi){
          areas_k[fi] = jtf::mesh::cal_face_area(tm.outside_face_(colon(),fi),m_.node_);
        }

      global_rotate_mesh(tm.outside_face_, areas_k);
      hj_func_ptr func = set_polycube_func_tet(original_node_, m_.node_, m_.node_(colon(),0), m_.mesh_, tm.outside_face_, adj_normal_w);

      log_.out(std::cerr, "# [info] add polycube surface normal L1 function, wegiht ", normal_align_w);
      log_.out(std::cerr, "# [info] L1_sqrt_eps = ", L1_sqrt_eps);

      vector<jtf_func_ptr> sum;
      vector<pair<jtf_func*, double> > wf;
      sum.push_back(jtf_func_ptr(jtf::function::least_square_warpper(func)));
      jtf_func_ptr arap_for_diagnose(jtf::function::least_square_warpper(*diagnose_));
      wf.push_back(make_pair(arap_for_diagnose.get(),1));
      if(normal_align_w > 0){
          sum.push_back(jtf_func_ptr(build_smooth_L1_area_normal(m_.node_, tm.outside_face_, areas_k, normal_align_w)));
          wf.push_back(make_pair(sum.back().get(), normal_align_w));
        }
      if(anti_flip_w > 0){
          matrix<double> weight;
          hj_func_ptr adj_normal_func(build_adj_normal_func(m_.node_, tm.outside_face_, weight));
          weight *= anti_flip_w;
          sum.push_back(jtf_func_ptr(jtf::function::neg_log_warpper(adj_normal_func, weight)));

          log_.out(std::cerr, "# [info] use anti-flip with weight: ", anti_flip_w);
        }
      jtf_func_ptr target(new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD>(sum));
      pair<jtf_func*,double> tmp(target.get(), 1+normal_align_w);
      wf.push_back(tmp);
      jtf_func_ptr constraint(new area_sum(m_.node_.size(2), tm.outside_face_, total_area));
      unnormalized_normal_quality_checker cb(m_.node_, m_.mesh_, tm.outside_face_, wf, *constraint);

      pt_.put("epsg.value", epsg*(1+normal_align_w));

      int rtn = -1;
      if(normal_align_w > 0){
          vector<jtf_func_ptr> constraint_vec;
          constraint_vec.push_back(constraint);
          rtn = jtf::optimize(*target, m_.node_, pt_, &constraint_vec, nullptr, &cb);
        }
      else
        rtn = jtf::optimize(*target, m_.node_, pt_, nullptr, nullptr, nullptr);

      // remove_degenerated_face_of_tet(m_.mesh_, original_node, m_.node_, tm.outside_face_, 5);

      translate_node();
      if(polycube_L1_area_quality(&m_.node_[0], m_.node_.size(2), tm.outside_face_) < 1e-3){
          log_.out(std::cout, "# [info] global converge.");
          break;
        }

      ostringstream vtk_path_pref;
      vtk_path_pref << adj_normal_w << "-" << i << ".vtk";
      draw_temp_model(vtk_path_pref.str());

      normal_align_w *= alpha_times;
      L1_sqrt_eps /= div_L;
      if(L1_sqrt_eps < 1e-2)
        L1_sqrt_eps = 1e-2;
    }

  log_.out(std::cerr, "# [info] finial polycube_quality", polycube_L1_area_quality(&m_.node_[0], m_.node_.size(2), tm.outside_face_));

  return 0;
}

int polycube_generator::generate_tri()
{
  if(load_input_model(pt_.get<string>("obj.value"))) return __LINE__;
  if(zjucad::has("init.value",pt_)) {
      if(set_init_node(pt_.get<string>("init.value")))
        return __LINE__;
    }
  original_node_ = m_.node_;
  matrix<double> area = zeros<double>(m_.mesh_.size(2),1);
  for(size_t fi = 0; fi < m_.mesh_.size(2); ++fi){
      area[fi] = jtf::mesh::cal_face_area(m_.mesh_(colon(),fi),m_.node_);
    }

  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  std::cerr << "# [info] total_area: " << total_area << std::endl;

  const size_t iter = pt_.get<size_t>("iter_w.value",1);
  for(size_t i = 0; i < iter; ++i){
      global_rotate_mesh(m_.mesh_, area);
    }
  return 0;
}


int polycube_generator::set_init_node(string path)
{
  size_t sufix_p = path.find_last_of('.');
  if(sufix_p == std::string::npos) return __LINE__;
  string sufix = path.substr(sufix_p);
  jtf::mesh::meshes init_mesh;
  if(sufix == "obj" || sufix == "OBJ"){
      if(jtf::mesh::load_obj(path.c_str(), init_mesh.mesh_, init_mesh.node_))
        return __LINE__;
    }else if(sufix == "tet" || sufix == "TET"){
      if(jtf::mesh::tet_mesh_read_from_zjumat(path.c_str(), &init_mesh.node_, &init_mesh.mesh_))
        return __LINE__;
    }else{
      std::cerr << "# [error] I don't know the type." << std::endl;
      return __LINE__;
    }

  if(max(abs(init_mesh.mesh_ - m_.mesh_)) != 0){
      std::cerr << "# [error] incompatible init value." << std::endl;
      return __LINE__;
    }
  m_.node_ = init_mesh.node_;

  return 0;
}

int polycube_generate(ptree &pt)
{
  polycube_generator pg(pt);

  if(zjucad::has("tet.value",pt)){
      pg.generate_tet();
    }else if(zjucad::has("obj.value",pt)){
      pg.generate_tri();
    }
  pg.output(pt.get<string>("output.value"));
  return 0;
}
