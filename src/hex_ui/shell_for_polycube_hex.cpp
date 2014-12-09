#include <zjucad/ptree/ptree.h>
#include <jtflib/optimizer/optimizer.h>
#include <hjlib/math/polar.h>
#include <hjlib/math_func/math_func.h>
#include <hjlib/math_func/operation.h>
#include "../tetmesh/tetmesh.h"
#include "../hex_frame_opt/optimal_sh_generator.h"
#include "../hexmesh/hexmesh.h"
#include "../hexmesh/util.h"
#include "../tetmesh/util.h"
#include "../common/vtk.h"
#include "../mesh_func/frame_func.h"

using namespace std;
using namespace zjucad::matrix;
using boost::property_tree::ptree;

class optimal_sh_generator_with_cons:  public optimal_sh_generator
{
public:
  optimal_sh_generator_with_cons(jtf::tet_mesh & tm,
                                 const matrix<matrix<double> > &full_face_cons)
    :optimal_sh_generator(tm), tm_(tm), full_face_cons_(full_face_cons){}
  virtual ~optimal_sh_generator_with_cons(){}
  virtual void opt(zjucad::matrix::matrix<double> &field,
                   boost::property_tree::ptree &pt){
    const int opt_num = pt.get<int>("input/opt_num.value",7);
    switch (opt_num) {
      case 1: opt1(field,pt); break;
      default: cerr << "# [error] wrong opt_num." << endl; break;
      }
  }
private:
  void opt1(zjucad::matrix::matrix<double> & field,
            boost::property_tree::ptree &pt);

  double calculate_dihedral_angle_degree(const matrix<double> &d1,
                                         const matrix<double> &d2)const;

  void init_feature_lines(vector<vector<size_t> > & feature_lines) const;
private:
  const jtf::tet_mesh &tm_;
  const matrix<matrix<double> > & full_face_cons_;
};

double optimal_sh_generator_with_cons::calculate_dihedral_angle_degree(
    const matrix<double> & n1, const matrix<double> & n2)const
{
  const double len1 = norm(n1);
  const double len2 = norm(n2);
  matrix<double> n1_norm = n1, n2_norm = n2;
  if(len1 > 1e-6) n1_norm /= len1;
  if(len2 > 1e-6) n2_norm /= len2;

  return jtf::math::cos2deg(dot(n1_norm, n2_norm));
}

void optimal_sh_generator_with_cons::init_feature_lines(vector<vector<size_t> > & feature_lines) const
{
  vector<pair<size_t,size_t> > edge_store;
  const double threshold = 45.;
  for(size_t ei = 0; ei < tm_.ea_outside_->edges_.size(); ++ei){
      const pair<size_t,size_t> & one_edge = tm_.ea_outside_->edges_[ei];
      const pair<size_t,size_t> & edge2cell = tm_.ea_outside_->edge2cell_[ei];
      if(tm_.ea_outside_->is_boundary_edge(edge2cell)) edge_store.push_back(one_edge);
      if(full_face_cons_[edge2cell.first].size(2) == 1 && full_face_cons_[edge2cell.second].size(2) == 1){
          const double dihedral_angle_degree =
              calculate_dihedral_angle_degree(tm_.outside_face_normal_(colon(), edge2cell.first),
                                              tm_.outside_face_normal_(colon(), edge2cell.second));
          if(dihedral_angle_degree > threshold){
              edge_store.push_back(one_edge);;
            }
        }
    }

  vector<deque<pair<size_t,size_t> > > feature_line_deque;
  jtf::util::extract_chain_from_edges(edge_store, feature_line_deque);
  feature_lines.clear();
  feature_lines.resize(feature_line_deque.size());
  for(size_t fi = 0; fi < feature_line_deque.size(); ++fi){
      feature_lines[fi].push_back(feature_line_deque[fi].front().first);
      for(size_t ei = 0; ei < feature_line_deque[fi].size(); ++ei){
          feature_lines[fi].push_back(feature_line_deque[fi][ei].second);
        }
    }
}

void optimal_sh_generator_with_cons::opt1(zjucad::matrix::matrix<double> &field,
                                          boost::property_tree::ptree &pt)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  init(pt);

  pt.put("input/smooth_scheme.desc", "input smooth scheme strategy [face/gradient]");
  pt.put("input/opt_type.desc", "init or zyz");
  pt.put("weight/feature.desc", "weight of feature align.");

  const double w_normal_align = pt.get<double>("weight/normal_align.value");
  const double w_feature = pt.get<double>("weight/feature.value");

  const string smooth_scheme = pt.get<string>("input/smooth_scheme.value");
  const string opt_type = pt.get<string>("input/opt_type.value");
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

  if(w_normal_align > 1e-8){
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_full_surface_cons_align(
                            tm_.tetmesh_.mesh_, *tm_.fa_, tm_.outside_face_idx_,
                            full_face_cons_, tm_.outside_face_area_,
                            w_normal_align, opt_type))));
      cerr << "# [info] add normal align function: " << w_normal_align << endl;
    }

  if(w_feature > 1e-8){
      vector<vector<size_t> > feature_lines;
      init_feature_lines(feature_lines);
      func->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_feature_line_align_func<double>(
                            tm_.tetmesh_.mesh_, tm_.tetmesh_.node_,
                            feature_lines,tm_.ortae_,w_feature, opt_type, true))));
      cerr << "# [info] add feature align function: " << w_feature << endl;
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

inline void generate_one_hex_frame(const matrix<size_t> & one_hex,
                                   const matrix<double> & hex_node,
                                   matrix<double> & frame)
{
  frame.resize(3,3);
  frame *= 0;

  size_t u_edge[] = {0,1,2,3,4,5,6,7};
  size_t v_edge[] = {0,2,1,3,4,6,5,7};
  size_t w_edge[] = {0,4,1,5,2,6,3,7};

  for(size_t pi = 0; pi < 4; ++pi){
      frame(colon(),0) += hex_node(colon(), one_hex[u_edge[2*pi+1]]) - hex_node(colon(), one_hex[u_edge[2*pi]]);
      frame(colon(),1) += hex_node(colon(), one_hex[v_edge[2*pi+1]]) - hex_node(colon(), one_hex[v_edge[2*pi]]);
      frame(colon(),2) += hex_node(colon(), one_hex[w_edge[2*pi+1]]) - hex_node(colon(), one_hex[w_edge[2*pi]]);
    }

  for(size_t di = 0; di < 3; ++di){
      frame(colon(), di) /= norm(frame(colon(),di));
    }

  hj::polar3d p;
  p(frame,2);
}

inline void generate_hex_frame_field(const jtf::hex_mesh &hm, matrix<matrix<double> > & frame)
{
  if(frame.size() != hm.hexmesh_.mesh_.size(2))
    frame.resize(hm.hexmesh_.mesh_.size(2),1);
  for(size_t hi = 0; hi < hm.hexmesh_.mesh_.size(2); ++hi){
      generate_one_hex_frame(hm.hexmesh_.mesh_(colon(), hi), hm.hexmesh_.node_, frame[hi]);
    }
}


inline void convert_quad_to_tri(const matrix<size_t> & quad,
                                matrix<size_t> & tri)
{
  tri.resize(3, quad.size(2) * 2);
  for(size_t qi = 0; qi < quad.size(2); ++qi){
      for(size_t pi = 0; pi < 3; ++pi){
          tri(pi, 2*qi+0) = quad(pi,qi);
          tri(pi, 2*qi+1) = quad((pi+2)%4,qi);
        }
    }
}

void generate_one_layer_tet(const jtf::hex_mesh &hm, jtf::mesh::meshes & tm,
                            map<size_t,size_t> & part_tri_face2hex_idx)
{
  matrix<size_t> tri_face;
  matrix<double> point_normal;
  convert_quad_to_tri(hm.outside_face_, tri_face);

  jtf::mesh::cal_point_normal(tri_face, hm.hexmesh_.node_, point_normal);

  jtf::tetmesh::extend_tetmesh(tri_face, hm.hexmesh_.node_, point_normal, tm);

  {
    ofstream ofs("tet_layer.vtk");
    tet2vtk(ofs, &tm.node_[0], tm.node_.size(2), &tm.mesh_[0], tm.mesh_.size(2));
  }

  shared_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  assert(fa.get());

  for(size_t fi = 0; fi < tri_face.size(2); ++fi){
      const size_t face_idx = fa->get_face_idx(&tri_face(0,fi));
      assert(face_idx != -1);
      const size_t qi = fi/2;
      const size_t h_face_idx = hm.fa_->get_face_idx(&hm.outside_face_(0,qi));
      assert(h_face_idx != -1);
      const pair<size_t,size_t> & hex_pair = hm.fa_->face2hex_[h_face_idx];
      const size_t other_hex_idx = (hex_pair.first==-1?hex_pair.second:hex_pair.first);
      part_tri_face2hex_idx[face_idx] = other_hex_idx;
    }
}

int shell_for_polycube_hex(ptree &pt)
{
  jtf::hex_mesh hm(pt.get<string>("input/hex.value").c_str());

  matrix<matrix<double> > frames;
  generate_hex_frame_field(hm, frames);

  map<size_t,size_t> tri_face2hex;
  jtf::mesh::meshes tm;
  generate_one_layer_tet(hm, tm, tri_face2hex);

  jtf::tet_mesh tm_(tm);
  matrix<matrix<double> > full_face_cons(tm_.outside_face_.size(2),1);

  for(size_t fi = 0; fi < tm_.outside_face_idx_.size(); ++fi){
      auto it = tri_face2hex.find(tm_.outside_face_idx_[fi]);
      if(it == tri_face2hex.end()) full_face_cons[fi] = tm_.outside_face_normal_(colon(),fi);
      else full_face_cons[fi] = frames[it->second];
    }

  matrix<double> zyz = zeros<double>(3,tm.mesh_.size(2));
  zjucad::matrix::matrix<double> sh(9, zyz.size(2));
  for(size_t ti = 0; ti < tm.mesh_.size(2); ++ti){
      sh(4,colon()) += sqrt(7.0);
      sh(8,colon()) += sqrt(5.0);
    }
  optimal_sh_generator_with_cons osg(tm_, full_face_cons);
  osg.opt(sh, pt);
  for(size_t i = 0; i < sh.size(2); ++i){
      sh(colon(),i) /= norm(sh(colon(),i));
    }
  sh *= sqrt(12.0);

  osg.write_zyz(pt.get<string>("output.value").c_str(),sh, true);
  return 0;
}
