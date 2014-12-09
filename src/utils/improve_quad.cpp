#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <hjlib/math_func/math_func.h>
#include "../common/K_near_points.h"
#include <zjucad/ptree/ptree.h>
#include <jtflib/optimizer/optimizer.h>
#include "../mesh_func/point_fix.h"
#include "../mesh_func/surface_smooth.h"
#include "../mesh_func/feature_line.h"
using namespace std;
using namespace zjucad::matrix;
using boost::property_tree::ptree;

class surf_mesh_improver{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;
public:
  surf_mesh_improver(jtf::mesh::meshes &goal, const jtf::mesh::meshes &reference)
    :goal_(goal), reference_(reference){}

  void run(const size_t iter = 5){
    init();
    for(size_t i = 0; i < iter; ++i){
        update_point_targets();
        build_equations();
        optimize();
      }
  }

  void run(const vector<deque<pair<size_t,size_t> > > & quad_fl,
           const vector<deque<pair<size_t,size_t> > > & tri_fl,
           const size_t iter = 5){
    init();
    for(size_t i = 0; i < iter; ++i){
        update_point_targets();
        build_equations(quad_fl, tri_fl);
        optimize();
      }
  }
  void output(jtf::mesh::meshes & new_goal) const{
    new_goal = goal_;
  }

  void run_iter(const vector<deque<pair<size_t,size_t> > > & quad_fl,
                const vector<deque<pair<size_t,size_t> > > & tri_fl,
                const size_t iter = 10){
    init();
    set<size_t> visited_set;
    for(size_t di = 0; di < quad_fl.size(); ++di){
        const deque<pair<size_t,size_t> > &one_fl = quad_fl[di];
        for(const auto & one_edge : one_fl){
            visited_set.insert(one_edge.first);
            visited_set.insert(one_edge.second);
          }
      }

    vector<bool> visited_node_flag(goal_.node_.size(2), false);
    for(const auto & one_idx : visited_set) visited_node_flag[one_idx] = true;

    for(size_t it = 0; it < iter ; ++it){
        update_points_along_feature_line(quad_fl);
        for(size_t j = 0; j < 1; ++j)
          update_other_points(visited_node_flag);
            update_points_along_feature_line(quad_fl);
        //            update_point_targets();
        //            for(size_t pi = 0;  pi < visited_node_flag.size(); ++pi){
        //                if(visited_node_flag[pi] == false)
        //                    goal_.node_(colon(), pi) = point_target_(colon(), pi);
        //            }
      }
  }

  void run_iter(const size_t iter = 5){
    init();

    vector<bool> visited_node_flag(goal_.node_.size(2), false);

    for(size_t it = 0; it < iter ; ++it){
        update_other_points(visited_node_flag);
        update_point_targets();
        goal_.node_ = point_target_;
      }
  }
protected:
  void init(){
    knp.reset(new K_near_points(reference_.node_));
    shared_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(reference_.mesh_));
    if(!ea.get()) throw std::logic_error("can not build edge2cell_adjacent.");
    orfap_.add_all_faces(reference_.mesh_, *ea);
    funcs_.reset(new vector<math_func_ptr>);
    shared_ptr<jtf::mesh::edge2cell_adjacent> goal_ea(jtf::mesh::edge2cell_adjacent::create(goal_.mesh_));
    if(!goal_ea.get()) throw std::logic_error("can not build edge2cell_adjacent.");
    orfap_goal_.add_all_faces(goal_.mesh_, *goal_ea);
  }

  void update_point_targets(){
    if(point_target_.size(2) != goal_.node_.size(2))
      point_target_.resize(3, goal_.node_.size(2));
    vector<int> kps;
    vector<double> dis;
    for(size_t i = 0; i < goal_.node_.size(2); ++i){
        if(i == 16490)
          cerr << endl;
        knp->query_k_near_points(goal_.node_(colon(), i), 1, kps, dis,0.0);
        assert(kps.size());
        auto it = orfap_.p2f_.find(kps.front());
        assert(it != orfap_.p2f_.end());
        point_target_(colon(),i) = get_nearest_geo_point(goal_.node_(colon(), i), kps.front(), it->second);
      }
  }

  void update_points_along_feature_line(const vector<deque<pair<size_t,size_t> > >  & fl)
  {
    const double alpha = 0.8;
    for(size_t di = 0; di < fl.size(); ++di){
        const deque<pair<size_t,size_t> > & one_line = fl[di];
        for(size_t i = 0; i + 1 < one_line.size(); ++i){
            goal_.node_(colon(), one_line[i].second) =
                ((1-alpha)*0.5*(goal_.node_(colon(), one_line[i].first)+
                                goal_.node_(colon(), one_line[i+1].second))
                + alpha*goal_.node_(colon(), one_line[i].second));
          }
      }
  }

  void update_other_points(const vector<bool> & visited_node_flag){
    matrix<double> center(3,1);
    const double alpha = 0.7;
    set<size_t> one_ring_points;
    for(const auto & p : orfap_goal_.p2f_){
        if(visited_node_flag[p.first]) continue;
        center *= 0;
        one_ring_points.clear();
        const vector<size_t> & one_ring_face = p.second;
        for(size_t fi = 0; fi < one_ring_face.size(); ++fi){
            if(one_ring_face[fi] != -1)
              one_ring_points.insert(
                    goal_.mesh_(colon(), one_ring_face[fi]).begin(), goal_.mesh_(colon(), one_ring_face[fi]).end());
          }
        for(const size_t &idx : one_ring_points){
            if(idx != p.first)
              center += goal_.node_(colon(), idx);
          }
        center /= one_ring_points.size()-1;
        goal_.node_(colon(), p.first) =
            (alpha*goal_.node_(colon(), p.first) + (1.0-alpha)*center);
      }
  }

  matrix<double> get_nearest_geo_point(const matrix<double> &p,
                                       const size_t center_idx,
                                       const vector<size_t> & faces)const{
    matrix<double> pt(3,1);
    for(size_t fi = 0; fi < faces.size(); ++fi){
        get_project_point(pt, p, reference_.mesh_(colon(), faces[fi]), reference_.node_);
        if(is_inside_triangle(reference_.mesh_(colon(), faces[fi]), reference_.node_, pt))
          return pt;
      }
    return reference_.node_(colon(), center_idx);
  }
  void build_equations(){
    assert(funcs_.get());
    double weight_point_fix = 1;
    double weight_smooth = 0.1;
    funcs_->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_surface_smooth_math_func(goal_.mesh_, goal_.node_, weight_smooth))));
    funcs_->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_point_fix_math_func(point_target_, weight_point_fix ))));
    //funcs_->push_back(build_orthogonal_func());
  }
  void build_equations(const vector<deque<pair<size_t,size_t> > > & quad_fl,
                       const vector<deque<pair<size_t,size_t> > > & tri_fl){
    assert(funcs_.get());
    double weight_point_fix = 1;
    double weight_smooth = 1;
    double weight_feature_line = 1;
    funcs_->push_back(math_func_ptr(
                        new hj::math_func::sumsqr<double,int32_t>(
                          build_surface_smooth_math_func(goal_.mesh_, goal_.node_, weight_smooth))));

    {
      set<size_t> feature_ends;
      for(size_t fi = 0; fi < quad_fl.size(); ++fi){
          feature_ends.insert(quad_fl[fi].front().first);
          feature_ends.insert(quad_fl[fi].back().second);
        }
      matrix<size_t> selectd_points(feature_ends.size(),1);
      std::copy(feature_ends.begin(), feature_ends.end(), selectd_points.begin());
      assert(selectd_points[selectd_points.size()-1] < point_target_.size(2));
      funcs_->push_back(math_func_ptr(
                          new hj::math_func::sumsqr<double,int32_t>(
                            build_point_fix_math_func(point_target_,selectd_points, point_target_.size(2), weight_point_fix ))));
    }
    //    funcs_->push_back(math_func_ptr(
    //                        new hj::math_func::sumsqr<double,int32_t>(
    //                          build_feature_line_math_func(quad_fl, point_target_, weight_feature_line))));
  }
  void optimize(){
    math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(funcs_));
    math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));
    ptree pt;
    pt.put("package.value", "jtf");
    pt.put("alg.value", "SQP");
    pt.put("iter.value",10);
    jtf::optimize(*obj, goal_.node_, pt, nullptr, nullptr, nullptr);
  }
  void get_project_point(zjucad::matrix::matrix<double> &pt,
                         const zjucad::matrix::matrix<double> & source,
                         const zjucad::matrix::matrix<size_t> & face,
                         const zjucad::matrix::matrix<double> & node)const{
    matrix<double> normal(3,1);
    jtf::mesh::cal_face_normal(face, node, normal);
    matrix<double> dir0 = node(colon(), face[0]) - source;
    if(dot(dir0, normal) < 0) normal *= -1;
    matrix<double> tet_node(3,4);
    tet_node(colon(), colon(0,2)) = node(colon(), face);
    tet_node(colon(),3) = source;
    const double height = 3*fabs(jtf::mesh::cal_tet_vol(tet_node))/fabs(jtf::mesh::cal_face_area(node(colon(), face)));
    pt = source + normal * height;
  }
  bool is_inside_triangle(const zjucad::matrix::matrix<size_t> & tri,
                          const zjucad::matrix::matrix<double> & node,
                          const zjucad::matrix::matrix<double> & pt)const{
    double edge_len[] = {
      norm(pt - node(colon(), tri[0])),
      norm(pt - node(colon(), tri[1])),
      norm(pt - node(colon(), tri[2])),
    };
    double edge_len_orig[] = {
      norm(node(colon(), tri[1]) - node(colon(), tri[2])),
      norm(node(colon(), tri[0]) - node(colon(), tri[2])),
      norm(node(colon(), tri[0]) - node(colon(), tri[1]))
    };
    double area_orig = jtf::mesh::cal_face_area(edge_len_orig[0], edge_len_orig[1], edge_len_orig[2]);
    double area_out[] = {jtf::mesh::cal_face_area(edge_len[0], edge_len[1], edge_len_orig[2]),
                         jtf::mesh::cal_face_area(edge_len[1], edge_len[2], edge_len_orig[0]),
                         jtf::mesh::cal_face_area(edge_len[2], edge_len[0], edge_len_orig[1])};
    const double err = fabs(std::accumulate(area_out, area_out+3, 0.0) - area_orig);
    if( err > 1e-5)
      return false;
    return true;
  }
private:
  shared_ptr<K_near_points> knp;
  shared_ptr<vector<math_func_ptr> > funcs_;
  jtf::mesh::one_ring_face_at_point orfap_;
  jtf::mesh::one_ring_face_at_point orfap_goal_;
  jtf::mesh::meshes &goal_;
  zjucad::matrix::matrix<double> point_target_;
  const jtf::mesh::meshes &reference_;

};

int improve_quad(int argc, char *argv[])
{
  if(argc != 3 && argc != 5){
      cerr << "# [usage] improve_quad quad tri [quad_fl] [tri_fl]" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes quad, tri;
  if(jtf::mesh::load_obj(argv[1], quad.mesh_, quad.node_) ||
     jtf::mesh::load_obj(argv[2], tri.mesh_, tri.node_)){
      cerr << "# [error] can not load quad or tri." << endl;
      return __LINE__;
    }

  surf_mesh_improver mi(quad, tri);

  if(argc == 3){
      mi.run_iter();
    }else if(argc == 5){
      vector<deque<pair<size_t,size_t> > > tri_fl, quad_fl;
      if(jtf::mesh::load_feature_line(argv[3], quad_fl) ||
         jtf::mesh::load_feature_line(argv[4], tri_fl)){
          cerr << "# [error] can not load feature line" << endl;
          return __LINE__;
        }
      mi.run_iter(quad_fl, tri_fl);
    }

  mi.output(quad);
  jtf::mesh::save_obj("opt.obj", quad.mesh_, quad.node_);
  return 0;
}
