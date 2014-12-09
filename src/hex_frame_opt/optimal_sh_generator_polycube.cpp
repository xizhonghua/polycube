#include "optimal_sh_generator.h"
#include "../numeric/util.h"
#include "../common/vtk.h"
#include "../common/visualize_tool.h"
#include "../hex_param/io.h"
#include "../common/zyz.h"
#include <stack>
#include <jtflib/mesh/util.h>
#include <jtflib/util/util.h>

using namespace std;
using namespace zjucad::matrix;

void optimal_sh_generator_polycube::opt(zjucad::matrix::matrix<double> &field,
                                        boost::property_tree::ptree &pt)
{
  const size_t opt_num = pt.get<size_t>("input/polycube_opt_num.value");
  switch (opt_num)
    {
    case 0: opt0(field,pt);break;
    case 1: opt1(field,pt);break;
    case 2: opt2(field,pt);break;
    default: throw std::invalid_argument("unknown opt num.");
    }
}

void optimal_sh_generator_polycube::update_normal(
    const set<size_t> & patch_to_use_real_noraml,
    const vector<vector<size_t> > & all_patches)
{
  jtf::mesh::cal_face_normal(orig_tm_.outside_face_, orig_tm_.tetmesh_.node_, orig_tm_.outside_face_normal_);

  matrix<double> avg_normal(3,1);
  for(size_t pi = 0; pi < all_patches.size(); ++pi){
      if(patch_to_use_real_noraml.find(pi) == patch_to_use_real_noraml.end()){
          cal_weighted_average_normal(all_patches[pi], orig_tm_.outside_face_normal_, orig_tm_.outside_face_area_, avg_normal);
          for(size_t fi = 0; fi < all_patches[pi].size(); ++fi){
              orig_tm_.outside_face_normal_(colon(), all_patches[pi][fi]) = avg_normal;
            }
        }
    }
}

void optimal_sh_generator_polycube::cal_weighted_average_normal(
    const std::vector<size_t> &one_patch,
    const zjucad::matrix::matrix<double> &normal,
    const zjucad::matrix::matrix<double> &area,
    zjucad::matrix::matrix<double> &avg_normal)
{
  avg_normal.resize(3,1);
  avg_normal *= 0;
  for(const auto & one_face : one_patch){
      avg_normal += fabs(area[one_face])*normal(colon(),one_face);
    }
  avg_normal /= norm(avg_normal);
}

bool optimal_sh_generator_polycube::valid_field(const zjucad::matrix::matrix<double> &field)
{
  static vector<deque<pair<size_t,size_t> > > chain_list;
  static vector<deque<size_t> > singularities_type;
  static set<size_t> surface_points(orig_tm_.outside_face_.begin(), orig_tm_.outside_face_.end());
  static matrix<double> point_normal;
  if(point_normal.size() == 0)
    jtf::mesh::cal_point_normal(orig_tm_.outside_face_, orig_tm_.tetmesh_.node_, point_normal);

  singularity_extractor se(orig_tm_);
  se.extract(field, chain_list,singularities_type);

  {
    static int test = 0;
    stringstream ss;
    ss << "singularity_" << test++ << ".vtk";
    dump_singularity_to_vtk(ss.str().c_str(), orig_tm_.tetmesh_.node_, chain_list);
  }
  for(size_t ci = 0; ci < chain_list.size(); ++ci){
      if(rule0(chain_list[ci], surface_points, point_normal) == false) {
          cerr << "# [info] against rule0" << endl;
          return false;
        }
    }
  cerr << "# [info] pass all rules." << endl;
  return true;
}

bool optimal_sh_generator_polycube::rule0(const deque<pair<size_t,size_t> > &one_chain,
                                          const set<size_t> &surface_points,
                                          const zjucad::matrix::matrix<double> &point_normal)
{// to judge a chain is othognal with surface normal
  if(one_chain.size() < 5) return true; // ignore small noisy singularity.

  set<size_t> points_touch_surface;
  if(surface_points.find(one_chain.front().first) != surface_points.end()){
      points_touch_surface.insert(one_chain.front().first);
    }
  if(surface_points.find(one_chain.back().second) != surface_points.end()){
      points_touch_surface.insert(one_chain.back().second);
    }

  matrix<double> avg_dir = zeros<double>(3,1);
  for(size_t ei = 0; ei < one_chain.size(); ++ei){
      const pair<size_t,size_t> & one_edge = one_chain[ei];
      avg_dir += orig_tm_.tetmesh_.node_(colon(), one_edge.second)
          - orig_tm_.tetmesh_.node_(colon(), one_edge.first);
    }
  avg_dir /= norm(avg_dir);

  const double angle_cos = cos(45.0*My_PI()/180.);
  for(const auto & one_point : points_touch_surface){
      if(fabs(dot(point_normal(colon(), one_point), avg_dir)) < angle_cos){
          cerr << "# [warning] chain <" << one_chain.front().first << " --> "
               << one_chain.back().second << "> with surface normal angle at "
               << one_point << " is " <<  jtf::math::safe_acos(dot(point_normal(colon(), one_point), avg_dir))*180./My_PI() << endl;
          return false;
        }
    }
  return true;
}

void draw_patches(
    const set<size_t> &patch_to_use_real_noraml,
    const vector<vector<size_t> > &all_patches,
    const zjucad::matrix::matrix<double> &node,
    const zjucad::matrix::matrix<size_t> & faces,
    const char * filename)
{
  set<size_t> face_selected;
  for(const auto & one_p : patch_to_use_real_noraml){
      face_selected.insert(all_patches[one_p].begin(), all_patches[one_p].end());
    }
  matrix<size_t> face_m(face_selected.size(),1);
  copy(face_selected.begin(), face_selected.end(), face_m.begin());

  matrix<size_t> selected_face = faces(colon(), face_m);
  ofstream ofs(filename);
  tri2vtk(ofs, &node[0], node.size(2), &selected_face[0], selected_face.size(2));
}

void optimal_sh_generator_polycube::opt0(zjucad::matrix::matrix<double> & field,
                                         boost::property_tree::ptree &pt)
{
  init(pt);
  jtf::mesh::patch_separater ps(orig_tm_.outside_face_, orig_tm_.ea_outside_);
  ps.separater(chain_);

  {
    ofstream ofs("patch.vtk");
    tri2vtk(ofs, &orig_tm_.tetmesh_.node_[0], orig_tm_.tetmesh_.node_.size(2),
        &orig_tm_.outside_face_[0], orig_tm_.outside_face_.size(2));
    cell_data(ofs, &ps.get_all_face2patch()[0], ps.get_all_face2patch().size(), "patch_idx");
  }
  set<size_t> patch_to_use_real_noraml;
  matrix<double> zyz;
  bool need_to_recal = false;
  for(size_t ci = 0; ci < chain_.size(); ++ci){
      cerr << "# [info] handling chain " << ci << "/" << chain_.size() << endl;
      if(is_sharp_chain_[ci]) continue;
      set<size_t> adj_patches = ps.get_chain_adj_patches(ci);
      set<size_t> patch_to_use_real_noraml_bkp = patch_to_use_real_noraml;
      patch_to_use_real_noraml.insert(adj_patches.begin(), adj_patches.end());

      update_normal(patch_to_use_real_noraml, ps.get_all_patches());
      optimal_sh_generator osg(orig_tm_);
      osg.opt(field, pt);
      write_zyz(field, zyz, true);
      if(valid_field(zyz) == false){
          patch_to_use_real_noraml = patch_to_use_real_noraml_bkp;
          need_to_recal = true;
        }
      else{
          stringstream ss;
          ss << "init_" << ci << ".zyz";
          osg.write_zyz(ss.str().c_str(),zyz, true);

          stringstream ss_patch;
          ss_patch << "patch_" << ci << ".vtk";
          draw_patches(patch_to_use_real_noraml, ps.get_all_patches(),
                       orig_tm_.tetmesh_.node_, orig_tm_.outside_face_ ,ss_patch.str().c_str());
        }
    }
  if(need_to_recal){
      update_normal(patch_to_use_real_noraml, ps.get_all_patches());
      optimal_sh_generator osg(orig_tm_);
      osg.opt(field, pt);
    }
}

int optimal_sh_generator_polycube::load_picked_patches(
    vector<size_t> & picked_patch_idx, boost::property_tree::ptree &pt)
{
  ifstream ifs(pt.get<string>("input/picked_patch.value").c_str());
  if(ifs.fail()){
      std::cerr << "can not open picked patch file." << endl;
      return __LINE__;
    }
  picked_patch_idx.clear();
  size_t num;
  ifs >> num;
  size_t temp;
  for(size_t i = 0; i < num; ++i){
      ifs >> temp;
      picked_patch_idx.push_back(temp);
    }
  cerr << "# [info] load picked_patches." << endl;
  return 0;
}
void optimal_sh_generator_polycube::opt1(zjucad::matrix::matrix<double> & field,
                                         boost::property_tree::ptree &pt)
{
  init(pt);
  jtf::mesh::patch_separater ps(orig_tm_.outside_face_, orig_tm_.ea_outside_);
  ps.separater(chain_);

  {
    ofstream ofs("patch.vtk");
    tri2vtk(ofs, &orig_tm_.tetmesh_.node_[0], orig_tm_.tetmesh_.node_.size(2),
        &orig_tm_.outside_face_[0], orig_tm_.outside_face_.size(2));
    cell_data(ofs, &ps.get_all_face2patch()[0], ps.get_all_face2patch().size(), "patch_idx");
  }

  matrix<matrix<double> > frame(orig_tm_.tetmesh_.mesh_.size(2),1);
  matrix<double> point_normal;
  if(point_normal.size() == 0)
    jtf::mesh::cal_point_normal(orig_tm_.outside_face_, orig_tm_.tetmesh_.node_, point_normal);

  int is_frame_ok = 0;
  size_t round = 0;
  jtf::tet_mesh try_tm = orig_tm_;
  while(is_frame_ok != 1){
      optimal_sh_generator osg(try_tm);
      osg.opt(field, pt);
      matrix<double> zyz;
      write_zyz(field, zyz, true);
      osg.zyz2frame(zyz, frame);
      {
        stringstream ss;
        ss << "init_" << round++ << ".zyz";
        osg.write_zyz(ss.str().c_str(),zyz, true);
      }

      is_frame_ok = adjust_surface_normal(try_tm.outside_face_normal_, frame, ps, polycube_tm_.outside_face_normal_);
      if(is_frame_ok == 1){
          cerr << "# [info] singularity is all ok." << endl;
        }else if(is_frame_ok == 0){
          cerr << "# [info] adjusted surface normal." << endl;
        }else{
          cerr << "# [info] all surface patches are adjusted, maybe not feasible." << endl;
          break;
        }
    }
  {
    jtf::mesh::write_matrix("test_normal", try_tm.outside_face_normal_);
  }
}

size_t optimal_sh_generator_polycube::get_point_convex_info(
    const size_t pi, const jtf::mesh::one_ring_face_at_point &orfap,
    const jtf::mesh::patch_separater & ps, const zjucad::matrix::matrix<double> & normal)
{
  cerr << "# [warning] you should ensure input one_ring_face_at_point is oriented in right hand order." << endl;
  auto it = orfap.p2f_.find(pi);
  if(it == orfap.p2f_.end()) {
      cerr << "# [error] can not find point in orfap." << endl;
      return -1;
    }
  const vector<size_t> & around_faces = it->second;
  assert(around_faces.front() == around_faces.back());

  vector<size_t> patch_idx;
  for(size_t i = 0 ; i < around_faces.size(); ++i){
      const size_t patch_i = ps.get_face2patch(around_faces[i]);
      if(patch_idx.empty()) patch_idx.push_back(patch_i);
      else{
          if(patch_idx.back() != patch_i)  patch_idx.push_back(patch_i);
        }
    }
  if(patch_idx.back() == patch_idx.front()) patch_idx.pop_back();

  if(patch_idx.size() < 3){
      cerr << "# [error] " << pi << "is not a patch corner." << endl;
      return -1;
    }

  matrix<double> patch_normal(3,3);
  for(size_t i = 0; i < 3; ++i){
      patch_normal(colon(), i) = normal(colon(), ps.get_all_patches()[patch_idx[i]].front());
    }
  if(dot(cross(patch_normal(colon(),0), patch_normal(colon(),1)), patch_normal(colon(),2)) < 0){
      return 1;
    }else return 0;
}

void optimal_sh_generator_polycube::opt2(zjucad::matrix::matrix<double> & field,
                                         boost::property_tree::ptree &pt)
{
  init(pt);
  jtf::mesh::patch_separater ps(orig_tm_.outside_face_, orig_tm_.ea_outside_);
  ps.separater(chain_);

  {
    ofstream ofs("patch.vtk");
    tri2vtk(ofs, &orig_tm_.tetmesh_.node_[0], orig_tm_.tetmesh_.node_.size(2),
        &orig_tm_.outside_face_[0], orig_tm_.outside_face_.size(2));
    cell_data(ofs, &ps.get_all_face2patch()[0], ps.get_all_face2patch().size(), "patch_idx");
  }

  matrix<double> point_normal;
  if(point_normal.size() == 0)
    jtf::mesh::cal_point_normal(orig_tm_.outside_face_, orig_tm_.tetmesh_.node_, point_normal);

  jtf::tet_mesh try_tm = orig_tm_;

  vector<size_t> picked_patch_idx;
  if(load_picked_patches(picked_patch_idx, pt)){
      throw std::invalid_argument("invalid picked patches file.");
    }

  for(size_t pi = 0; pi < picked_patch_idx.size(); ++pi){
      const size_t patch_idx = picked_patch_idx[pi];
      const vector<size_t> & one_patch = ps.get_patch(patch_idx);
      for(size_t fi = 0; fi < one_patch.size(); ++fi){
          try_tm.outside_face_normal_(colon(), one_patch[fi]) = polycube_tm_.outside_face_normal_(colon(), one_patch[fi]);
        }
    }

  {
    ofstream ofs("selected_patch.vtk");
    vector<size_t> selected_patch_faces;
    for(size_t pi = 0; pi < picked_patch_idx.size(); ++pi){
        selected_patch_faces.insert(selected_patch_faces.end(), ps.get_patch(picked_patch_idx[pi]).begin(),
                                    ps.get_patch(picked_patch_idx[pi]).end());
      }
    itr_matrix<const size_t*> selected_face(selected_patch_faces.size(),1, &selected_patch_faces[0]);

    matrix<size_t> selected_faces = try_tm.outside_face_(colon(), selected_face);
    const vector<size_t> &f2p = ps.get_all_face2patch();
    itr_matrix<const size_t*> f2p_m(f2p.size(),1,&f2p[0]);
    matrix<size_t> selected_face_type = f2p_m(selected_face,0);

    tri2vtk(ofs, &orig_tm_.tetmesh_.node_[0], orig_tm_.tetmesh_.node_.size(2),
        &selected_faces[0], selected_faces.size(2));
    cell_data(ofs, &selected_face_type[0], selected_face_type.size(), "patch_idx");
  }

  optimal_sh_generator osg(try_tm);
  osg.opt(field, pt);
  matrix<double> zyz;
  write_zyz(field, zyz, true);

  valid_field(zyz);

  {
    jtf::mesh::write_matrix("test_normal", try_tm.outside_face_normal_);
  }
}

void optimal_sh_generator_polycube::search_nearest_k_chain(
    const deque<pair<size_t,size_t> > &one_chain, const size_t k_nearest, set<size_t> &nearest_chains)
{
  set<size_t> points_of_this_chain;
  for(const auto & one_edge : one_chain) {
      points_of_this_chain.insert(one_edge.first);
      points_of_this_chain.insert(one_edge.second);
    }

  vector<pair<double,size_t> > dis;
  for(size_t ci = 0; ci < chain_.size(); ++ci){
      dis.push_back(make_pair(dis_of_two_chains(one_chain, chain_[ci], orig_tm_.tetmesh_.node_), ci));
    }
  sort(dis.begin(), dis.end());
  nearest_chains.clear();
  for(size_t i = 0; i < k_nearest; ++i)
    nearest_chains.insert(dis[i].second);
}

double optimal_sh_generator_polycube::dis_of_two_chains(
    const deque<pair<size_t,size_t> > & chain0, const deque<pair<size_t,size_t> > & chain1,
    const matrix<double> & node)
{
  vector<double> dis(chain0.size()+1,std::numeric_limits<double>::max());
  for(size_t ei = 0; ei < chain0.size(); ++ei) {
      for(size_t ej = 0; ej < chain1.size(); ++ej){
          const double len = norm(node(colon(), chain1[ej].first) - node(colon(), chain0[ei].first));
          if(len < dis[ei]) dis[ei] = len;
        }
      const double len = norm(node(colon(), chain1.back().second) - node(colon(), chain0[ei].first));
      if(len < dis[ei]) dis[ei] = len;
    }

  for(size_t ej = 0; ej < chain1.size(); ++ej){
      const double len = norm(node(colon(), chain1[ej].first) - node(colon(), chain0.back().second));
      if(len < dis.back()) dis.back() = len;
    }
  const double len = norm(node(colon(), chain1.back().second) - node(colon(), chain0.back().second));
  if(len < dis.back()) dis.back() = len;
  return std::accumulate(dis.begin(), dis.end(),0.0)/dis.size();
}

int optimal_sh_generator_polycube::adjust_surface_normal(
    zjucad::matrix::matrix<double> & surface_normal,
    const zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & frame,
    const jtf::mesh::patch_separater &ps,
    const zjucad::matrix::matrix<double> &polycube_face_normal)
{
  static vector<deque<pair<size_t,size_t> > > chain_list;
  static vector<deque<size_t> > singularities_type;
  static set<size_t> surface_points(orig_tm_.outside_face_.begin(), orig_tm_.outside_face_.end());
  static matrix<double> point_normal;
  if(point_normal.size() == 0)
    jtf::mesh::cal_point_normal(orig_tm_.outside_face_, orig_tm_.tetmesh_.node_, point_normal);

  singularity_extractor se(orig_tm_);
  se.extract(frame, chain_list, singularities_type);

  {
    static int count = 0;
    stringstream ss;
    ss << "singularity_" << count++ << ".vtk";
    dump_singularity_to_vtk(ss.str().c_str(),orig_tm_.tetmesh_.node_, chain_list);
  }

  static set<size_t> prev_patches_handled;
  const size_t k = 1;
  set<size_t> patches_to_use_fixed_normal;
  set<size_t> nearest_chains;
  matrix<double> avg_normal(3,1);
  int is_singularity_ok = 1;
  for(size_t ci = 0; ci < chain_list.size(); ++ci){
      nearest_chains.clear();
      if(rule0(chain_list[ci], surface_points, point_normal) == false){
          search_nearest_k_chain(chain_list[ci], k, nearest_chains);
          for(const auto & chain_idx : nearest_chains){
              set<size_t> adj_patches = ps.get_chain_adj_patches(chain_idx);
              patches_to_use_fixed_normal.insert(adj_patches.begin(), adj_patches.end());
            }
          /// update the patch normal
          for(const auto & patch_idx : patches_to_use_fixed_normal){
              cal_weighted_average_normal(ps.get_patch(patch_idx), polycube_face_normal, polycube_tm_.outside_face_area_, avg_normal);
              for(const auto & one_face : ps.get_patch(patch_idx)){
                  surface_normal(colon(), one_face) = avg_normal;
                }
            }
          is_singularity_ok = 0;
        }
    }
  if(prev_patches_handled == patches_to_use_fixed_normal){
      return -1;
    }else{
      prev_patches_handled = patches_to_use_fixed_normal;
    }
  return is_singularity_ok;
}

void optimal_sh_generator_polycube::init(boost::property_tree::ptree &pt)
{
  if(norm(orig_tm_.tetmesh_.mesh_ - polycube_tm_.tetmesh_.mesh_) > 1e-6){
      throw std::invalid_argument("orig tetmesh is not compatible with polycube tet.");
    }

  const double angle_cos = cos(45.*My_PI()/180.);

  if(zjucad::has("input/surface_type.value",pt)){
      boost::unordered_map<size_t,size_t> surface_type;
      if(load_surface_restricted_type(pt.get<string>("input/surface_type.value").c_str(),
                                      surface_type))
        throw std::invalid_argument("can not open surface type file.");
      matrix<size_t> surface_face_type = ones<size_t>(orig_tm_.outside_face_idx_.size(),1)*-1;
      for(size_t fi = 0; fi < orig_tm_.outside_face_idx_.size(); ++fi){
          auto it = surface_type.find(orig_tm_.outside_face_idx_[fi]);
          if(it != surface_type.end())
            surface_face_type[fi] = it->second;
          else{
              cerr << "# [error] strange, can not find surface_type." << endl;
            }
        }
      vector<pair<size_t,size_t> > all_edges;
      for(size_t ei = 0; ei < orig_tm_.ea_outside_->edge2cell_.size(); ++ei){
          const pair<size_t,size_t> & one_edge = orig_tm_.ea_outside_->edges_[ei];
          const pair<size_t,size_t> & tri_pair = orig_tm_.ea_outside_->edge2cell_[ei];
          assert(orig_tm_.ea_outside_->is_boundary_edge(tri_pair) != true);
          if(surface_face_type[tri_pair.first] != surface_face_type[tri_pair.second])
            all_edges.push_back(one_edge);
        }
      jtf::util::extract_chain_from_edges(all_edges, chain_);
      cerr << "# [info] use surface type info." << endl;
    }else{
      // find polycube_sharp_edges
      zjucad::matrix::matrix<size_t> polycube_feature_line;
      jtf::mesh::extract_mesh_feature_line(
            polycube_tm_.outside_face_, polycube_tm_.tetmesh_.node_, polycube_feature_line, angle_cos);

      vector<pair<size_t,size_t> > all_edges;
      for(size_t ei = 0; ei < polycube_feature_line.size(2); ++ei){
          all_edges.push_back(make_pair(polycube_feature_line(0,ei), polycube_feature_line(1,ei)));
        }
      jtf::util::extract_chain_from_edges(all_edges, chain_);
    }

  is_sharp_chain_.resize(chain_.size(), false);
  for(size_t ci = 0; ci < chain_.size(); ++ci){
      const deque<pair<size_t,size_t> > & one_chain = chain_[ci];
      for(size_t ei = 0; ei < one_chain.size(); ++ei){
          const size_t edge_idx = orig_tm_.ea_outside_->get_edge_idx(one_chain[ei].first, one_chain[ei].second);
          const pair<size_t,size_t> & one_face_pair = orig_tm_.ea_outside_->edge2cell_[edge_idx];
          if(fabs(dot(orig_tm_.outside_face_normal_(colon(), one_face_pair.first),
                      orig_tm_.outside_face_normal_(colon(), one_face_pair.second))) < angle_cos){
              is_sharp_chain_[ci] = true;
              break;
            }
        }
    }
  {
    dump_singularity_to_vtk("chain.vtk", orig_tm_.tetmesh_.node_, chain_);

    vector<deque<pair<size_t,size_t> > > sharp_chain;
    for(size_t ci = 0; ci < chain_.size(); ++ci){
        if(is_sharp_chain_[ci]) sharp_chain.push_back(chain_[ci]);
      }
    dump_singularity_to_vtk("sharp_chain.vtk", orig_tm_.tetmesh_.node_, sharp_chain);
  }
}
