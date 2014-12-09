/* detail control param
/usr/bin/time ../../bin/polycube prog=polycube3 linear_solver/type=direct linear_solver/name=cholmod iter_w=20 output="a.tet" iter=100 tet=../../dat/kitty-4.8k.mesh-split-surface-tet.mesh normal_align_w=5e-2 L1_sqrt_eps=5e-1 epsg=1e-3 adj_normal_w=2e1
 */
#include <boost/property_tree/ptree.hpp>

#include <string>
#include <fstream>
#include <numeric>
#include <hjlib/function/func_aux.h>
#include <hjlib/math/polar.h>
#include <hjlib/sparse/sparse.h>

#include <zjucad/optimizer/optimizer.h>
#include <zjucad/ptree/ptree.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"
#include "../common/util.h"
#include "../common/IO.h"
#include <jtflib/mesh/mesh.h>
#include <jtflib/optimizer/optimizer.h>

#include "io.h"
#include "polycube_surface_func.h"
#include "smooth_L1.h"
#include "../mesh_func/tri-area-normal-func.h"

#include "polycube_surface_func.h"

#include "tet_func.h"
#include "../mesh_func/tri-area.h"
#include "../mesh_func/tri-normal.h"
#include "../mesh_func/tri-area-normal.h"
#include "util.h"

#include "quality.h"

using namespace std;
using boost::property_tree::ptree;
using namespace hj::function;
using namespace zjucad::matrix;

extern double L1_sqrt_eps;
double arap_w_control;
double fix_zero_w_control;
double adj_normal_w_control;
double fl_w_control;
double patch_w_control;

// move surface point to the center of one-ring neighborhood by a little
void pre_smooth_control(matrix<double> &node, const matrix<double> &faces)
{
  matrix<size_t> nb_num = zeros<size_t>(node.size(2), 1);
  matrix<double> ct = zeros<double>(3, node.size(2));
  for(size_t fi = 0; fi < faces.size(2); ++fi) {
      nb_num(faces(colon(), fi)) += 1;
      for(int i = 0; i < 3; ++i) {
          ct(colon(), faces(i, fi)) += node(colon(), faces((i+1)%3, fi)) + node(colon(), faces((i+2)%3, fi));
        }
    }
  const double alpha = 0.3;
  for(size_t ni = 0; ni < node.size(2); ++ni) {
      if(nb_num[ni] == 0) continue;
      node(colon(), ni) = (1-alpha)*node(colon(), ni) + alpha*ct(colon(), ni) / nb_num[ni] / 2;
    }
}

int add_control(const matrix<size_t> & tet,
                boost::property_tree::ptree &pt,
                matrix<double>  &point_weight_diffused,
                boost::unordered_set<pair<size_t,size_t> > &feature_edgs,
                vector<matrix<size_t> > &patch_face)
{
  if(zjucad::has("point_diffused.value",pt)){
      if(load_point_diffused_value(pt.get<string>("point_diffused.value").c_str(),
                                   point_weight_diffused))
        return __LINE__;
      cerr << "# [info] control: use point_diffused value " << endl;
    }

  if(zjucad::has("feature_line.value",pt) && zjucad::has("s2v.value",pt)){
      vector<vector<size_t> > feature_line;
      if(jtf::mesh::load_feature_line(pt.get<string>("feature_line.value").c_str(),
                                      pt.get<string>("s2v.value").c_str(),
                                      feature_line))
        return __LINE__;

      pt.put("fl_w.desc", "weighting align given feature line");
      fl_w_control = pt.get<double>("fl_w.value", 0.1);

      cerr << "# [info] control: use feature line " << endl;

      pair<size_t,size_t> one_edge;

      for(size_t li = 0; li < feature_line.size(); ++li){
          const vector<size_t> & one_line = feature_line[li];
          for(size_t pi = 0; pi < one_line.size() -1; ++pi){
              one_edge.first = one_line[pi];
              one_edge.second = one_line[pi+1];
              if(one_edge.first > one_edge.second)
                swap(one_edge.first, one_edge.second);
              feature_edgs.insert(one_edge);
            }
        }
#if 0 // visual
      {
        ofstream ofs("feature_line.vtk");
        vector<size_t> edges;

        for(size_t li = 0; li < feature_line.size(); ++li){
            const vector<size_t> & one_line = feature_line[li];
            for(size_t pi = 0; pi < one_line.size()-1; ++pi){
                edges.push_back(one_line[pi]);
                edges.push_back(one_line[pi+1]);
              }
          }
        line2vtk(ofs, &node[0], node.size(2), &edges[0], edges.size()/2);
      }
#endif
    }

  if(zjucad::has("patch.value",pt) && zjucad::has("s2v.value",pt)
     && zjucad::has("obj.value",pt)){
      if(load_surface_patch(pt.get<string>("patch.value").c_str(),
                            pt.get<string>("obj.value").c_str(),
                            pt.get<string>("s2v.value").c_str(),
                            patch_face))
        return __LINE__;

      pt.put("patch_w.desc", "weighting smooth given patch surface");
      patch_w_control = pt.get<double>("patch_w.value", 10);

      cerr << "# [info] control: use patch smooth." << endl;
#if 0 // visual
      {
        ofstream ofs("patch.vtk");
        vector<size_t> patch_faces_vis;
        vector<size_t> patch_type;
        for(size_t pi = 0; pi < patch_faces.size(); ++pi){
            patch_faces_vis.insert(patch_faces_vis.end(), patch_faces[pi].begin(),
                                   patch_faces[pi].end());
            patch_type.insert(patch_type.end(), patch_faces[pi].size(2), pi);
          }
        tri2vtk(ofs, &node[0], node.size(2), &patch_faces_vis[0], patch_faces_vis.size()/3);
        cell_data(ofs, &patch_type[0], patch_type.size(), "patch_idx");
      }
#endif
    }

  boost::unordered_map<size_t,size_t> surface_type;
  if(zjucad::has("surface_type.value",pt)){
      if(load_surface_restricted_type_static(
           pt.get<string>("surface_type.value").c_str(), surface_type))
        return __LINE__;
      cerr << "# [info] use surface type as patches." << endl;
      convert_surface_type_to_surface_patches(tet, surface_type,patch_face,
                                              feature_edgs);
    }
  return 0;
}

int polycube3_control(ptree &pt)
{

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm.node_, &tm.mesh_))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
      cerr  << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }
  matrix<size_t> faces;
  get_outside_face(*fa, faces);

  cerr << "# read in tet " << tm.node_.size(2) << " " << tm.mesh_.size(2) << " " << faces.size(2) << endl;

  pt.put("L1_sqrt_eps.desc", "L1 sqrt eps, default is 1e-1.");
  pt.put("normal_align_w.desc", "surface L1 weight, default is 1e2.");
  pt.put("adj_normal_w.desc", "surface smoothness weight, default is 0.");
  pt.put("arap_w.desc", "arap weight, default is 1.");
  pt.put("iter_w.desc", "iterative weight, default is 1");
  pt.put("fix_zero_w.desc", "fix first node to zero, default is 1");
  pt.put("point_diffused.desc", "point_diffuesd weighting file.");
  pt.put("div_L.desc", "L1_sqrt_eps div L, default is 1.0");

  L1_sqrt_eps = pt.get<double>("L1_sqrt_eps.value", 5e-1);
  arap_w_control = pt.get<double>("arap_w.value",1);
  double normal_align_w = pt.get<double>("normal_align_w.value", 5e-2);
  const double div_L = pt.get<double>("div_L.value",1.0);
  adj_normal_w_control = pt.get<double>("adj_normal_w.value", 0.0);
  fix_zero_w_control = pt.get<double>("fix_zero_w.value", 1);
  const size_t iter_w = pt.get<size_t>("iter_w.value",1);
  double anti_flip_w = pt.get<double>("anti_flip_w.value", 1e-2);

  const double epsg = pt.get<double>("epsg.value", 1e-6);
  const double epsf = pt.get<double>("epsf.value", 1e-3);
  matrix<double> areas = zeros<double>(faces.size(2),1);
  for(size_t fi = 0; fi < faces.size(2); ++fi){
      areas[fi] = jtf::mesh::cal_face_area(faces(colon(),fi), tm.node_);
    }

  const double total_area = std::accumulate(areas.begin(), areas.end(), 0.0);
  cerr << "total_area: " << total_area << endl;
  matrix<double> R = eye<double>(3);
  ostringstream vtk_path_pref;
  vtk_path_pref << normal_align_w << "-" << L1_sqrt_eps
                << "-" << adj_normal_w_control << "-";

  matrix<double> node = tm.node_;
  if(pt.get_child_optional("init.value")) {
      jtf::mesh::meshes init;
      if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("init.value").c_str(), &init.node_, &init.mesh_))
        return __LINE__;
      if(max(abs(init.mesh_ - tm.mesh_)) != 0) {
          cerr << "#incompatible init value" << endl;
          return __LINE__;
        }
      node = init.node_;
      cerr << "# use init value." << endl;
      cerr << tm.node_(colon(), colon(0, 5))
           << node(colon(), colon(0, 5)) << endl;
      cerr << "# pre-smooth" << endl;
      pre_smooth_control(node, faces);
    }


  ///////////////// add control ///////////////////
  matrix<double> point_weight_diffused = ones<double>(node.size(2),1);

  boost::unordered_set<pair<size_t,size_t> > feature_edgs;
  vector<matrix<size_t> > patch_faces;

  add_control(tm.mesh_, pt, point_weight_diffused, feature_edgs, patch_faces);

  {
    ostringstream patch_face_os;
    patch_face_os << "orig-patch_face.vtk";
    ofstream ofs(patch_face_os.str().c_str());
    vector<size_t> patch_face_vec;
    for(size_t pi = 0; pi < patch_faces.size(); ++pi){
        patch_face_vec.insert(patch_face_vec.end(), patch_faces[pi].begin(), patch_faces[pi].end());
      }
    tri2vtk(ofs, &tm.node_[0], tm.node_.size(2), &patch_face_vec[0], patch_face_vec.size()/3);
  }

  pt.put("use_remesh.desc", "whether use remesh or not, y/n[1/0]");
  const size_t use_remesh = pt.get<size_t>("use_remesh.value", 0);
  const matrix<double> ct = tm.node_*ones<double>(tm.node_.size(2), 1)/tm.node_.size(2);
  double opt_time = 0;
  double remesh_time = 0;

  for(size_t i = 0; i < iter_w; ++i){
      clock_t beg = clock();

      matrix<double> areas_k = zeros<double>(faces.size(2),1);
      for(size_t fi = 0; fi < faces.size(2); ++fi){
          areas_k[fi] = jtf::mesh::cal_face_area(faces(colon(),fi),node);
        }

      if(1) { // opt global R
          R = eye<double>(3);

          unique_ptr<jtf::function::functionN1_t<double,int32_t> > func(build_polycube_rot_func2(node, faces,areas_k));
          jtf::optimize(*func,R,pt, nullptr, nullptr, nullptr);

          cout << R << endl;
          hj::polar3d p;
          p(R);
          cout << R << endl;
          node = temp(R*node);
        }

      //        if(use_remesh == 1) // local remesh
      //        {
      //            cerr << "# +++++++++++++++ begin local remesh ++++++++++++++++" << endl;
      //            cerr << "# [before remesh] tet number " << tm.mesh.size(2)
      //                 << " node number " << tm.node.size(2) << endl;

      //            local_remesh_tet(tm.mesh, tm.node, node);

      //            cerr << "# [after remesh] tet number " << tm.mesh.size(2)
      //                 << " node number " << tm.node.size(2) << endl;
      //            cerr << "# +++++++++++++++ end local remesh ++++++++++++++++\n" << endl ;
      //        }else{
      //            cerr << "# [info] no local remesh." << endl;
      //        }

      hj::function::function_t<double,int32_t> *diagnose;
      shared_ptr<function_t<double, int32_t> > func(
            build_polycube_function(
              tm.node_, node, node(colon(),0), tm.mesh_, faces, diagnose, &feature_edgs,
              &patch_faces, adj_normal_w_control, 0,0));

      cerr << "# [info] add polycube surface normal L1 function, wegiht "
           << normal_align_w << endl;
      cerr << "# [info] L1_sqrt_eps = " << L1_sqrt_eps << endl;

      vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > sum;
      vector<pair<jtf::function::functionN1_t<double,int32_t> *, double> > wf;
      sum.push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(jtf::function::least_square_warpper(func)));
      wf.push_back(make_pair(sum[0].get(), 1));
      if(normal_align_w > 0) {
          matrix<double> areas_k = zeros<double>(faces.size(2),1);
          for(size_t fi = 0; fi < faces.size(2); ++fi){
              areas_k[fi] = jtf::mesh::cal_face_area(faces(colon(),fi),node);
            }
          sum.push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                          build_smooth_L1_area_normal(node, faces, areas_k, normal_align_w)));
          wf.push_back(make_pair(sum[1].get(), normal_align_w));
        }
      if(anti_flip_w > 0) {// use unflip
          matrix<double> weight;
          shared_ptr<function_t<double, int32_t> > adj_normal
              (build_adj_normal_func(node, faces, weight));
          weight *= anti_flip_w;
          sum.push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(jtf::function::neg_log_warpper(adj_normal, weight)));
          cerr << "# use anti-flip with weight: " << anti_flip_w << endl;
        }
      shared_ptr<jtf::function::functionN1_t<double,int32_t> > target(new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD>(sum));
      wf.push_back(make_pair(target.get(), (1+normal_align_w)));

      shared_ptr<jtf::function::functionN1_t<double,int32_t> > constraint(new area_sum(tm.node_.size(2), faces, total_area));
      unnormalized_normal_quality_checker cb(node, tm.mesh_, faces, wf, *constraint);

      pt.put("epsg.value", epsg*(1+normal_align_w));

      int rtn = 0;
      if(normal_align_w > 0){
          vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > constraint_vec;
          constraint_vec.push_back(constraint);
          rtn = jtf::optimize(*target, node, pt, &constraint_vec, nullptr,&cb);
        }else
        rtn = jtf::optimize(*target, node, pt, nullptr, nullptr, nullptr);

      double opt_time_once = double(clock()-beg)/CLOCKS_PER_SEC ;
      cout << "# a loop: " << opt_time_once << endl;
      opt_time += opt_time_once;

      //        clock_t remesh_beg = clock();
      //        remove_degenerated_face_of_tet(tm.mesh, tm.node, node, faces, 5);
      //        double remesh_time_once = double(clock()-remesh_beg)/CLOCKS_PER_SEC ;
      //        remesh_time += remesh_time_once;

      const matrix<double> cur_ct = node*ones<double>(tm.node_.size(2), 1)/tm.node_.size(2);
      node += (ct-cur_ct)*ones<double>(1, tm.node_.size(2));

      //        if(rtn == -1){ // interior fail, need to reoptimize the shape with current weighting
      //            static vector<size_t> vis_log;

      //            if(vis_log.empty() || vis_log.back() != i){ // meet interior point at new step
      //                vis_log.clear();
      //                vis_log.push_back(i);
      //                --i;
      //                continue;
      //            }else{ // meet interior point many times,
      //                vis_log.push_back(i);
      //                if(vis_log.size() < 3){ // if stuck in this step for more than three times, skip on
      //                    --i;
      //                    continue;
      //                }
      //            }
      //        }

      { // visualize
        ostringstream vtk_path;
        vtk_path << vtk_path_pref.str() << i << ".vtk";
        ofstream ofs(vtk_path.str().c_str());
        tet2vtk(ofs, &node[0], node.size(2), &tm.mesh_[0], tm.mesh_.size(2));
      }

      { // visualize
        ostringstream vtk_path;
        vtk_path << vtk_path_pref.str() << i << ".orig.vtk";
        ofstream ofs(vtk_path.str().c_str());
        tet2vtk(ofs, &tm.node_[0], tm.node_.size(2), &tm.mesh_[0], tm.mesh_.size(2));
      }

      {
        ostringstream patch_face_os;
        patch_face_os << "patch_face-" << vtk_path_pref.str() << i << ".vtk";
        ofstream ofs(patch_face_os.str().c_str());
        vector<size_t> patch_face_vec;
        for(size_t pi = 0; pi < patch_faces.size(); ++pi){
            patch_face_vec.insert(patch_face_vec.end(), patch_faces[pi].begin(), patch_faces[pi].end());
          }
        tri2vtk(ofs, &node[0], node.size(2), &patch_face_vec[0], patch_face_vec.size()/3);
      }

      {
        ostringstream feature_line_os;
        feature_line_os << "feature_line-" << vtk_path_pref.str() << i << ".vtk";
        ofstream ofs(feature_line_os.str().c_str());
        vector<size_t> feature_edge_vec;
        for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit =
            feature_edgs.begin(); cit != feature_edgs.end(); ++cit){
            feature_edge_vec.push_back(cit->first);
            feature_edge_vec.push_back(cit->second);
          }
        line2vtk(ofs, &node[0], node.size(2), &feature_edge_vec[0], feature_edge_vec.size()/2);
      }

      if(polycube_L1_area_quality(&node[0], node.size(2), faces) < 1e-3) { // the global condition
          cout << "global converge." << endl;
          break;
        }

      normal_align_w *= 2.0;
      L1_sqrt_eps /= sqrt(2.0);
      if(L1_sqrt_eps < 1e-2)
        L1_sqrt_eps = 1e-2;
    }

  jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("output.value").c_str(), &node, &tm.mesh_);
  jtf::mesh::tet_mesh_write_to_zjumat("after_collapse_output.tet", &tm.node_, &tm.mesh_);

  cerr << "# [info] total opt time " << opt_time << endl;
  cerr << "# [info] total remesh time " << remesh_time << endl;
  cerr << "success." << endl;
  return 0;
}
