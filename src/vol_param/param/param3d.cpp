#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include <zjucad/ptree/ptree.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <zjucad/matrix/io.h>

#include "../descriptor/descriptor_vol.h"
#include "../descriptor/func_terms/linear_equation.h"
#include "../descriptor/func_terms/quality_check.h"
#include "../descriptor/func_terms/l1-normal.h"
#include "../solver/solver_ipopt.h"
#include "../solver/solver_jtf.h"

#include <boost/algorithm/string.hpp>

#include "../../hex_param/io.h"
#include "../common/util.h"
#include "../../common/vtk.h"
#include "../../common/zyz.h"
#include "util.h"
#include "../solver/util.h"
#include "../common/io.h"
#include <jtflib/math/null.h>

using namespace std;
using namespace jtf::mesh;
using namespace zjucad::matrix;
using boost::property_tree::ptree;

static void pt_description(ptree &pt)
{
  pt.put("input/tet.desc", "insput tet file.");
  pt.put("object/type.desc", "objetive type: arap/frame");
  pt.put("cons/weight.desc", "weight for each constraint, separated by ',' ");
  pt.put("solver.desc", "zjucad/ipopt");
  pt.put("cons/type.desc", "constraint type: equation");
  pt.put("cons/type.is_required", "n");
  pt.put("input/equation.desc", "input equation file");
  pt.put("input/node_mapping.desc","input node_mapping file");
  pt.put("input/node_mapping.is_required","n");
  pt.put("input/group.desc", "input variable group file.");
  pt.put("weight/linear_eqn_fit.desc", "weight to fit linear eqn");
  pt.put("input/iter_solver.desc", "iter number to dynamically adjust linear eqn/l1 normal align.");
  pt.put("weight/normal-diff.desc", "weight to control normal difference between two normals");
  pt.put("output/tet.desc", "output tet after parameterization");
  pt.put("output/new_uncut_tet.desc", "output new uncut_tet after local remesh.");
  pt.put("output/new_cut_tet.desc", "output new cut tet after local remesh");
  pt.put("output/inner_type.desc", "output inner face type after local remesh");
  pt.put("output/frame.desc", "output frame after local remesh");
  pt.put("output/equation.desc", "output equation after local remesh");
  pt.put("input/need_node_compression.desc", "need node compression, [y/n]");
}

template <typename OS>
void ouput(OS & os, const matrix<size_t> & tet,
           const matrix<double> & node,
           const matrix<size_t> * node_mapping = 0)
{
  tet2vtk(os, &node[0], node.size(2), &tet[0], tet.size(2));
  //    } else {
  //      //      matrix<double> assemble_node(3, *(max_element(tet.begin(), tet.end()))+1);
  //      //      for(size_t vi = 0; vi < assemble_node.size(); ++vi)
  //      //        assemble_node[vi] = node[(*node_mapping)[vi]];
  //      //      tet2vtk(os, &assemble_node[0], assemble_node.size(2), &tet[0], tet.size(2));
  //      matrix<double> assemble_node(node_mapping->size(),1);
  //      for(size_t i = 0; i < assemble_node.size(); ++i)
  //        assemble_node[i] = node[(*node_mapping)[i]];
  //      tet2vtk(os, &assemble_node[0], assemble_node.size()/3, &tet[0], tet.size(2));
  //    }
}

void output_remesh_configuration(
    const shared_ptr<const descriptor_vol> des,
    const boost::property_tree::ptree & pt,
    const matrix<size_t> & mesh,
    const matrix<double> & node)
{
  if(des->bi_.uncut_tet.size()){
      const string new_uncut_tet_str = pt.get<string>("output/new_uncut_tet.value");
      if(jtf::mesh::tet_mesh_write_to_zjumat(
           new_uncut_tet_str.c_str(),
           &des->bi_.uncut_node, &des->bi_.uncut_tet))
        return ;
    }
  {// cut tet
    const string new_cut_tet_str = pt.get<string>("output/new_cut_tet.value");
    if(jtf::mesh::tet_mesh_write_to_zjumat(new_cut_tet_str.c_str(),
                                           &node, &mesh))
      return ;
  }
  {// inner_type
    const string inner_type_str = pt.get<string>("output/inner_type.value");
    if(dump_inner_face_jump_type(inner_type_str.c_str(), des->bi_.inner_type))
      return ;
  }
  {// frame
    if(des->bi_.frames.size()){
        const string frame_str = pt.get<string>("output/frame.value");
        matrix<double> new_zyz(3, des->bi_.frames.size());
        for(size_t ti = 0; ti < des->bi_.frames.size(); ++ti){
            rotation_matrix_2_zyz_angle(&des->bi_.frames[ti][0], &new_zyz(0,ti),0);
          }
        if(jtf::mesh::write_matrix(frame_str.c_str(), new_zyz))
          return ;
      }
  }
}

int param3d(ptree &pt)
{
  pt_description(pt);

  meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(
       pt.get<string>("input/tet.value").c_str(),
       &tm.node_, &tm.mesh_))
    return __LINE__;

  shared_ptr<descriptor_vol> des(new descriptor_vol);;

  des->init(tm.mesh_, tm.node_, pt);
  des->set_objective(tm.mesh_, tm.node_, pt);
  des->set_constraint(tm.mesh_, tm.node_, pt);

  shared_ptr<solver_base> sv;//(new solver_zjucad);

  switch(str2int(pt.get<string>("package.value").c_str())) {
    case str2int("ipopt"): {
        sv.reset(new solver_ipopt); break;
      }
    case str2int("jtf"): {
        sv.reset(new solver_jtf); break;
      }
dafault:
      throw std::invalid_argument("unknown solver.");
    }

  matrix<double> init_node = tm.node_;
  if(zjucad::has("input/init_tet.value",pt)) {
      cerr << "# [info] use init_tet" << endl;
      matrix<size_t> temp_tet;
      if(tet_mesh_read_from_zjumat(pt.get<string>("input/init_tet.value").c_str(),
                                   &init_node, &temp_tet))
        return __LINE__;
    }

  double w_l1 = pt.get<double>("weight/l1_normal_align.value");
  size_t iter_solver = pt.get<size_t>("input/iter_solver.value",1);
  double eps = pt.get<double>("weight/l1_eps.value", 2.0);

  if(des->bi_.inner_type.empty()){
      if(load_inner_face_jump_type(pt.get<string>("input/inner_type.value").c_str(),
                                   des->bi_.inner_type))
        return __LINE__;
    }

  const string need_remesh = pt.get<string>("input/need_remesh.value", "n");
  bool need_remesh_bool ;
  string remesh_strategy = "surface";
  if(need_remesh == "y" || need_remesh == "Y"
     || need_remesh == "yes" || need_remesh == "YES"){
      need_remesh_bool = true;
      remesh_strategy = pt.get<string>("input/remesh_strategy.value");
    }
  else
    need_remesh_bool = false;

  size_t remesh_i = 0;
  vector<size_t> vis_log;
  for(size_t i = 0 ; i < iter_solver; ++i){

      w_l1 *= 2;
      eps /= 2;
      if(eps < 0.01) eps = 0.01;

      cerr << "# [info] smooth_L1::scalar_ "
           << smooth_L1::scalar_ << endl;

      set_l1_normal_align_weight(w_l1);
      set_l1_eps_weight(eps);

      if(des->get_callback()){
          sv->set_callback(des->get_callback());
        }
      const int solver_rtn = sv->solve(init_node, des, pt);

      //      for(size_t k = 0; k < local_stiffen_iter; ++k){ // local stiffening
      //          dv_ptr->frame_local_stiffening(tm.mesh_, tm.node_, init_node);
      //          sv->solve(init_node, des, pt);
      //        }

      if(need_remesh_bool && solver_rtn == -1){
          int rtn = remove_degenerated_face_of_tet(
                des->bi_.orig_face_in_cut, *(des->bi_.fa_uncut), des->bi_.uncut_tet, des->bi_.uncut_node,
                tm.mesh_, tm.node_, init_node,des->bi_.inner_type, des->bi_.frames, *(des->bi_.ea_orig_in_cut) ,remesh_strategy);
          if(rtn > 0){// removed degenerated faces
              if(vis_log.size() > 0 && vis_log.back() != i) {
                  vis_log.clear();
                  vis_log.push_back(i);

                  des->reset(des->bi_.uncut_tet, des->bi_.uncut_node,
                             tm.mesh_, tm.node_, pt);
                  --i;
                  if(i == -1) i = 0;

                  w_l1 /= 2;
                  eps *= 2;
                }else{
                  if(vis_log.size() < 3){ // only run 3 times
                      vis_log.push_back(i);
                      des->reset(des->bi_.uncut_tet, des->bi_.uncut_node,
                                 tm.mesh_, tm.node_, pt);
                      --i;
                      if(i == -1) i = 0;
                      w_l1 /= 2;
                      eps *= 2;
                    }
                }
              const size_t smooth_iter = 2;
              smooth_mesh_surface(des->bi_.orig_face_in_cut, tm.node_, smooth_iter);
              smooth_mesh_surface(des->bi_.outside_face_uncut,des->bi_.uncut_node,smooth_iter);
            }
          stringstream ss,ss2;
          ss << "remesh_cut_" << remesh_i << ".vtk";
          ofstream ofs(ss.str().c_str());
          tet2vtk(ofs, &init_node[0], init_node.size(2), &tm.mesh_[0], tm.mesh_.size(2));

          ss2 << "remesh_uncut_" << remesh_i++ << ".vtk";
          ofstream ofs_uncut(ss.str().c_str());
          tet2vtk(ofs_uncut, &des->bi_.uncut_node[0], des->bi_.uncut_node.size(2), &des->bi_.uncut_tet[0], des->bi_.uncut_tet.size(2));

        }

      stringstream ss;
      ss << pt.get<string>("output/tet.value") << "_" << i << ".vtk";
      ofstream os(ss.str().c_str());
      ouput(os, tm.mesh_, init_node, 0);

      if(des->get_objective("l1-normal") == nullptr) continue;

      if(polycube_L1_area_quality(&init_node[0], init_node.size(2),
                                  des->bi_.orig_face_in_cut) < 1e-3){
          cout << "global convege." << endl;
          break;
        }
    }
  string out_tet = pt.get<string>("output/tet.value");

  jtf::mesh::tet_mesh_write_to_zjumat(out_tet.c_str(), &init_node, &tm.mesh_);

  output_remesh_configuration(des,pt,tm.mesh_, tm.node_);

  return 0;
}


int integer_rounding(std::shared_ptr<descriptor_base> des,
                     const matrix<size_t> &mesh,
                     matrix<double> &node,
                     boost::property_tree::ptree &pt,
                     vector<vector<size_t> > & integer_variants,
                     vector<vector<pair<size_t,size_t> > > & uvw_restricted_edges);

int param3d_integer(ptree &pt)
{
  meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(
       pt.get<string>("input/tet.value").c_str(),
       &tm.node_, &tm.mesh_))
    return __LINE__;

  shared_ptr<solver_base> sv;//(new solver_zjucad);

  switch(str2int(pt.get<string>("package.value").c_str()))
    {
    case str2int("ipopt"):{
        sv.reset(new solver_ipopt); break;
      }
    case str2int("jtf"): {
        sv.reset(new solver_jtf); break;
      }
defalut:
      throw std::invalid_argument("unknown solver.");
    }

  matrix<double> init_node = tm.node_;
  if(zjucad::has("input/init_tet.value",pt))
    {
      matrix<size_t> temp_tet;
      if(tet_mesh_read_from_zjumat(pt.get<string>("input/init_tet.value").c_str(),
                                   &init_node, &temp_tet))
        return __LINE__;
      cerr << "# [info] use init node " << pt.get<string>("input/init_tet.value") << endl;
    }

  vector<vector<size_t> > integer_variants;
  if(load_integer_groups(pt.get<string>("input/integer_group.value").c_str(),
                         integer_variants))
    return __LINE__;

  vector<vector<pair<size_t,size_t> > > uvw_restricted_edges;
  if(load_restricted_path(pt.get<string>("input/restricted_edge.value").c_str(),
                          uvw_restricted_edges))
    return __LINE__;

  cerr << "# [info] restricted edges number " << uvw_restricted_edges[0].size()
          + uvw_restricted_edges[1].size() + uvw_restricted_edges[2].size() << endl;

  double bb_size = calc_bounding_sphere_size(init_node);
  const double min_dis = get_minimal_distance(integer_variants, init_node, bb_size/100);
  double scale = 1.0;
  if(1){
      scale = 1.0/min_dis;
    }
  if(zjucad::has("weight/scale.value",pt)) {
      scale = pt.get<double>("weight/scale.value");
    }

  cerr << "# [info] bounding_sphere_size " << bb_size << endl;
  cerr << "# [info] minimal distance of integer varibles is " << min_dis << endl;
  cerr << "# [info] global scale " << scale << endl;
  init_node *= scale;
  tm.node_ *= scale;

  double xyz_scale[] = {1,1,1};
  if(zjucad::has("weight/x_scale.value",pt)){
      xyz_scale[0] = pt.get<double>("weight/x_scale.value");
    }
  if(zjucad::has("weight/y_scale.value",pt)){
      xyz_scale[1] = pt.get<double>("weight/y_scale.value");
    }
  if(zjucad::has("weight/z_scale.value",pt)){
      xyz_scale[2] = pt.get<double>("weight/z_scale.value");
    }
  for(size_t i = 0; i < 3; ++i){
      init_node(i,colon()) *= xyz_scale[i];
      tm.node_(i, colon()) *= xyz_scale[i];
    }


  shared_ptr<descriptor_base> des(new descriptor_vol);

  des->init(tm.mesh_, tm.node_, pt);
  des->set_objective(tm.mesh_, tm.node_, pt);
  des->set_constraint(tm.mesh_, tm.node_, pt);

  integer_rounding(des, tm.mesh_, init_node, pt, integer_variants, uvw_restricted_edges);

  string out_tet = pt.get<string>("output/tet.value");
  if(tet_mesh_write_to_zjumat(out_tet.c_str(), &init_node, &tm.mesh_))
    return __LINE__;

  return 0;
}
