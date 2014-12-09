#include <fstream>
#include <memory>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <hjlib/function/operation.h>
#include <zjucad/optimizer/optimizer.h>
#include <zjucad/ptree/ptree.h>

#include "hex_process.h"
#include "hex_deform_function.h"
#include "../common/vtk.h"
#include <jtflib/util/container_operation.h>

#include "../quadmesh/util.h"
#include "../hexmesh/util.h"
#include "../tetmesh_refine/tetmesh_refine.h"
#include "../tetmesh_refine/relax_function.h"
#include "../common/cell_quality.h"

using namespace std;
using namespace zjucad::matrix;
using namespace hj::function;
//using namespace  jtf::hexmesh;

int deform_hex_to_original_tet(matrixst & hex,
                               matrixd & hex_node,
                               const matrixst & tet,
                               const matrixd & tet_node,
                               boost::property_tree::ptree &pt)
{
  //create_map_from_hex_to_tet(tet_node, );
  matrixst tet_surface, hex_surface;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa_tet(jtf::mesh::face2tet_adjacent::create(tet));
  unique_ptr<jtf::mesh::face2hex_adjacent> fa_hex(jtf::mesh::face2hex_adjacent::create(hex));
  if(!fa_tet.get() || !fa_hex.get())
    return __LINE__;

  //cerr << "# [info] glue hex on surface " << endl;
  //glue_hex_on_surface(hex, *fa_hex);
  jtf::mesh::get_outside_face(*fa_tet, tet_surface);
  jtf::mesh::get_outside_face(*fa_hex, hex_surface);


  ///////////////// pre-process //////////////////////////
  matrixst hex_surface_points, tet_surface_points;
  boost::unordered_map<size_t, boost::unordered_set<size_t> > point_adj_points;
  std::unique_ptr<jtf::mesh::edge2cell_adjacent> tet_surface_edges(
        jtf::mesh::edge2cell_adjacent::create(tet_surface));
  //get one ring face at point
  jtf::mesh::one_ring_face_at_point ring_faces;
  ring_faces.add_all_faces(tet_surface, *tet_surface_edges);
  boost::unordered_set<size_t> hex_surface_point_set(hex_surface.begin(),
                                                     hex_surface.end());
  { // get outside face points
    boost::unordered_set<size_t> tet_surface_point_set(tet_surface.begin(),
                                                       tet_surface.end());
    tet_surface_points.resize(tet_surface_point_set.size(),1);
    copy(tet_surface_point_set.begin(), tet_surface_point_set.end(),
         tet_surface_points.begin());

    boost::unordered_set<size_t> hex_surface_point_set(hex_surface.begin(),
                                                       hex_surface.end());
    hex_surface_points.resize(hex_surface_point_set.size(),1);
    copy(hex_surface_point_set.begin(), hex_surface_point_set.end(),
         hex_surface_points.begin());


    matrixst face_of_one_hex;

    // gather adj_points for each point p
    // if p is surface point, it only has surface points as its adjacents
    for(size_t hi = 0; hi < hex.size(2); ++hi){
      //      for(size_t pi = 0; pi < hex.size(1); ++pi){
      //        for(size_t npi = 1; npi < hex.size(1); ++npi){
      //          point_adj_points[hex(pi, hi)].insert(hex((pi + npi)%hex.size(1),hi));
      //        }
      //      }

      jtf::mesh::get_faces_for_one_hex(hex(colon(),hi),face_of_one_hex);
      for(size_t fi = 0; fi < face_of_one_hex.size(2); ++fi){
        for(size_t pi = 0; pi < face_of_one_hex.size(1); ++pi){
          point_adj_points[face_of_one_hex[pi]].insert(
                face_of_one_hex[(pi+1)%face_of_one_hex.size(1)]);
          point_adj_points[face_of_one_hex[(pi+1)%face_of_one_hex.size(1)]].insert(
                face_of_one_hex[pi]);
        }
      }
    }
    typedef boost::unordered_map<size_t, boost::unordered_set<size_t> >::iterator
        bucit;

    for(bucit bucit_ = point_adj_points.begin(); bucit_ !=  point_adj_points.end();
        ++bucit_){
      boost::unordered_set<size_t>::const_iterator buscit_p =
          hex_surface_point_set.find(bucit_->first);

      if(buscit_p != hex_surface_point_set.end()){// surface node
        boost::unordered_set<size_t> & adj_points_set = bucit_->second;
        for(boost::unordered_set<size_t>::const_iterator buscit =
            adj_points_set.begin(); buscit != adj_points_set.end();){
          if(hex_surface_point_set.find(*buscit) == hex_surface_point_set.end())
            adj_points_set.erase(buscit++);
          else
            ++buscit;
        }
      }
    }

#if 1
    {
      for(size_t fi = 0; fi < hex_surface.size(2); ++fi){
        for(size_t p = 0; p < hex_surface.size(1); ++p){
          point_adj_points[hex_surface(p,fi)].insert(hex_surface((p+1)%4, fi));
          point_adj_points[hex_surface(p,fi)].insert(hex_surface((p+2)%4, fi));
          point_adj_points[hex_surface(p,fi)].insert(hex_surface((p+3)%4, fi));
        }
      }
    }
#endif

  }

  matrixd tri_face_normal;
  matrixst tet_surface_idx;
  get_outside_face_idx(*fa_tet, tet_surface_idx);
  jtf::mesh::cal_face_normal(tet_surface, tet_node, tri_face_normal);
  jtf::tetmesh::orient_face_normal_outside_tetmesh(
        tet, tet_node, tet_surface,
        tet_surface_idx,*fa_tet, tri_face_normal);

  /////////////////////////////////////////////////////////////////

  pt.put("iter_num.desc", "iter number");
  const size_t iter_num = pt.get<size_t>("iter_num.value");
  for(size_t t = 0; t < iter_num; ++t){
    one_deformation(hex, hex_node, tet, tet_node, tet_surface,
                    hex_surface, tet_surface_points, hex_surface_points,
                    hex_surface_point_set, *fa_hex, *tet_surface_edges,
                    ring_faces, point_adj_points,tri_face_normal,
                    pt);
  }
  return 0;
}

hj::function::function_t<double,int32_t> *
build_hex_deform_func(
    const matrixst & hex,
    const matrixd & hex_node,
    const boost::unordered_set<size_t> & hex_surface_point_set,
    const boost::unordered_map<size_t,boost::unordered_set<size_t> > &point_adj_points,
    const string smooth_strategy)
{
  std::shared_ptr<vector<std::shared_ptr<function_t<double,int32_t> > > > deformation_func(
        new vector<std::shared_ptr<function_t<double,int32_t> > >());
  //#if use_laplacian
  if(smooth_strategy == "laplacian"){
    typedef boost::unordered_map<size_t, boost::unordered_set<size_t> >
        ::const_iterator bucit;
    vector<size_t> adj_points;
    for(bucit bucit_ = point_adj_points.begin(); bucit_ !=  point_adj_points.end();
        ++bucit_){
      boost::unordered_set<size_t>::const_iterator buscit_p =
          hex_surface_point_set.find(bucit_->first);

      if(buscit_p == hex_surface_point_set.end()){// surface node
        const boost::unordered_set<size_t> & adj_points_set = bucit_->second;
        adj_points.resize(adj_points_set.size());
        copy(adj_points_set.begin(), adj_points_set.end(), adj_points.begin());
        std::shared_ptr<laplacian_smooth_function> lsf(
              new laplacian_smooth_function(bucit_->first, hex_node.size(2),
                                            adj_points));
        deformation_func->push_back(lsf);
      }
    }
  }else if(smooth_strategy == "meshless") {
    const matrixd weight = ones<double>(8,1);
    matrix<double> hex_vol = zeros<double>(hex.size(2),1);
    for(size_t hi = 0; hi < hex.size(2); ++hi){
      hex_vol[hi] = fabs(hex_volume(hex_node(colon(),hex(colon(),hi))));
    }

    hex_vol /= std::accumulate(hex_vol.begin(), hex_vol.end(),0.0);

    for(size_t hi = 0; hi < hex.size(2); ++hi){
      const double hex_jac = hex_scaled_jacobian(hex_node(colon(), hex(colon(),hi)));

      const double w_ = (hex_jac>0?1e-3*hex_vol[hi]: hex_vol[hi]);

      for(size_t pi = 0; pi < 8; ++pi){
        std::shared_ptr<rigid_matching_function_p> rmf(
              new rigid_matching_function_p(
                hex(colon(), hi), hex_node, weight, pi, w_));
        deformation_func->push_back(rmf);
      }
    }
  }else if(smooth_strategy == "arap")
  {
    cerr << "# [error] hex arap is not well defined." << endl;
    return 0;
    //    matrixst tet_from_subdivied_hex;
    //    create_tet_from_hex_cut_raw(hex, hex_node, tet_from_subdivied_hex);

    //    for(size_t t = 0; t < tet_from_subdivied_hex.size(2); ++t){
    //      unique_ptr<arap_func> ddf(
    //            new arap_func(tet_from_subdivied_hex(colon(),t),hex_node,t));
    //      deformation_func->push_back(func_ptr(ddf.release()));
    //    }
  }else {
    cerr << "# [error] error smooth stratrgy." << endl;
    return 0;
  }
  cerr << "# [info] deformation func num: " << deformation_func->size() << endl;
  return new_catenated_function<double,int32_t>(deformation_func);
}

hj::function::function_t<double,int32_t> *
build_boundary_align_func(const matrixd & hex_node,
                          const matrixst & hex_surface_points,
                          const matrixd & nearest_node_of_hex_surface,
                          const boost::unordered_map<size_t, matrixd > &q2t)
{
  std::shared_ptr<vector<std::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<std::shared_ptr<function_t<double, int32_t> > >());
  for(size_t pi = 0; pi < hex_surface_points.size(); ++pi){
    if(q2t.find(hex_surface_points[pi]) == q2t.end()){
      std::shared_ptr<function_t<double,int32_t> > baf(
            new point_fix_function(hex_node.size(2),hex_surface_points[pi],
                                   nearest_node_of_hex_surface(colon(), pi)));
      funcs->push_back(baf);
    }
  }
  cerr << "# [info] boundary align func num: " << funcs->size() << endl;
  return new_catenated_function<double,int32_t>(funcs);
}

hj::function::function_t<double,int32_t> *
build_feature_align_func(const matrixd &hex_node,
                         const boost::unordered_map<size_t, matrixd >  &q2t)
{
  std::shared_ptr<vector<std::shared_ptr<function_t<double,int32_t> > > >
      funcs(new vector<std::shared_ptr<function_t<double,int32_t> > >());
  for(boost::unordered_map<size_t, matrixd >::const_iterator bmit = q2t.begin();
      bmit != q2t.end(); ++bmit){
    std::shared_ptr<point_fix_function> pff(
          new point_fix_function(hex_node.size(2),bmit->first,
                                 bmit->second));
    funcs->push_back(pff);
  }
  cerr << "# [info] feature align func num: " << funcs->size() << endl;
  return new_catenated_function<double,int32_t>(funcs);
}

hj::function::function_t<double,int32_t> *
build_tangent_moving_func(
    const matrixd &hex_node,
    const matrixst & hex_surface_points,
    const matrixd &nearest_node_of_hex_surface,
    const matrixst &nearest_node_of_hex_surface_belong_faces,
    const matrixd &nearest_node_of_hex_surface_weights,
    const matrixd & tet_node,
    const matrixst & tet_surface,
    const matrixd & tri_face_normal,
    const jtf::mesh::one_ring_face_at_point & orfap)
{
  std::shared_ptr<vector<std::shared_ptr<function_t<double,int32_t> > > >
      func(new vector<std::shared_ptr<function_t<double,int32_t> > >());
  // move points on tangent plane
  matrixd normal = zeros<double>(3,1);
  // matrixd advance_nearset_nodes = nearest_node_of_hex_surface;

  for(size_t pi = 0; pi < hex_surface_points.size(); ++pi){
    calc_point_normal(nearest_node_of_hex_surface_belong_faces[pi],
                      nearest_node_of_hex_surface_weights(colon(),pi),
                      orfap, tet_surface, tet_node, tri_face_normal,normal);
    //advance_nearset_nodes(colon(), pi) += 0.1 * normal;
    std::shared_ptr<function_t<double,int32_t> > tmf(
          new tangent_moveing_function(
            hex_surface_points[pi],
            hex_node.size(2),
            nearest_node_of_hex_surface(colon(), pi),
            normal));
    func->push_back(tmf);
  }
  cerr << "# [info] tangent moving func num: "
       << func->size() << endl;
  return new_catenated_function<double,int32_t>(func);
}

hj::function::function_t<double,int32_t>  *
build_surface_smooth_func(
    const matrixd &hex_node,
    const boost::unordered_set<size_t> & hex_surface_point_set,
    const boost::unordered_map<size_t,boost::unordered_set<size_t> > &point_adj_points,
    const boost::unordered_map<size_t, matrixd >  &q2t)
{
  std::shared_ptr<vector<std::shared_ptr<function_t<double,int32_t> > > >
      funcs(new vector<std::shared_ptr<function_t<double,int32_t> > >());

  typedef boost::unordered_map<size_t, boost::unordered_set<size_t> >
      ::const_iterator bucit;

  vector<size_t> adj_points;
  for(bucit bucit_ = point_adj_points.begin(); bucit_ !=  point_adj_points.end();
      ++bucit_){
    boost::unordered_set<size_t>::const_iterator buscit_p =
        hex_surface_point_set.find(bucit_->first);

    if(buscit_p != hex_surface_point_set.end()){// surface node

      const boost::unordered_set<size_t> & adj_points_set = bucit_->second;
      adj_points.resize(adj_points_set.size());
      copy(adj_points_set.begin(), adj_points_set.end(), adj_points.begin());

      std::shared_ptr<laplacian_smooth_function> lsf(
            new laplacian_smooth_function(bucit_->first, hex_node.size(2),
                                          adj_points));
      funcs->push_back( lsf);
    }
  }
  cerr << "# [info] smooth surface func num: " << funcs->size() << endl;
  return new_catenated_function<double,int32_t>(funcs);
}


int one_deformation(
    const matrixst & hex,
    matrixd & hex_node,
    const matrixst & tet,
    const matrixd & tet_node,
    const matrixst & tet_surface,
    const matrixst & hex_surface,
    const matrixst & tet_surface_points,
    const matrixst & hex_surface_points,
    const boost::unordered_set<size_t> & hex_surface_point_set,
    const jtf::mesh::face2hex_adjacent &fa_hex,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::one_ring_face_at_point & orfap,
    const boost::unordered_map<size_t,boost::unordered_set<size_t> > &point_adj_points,
    const matrixd & tri_face_normal,
    boost::property_tree::ptree &pt)
{
  matrixd nearest_node_of_hex_surface;
  matrixst nearest_node_of_hex_surface_belong_faces;
  matrixd nearest_node_of_hex_surface_weights;
  //  // find the corresponding node of hex surface
  create_map_from_hex_to_tet(tet_node,tet_surface, hex_node, hex_surface, hex,
                             tet_surface_points,fa_hex,ea,orfap,hex_surface_points,
                             nearest_node_of_hex_surface,
                             nearest_node_of_hex_surface_belong_faces,
                             nearest_node_of_hex_surface_weights);

  static matrixst tfl;
  matrixst qfl;
  boost::unordered_map<size_t, matrixd > q2t;
  {
    static size_t execute_number = 0;
    if(execute_number++ == 0){
      jtf::mesh::extract_mesh_feature_line(tet_surface,tet_node,tfl);
      ofstream ofs("tri_feature_line.vtk");
      line2vtk(ofs, &tet_node[0], tet_node.size(2), &tfl[0], tfl.size(2));
    }
  }
//  if(zjucad::has("quad_fl.value", pt)){
//    vector<deque<pair<size_t,size_t> > > feature_line;
//    if(load_feature_lines(pt.get<string>("quad_fl.value").c_str(), feature_line))
//      return __LINE__;
//    vector<size_t> temp_fl;
//    for(size_t di = 0; di < feature_line.size(); ++di){
//      const deque<pair<size_t,size_t> > & one_deq = feature_line[di];
//      for(size_t ei = 0; ei < one_deq.size(); ++ei){
//        temp_fl.push_back(one_deq[ei].first);
//        temp_fl.push_back(one_deq[ei].second);
//      }
//    }
//    qfl.resize(2, temp_fl.size()/2);
//    std::copy(temp_fl.begin(), temp_fl.end(), qfl.begin());
//  }else
//    jtf::mesh::extract_mesh_feature_line(hex_surface,hex_node,qfl);
//  {
//    ofstream ofs("quad_feature_line.vtk");
//    line2vtk(ofs, &hex_node[0], hex_node.size(2), &qfl[0], qfl.size(2));
//  }
//  project_quad_feature_to_tri_feature_new(qfl, hex_node, tfl, tet_node, q2t);

  {
    matrix<double> project_nodes(3, q2t.size()) ;
    matrixst project_nodes_idx(q2t.size()) ;
    matrix<double> orig_nodes(3, q2t.size()) ;
    matrixst orig_nodes_idx(q2t.size()) ;
    size_t idx = 0;
    for(boost::unordered_map<size_t, matrixd >::const_iterator cit = q2t.begin();
        cit != q2t.end(); ++cit, ++idx){
      project_nodes(colon(), idx) = cit->second;
      orig_nodes(colon(), idx) = hex_node(colon(),cit->first);
      project_nodes_idx[idx] = idx;
      orig_nodes_idx[idx] = idx;
    }

    ofstream ofs("point_project.vtk");
    point2vtk(ofs, &project_nodes[0], project_nodes.size(2),
              &project_nodes_idx[0], project_nodes_idx.size());
    cell_data(ofs, &project_nodes_idx[0], project_nodes_idx.size(), "idx");

    ofstream ofs_orig("orig_project.vtk");
    point2vtk(ofs_orig, &orig_nodes[0], orig_nodes.size(2),
              &orig_nodes_idx[0], orig_nodes_idx.size());
    cell_data(ofs_orig, &orig_nodes_idx[0], orig_nodes_idx.size(), "idx");
  }

#if 0
  {
    static size_t enter_num = 0;
    if(enter_num++ == 0){
      boost::unordered_map<size_t,double> error_on_each_point;
      for(size_t t = 0; t < hex_surface_points.size(); ++t){
        error_on_each_point[hex_surface_points[t]] =
            norm(hex_node(colon(), hex_surface_points[t])
                 - nearest_node_of_hex_surface(colon(),t));
      }
      matrixd hex_surface_error = zeros<double>(hex_surface.size(2), 1);
      for(size_t fi = 0; fi < hex_surface_error.size(); ++fi){
        for(size_t pi = 0; pi < hex_surface.size(1); ++pi)
          hex_surface_error[fi] += error_on_each_point[hex_surface(pi,fi)];
      }
      hex_surface_error /= 4.0;
      ofstream ofs("surface_distance_error.vtk");
      quad2vtk(ofs, &hex_node[0], hex_node.size(2), &hex_surface[0], hex_surface.size(2));
      cell_data(ofs, &hex_surface_error[0], hex_surface_error.size(),
                "distance");
    }
  }
#endif

  cerr << "# [info] finish nearest point calculation." << endl;
  ///////////////////////////////////////////////

  pt.put("smooth_strategy.desc", "smooth stratrgy,[laplacian/meshless/arap]");
  const string smooth_strategy = pt.get<string>("smooth_strategy.value");


  std::shared_ptr<function_t<double,int32_t> > deform_func(
        build_hex_deform_func(hex, hex_node, hex_surface_point_set,
                              point_adj_points, smooth_strategy));

  std::shared_ptr<function_t<double,int32_t> > boundary_align_func(
        build_boundary_align_func(hex_node, hex_surface_points,
                                  nearest_node_of_hex_surface, q2t));

//  std::shared_ptr<function_t<double,int32_t> > feature_align_func(
//        build_feature_align_func(hex_node, q2t));

  std::shared_ptr<function_t<double,int32_t> > surface_smooth_func(
        build_surface_smooth_func(hex_node, hex_surface_point_set, point_adj_points, q2t));

  std::shared_ptr<function_t<double,int32_t> > tangent_moving_func(
        build_tangent_moving_func(
          hex_node, hex_surface_points,nearest_node_of_hex_surface,
          nearest_node_of_hex_surface_belong_faces,
          nearest_node_of_hex_surface_weights,tet_node, tet_surface,
          tri_face_normal, orfap));

  pt.put("deformation_w.desc","the weight to deform the inner tets");
  pt.put("boundary_align_w.desc","weight to control boundary alignment");
  pt.put("tangent_move_w.desc", "weight of tangent moving");
  pt.put("smooth_surface_w.desc", "smooth surface weight");
  pt.put("feature_align_w.desc", "feature align weight");

  const double ssw = pt.get<double>("smooth_surface_w.value");
  const double dfw = pt.get<double>("deformation_w.value");
  const double baw = pt.get<double>("boundary_align_w.value");
  const double tmw = pt.get<double>("tangent_move_w.value");
  const double faw = pt.get<double>("feature_align_w.value");

  ///// optimization //////
  vec_func_ptr all_funcs(new vec_func);
  if(boundary_align_func.get() && baw > 0)
    all_funcs->push_back(func_ptr(func_ptr(boundary_align_func) *sqrt(fabs(baw))));

  if(deform_func.get() && dfw > 0)
    all_funcs->push_back(
          func_ptr(func_ptr(deform_func) * sqrt(fabs(dfw))));

  if(surface_smooth_func.get() && ssw > 0)
    all_funcs->push_back(
          func_ptr(func_ptr(surface_smooth_func) * sqrt(fabs(ssw))));

  if(tangent_moving_func.get() && tmw > 0)
    all_funcs->push_back(
          func_ptr(func_ptr(tangent_moving_func) * sqrt(fabs(tmw))));

//  if(feature_align_func.get() && faw > 0)
//    all_funcs->push_back(
//          func_ptr(func_ptr(feature_align_func) * sqrt(fabs(faw))));

  if(all_funcs->empty()){
    cerr << "# [error] empty functions." << endl;
    return __LINE__;
  }

  func_ptr func(new_catenated_function<double, int32_t>(all_funcs));
  matrixd residual(func->dim_of_f());
  zjucad::optimize(*func, hex_node, residual, pt);
  return 0;
}
