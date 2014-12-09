#include "tetmesh_refine.h"
#include "relax_function.h"
#include <vector>
#include <fstream>
#include "../common/util.h"
#include "../common/transition_type.h"
#include "../hex_param/topology_analysis.h"
#include "../common/visualize_tool.h"
#include "../common/vtk.h"
#include <jtflib/util/container_operation.h>
#include <jtflib/mesh/mesh.h>
#include "zjucad/optimizer/optimizer.h"
using namespace std;
using namespace zjucad::matrix;
using namespace hj::function;

int arap_with_degenerated_points_driven(
    vec_func_ptr func,
    const matrixst & tet,
    const std::vector<std::vector<size_t> > & degenerated_points,
    matrixd & node,
    boost::property_tree::ptree &pt)
{
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  matrixst outside_face;
  get_outside_face(*fa,outside_face);

  // there are three mainly term in function

  vec_func_ptr arap_deformation_func(new vec_func);
#if 1
  { // arap deformation
    for(size_t t = 0; t < tet.size(2); ++t){
      unique_ptr<arap_func> ddf(new arap_func(tet(colon(),t),node,t));
      arap_deformation_func->push_back(func_ptr(ddf.release()));
      //unique_ptr<deformation_dest_func> ddf(new deformation_dest_func(tet(colon(),t),
      //                                                               node,t));
      //arap_deformation_func->push_back(func_ptr(ddf.release()));
    }
  }
#endif

  vec_func_ptr boundary_align_func(new vec_func);
  { // bonundary align
    set<size_t> boundary_points(outside_face.begin(),outside_face.end());
    for(set<size_t>::const_iterator scit = boundary_points.begin();
        scit != boundary_points.end(); ++scit){
      unique_ptr<point_fix_function> baf(
            new point_fix_function(node.size(2),*scit,node(colon(),*scit)));
      boundary_align_func->push_back(func_ptr(baf.release()));
    }
  }

  vec_func_ptr degenerated_points_concentrate_func(new vec_func);
  { // degenerated points concentrated
    matrixd average_node = zeros<double>(3,1);
    for(size_t t = 0; t < degenerated_points.size(); ++t){
      const vector<size_t> & one_group = degenerated_points[t];
      assert(one_group.size() > 1);
      for(size_t i = 0; i < one_group.size(); ++i){
        for(size_t j = i + 1; j < one_group.size(); ++j){
          unique_ptr<equal_function> ef(
                new equal_function(node.size(2),one_group[i],one_group[j]));
          degenerated_points_concentrate_func->push_back(func_ptr(ef.release()));
        }
      }
    }
  }

  pt.put("degenerated_cons_w.desc","the weight to concentrated degenerated points");
  const double dcw = pt.get<double>("degenerated_cons_w.value");
  pt.put("boundary_align_w.desc","weight to control boundary alignment");
  const double baw = pt.get<double>("boundary_align_w.value");

  cerr << "# [info] degenerated points concentrated func num: "
       << degenerated_points_concentrate_func->size() << endl;
  cerr << "# [info] boundary align func num: "
       << boundary_align_func->size() << endl;
  cerr << "# [info] arap deformation func num: "
       << arap_deformation_func->size() << endl;

  const double bounding_sphere_radius =
      calc_bounding_sphere_size(node)/2.0;
  if(!degenerated_points_concentrate_func->empty() && dcw > 0)
    func->push_back(
          func_ptr(func_ptr(new_catenated_function<double,int32_t>(
                              degenerated_points_concentrate_func))
                   *sqrt(fabs(dcw)/bounding_sphere_radius)));

  if(!boundary_align_func->empty() && baw > 0)
    func->push_back(
          func_ptr(func_ptr(new_catenated_function<double,int32_t>(
                              boundary_align_func))
                   *sqrt(fabs(baw))));

  if(!arap_deformation_func->empty())
    func->push_back(
          func_ptr(func_ptr(new_catenated_function<double,int32_t>(
                              arap_deformation_func))));
  return 0;
}


int arap_to_concentrate_points(
    const matrixst & tet,
    matrixd & node,
    boost::property_tree::ptree &pt)
{
  vector<vector<size_t> > degenerated_points;
  pt.put("strategy.desc",
         "strategy to choose handle original black edge or integer conflicts"
         "[black/integer]");
  const string strategy = pt.get<string>("strategy.value");
  if(strategy == "integer")
  {
    pt.put("singularity_degenerated_points.desc",
           "file store the singularity degenerated points");
    if(load_singularity_degenerated_points(
         pt.get<string>("singularity_degenerated_points.value").c_str(),
         degenerated_points))
      return __LINE__;
  }else if(strategy == "black"){
    pt.put("singularity_chain.desc","singularity chain file");
    vector<deque<pair<size_t,size_t> > > singularity_chain;
    vector<deque<size_t> > singularity_type;
    if(load_singularity_chain_new(
         pt.get<string>("singularity_chain.value").c_str(),
         singularity_chain,singularity_type))
      return __LINE__;

    for(size_t t = 0; t < singularity_chain.size(); ++t){
      const deque<pair<size_t,size_t> > & one_chain = singularity_chain[t];
      vector<size_t> degenerated_points_vec;
      if(is_black_line_new(singularity_type[t].front()))
      {
        for(size_t i = 0; i < one_chain.size(); ++i){
          degenerated_points_vec.push_back(one_chain[i].first);
        }
        if(one_chain.back().second != one_chain.front().first)
          degenerated_points_vec.push_back(one_chain.back().second);
        degenerated_points.push_back(degenerated_points_vec);
        continue;
      }
    }
    cerr << "# [info] degenerated points group: " << degenerated_points.size() << endl;
    for(size_t t = 0;t < degenerated_points.size(); ++t){
      cerr << "# [info] gropu " << t << ": ";
      copy(degenerated_points[t].begin(),degenerated_points[t].end(),
           ostream_iterator<size_t>(cerr, " "));
      cerr << endl;
    }
  }else{
    cerr << "# [error] can not handle such strategy." << endl;
    return __LINE__;
  }

  vec_func_ptr all_funcs(new vec_func);
  arap_with_degenerated_points_driven(all_funcs,tet,degenerated_points,node,pt);
  func_ptr func(new_catenated_function<double, int32_t>(all_funcs));
  matrixd residual(func->dim_of_f());

#if 0
  {
    //adjust the initial node
    matrixd average_node = zeros<double>(3,1);
    for(size_t t = 0; t < degenerated_points.size(); ++t){
      const vector<size_t> & one_group = degenerated_points[t];
      assert(one_group.size() > 1);
      average_node = zeros<double>(3,1);
      for(size_t i = 0; i < one_group.size(); ++i)
        average_node += node(colon(),one_group[i]);
      average_node /= one_group.size();
      for(size_t i  = 0; i < one_group.size(); ++i){
        node(colon(),one_group[i]) = average_node;
      }
    }
    jtf::mesh::tet_mesh_write_to_zjumat("init_fandisk.tet",
                             &node,&tet);
  }
#endif

#ifdef fun // just for fun
  for(size_t t = 0; t < 100; ++t){
    ostringstream os;
    os << t;
    ofstream ofs((os.str() + "_tet.vtk").c_str());
    tet2vtk(ofs,&node[0],node.size(2),&tet[0],tet.size(2));
    pt.put("iter.value",50);
#endif
    zjucad::optimize(*func, node, residual, pt);
#ifdef fun
  }
#endif

  cerr << "# [info] final node value " << node(colon(),colon(0,4)) << endl;
  return 0;
}


int tet_mesh_inflation_for_ball(
    const matrixst & tet,
    matrixd & node,
    boost::property_tree::ptree &pt)
{
  const double radius = calc_bounding_sphere_size(node)/2.0;
  matrixd center_node;
  cal_average_node(node, center_node);
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));

  matrixst outside_face;
  get_outside_face(*fa,outside_face);


  vec_func_ptr arap_deformation_func(new vec_func);
  { // arap deformation
    for(size_t t = 0; t < tet.size(2); ++t){
      unique_ptr<arap_func> ddf(new arap_func(tet(colon(),t),node,t));
      arap_deformation_func->push_back(func_ptr(ddf.release()));
    }
  }

  vec_func_ptr boundary_to_sphere_func(new vec_func);
  { // bonundary align
    set<size_t> boundary_points(outside_face.begin(),outside_face.end());
    for(set<size_t>::const_iterator scit = boundary_points.begin();
        scit != boundary_points.end(); ++scit){
      unique_ptr<point_to_sphere_function> baf(
            new point_to_sphere_function(node.size(2),*scit,radius,center_node));
      boundary_to_sphere_func->push_back(func_ptr(baf.release()));
    }
  }

  pt.put("boundary_to_sphere_w.desc","weight to control boundary to sphere");
  const double btsw = pt.get<double>("boundary_to_sphere_w.value");

  cerr << "# [info] boundary align func num: "
       << boundary_to_sphere_func->size() << endl;
  cerr << "# [info] arap deformation func num: "
       << arap_deformation_func->size() << endl;

  vec_func_ptr func(new vec_func);
  if(!boundary_to_sphere_func->empty() && btsw > 0)
    func->push_back(
          func_ptr(func_ptr(new_catenated_function<double,int32_t>(
                              boundary_to_sphere_func))
                   *sqrt(fabs(btsw))));

  if(!arap_deformation_func->empty())
    func->push_back(
          func_ptr(func_ptr(new_catenated_function<double,int32_t>(
                              arap_deformation_func))));

  func_ptr all_func(new_catenated_function<double, int32_t>(func));
  matrixd residual(all_func->dim_of_f());
  zjucad::optimize(*all_func, node, residual, pt);

  return 0;
}
