#include "tetmesh_refine.h"
#include "relax_function.h"
#include <boost/scoped_ptr.hpp>
#include <set>
#include "../common/util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace boost::property_tree;
using namespace boost;
using namespace hj::function;

int construct_singularity_relaxing_func(
    vec_func_ptr func,
    const std::map<size_t,matrixd > &desired_vertex_coord,
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    ptree &pt)
{
  vec_func_ptr singularity_align_func(new vec_func);
  typedef map<size_t,matrixd >::const_iterator mmcit;
#if 1
  {// add singularity_alignment functions
    singularity_align_func->reserve(desired_vertex_coord.size());
    for(mmcit mcit = desired_vertex_coord.begin();
        mcit != desired_vertex_coord.end(); ++mcit){
      unique_ptr<point_fix_function> saf(
            new point_fix_function(node.size(2),
                                   mcit->first,
                                   mcit->second));
      singularity_align_func->push_back(func_ptr(saf.release()));
    }
  }
#endif

  vec_func_ptr boundary_align_func(new vec_func);
#if 1
  {// add boundary alignment, we hope the boundary point still remain on surface
    std::set<size_t> outside_point(outside_face.begin(),outside_face.end());
    boundary_align_func->reserve(outside_point.size());
    for(set<size_t>::const_iterator scit = outside_point.begin();
        scit != outside_point.end(); ++scit){
      assert(*scit < node.size(2));
      unique_ptr<point_fix_function> baf(new point_fix_function(
                                         node.size(2),
                                         *scit,
                                         node(colon(),*scit)));
      boundary_align_func->push_back(func_ptr(baf.release()));
    }
  }
#endif

  vec_func_ptr deformation_func(new vec_func);
  deformation_func->reserve(tet.size(2));

#if 1 // if original model is scaled as A, the this term will be scaled as A^2,
  // so it's necessary to choose a good weight, for example divide with bounding box length
  {
    for(size_t t = 0; t < tet.size(2); ++t){
      unique_ptr<deformation_dest_func> ddf(new deformation_dest_func(
                                            tet(colon(),t),
                                            node,t));
      deformation_func->push_back(func_ptr(ddf.release()));
    }
  }
#endif

  pt.put("singularity_align_w.desc","the weight to control singularity alignment");
  const double saw = pt.get<double>("singularity_align_w.value");
  pt.put("boundary_align_w.desc","weight to control boundary alignment");
  const double baw = pt.get<double>("boundary_align_w.value");

  cerr << "# [info] singularity align func num: " << singularity_align_func->size() << endl;
  cerr << "# [info] boundary align func num: " << boundary_align_func->size() << endl;
  cerr << "# [info] deformation destorition func num: " << deformation_func->size() << endl;


  const double bounding_sphere_radius = calc_bounding_sphere_size(node)/2.0;
  if(!singularity_align_func->empty() && saw > 0)
    func->push_back(func_ptr(func_ptr(new_catenated_function<double,int32_t>(singularity_align_func))
                             *sqrt(fabs(saw)/bounding_sphere_radius)));

  if(!boundary_align_func->empty() && baw > 0)
    func->push_back(func_ptr(func_ptr(new_catenated_function<double,int32_t>(boundary_align_func))
                             *sqrt(fabs(baw)/bounding_sphere_radius)));

  if(!deformation_func->empty())
    func->push_back(func_ptr(func_ptr(new_catenated_function<double,int32_t>(deformation_func))
                             *sqrt(1.0/(bounding_sphere_radius*bounding_sphere_radius))));

  return 0;
}

