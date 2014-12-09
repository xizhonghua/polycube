#include "tetmesh_refine.h"
#include "zjucad/optimizer/optimizer.h"
#include "../common/visualize_tool.h"
using namespace std;
using namespace boost;
using namespace zjucad::matrix;
using namespace hj::function;
using namespace boost::property_tree;

typedef function_t<double, int32_t> function_di;
typedef std::shared_ptr<const function_di> func_ptr;
typedef vector<func_ptr> vec_func;
typedef std::shared_ptr<vector<func_ptr> > vec_func_ptr;

int design_singularities_coord(const matrixd &node,
                               const vector<deque<pair<size_t,size_t> > > &singularities,
                               map<size_t,matrixd > &design_singularities_coord,
                               ptree &pt)
{
  typedef map<size_t,matrixd >::const_iterator msmcit;
  pt.put("singularity_refine.desc","the method to get desired singularities,[laplace]");
  const string method = pt.get<string>("singularity_refine.value","laplace");
  if(method == "laplace"){
    for(size_t t = 0; t < singularities.size(); ++t){
      const deque<pair<size_t,size_t> > & chain = singularities[t];
      if(chain.size() > 1){ // only chain.size() > 1, need to handle this
        if(chain.front().first == chain.back().second){ // it's a loop
          for(size_t i = 0; i < chain.size(); ++i){
            msmcit v0 = design_singularities_coord.find(chain[i].second);
            msmcit v1 = design_singularities_coord.find(chain[(chain.size() + i - 1)% chain.size()].first);

            design_singularities_coord[chain[i].first] =
                0.5 * (((v0 == design_singularities_coord.end())?node(colon(),chain[i].second):v0->second)
                + ((v1 == design_singularities_coord.end())?node(colon(),chain[(chain.size() + i - 1)% chain.size()].first):v1->second));
          }
        }else{ // open chain
          assert(chain.size() > 1);
          for(size_t i = 1; i < chain.size(); ++i){
            msmcit v0 = design_singularities_coord.find(chain[i].second);
            msmcit v1 = design_singularities_coord.find(chain[(chain.size() + i - 1)% chain.size()].first);

            design_singularities_coord[chain[i].first] =
                0.5 * (((v0 == design_singularities_coord.end())?node(colon(),chain[i].second):v0->second)
                + ((v1 == design_singularities_coord.end())?node(colon(),chain[(chain.size() + i - 1)% chain.size()].first):v1->second));
          }
        }// end else
      }// end size > 1
    } // end one chain
  }// end method == laplace
  return 0;
}

int relax_singularities(const matrixst & tet,
                        matrixd &node,
                        const vector<deque<pair<size_t,size_t> > > &singularities,
                        const jtf::mesh::face2tet_adjacent &fa,
                        const matrixst &outside_face,
                        ptree &pt)
{
  map<size_t,matrixd > desired_vertex_coord;

  design_singularities_coord(node,singularities,desired_vertex_coord,pt);
  cerr << "# [info] finished smoothed singularities coord calculation." << endl;
  matrixd temp_node = node;
  { // dump out the desired singularity chain
    for(map<size_t,matrixd >::const_iterator mcit = desired_vertex_coord.begin();
        mcit != desired_vertex_coord.end(); ++mcit){
      temp_node(colon(),mcit->first) = mcit->second;
    }
    dump_singularity_to_vtk("desired_singularity.vtk",temp_node,singularities);
  }

  {// construct the relaxing functions
    vec_func_ptr all_funcs(new vec_func);
    construct_singularity_relaxing_func(all_funcs,desired_vertex_coord,tet,
                                        node,outside_face,pt);
    func_ptr func(new_catenated_function<double, int32_t>(all_funcs));
    matrixd residual(func->dim_of_f());
    cerr << "# [info] start global relaxing." << endl;
   // node += 0.1 * rand<double>(3,node.size(2));
    cerr << "# [info] init node value " << node(colon(),colon(0,4)) << endl;
    zjucad::optimize(*func, temp_node, residual, pt);
    node = temp_node;
    cerr << "# [info] final node value " << node(colon(),colon(0,4)) << endl;
  }
  return 0;
}
