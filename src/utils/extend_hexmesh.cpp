#include <set>
#include <fstream>
#include <iostream>
//#include <localsolver.h>
#include <zjucad/matrix/io.h>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/util/util.h>
#include "../hexmesh/io.h"
#include "../hexmesh/util.h"
#include "../common/vtk.h"
#include "../hexmesh/hexmesh.h"
using namespace std;
using namespace jtf::hexmesh;
using namespace zjucad::matrix;

int extend_hexmesh(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] extend_hexmesh input_mesh output_mesh extend_length(0-1)"  << endl;
      return __LINE__;
    }

  jtf::hex_mesh hm(argv[1]);

  matrix<size_t> select_faces;
  matrix<double> select_normal;

  select_faces = hm.outside_face_;
  select_normal = hm.outside_face_normal_;

  jtf::mesh::meshes new_hexmesh;

  extend_hexmesh_based_on_faces(select_faces, select_normal, hm.hexmesh_, new_hexmesh, atof(argv[3]));

  ofstream ofs("extend_hex.vtk");
  hex2vtk(ofs, &new_hexmesh.node_[0], new_hexmesh.node_.size(2), &new_hexmesh.mesh_[0], new_hexmesh.mesh_.size(2));

  hex_mesh_write_to_wyz(argv[2], new_hexmesh.mesh_, new_hexmesh.node_);

  return 0;
}


int pick_flatten_singularity(
    const vector<deque<pair<size_t,size_t> > > &singularity_chains,
    const vector<size_t> &chain_degree,
    const jtf::hex_mesh &hm, vector<size_t> &select_flatten_singularity,
    const double threshold_degree)
{
  select_flatten_singularity.clear();
  for(size_t ci = 0; ci < singularity_chains.size(); ++ci){
      if(chain_degree[ci] != 1) continue;
      const deque<pair<size_t,size_t> > & one_chain = singularity_chains[ci];
      size_t segments_above_thresold = 0;
      for(size_t ei = 0; ei < one_chain.size(); ++ei){
          const size_t edge_idx = hm.ea_outside_->get_edge_idx(one_chain[ei].first, one_chain[ei].second);
          const pair<size_t,size_t> & face_pair = hm.ea_outside_->edge2cell_[edge_idx];
          assert(hm.ea_outside_->is_boundary_edge(face_pair) == false);
          // assert normal is normalized
          const double angle = jtf::math::cos2deg(dot(hm.outside_face_normal_(colon(), face_pair.first),
                                                      hm.outside_face_normal_(colon(), face_pair.second)));
          if(angle < threshold_degree) ++segments_above_thresold;
        }
      if(segments_above_thresold*1.0 > 0.5*one_chain.size())
        select_flatten_singularity.push_back(ci);
    }
  return 0;
}

inline int edge_degree2chain_degree(
    const matrix<size_t> & surface_singularity_degree,
    const matrix<size_t> & surface_singularity,
    const vector<deque<pair<size_t,size_t> > > & singularity_chains,
    vector<size_t> & chain_degree)
{
  map<pair<size_t,size_t>,size_t> edge2degree;
  pair<size_t,size_t> one_edge;
  for(size_t ei = 0; ei < surface_singularity.size(2); ++ei){
      one_edge.first = surface_singularity(0,ei);
      one_edge.second = surface_singularity(1,ei);
      if(one_edge.first > one_edge.second) swap(one_edge.first, one_edge.second);
      edge2degree[one_edge] = surface_singularity_degree[ei];
    }
  for(size_t ci = 0; ci < singularity_chains.size(); ++ci){
      const deque<pair<size_t,size_t> > & one_chain = singularity_chains[ci];
      for(const auto & one_edge : one_chain){
          map<pair<size_t,size_t>,size_t>::const_iterator it;
          if(one_edge.first > one_edge.second)
            it = edge2degree.find(make_pair(one_edge.second, one_edge.first));
          else it = edge2degree.find(one_edge);
          if(it == edge2degree.end()){
              cerr << "# [error] can not find edge2degree of " << one_edge.first << " " << one_edge.second << endl;
            }
          if(chain_degree[ci] == 0) chain_degree[ci] = it->second;
          else if(chain_degree[ci] != it->second){
              cerr << "# [error] this chain contains different degree edges " << ci << endl;
              return __LINE__;
            }
        }
    }
  return 0;
}

int extend_hexmesh2(int argc, char *argv[])
{
  if(argc != 4){
      cerr << "# [usage] extend_hexmesh2 input_mesh output_mesh extend_length(0-1)" << endl;
      return __LINE__;
    }

  jtf::hex_mesh hm(argv[1]);
  jtf::hexmesh::hex_singularity_extractor hse(hm);

  matrix<size_t> surface_singularity = hse.get_surface_singularity_edges();
  matrix<size_t> surface_singularity_degree = hse.get_surface_singularity_edges_degree();

  vector<deque<pair<size_t,size_t> > > singularity_chains;
  {
    vector<pair<size_t,size_t> > singularity_edges(surface_singularity.size(2));
    for(size_t ei = 0; ei < surface_singularity.size(2); ++ei)
      singularity_edges[ei] = make_pair(surface_singularity(0,ei), surface_singularity(1,ei));

    jtf::util::extract_chain_from_edges(singularity_edges, singularity_chains);
  }

  vector<size_t> chain_degree(singularity_chains.size(),0);
  jtf::mesh::patch_separater ps(hm.outside_face_);
  ps.separater(singularity_chains);
  cerr << "# [info] finish patch separater" << endl;

  edge_degree2chain_degree(surface_singularity_degree, surface_singularity, singularity_chains, chain_degree);

  vector<size_t> select_flatten_singularity;
  pick_flatten_singularity(singularity_chains, chain_degree, hm, select_flatten_singularity, 45.0);
  cerr << "# [info] pick " << select_flatten_singularity.size() << " singularity chains." << endl;

  {
    ofstream ofs("flatten_singularity.vtk");
    vector<size_t> select_edges;
    for(size_t ci =0; ci < select_flatten_singularity.size(); ++ci){
        const deque<pair<size_t,size_t> > & one_chain = singularity_chains[select_flatten_singularity[ci]];
        for(size_t ei = 0; ei < one_chain.size(); ++ei){
            select_edges.push_back(one_chain[ei].first);
            select_edges.push_back(one_chain[ei].second);
          }
      }
    line2vtk(ofs, &hm.hexmesh_.node_[0], hm.hexmesh_.node_.size(2), &select_edges[0], select_edges.size()/2);
  }

  set<size_t> select_patches;
  vector<bool> visit_patch (ps.get_all_patches().size(),false);
  deque<size_t> pd;
  for(size_t ci = 0; ci < select_flatten_singularity.size(); ++ci){
      set<size_t> patch_idx = ps.get_chain_adj_patches(select_flatten_singularity[ci]);
      select_patches.insert(patch_idx.begin(), patch_idx.end());
    }

  pd.insert(pd.end(), select_patches.begin(), select_patches.end());

  while(!pd.empty()){
      const size_t one_p = pd.front();
      pd.pop_front();
      if(visit_patch[one_p]) continue;
      visit_patch[one_p] = true;
      select_patches.insert(one_p);
      set<size_t> chain_idx = ps.get_chain_of_patch(one_p);
      for(const auto & one_chain_idx : chain_idx){
          if(chain_degree[one_chain_idx] == 3){
              set<size_t> adj_patches = ps.get_chain_adj_patches(one_chain_idx);
              for(const auto & one_adj_p : adj_patches){
                  if(visit_patch[one_adj_p]) continue;
                  pd.push_back(one_adj_p);
                }
            }
        }
    }

  cerr << "# [info] find " << select_patches.size() << "patches to be extended" << endl;

  vector<size_t> faces;
  const vector<vector<size_t> > & patches_face = ps.get_all_patches();
  for(const auto & one_p : select_patches){
      faces.insert(faces.end(), patches_face[one_p].begin(), patches_face[one_p].end());
    }
  itr_matrix<const size_t*> face_m(faces.size(),1, &faces[0]);
  matrix<size_t> select_faces = hm.outside_face_(colon(), face_m);

  matrix<double> select_normal = hm.outside_face_normal_(colon(), face_m);
  jtf::mesh::meshes new_hexmesh;
  extend_hexmesh_based_on_faces(select_faces, select_normal, hm.hexmesh_, new_hexmesh, atof(argv[3]));

  jtf::hexmesh::hex_mesh_write_to_wyz(argv[2], new_hexmesh.mesh_, new_hexmesh.node_);
  return 0;
}

double get_theta_angle(const pair<size_t,size_t> & one_edge, const jtf::hex_mesh &hm)
{
  const size_t edge_idx = hm.ea_outside_->get_edge_idx(one_edge.first, one_edge.second);
  assert(edge_idx != -1);
  const pair<size_t,size_t> & cell_pair = hm.ea_outside_->edge2cell_[edge_idx];
  const matrix<double> n1 = -1*hm.outside_face_normal_(colon(), cell_pair.first);
  const matrix<double> n2 = hm.outside_face_normal_(colon(), cell_pair.second);
  double theta = jtf::math::safe_acos(dot(n1,n2));
  matrix<double> edge_dir = hm.hexmesh_.node_(colon(), one_edge.second) - hm.hexmesh_.node_(colon(), one_edge.first);
  const double a = dot(cross(n2,n1),edge_dir);
  //  matrix<double> n2crossn1 = cross(n2,n1);
  //  cerr << hm.outside_face_(colon(), cell_pair.first) << endl;
  //  cerr << n1 << n2 << endl;
  //  cerr << n2crossn1 << endl;
  //  cerr << edge_dir/norm(edge_dir) << endl;
  if(a > 0) theta *= -1;
  if(theta < 0) theta+= 2*jtf::math::My_PI();
  return theta;
}


int opt_to_select_faces(
    const jtf::mesh::patch_separater & ps, const vector<size_t> &chain_degree,
    const vector<deque<pair<size_t,size_t> > > & singularity_chains,
    const jtf::hex_mesh &hm, vector<size_t> &select_faces_idx)
{
    throw std::logic_error("not finished.");
//  using namespace localsolver;

//  LocalSolver localsolver;
//  LSModel model = localsolver.getModel();

//  const std::vector<std::vector<size_t> > & all_patches = ps.get_all_patches();
//  // 0-1 decisions
//  vector<LSExpression> x(all_patches.size(), model.createExpression(O_Bool));

//  LSExpression obj = model.createExpression(O_Sum);
//  for(size_t ci = 0; ci < singularity_chains.size(); ++ci){
//      const size_t degree = chain_degree[ci];
//      for(const auto & one_edge : singularity_chains[ci]){
//          const double w_e = 1.0;
//          const double angle = 2*jtf::math::My_PI()-get_theta_angle(one_edge, hm);
//          const set<size_t> & adj_patches = ps.get_edge_adj_patches(one_edge);
//          if(adj_patches.size() != 2) continue;
//          auto it = adj_patches.begin();
//          const size_t bi_idx = *(it++);
//          const size_t bj_idx = *it;

//          const double piD_2 = 0.5*jtf::math::My_PI()*degree;
//          LSExpression one_item = model.createExpression(O_Sum);
//          // E+= angle-pid/2
//          one_item.addOperand(angle-piD_2);
//          // E+= (xi or xj)*pid/2
//          one_item.addOperand(model.createExpression(O_Prod, model.createExpression(O_Or, x[bi_idx], x[bj_idx]),piD_2));
//          // E-= (xi+xj)*pi/2
//          one_item.addOperand(model.createExpression(O_Prod, model.createExpression(O_Sum, x[bi_idx], x[bj_idx]), -0.5*jtf::math::My_PI()));
//          LSExpression square_one_item = model.createExpression(O_Prod, model.createExpression(O_Pow, one_item, static_cast<lsint>(2)), w_e);
//          obj.addOperand(square_one_item);
//        }
//    }

//  // maximize knapsackValue;
//  model.addObjective(obj, OD_Minimize);

//  // close model, then solve
//  model.close();
//  LSPhase phase = localsolver.createPhase();
//  phase.setTimeLimit(1);
//  localsolver.solve();

////  for(size_t i = 0; i < x.size()-2; ++i)
////    x[i].setIntValue(lsint(1));

////  x[4].setIntValue(lsint(0));
////  x[5].setIntValue(lsint(0));

//  for(size_t pi = 0; pi < x.size(); ++pi){
//      cerr << x[pi].getIntValue() << endl;
//    }

//  cerr << obj.getDoubleValue() << endl;
}

int extend_hexmesh3(int argc, char *argv[])
{
  if(argc != 3){
      cerr << "# [usage] extend_hexmesh3 input_hex out_hex" << endl;
      return __LINE__;
    }

  try{
    jtf::hex_mesh hm(argv[1]);

    jtf::hexmesh::hex_singularity_extractor hse(hm);

    matrix<size_t> surface_singularity = hse.get_surface_singularity_edges();
    matrix<size_t> surface_singularity_degree = hse.get_surface_singularity_edges_degree();

    vector<deque<pair<size_t,size_t> > > singularity_chains;
    {
      vector<pair<size_t,size_t> > singularity_edges(surface_singularity.size(2));
      for(size_t ei = 0; ei < surface_singularity.size(2); ++ei)
        singularity_edges[ei] = make_pair(surface_singularity(0,ei), surface_singularity(1,ei));

      jtf::util::extract_chain_from_edges(singularity_edges, singularity_chains);
    }

    vector<size_t> chain_degree(singularity_chains.size(),0);
    jtf::mesh::patch_separater ps(hm.outside_face_);
    ps.separater(singularity_chains);

    {
      const std::vector<std::vector<size_t> > & all_patches = ps.get_all_patches();
      matrix<size_t> face_type(hm.outside_face_.size(2),1);
      for(size_t pi = 0; pi < all_patches.size(); ++pi){
          const vector<size_t> & one_patch = all_patches[pi];
          for(size_t fi = 0; fi < one_patch.size(); ++fi){
              face_type[one_patch[fi]] = pi;
            }
        }
      ofstream ofs("surface_patch.vtk");
      quad2vtk(ofs, &hm.outside_face_[0], hm.outside_face_.size(2),
          &hm.hexmesh_.node_[0], hm.hexmesh_.node_.size(2));
      cell_data(ofs, &face_type[0], face_type.size(), "patch_idx");
    }
    edge_degree2chain_degree(surface_singularity_degree, surface_singularity,
                             singularity_chains,  chain_degree);

    vector<size_t> select_faces_idx;
    opt_to_select_faces(ps, chain_degree, singularity_chains, hm, select_faces_idx);

    itr_matrix<const size_t*> select_faces_idx_m(select_faces_idx.size(), 1, &select_faces_idx[0]);
    matrix<size_t> select_faces = hm.outside_face_(colon(), select_faces_idx_m);
    matrix<double> select_normal = hm.outside_face_normal_(colon(), select_faces_idx_m);
    jtf::mesh::meshes new_hexmesh;
    extend_hexmesh_based_on_faces(select_faces, select_normal, hm.hexmesh_, new_hexmesh, 0.01);
    jtf::hexmesh::hex_mesh_write_to_wyz(argv[2], new_hexmesh.mesh_, new_hexmesh.node_);
  }catch(std::exception &e){
    cerr << e.what() << endl;
  }

  return 0;
}
