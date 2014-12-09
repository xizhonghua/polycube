#include <zjucad/ptree/ptree.h>
#include <jtflib/util/util.h>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/io.h>
#include <jtflib/optimizer/optimizer.h>
#include <boost/unordered_map.hpp>
#include <hjlib/math_func/operation.h>
#include <hjlib/math_func/math_func.h>

#include <string>
#include "../hexmesh/io.h"
#include "../common/vtk.h"
#include "../tetmesh/tetmesh.h"
#include "../hex_param/io.h"
#include "../common/visualize_tool.h"
#include "../common/util.h"
#include "../equation_graph/equation_graph.h"
#include "../mesh_func/arap.h"
#include "../mesh_func/variable_fix.h"

#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include "../common/K_near_points.h"

using namespace std;
using namespace zjucad::matrix;
using boost::property_tree::ptree;

class meta_mesh_extractor
{
public:
  meta_mesh_extractor(const jtf::tet_mesh &tm,
                      const boost::unordered_map<size_t,size_t> & surface_type)
    :tm_(tm), surface_type_(surface_type){}
  void run(boost::property_tree::ptree & pt){
    init();
    const string need_coarsening = pt.get<string>("input/need_coarsening.value","Y");
    if(need_coarsening == "y" || need_coarsening == "Y" ||
       need_coarsening == "yes" || need_coarsening == "YES"){
        coarsening();
        update_patch_list();
      }
    extract();
  }
  void get_meta_mesh(jtf::mesh::meshes & meta_mesh){ meta_mesh = meta_mesh_;}
protected:
  void init(){
    tri_type_ = tm_.outside_face_idx_;

    for(size_t fi = 0; fi < tri_type_.size(); ++fi){
        auto it = surface_type_.find(tri_type_[fi]);
        assert(it != surface_type_.end());
        tri_type_[fi] = it->second;
      }

    vector<deque<pair<size_t,size_t> > > chain;;
    ps_.reset(new jtf::mesh::patch_separater(tm_.outside_face_));
    ps_->separater(tri_type_, chain);
    patches_ = ps_->get_all_patches();
    patch_ea_.resize(patches_.size());
    update_patch_list();
  }

  void update_patch_list(){
    sorted_list_.clear();
    sorted_list_.resize(3);
    set<size_t> patch_points;
    for(size_t pi = 0; pi < patches_.size(); ++pi){
        const size_t one_type = tri_type_[patches_[pi].front()];
        double v = tm_.tetmesh_.node_(one_type, tm_.outside_face_(0,patches_[pi].front()));
        patch_points.clear();
        itr_matrix<const size_t*> pfm(patches_[pi].size(),1, &patches_[pi][0]);
        matrix<size_t> select_faces = tm_.outside_face_(colon(), pfm);
        patch_points.insert(select_faces.begin(), select_faces.end());
        for(const auto & one_p : patch_points) tm_.tetmesh_.node_(one_type, one_p) = v;
        sorted_list_[one_type].push_back(make_pair(v,pi));
      }

    patch_bound_.resize(patches_.size());
    for(size_t pi = 0; pi < patches_.size(); ++pi){
        patch_bound_[pi].resize(3,2);
        cal_patch_xyz_bound(patches_[pi],patch_bound_[pi]);
      }

    sorted_list_.resize(3);
    cal_xyz_sorted_list();

    average_list_interval_.resize(3,1);
    average_list_interval_ *= 0;
    for(size_t di = 0; di < sorted_list_.size(); ++di){
        const vector<pair<double,size_t> >  & one_list = sorted_list_[di];
        for(size_t pi = 0; pi != one_list.size()-1; ++pi){
            average_list_interval_[di] += fabs(one_list[pi+1].first - one_list[pi].first);
          }
        average_list_interval_[di] /= one_list.size()-1;
      }
  }

  void coarsening(){
    equation_graph eg(tm_.tetmesh_.node_.size());
    eg.bind_node(&tm_.tetmesh_.node_);

    matrix<size_t> cut_tet = tm_.tetmesh_.mesh_;
    matrix<size_t> cut_tet2tet(max(cut_tet)+1);
    cut_tet2tet(cut_tet) = tm_.tetmesh_.mesh_(colon());

    vector<pair<size_t,size_t> > edges;
    get_edges_from_tet(edges);

    bool is_restricted_type = true;
    boost::unordered_map<pair<size_t,size_t>,size_t> fake_inner_type;

    unique_ptr<transition_elimination> tre(
          transition_elimination::create(
            tm_.tetmesh_.mesh_, cut_tet, cut_tet2tet, tm_.tetmesh_.node_,
            fake_inner_type, surface_type_,is_restricted_type ));
    if(!tre.get()){
        cerr << "# [error] can not get transition_elimination." << endl;
        return ;
      }

    vector<group<size_t> > node_group = tre->out();
    vector<size_t> integer_group_idx = tre->get_integer_group_idx();
    vector<vector<vector<size_t> > > node_group_vec(3);

    sort_node_groups(node_group, integer_group_idx, node_group_vec);

    for(size_t di = 0; di < node_group_vec.size(); ++di){
        vector<vector<size_t> >  one_group_dim_bkp = node_group_vec[di];
        vector<vector<size_t> > &one_group_dim = node_group_vec[di];
        for(size_t gi = 0; gi != one_group_dim.size()-1; ++gi){
            one_group_dim[gi].insert(one_group_dim[gi].end(), one_group_dim[gi+1].begin(), one_group_dim[gi+1].end());
            one_group_dim.erase(one_group_dim.begin()+gi+1, one_group_dim.begin()+gi+2);
            eg.assemble_equation_graph_generally(node_group_vec, edges);

            bool rtn = eg.check_valid(false);
            if(rtn == false){
                one_group_dim = one_group_dim_bkp;
              }else{
                one_group_dim_bkp = one_group_dim;
                --gi;
              }
            eg.reset(tm_.tetmesh_.node_.size());
          }

      }

    {
      matrix<int> node_type = ones<size_t>(3,tm_.tetmesh_.node_.size(2))*-1;
      size_t accumulate_groups_numbers[] =
      {0, node_group_vec[0].size(), node_group_vec[0].size() + node_group_vec[1].size()};
      for(size_t di = 0; di < node_group_vec.size(); ++di){
          const vector<vector<size_t> > & one_dim_groups = node_group_vec[di];
          for(size_t gi = 0; gi < one_dim_groups.size(); ++gi){
              const vector<size_t> & one_group = one_dim_groups[gi];
              for(size_t vi = 0; vi < one_group.size(); ++vi){
                  const size_t idx = one_group[vi]/3;
                  node_type(di, idx) = gi + accumulate_groups_numbers[di];
                }
            }
        }
      matrix<size_t> tri_type(tm_.outside_face_.size(2),1);
      matrix<int> one_face_delta1(3,1), one_face_delta2(3,1);

      for(size_t fi = 0; fi < tm_.outside_face_.size(2); ++fi){

          one_face_delta1 = 2*node_type(colon(), tm_.outside_face_(0,fi)) -
              (node_type(colon(), tm_.outside_face_(1,fi)) + node_type(colon(), tm_.outside_face_(2,fi)));
          one_face_delta2 = 2*node_type(colon(), tm_.outside_face_(1,fi)) -
              (node_type(colon(), tm_.outside_face_(0,fi)) + node_type(colon(), tm_.outside_face_(2,fi)));

          for(int i = 0 ; i < 3; ++i){
              if(one_face_delta1[i] == 0 && one_face_delta2[i] == 0 && node_type(i, tm_.outside_face_(0,fi)) != -1){
                  tri_type[fi] = node_type(i, tm_.outside_face_(0,fi));
                  break;
                }
            }
        }
      ofstream ofs("surface_group.vtk");
      tri2vtk(ofs, &tm_.tetmesh_.node_[0], tm_.tetmesh_.node_.size(2), &tm_.outside_face_[0], tm_.outside_face_.size(2));
      cell_data(ofs, &tri_type[0], tri_type.size(), "group");
    }

    optimize_tet(node_group_vec);
  }

  void extract(){
    vector<pair<double,size_t> > & z_list = sorted_list_[2];

    vector<size_t> raw_hex_mesh;
    vector<double> raw_hex_node;
    for(size_t zi = 0; zi != z_list.size()-1; ++zi){
        if(fabs(z_list[zi].first - z_list[zi+1].first) < 1e-6) continue;
        extrac_hex_between_z_layers(z_list[zi], z_list[zi+1], raw_hex_mesh, raw_hex_node);
      }

    itr_matrix<const size_t*> rm(8, raw_hex_mesh.size()/8, &raw_hex_mesh[0]);
    itr_matrix<const double*> rn(3, raw_hex_node.size()/3, &raw_hex_node[0]);
    meta_mesh_.mesh_ = rm;
    meta_mesh_.node_ = rn;

    cerr << "# [info] meta mesh node " << meta_mesh_.node_.size(2) << endl;
    remove_duplicated_node(meta_mesh_.mesh_, meta_mesh_.node_,1e-4);
    cerr << "# [info] meta mesh node " << meta_mesh_.node_.size(2) << endl;
  }

  int extrac_hex_between_z_layers(
      const pair<double,size_t>  &low_zyz_layer,
      const pair<double, size_t> &high_zyz_layer,
      vector<size_t> &raw_hex_mesh,
      vector<double> &raw_hex_node)
  {
    assert(low_zyz_layer.first < high_zyz_layer.first);

    const vector<pair<double,size_t> > & x_list = sorted_list_[0];
    const vector<pair<double,size_t> > & y_list = sorted_list_[1];

    static matrix<double> center_xyz(3,1);
    for(size_t yi = 0; yi != y_list.size()-1; ++yi){
        if(fabs(y_list[yi].first-y_list[yi+1].first) < 0.01*average_list_interval_[1]) continue;
        for(size_t xi = 0; xi != x_list.size()-1; ++xi){
            if(fabs(x_list[xi].first-x_list[xi+1].first) < 0.01*average_list_interval_[0]) continue;
            center_xyz[0] = (x_list[xi].first+x_list[xi+1].first)*0.5;
            center_xyz[1] = (y_list[yi].first+y_list[yi+1].first)*0.5;
            center_xyz[2] = (low_zyz_layer.first+high_zyz_layer.first)*0.5;

            if(is_point_inside(center_xyz)){
                add_one_hex(
                      x_list[xi].first, x_list[xi+1].first,
                    y_list[yi].first, y_list[yi+1].first,
                    low_zyz_layer.first, high_zyz_layer.first,
                    raw_hex_mesh, raw_hex_node);
              }
          }
      }
    return 0;
  }

  bool is_point_inside(const matrix<double> & point_xyz)
  {
    const vector<pair<double,size_t> > & x_list = sorted_list_[0];

    size_t xi = 0;
    for(; xi < x_list.size(); ++xi){
        if(x_list[xi].first > point_xyz[0]) break;
      }

    int left_x_walls = 1;
    const double y0 = patch_bound_[x_list[xi].second](1,0);
    const double y1 = patch_bound_[x_list[xi].second](1,1);
    const double z0 = patch_bound_[x_list[xi].second](2,0);
    const double z1 = patch_bound_[x_list[xi].second](2,1);
    if(y0 > point_xyz[1] || y1 < point_xyz[1]
       || z0 > point_xyz[2] || z1 < point_xyz[2])
      --left_x_walls;
    else{
        if(!is_point_inside_patch(patch_bound_[x_list[xi].second](0,0),
                                  point_xyz[1], point_xyz[2], 0,
                                  x_list[xi].second))
          --left_x_walls;
      }
    for(size_t c = xi+1; c < x_list.size(); ++c){
        //if(fabs(x_list[c].first - x_list[c-1].first) < 1e-6) continue;
        ++left_x_walls;
        const double yy0 = patch_bound_[x_list[c].second](1,0);
        const double yy1 = patch_bound_[x_list[c].second](1,1);
        const double zz0 = patch_bound_[x_list[c].second](2,0);
        const double zz1 = patch_bound_[x_list[c].second](2,1);
        if(yy0 > point_xyz[1] || yy1 < point_xyz[1]
           || zz0 > point_xyz[2] || zz1 < point_xyz[2])
          --left_x_walls;
        else{
            if(!is_point_inside_patch(patch_bound_[x_list[c].second](0,0),
                                      point_xyz[1], point_xyz[2],0,
                                      x_list[c].second))
              --left_x_walls;
          }
      }
    if(left_x_walls < 0)
      throw std::logic_error("strange, can not tell whether this point is inside or outside");
    if(left_x_walls % 2 == 0) return false;
    return true;
  }

  void add_one_hex(const double x0, const double x1,
                   const double y0, const double y1,
                   const double z0, const double z1,
                   vector<size_t> & mesh, vector<double> & node)
  {
    double hex_node [] = {
      x0, y0, z0, // p0
      x1, y0, z0, // p1
      x0, y1, z0, // p2
      x1, y1, z0, // p3
      x0, y0, z1, // p4
      x1, y0, z1, // p5
      x0, y1, z1, // p6
      x1, y1, z1 // p7
    };

    size_t node_idx = node.size()/3;
    node.insert(node.end(), hex_node, hex_node+8*3);
    for(size_t i = 0; i < 8; ++i)  mesh.push_back(node_idx+i);
  }

  bool is_point_inside_patch(const double x, const double y, const double z,
                             const size_t same_axis, const size_t patch_idx)
  {
    double point_xyz[] = {x,y,z};
    itr_matrix<const size_t*> pm(patches_[patch_idx].size(),1, &patches_[patch_idx][0]);
    if(!patch_ea_[patch_idx].get())
      patch_ea_[patch_idx].reset(jtf::mesh::edge2cell_adjacent::create(tm_.outside_face_(colon(), pm)));

    matrix<size_t> boundary_edges;
    jtf::mesh::get_boundary_edge(*patch_ea_[patch_idx], boundary_edges);

    size_t hit_edges = 0;
    for(size_t ei = 0; ei < boundary_edges.size(2); ++ei){
        if(((point_xyz[(same_axis+1)%3]-tm_.tetmesh_.node_((same_axis+1)%3, boundary_edges(0,ei)))
            *(point_xyz[(same_axis+1)%3]-tm_.tetmesh_.node_((same_axis+1)%3, boundary_edges(1,ei))) < 0)
           &&((point_xyz[(same_axis+2)%3] < tm_.tetmesh_.node_((same_axis+2)%3, boundary_edges(0,ei)))
              ||(point_xyz[(same_axis+2)%3] < tm_.tetmesh_.node_((same_axis+2)%3, boundary_edges(1,ei)))) )
          ++hit_edges;
      }
    if(hit_edges%2 == 0) return false;
    return true;
  }

  void cal_patch_xyz_bound(const vector<size_t> & one_patch,
                           matrix<double> & xyz_bound)
  {
    if(xyz_bound.size(1) != 3 || xyz_bound.size(2) != 2) {
        xyz_bound.resize(3,2);
      }

    xyz_bound(colon(),0) = ones<double>(3,1)*std::numeric_limits<double>::max();
    xyz_bound(colon(),1) = -1*ones<double>(3,1)*std::numeric_limits<double>::max();

    itr_matrix<const size_t*> pfm(one_patch.size(),1, &one_patch[0]);
    matrix<size_t> select_faces = tm_.outside_face_(colon(), pfm);
    set<size_t> points;
    points.insert(select_faces.begin(), select_faces.end());
    matrix<size_t> point_m(points.size(),1);
    std::copy(points.begin(), points.end(), point_m.begin());

    calc_bounding_box(tm_.tetmesh_.node_(colon(), point_m), &xyz_bound[0]);

    for(size_t i = 0; i < tm_.tetmesh_.node_.size(1); ++i){
        if(fabs(xyz_bound(i,0) - xyz_bound(i,1)) < 1e-6){
            xyz_bound(i,0) = xyz_bound(i,1);
          }
      }
  }

  void cal_xyz_sorted_list()
  {
    if(sorted_list_.size() != 3) sorted_list_.resize(3);
    for(size_t ci = 0; ci < sorted_list_.size(); ++ci) sorted_list_[ci].clear();

    for(size_t ci = 0 ;ci < patch_bound_.size(); ++ci){
        for(size_t di = 0; di < 3; ++di){
            if(fabs(patch_bound_[ci](di,0)-patch_bound_[ci](di,1)) < 1e-6){
                sorted_list_[di].push_back(make_pair(patch_bound_[ci](di,0), ci));
              }
          }
      }
    for(size_t ci = 0; ci < sorted_list_.size(); ++ci)
      sort(sorted_list_[ci].begin(), sorted_list_[ci].end());
  }

  void average_node(const deque<pair<size_t,size_t> > & one_chain,
                    matrix<double> & node,
                    const size_t dim_type,
                    const double v)
  {
    set<size_t> point;
    for(size_t ei = 0; ei < one_chain.size(); ++ei){
        point.insert(one_chain[ei].first);
        point.insert(one_chain[ei].second);
      }

    for(const auto & p : point)
      node(dim_type, p) = v;
  }


  void get_edges_from_tet(vector<pair<size_t,size_t> > & edges)
  {
    edges.clear();
    edges.reserve(tm_.ortae_.e2t_.size());
    pair<size_t,size_t> edge_point_pair;
    for(const auto & one_edge : tm_.ortae_.e2t_){
        edge_point_pair = one_edge.first;
        if(edge_point_pair.first > edge_point_pair.second)
          swap(edge_point_pair.first, edge_point_pair.second);
        edges.push_back(edge_point_pair);
      }
  }

  void optimize_tet(const vector<vector<vector<size_t> > > & variable_group){
    typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
    typedef std::shared_ptr<math_func_type> math_func_ptr;

    for(size_t i = 1; i < 10; ++i){
        shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);

        const double weight = pow(3,i*1.0);

        func->push_back(math_func_ptr(
                          new hj::math_func::sumsqr<double,int32_t>(
                            build_arap_math_func(tm_.tetmesh_.mesh_, tm_.tetmesh_.node_, 1.0))));
        for(size_t di = 0; di < variable_group.size(); ++di)
          for(size_t gi = 0; gi < variable_group[di].size(); ++gi){
              double target = 0;
              for(auto it = variable_group[di][gi].begin(); it != variable_group[di][gi].end(); ++it)
                target += tm_.tetmesh_.node_[*it];
              target /= variable_group[di][gi].size();

              func->push_back(math_func_ptr(
                                new hj::math_func::sumsqr<double,int32_t>(
                                  build_variable_fix_math_func(
                                    tm_.tetmesh_.node_.size(), variable_group[di][gi].begin(),
                                    variable_group[di][gi].end(),target, weight/variable_group[di][gi].size()))));
            }

        math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(func));
        math_func_ptr obj(new hj::math_func::sum<double,int32_t>(fun_cat));

        boost::property_tree::ptree pt;
        pt.put("package.value", "jtf");
        pt.put("alg.value","SQP");
        pt.put("iter.value",50);
        pt.put("epsg.value",1e-9);
        pt.put("linear_solver/type.value","direct");
        pt.put("linear_solver/name.value","cholmod");

        jtf::optimize(*obj,tm_.tetmesh_.node_,pt,nullptr,nullptr,nullptr);
      }

    {
      ofstream ofs("after_deform.vtk");
      tet2vtk(ofs, &tm_.tetmesh_.node_[0], tm_.tetmesh_.node_.size(2), &tm_.tetmesh_.mesh_[0], tm_.tetmesh_.mesh_.size(2));
    }
  }

  void sort_node_groups(const vector<group<size_t> > &node_group,
                        const vector<size_t> & integer_group_idx,
                        vector<vector<vector<size_t> > > &node_group_vec) const {
    node_group_vec.clear();
    node_group_vec.resize(3);
    vector<size_t> one_group_vec;
    for(size_t ii = 0; ii < integer_group_idx.size(); ++ii){
        const auto & one_group = node_group[integer_group_idx[ii]];
        if(one_group.empty()) continue;
        one_group_vec.resize(one_group.size());
        copy(one_group.begin(), one_group.end(), one_group_vec.begin());
        const size_t dim = one_group_vec.front()%3;
        node_group_vec[dim].push_back(one_group_vec);
      }

    vector<pair<double,size_t> > sorted_groups;

    for(size_t dim = 0; dim < node_group_vec.size(); ++dim){
        sorted_groups.clear();
        sorted_groups.reserve(node_group_vec[dim].size());
        for(size_t gi = 0; gi < node_group_vec[dim].size(); ++gi){
            vector<size_t> & one_group = node_group_vec[dim][gi];
            const size_t point_idx = one_group.front()/3;
            sorted_groups.push_back(make_pair(tm_.tetmesh_.node_(dim, point_idx), gi));
          }
        sort(sorted_groups.begin(), sorted_groups.end());

        for(size_t i = 0; i < sorted_groups.size(); ++i){
            node_group_vec[dim].push_back(node_group_vec[dim][sorted_groups[i].second]);
          }
        node_group_vec[dim].erase(node_group_vec[dim].begin(), node_group_vec[dim].begin()+sorted_groups.size());
      }
  }
private:
  jtf::tet_mesh tm_;
  jtf::mesh::meshes meta_mesh_;
  matrix<size_t> tri_type_;
  vector<vector<size_t> > patches_;
  vector<matrix<double> > patch_bound_;
  shared_ptr<jtf::mesh::patch_separater> ps_;
  vector<shared_ptr<jtf::mesh::edge2cell_adjacent> > patch_ea_;
  vector<vector<pair<double,size_t> > > sorted_list_;
  matrix<double> average_list_interval_;
  const boost::unordered_map<size_t,size_t> & surface_type_;
};




void cal_tet_coordinate1(const zjucad::matrix::matrix<double> & tet_node,
                         const zjucad::matrix::matrix<double> & one_p,
                         zjucad::matrix::matrix<double> & coordinate)
{
  matrix<double> edge(3,3);
  for(size_t pi = 1; pi < tet_node.size(2); ++pi){
      edge(colon(), pi-1) = tet_node(colon(), pi) - tet_node(colon(), 0);
    }
  inv(edge);
  coordinate = edge*(one_p-tet_node(colon(), 0));
}

bool is_inside_tet(const zjucad::matrix::matrix<double>& coordinate){
  //    for(size_t i = 0; i < coordinate.size(); ++i) {
  //        if(fabs(coordinate[i]) > 5e-2 && coordinate[i] < 0) return false;
  //        //if(coordinate[i] > 1 && fabs(coordinate[i]-1) > 5e-2) return false;
  //      }
  return true;
}

void cal_tet_coordinate(const zjucad::matrix::matrix<double> & one_p,
                        const size_t npt, const jtf::tet_mesh &tm,
                        tuple<size_t,double,double,double> & one_tet_coordinate)
{
  jtf::mesh::one_ring_tet_at_point ortap;
  ortap.add_tets(tm.tetmesh_.mesh_);
  auto it = ortap.p2t_.find(npt);
  assert(it != ortap.p2t_.end());
  const vector<size_t> & one_ring_tets = it->second;

  matrix<double> coordinate;
  for(size_t ti = 0; ti < one_ring_tets.size(); ++ti){
      cal_tet_coordinate1(tm.tetmesh_.node_(colon(), tm.tetmesh_.mesh_(colon(), one_ring_tets[ti])),
                          one_p, coordinate);
      if(is_inside_tet(coordinate)) {
          get<0>(one_tet_coordinate) = one_ring_tets[ti];
          get<1>(one_tet_coordinate) = coordinate[0];
          get<2>(one_tet_coordinate) = coordinate[1];
          get<3>(one_tet_coordinate) = coordinate[2];
          return ;
        }
    }
  cerr << "# [error] strange, do not find any tets contains given point." << endl;
  return;
}

void update_node(zjucad::matrix::matrix<double> &node,
                 const vector<tuple<size_t,double,double,double> > & coordinate,
                 const jtf::tet_mesh & tm_orig)
{
  matrix<double> edge(3,3);
  matrix<double> lambda(3,1);
  matrix<size_t> one_tet;
  for(size_t pi = 0; pi < coordinate.size(); ++pi){
      one_tet = tm_orig.tetmesh_.mesh_(colon(), get<0>(coordinate[pi]));
      for(size_t di = 0; di < 3; ++di){
          edge(colon(), di) = tm_orig.tetmesh_.node_(colon(), one_tet[di + 1]) - tm_orig.tetmesh_.node_(colon(), one_tet[0]);
        }
      lambda[0] = get<1>(coordinate[pi]);
      lambda[1] = get<2>(coordinate[pi]);
      lambda[2] = get<3>(coordinate[pi]);
      node(colon(), pi) = edge*lambda + tm_orig.tetmesh_.node_(colon(), one_tet[0]);
    }
}

void map_meta_to_orig(const jtf::tet_mesh &tm,
                      const jtf::tet_mesh &tm_orig,
                      jtf::mesh::meshes & meta)
{
  K_near_points knp(tm.tetmesh_.node_);
  vector<int> kpts;
  vector<double> dis;
  // tet_idx, lambda of <p0p1>, <p0p2>, <p0p3>
  vector<tuple<size_t,double,double,double> > tet_coordinate(meta.node_.size(2));
  matrix<double> new_meta_node = meta.node_;
  for(size_t i = 0; i < meta.node_.size(2); ++i) {
      knp.query_k_near_points(meta.node_(colon(),i), 1, kpts, dis);
      cal_tet_coordinate(meta.node_(colon(), i), kpts.front(), tm, tet_coordinate[i]);
    }
  update_node(new_meta_node, tet_coordinate, tm_orig);

  meta.node_ = new_meta_node;
}
int meta_mesh_extraction(ptree &pt)
{
  jtf::tet_mesh tm(pt.get<string>("input/polycube_tet.value").c_str());

  boost::unordered_map<size_t,size_t> surface_type;
  if(load_surface_type(pt.get<string>("input/surface_type.value").c_str(),
                       surface_type, tm.fa_.get())){
      cerr << "# [error] load surface type error." << endl;
      return __LINE__;
    }

  meta_mesh_extractor mme(tm, surface_type);
  mme.run(pt);

  jtf::mesh::meshes meta;
  mme.get_meta_mesh(meta);

  ofstream ofs("meta.vtk");
  hex2vtk(ofs, &meta.node_[0], meta.node_.size(2), &meta.mesh_[0], meta.mesh_.size(2));

  if(zjucad::has("input/orig_tet.value", pt)){
      jtf::tet_mesh orig_tm(pt.get<string>("input/orig_tet.value").c_str());
      if(norm(orig_tm.tetmesh_.mesh_ - tm.tetmesh_.mesh_) > 0) {
          cerr << "# [error] wrong tetmesh." << endl;
          return __LINE__;
        }
      map_meta_to_orig(tm, orig_tm, meta);

      ofstream ofs("meta_orig.vtk");
      hex2vtk(ofs, &meta.node_[0], meta.node_.size(2), &meta.mesh_[0], meta.mesh_.size(2));
    }

  return 0;
}


class sheet{
public:
  sheet(const zjucad::matrix::matrix<size_t> & hex,
        const zjucad::matrix::matrix<double> & node)
    :hex_(hex), node_(node){
    assert(node.size(1) == 3);
    assert(hex.size(1) == 8);
  }
  void run(){
    generate_node();
    generate_quads();
    generate_sheets();
  }
  void get_sheet_mesh(zjucad::matrix::matrix<size_t> & mesh,
                      zjucad::matrix::matrix<double> & node) const{
    mesh = sheet_mesh_;
    node = sheet_node_;
  }
  vector<vector<size_t> > get_sheet_patch()const{
    return sheet_patch_;
  }
  matrix<size_t> get_hex2sheets()const{
    return hex2sheets_;
  }
  matrix<size_t> get_sheet2hex()const{
    return sheet2hex_;
  }
private:
  void generate_node() {
    unique_ptr<jtf::mesh::face2hex_adjacent> fa(jtf::mesh::face2hex_adjacent::create(hex_));
    if(!fa.get()) {
        throw std::invalid_argument("# [error] can not build face2hex_adjacent.");
      }
    jtf::mesh::one_ring_hex_at_edge orhae_;
    orhae_.add_hex(hex_, node_, *fa);
    sheet_node_.resize(3, orhae_.e2h_.size());
    size_t ei = 0;
    for(const auto & one_edge : orhae_.e2h_){
        const pair<size_t,size_t> & one_edge_pair = one_edge.first;
        sheet_node_(colon(), ei++) = 0.5*(node_(colon(), one_edge_pair.first) + node_(colon(), one_edge_pair.second));
        if(one_edge_pair.first > one_edge_pair.second)
          hex_edge2point_[make_pair(one_edge_pair.second, one_edge_pair.first)] = ei-1;
        else
          hex_edge2point_[one_edge_pair] = ei - 1;
      }
  }
  void generate_quads(){
    sheet_mesh_.resize(4, 3*hex_.size(2));
    sheet2hex_.resize(3*hex_.size(2));
    hex2sheets_.resize(3, hex_.size(2));
    size_t mid_points_idx[12];
    for(size_t hi = 0; hi < hex_.size(2); ++hi){
        mid_points_idx[0] = get_edge_mid_idx(make_pair(hex_(0,hi), hex_(1, hi)));
        mid_points_idx[1] = get_edge_mid_idx(make_pair(hex_(1,hi), hex_(3, hi)));
        mid_points_idx[2] = get_edge_mid_idx(make_pair(hex_(3,hi), hex_(2, hi)));
        mid_points_idx[3] = get_edge_mid_idx(make_pair(hex_(2,hi), hex_(0, hi)));

        mid_points_idx[4] = get_edge_mid_idx(make_pair(hex_(4,hi), hex_(5, hi)));
        mid_points_idx[5] = get_edge_mid_idx(make_pair(hex_(5,hi), hex_(7, hi)));
        mid_points_idx[6] = get_edge_mid_idx(make_pair(hex_(7,hi), hex_(6, hi)));
        mid_points_idx[7] = get_edge_mid_idx(make_pair(hex_(6,hi), hex_(4, hi)));

        mid_points_idx[8] = get_edge_mid_idx(make_pair(hex_(0,hi), hex_(4, hi)));
        mid_points_idx[9] = get_edge_mid_idx(make_pair(hex_(1,hi), hex_(5, hi)));
        mid_points_idx[10] = get_edge_mid_idx(make_pair(hex_(7,hi), hex_(3, hi)));
        mid_points_idx[11] = get_edge_mid_idx(make_pair(hex_(6,hi), hex_(2, hi)));

        sheet_mesh_(0,3*hi+0) = mid_points_idx[0];
        sheet_mesh_(1,3*hi+0) = mid_points_idx[2];
        sheet_mesh_(2,3*hi+0) = mid_points_idx[6];
        sheet_mesh_(3,3*hi+0) = mid_points_idx[4];

        sheet_mesh_(0,3*hi+1) = mid_points_idx[1];
        sheet_mesh_(1,3*hi+1) = mid_points_idx[3];
        sheet_mesh_(2,3*hi+1) = mid_points_idx[7];
        sheet_mesh_(3,3*hi+1) = mid_points_idx[5];

        sheet_mesh_(0,3*hi+2) = mid_points_idx[8];
        sheet_mesh_(1,3*hi+2) = mid_points_idx[9];
        sheet_mesh_(2,3*hi+2) = mid_points_idx[10];
        sheet_mesh_(3,3*hi+2) = mid_points_idx[11];

        for(size_t i = 0; i < 3; ++i){
            sheet2hex_[3*hi+i] = hi;
            hex2sheets_(i,hi) = 3*hi+i;
          }
      }

  }
  void generate_sheets(){
    jtf::mesh::patch_separater ps(sheet_mesh_);
    ps.separater();
    sheet_patch_ = ps.get_all_patches();
  }

  size_t get_edge_mid_idx(const pair<size_t,size_t> & one_edge)const{
    map<pair<size_t,size_t>,size_t>::const_iterator it;
    if(one_edge.first < one_edge.second){
        it = hex_edge2point_.find(one_edge);
      }else
      it = hex_edge2point_.find(make_pair(one_edge.second, one_edge.first));
    if(it == hex_edge2point_.end())
      throw std::logic_error("can not find edge.") ;
    return it->second;
  }


private:
  sheet();
  const sheet& operator =(const sheet&);
private:
  const zjucad::matrix::matrix<size_t> & hex_;
  const zjucad::matrix::matrix<double> & node_;

  map<pair<size_t,size_t>,size_t> hex_edge2point_;
  vector<vector<size_t> > sheet_patch_;
  zjucad::matrix::matrix<size_t> sheet_mesh_;
  zjucad::matrix::matrix<size_t> sheet2hex_;
  zjucad::matrix::matrix<size_t> hex2sheets_;
  zjucad::matrix::matrix<double> sheet_node_;
};

int extract_sheets(ptree &pt)
{
  zjucad::matrix::matrix<size_t> hex;
  zjucad::matrix::matrix<double> node;
  if(jtf::hexmesh::hex_mesh_read_from_wyz(pt.get<string>("input/hex.value").c_str(),
                                          hex, node,1)){
      cerr << "# [error] can not load hexmesh." << endl;
      return __LINE__;
    }
  cerr << "# [info] hex number " << hex.size(2) << " node number "  << node.size(2) << endl;

  sheet st(hex, node);
  st.run();
  matrix<size_t> sheet_quad;
  matrix<double> sheet_node;
  st.get_sheet_mesh(sheet_quad, sheet_node);

  {
    ofstream ofs("sheet.vtk");
    quad2vtk(ofs, &sheet_node[0], sheet_node.size(2), &sheet_quad[0], sheet_quad.size(2));
  }

    {
      vector<vector<size_t> > sheet_patch = st.get_sheet_patch();
      for(size_t pi = 0; pi < sheet_patch.size(); ++pi){
          ostringstream os;
          os << "sheet_patch_" << pi << ".vtk";
          ofstream ofs(os.str().c_str());
          itr_matrix<const size_t*> select_face_idx(sheet_patch[pi].size(),1, &sheet_patch[pi][0]);
          zjucad::matrix::matrix<size_t> select_face = sheet_quad(colon(), select_face_idx);
          quad2vtk(ofs, &sheet_node[0], sheet_node.size(2), &select_face[0], select_face.size(2));
        }
    }


  matrix<size_t> hex2sheet = st.get_hex2sheets();
  vector<size_t> non_orthognal_hex;
  matrix<double> normal(3, sheet_quad.size(2));
  jtf::mesh::cal_face_normal(sheet_quad, sheet_node, normal);
  matrix<double> dir3(3,3);
//  vector<double>
//  for(size_t hi = 0; hi < hex2sheet.size(2); ++hi){
//      for(size_t di = 0; di < 3; ++di)
//        dir3(colon(),di) = normal(colon(), hex2sheet(di, hi));

//    }
  return 0;
}
