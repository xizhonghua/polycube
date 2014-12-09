#include <fstream>
#include <numeric>
#include <hjlib/function/function.h>
#include <hjlib/function/func_aux.h>
#include <zjucad/optimizer/optimizer.h>
#include "../tetmesh/tetmesh.h"
#include <jtflib/mesh/mesh.h>
#include "../common/vtk.h"
#include "../common/IO.h"
using namespace std;
using namespace hj::function;
using namespace zjucad::matrix;
////////////////////////////////////////////////////
// io:
////////////////////////////////////////////////////

int load_weight_constraint_file(const char * file,
                                boost::unordered_map<size_t,double> & p2w)
{
  ifstream ifs(file);
  if(ifs.fail()){
    cerr << "# [error] can not open file." << endl;
    return __LINE__;
  }

  p2w.clear();
  size_t p0, p1,p2;
  double val;
  while(!ifs.eof()){
    ifs >> p0 >> p1 >> p2 >> val;
    if(ifs.eof())
      break;

    p2w[p0] = val;
    p2w[p1] = val;
    p2w[p2] = val;
  }

  return 0;
}

class cot_lap_function : public hj::function::function_t<double,int32_t>
{
public:
  cot_lap_function(const matrix<double> & node,
                   const size_t p0,
                   const vector<size_t> & adj_point,
                   const vector<double> & adj_weight,
                   const double weight = 1.0)
    :node_num_(node.size(2)), p0_(p0) , adj_point_(adj_point),
      adj_weight_(adj_weight), w_(weight){
    total_weight_ = std::accumulate(adj_weight.begin(), adj_weight.end(),0.0);
  }

  virtual size_t dim_of_x(void) const {
    return node_num_;
  }
  virtual size_t dim_of_f(void) const {
    return 1;
  }
  virtual int val(const double*x, double *f, hj::function::func_ctx *ctx = 0)const{
    *f = x[p0_] ;//total_weight_;
    for(size_t pi = 0; pi < adj_point_.size(); ++pi){
      *f -= adj_weight_[pi] * x[adj_point_[pi]] / total_weight_;
    }
    *f *= w_;
    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {

    ptr[1] = ptr[0] + (1 + adj_point_.size());
    idx[ptr[0] + 0] = p0_ ;
    val[ptr[0] + 0] =  w_;

    for(size_t pi = 0; pi < adj_point_.size(); ++pi){
      idx[ptr[0] + pi + 1] = adj_point_[pi] ;
      val[ptr[0] + pi + 1] = -1.0 * w_ * adj_weight_[pi]/total_weight_;
    }

    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return (1 + adj_point_.size()) * dim_of_f();
  }
private:
  const size_t node_num_;
  const size_t p0_;
  const double w_;
  const vector<size_t> adj_point_;
  double total_weight_;
  const vector<double> adj_weight_;
};

double calc_cot_weight(
    const size_t & p0, const size_t &p1,
    const matrix<size_t> & faces,
    const matrix<double> & node,
    const jtf::mesh::edge2cell_adjacent & ea)
{
  const size_t edge_idx = ea.get_edge_idx(p0, p1);
  if(edge_idx == -1) {
    cerr << "# [error] can not find edge " << p0 << " " << p1 << endl;
    return 0;
  }
  const pair<size_t,size_t> & tri_pair = ea.edge2cell_[edge_idx];
  vector<size_t> tri_pair_vec;
  tri_pair_vec.push_back(tri_pair.first);
  tri_pair_vec.push_back(tri_pair.second);
  matrix<double> e(3,2);
  vector<double> cot_val_vec;
  for(size_t i = 0; i < tri_pair_vec.size(); ++i){
    if(tri_pair_vec[i] == -1) continue;
    const size_t other_idx =
        std::accumulate(faces(colon(),tri_pair.first).begin(),
                        faces(colon(),tri_pair.first).end(),0)
        - p0 - p1;
    e(colon(),0) = node(colon(),p0) - node(colon(),other_idx);
    e(colon(),1) = node(colon(),p1) - node(colon(),other_idx);
    double len_0 = norm(e(colon(),0));
    if(len_0 < 1e-6) len_0 = 1.0;
    double len_1 = norm(e(colon(),1));
    if(len_1 < 1e-6) len_1 = 1.0;
    e(colon(),0) /= len_0;
    e(colon(),1) /= len_1;
    const double cos_val = dot(e(colon(),0), e(colon(),1));
    const double sin_val = sqrt(1 - cos_val * cos_val);
    const double cot_val = cos_val / sin_val;
    cot_val_vec.push_back(cot_val);
  }

  return std::accumulate(cot_val_vec.begin(), cot_val_vec.end(), 0.0) / cot_val_vec.size();
}

hj::function::function_t<double, int32_t> *
build_laplacian_function(const matrix<double> & node,
                         const matrix<size_t> & face,
                         const jtf::mesh::edge2cell_adjacent &ea,
                         const boost::unordered_map<size_t,double> & constraint_compact,
                         const double weight = 1.0)
{
  unique_ptr<jtf::mesh::one_ring_point_at_point> orpap(
        jtf::mesh::one_ring_point_at_point::create(face));
  if(!orpap.get()){
    cerr << "# [error] can not build one_ring_point_at_point." << endl;
    return 0;
  }

  vector<double> cot_weight;
  std::shared_ptr<vector<std::shared_ptr<function_t<double,int32_t> > > >
      funcs(new vector<std::shared_ptr<function_t<double,int32_t> > >);

  for(boost::unordered_map<size_t, vector<size_t> >::const_iterator
      cit = orpap->p2p_.begin(); cit != orpap->p2p_.end(); ++cit){
//    if(constraint_compact.find(cit->first) != constraint_compact.end())
//      continue;
    const vector<size_t> & arounding_points = cit->second;
    cot_weight.resize(arounding_points.size());
    for(size_t p = 0; p < arounding_points.size(); ++p){
      cot_weight[p] = 1.0;//calc_cot_weight(cit->first, arounding_points[p],face, node, ea);
    }

    std::shared_ptr<function_t<double, int32_t> > pf(
          new cot_lap_function(node, cit->first, arounding_points,
                               cot_weight,sqrt(weight)));
    funcs->push_back(pf);

#if 0 // validate jac
    {
      const double err = jac_err(*pf, &node[0]);
      if(err > 1e-5) {
        cerr << "large error in jac: " << err << endl;
      }
    }
#endif
  }
  return new_catenated_function<double,int32_t>(funcs);
}

class fix_value_func : public hj::function::function_t<double, int32_t>
{
public:
  fix_value_func(const zjucad::matrix::matrix<double> &node,
                 const size_t p,
                 const double val,
                 const double weight)
    :p_(p), val_(val), node_num_(node.size(2)), weight_(weight){}
  virtual size_t dim_of_x(void) const {
    return node_num_;
  }
  virtual size_t dim_of_f(void) const {
    return 1;
  }
  virtual int val(const double *x, double *f,
                  hj::function::func_ctx *ctx = 0) const {
    *f = (x[p_] - val_) * weight_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    ptr[1] = ptr[0]+1;
    idx[ptr[0]] = p_;
    val[ptr[0]] = weight_;
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 1*dim_of_f();
  }
protected:
  const size_t p_;
  const double val_;
  size_t node_num_;
  const double weight_;
};

hj::function::function_t<double, int32_t> *
build_equation_function(
    const matrix<double> &compact_node_mat,
    const boost::unordered_map<size_t,double> & constraint_compact,
    const double weight)
{
  std::shared_ptr<vector<std::shared_ptr<function_t<double,int32_t> > > >
      funcs(new vector<std::shared_ptr<function_t<double, int32_t> > >);
  for(boost::unordered_map<size_t, double>::const_iterator cit =
      constraint_compact.begin(); cit != constraint_compact.end(); ++cit){
    funcs->push_back(std::shared_ptr<function_t<double,int32_t> >
                     (new fix_value_func(
                        compact_node_mat, cit->first, cit->second,
                        sqrt(weight))));
  }
  return new_catenated_function<double, int32_t>(funcs);
}

int sovle_weight(const matrix<size_t> & outside_face,
                 const matrix<double> & node,
                 matrix<double> & node_weight,
                 const boost::unordered_map<size_t,double> & constraint,
                 boost::property_tree::ptree &pt)
{
  boost::unordered_map<size_t,size_t> p2cp; // read_idx_to_compact_idx
  set<size_t> compact_node(outside_face.begin(), outside_face.end());
  vector<size_t> compact_node_idx(compact_node.begin(), compact_node.end());
  for(size_t pi  = 0; pi < compact_node_idx.size(); ++pi){
    p2cp[compact_node_idx[pi]] = pi;
  }

  matrix<size_t> compact_face = outside_face;
  matrix<double> compact_node_mat = zeros<double>(3, compact_node.size());
  for(size_t p = 0; p< compact_face.size(); ++p)
    compact_face[p] = p2cp[compact_face[p]];

  for(size_t p = 0; p < compact_node_idx.size(); ++p){
    compact_node_mat(colon(), p) = node(colon(), compact_node_idx[p]);
  }

  boost::unordered_map<size_t,double> constraint_compact;
  {
    for(boost::unordered_map<size_t,double>::const_iterator cit
        = constraint.begin(); cit != constraint.end(); ++cit){
      constraint_compact[p2cp[cit->first]] = cit->second;
    }
  }
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(compact_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2tria_adjacent." << endl;
    return __LINE__;
  }

  std::shared_ptr<vector<std::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<std::shared_ptr<function_t<double, int32_t> > >);
  std::shared_ptr<function_t<double, int32_t> > lap_function(
        build_laplacian_function(compact_node_mat, compact_face, *ea, constraint_compact));

  std::shared_ptr<function_t<double, int32_t> > constraint_function(
        build_equation_function(compact_node_mat, constraint_compact,1.0));

  if(lap_function.get()){
    funcs->push_back(lap_function);
    cerr << "# [add] lap function." << endl;
  }

  if(constraint_function.get()){
    funcs->push_back(constraint_function);
    cerr << "# [add] constraint_function" << endl;
  }

  unique_ptr<function_t<double,int32_t> > total_func(
        new_catenated_function<double,int32_t>(funcs));

  matrix<double> node_temp = ones<double>(1,total_func->dim_of_x());
  for(boost::unordered_map<size_t,double>::const_iterator cit =
      constraint_compact.begin(); cit != constraint_compact.end(); ++cit){
    node_temp[cit->first] = cit->second;
  }
  cerr << "# [info] before opt "  << norm(node_temp) << endl;

  matrix<double> residual = zeros<double>(total_func->dim_of_f(),1);
  zjucad::optimize(*total_func, node_temp, residual, pt);

  cerr << "# [info] after opt "  << norm(node_temp) << endl;

  node_weight = zeros<double>(1, node.size(2));
  for(boost::unordered_map<size_t,size_t>::const_iterator cit = p2cp.begin();
      cit != p2cp.end(); ++cit){
    node_weight[cit->first] = node_temp[cit->second];
  }
  return 0;
}


int diffustion_weight(int argc, char * argv[])
{
  if(argc != 4){
    cerr << "# [usage] diffustion_weight tet detail_file no_detail_file" << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  matrix<size_t> outside_face;
  {
    unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
    if(!fa.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }
    get_outside_face(*fa, outside_face);
  }

  boost::unordered_map<size_t, double> constraint;
//  if(load_weight_constraint_file(argv[2], constraint))
//    return __LINE__;
  // temp load detal_file and no_detal_file
  {
    ifstream ifs_detail(argv[2]);
    ifstream ifs_no_detail(argv[3]);
    if(ifs_detail.fail() || ifs_no_detail.fail()){
      cerr << "# [error] can not load detail control file." << endl;
      return __LINE__;
    }
    size_t idx;
    while(!ifs_detail.eof()){
      ifs_detail >> idx;
      if(ifs_detail.eof()) break;
      constraint[idx] = 0.0;
    }

    const size_t detail_remain_points = constraint.size();

    while(!ifs_no_detail.eof()){
      ifs_no_detail >> idx;
      if(ifs_no_detail.eof()) break;
      constraint[idx] = 1.0;
    }

    const size_t smooth_points = constraint.size() - detail_remain_points;

    cerr << "# [info] detail remain points: " << detail_remain_points << endl;
    cerr << "# [info] detail smooth points: " << smooth_points << endl;

  }

  matrix<double> node_weight = zeros(tm.node_.size(2), 1);

  boost::property_tree::ptree pt;
  pt.put("package.value","hj");
  pt.put("alg.value","More");
  pt.put("iter.value","100");
  sovle_weight(outside_face, tm.node_, node_weight, constraint,pt);

  double max_e = *max_element(node_weight.begin(), node_weight.end());
  double min_e = *min_element(node_weight.begin(), node_weight.end());
  if(fabs(max_e - min_e) < 1e-6) node_weight -= min_e;

  node_weight /= (max_e - min_e);
  node_weight -= min_e/(max_e-min_e);

  ofstream ofs("point_val.vtk");
  vector<size_t> node_idx(tm.node_.size(2));
  for(size_t pi = 0; pi < tm.node_.size(2); ++pi) node_idx[pi] = pi;
  point2vtk(ofs, &tm.node_[0], tm.node_.size(2), &node_idx[0], node_idx.size());
  cell_data(ofs, &node_weight[0], node_weight.size(), "val");

  {
    ofstream ofs_tri("face_val.vtk");
    matrix<double> tri_val = zeros<double>(outside_face.size(2),1);
    tri2vtk(ofs_tri, &tm.node_[0], tm.node_.size(2), &outside_face[0], outside_face.size(2));
    for(size_t fi = 0; fi < outside_face.size(2); ++fi){
      for(size_t pi = 0; pi < outside_face.size(1); ++pi)
        tri_val[fi] += node_weight[outside_face(pi,fi)];
      tri_val[fi] /= 3.0;
    }
    cell_data(ofs_tri, &tri_val[0], tri_val.size(), "diffused_val");
  }

  jtf::mesh::write_matrix("diffustion_value", node_weight);
  return 0;
}
