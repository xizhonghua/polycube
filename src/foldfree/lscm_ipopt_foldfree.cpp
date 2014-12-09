#include "lscm_ipopt_foldfree.h"

#include <set>
#include <fstream>
#include <Ipopt/coin/IpIpoptApplication.hpp>

namespace jy{
using namespace Ipopt;
typedef std::set<size_t>::const_iterator set_iter;

lscm_ipopt_foldfree::lscm_ipopt_foldfree(const char *file_mesh,
                                         const char *file_fix_vert,
                                         const char *output_uv)
  :reverse_(true), output_uv_(output_uv)
{
  if(!init(file_mesh,file_fix_vert)){
    std::cout<<"initilize failure!\n";
  }
}

lscm_ipopt_foldfree::~lscm_ipopt_foldfree()
{
}

bool lscm_ipopt_foldfree::init(const char *file_mesh,
                               const char *file_fix_vert)
{

  std::set<size_t> verts;

  mesh_.reset(new tri_mesh);

  if(!read_obj_mesh(file_mesh,*mesh_.get())){
    std::cerr<<"read mesh file error!"<<std::endl;
    return false;
  }

  //set_features(fpts);
  if(!read_fix_and_init_vertices(file_fix_vert,
                                 fix_ids_,fix_coords_,init_x_)){
    std::cout<<"read error!"<<std::endl;
    return false;
  }

  lscm_weight_.reserve(6*mesh_->num_faces_);
  int_weight_.reserve(mesh_->num_faces_);
  size_t num_fold = 0;
  avg_area_ = 0.0;
  for(size_t i = 0;i < mesh_->num_faces_;i++){
    size_t *f =  mesh_->faces_[i];
    double *v0 = mesh_->verts_[f[0]];
    double *v1 = mesh_->verts_[f[1]];
    double *v2 = mesh_->verts_[f[2]];

    size_t fv[3] = {2*f[0],2*f[1],2*f[2]};

    double area_ =
        (init_x_[fv[1]]-init_x_[fv[0]])*(init_x_[fv[2]+1]-init_x_[fv[0]+1])
        - (init_x_[fv[2]]-init_x_[fv[0]])*(init_x_[fv[1]+1]-init_x_[fv[0]+1]);
    if( area_ <= 0){
      ++num_fold;
    }
    avg_area_ += fabs(area_);

    verts.insert(f[0]);
    verts.insert(f[1]);
    verts.insert(f[2]);

    double aa = dist2(v0,v1);
    double bb = dist2(v1,v2);
    double cc = dist2(v2,v0);
    double a  = sqrt(aa);
    double c  = sqrt(cc);

    double angle = acos((aa+cc-bb)/(2*a*c));

    double x2 = c * cos(angle);
    double y2 = c * sin(angle);

    double area = a * y2;
    int_weight_.push_back(area);


    lscm_weight_.push_back((x2 - a)/area);
    lscm_weight_.push_back(y2 / area);
    lscm_weight_.push_back(-x2 / area);
    lscm_weight_.push_back(-y2 / area);
    lscm_weight_.push_back(a / area);
    lscm_weight_.push_back(0);
  }

  reverse_ = (num_fold > mesh_->num_faces_ - num_fold) ? false : true;

  front_map_.reserve(verts.size());
  back_map_.resize(mesh_->num_verts_,-1);

  size_t i = 0;
  for(set_iter begin = verts.begin(),end = verts.end();
      begin != end;++begin,++i){
    front_map_.push_back(*begin);
    back_map_[*begin] = i;
  }

  hess_.resize(2*front_map_.size());
  for(size_t i = 0;i < mesh_->num_faces_;++i){
    const double *lscm_weight = &lscm_weight_[6*i];
    size_t *face = mesh_->faces_[i];
    //size_t f[3] = {2*face[0],2*face[1],2*face[2]};
    size_t f[3] = {2*back_map_[face[0]],2*back_map_[face[1]],2*back_map_[face[2]]};

    double factor = 2.0 * int_weight_[i];
    for(size_t k,kk,j = 0,jj = 0;j < 3;++j,jj+=2){
      size_t row    = f[j];
      double ww1 = lscm_weight[jj];
      double ww2 = lscm_weight[jj+1];

      for(k = kk = 0;k < 3;++k,kk += 2){
        int col = f[k];
        if(row + 1 <= col){
          hess_[row]  [col] += factor*( ww1*lscm_weight[kk]   + ww2*lscm_weight[kk+1]);
          hess_[row][col+1] += factor*(-ww1*lscm_weight[kk+1] + ww2*lscm_weight[kk]);

          hess_[row+1][col]   += factor*(-ww2*lscm_weight[kk]  + ww1*lscm_weight[kk+1]);
          hess_[row+1][col+1] += factor*( ww2*lscm_weight[kk+1]+ ww1*lscm_weight[kk]);
        }else if(row <= col){
          hess_[row]  [col]   += factor*(ww1*lscm_weight[kk] + ww2*lscm_weight[kk+1]);

          hess_[row][col+1] += factor*(-ww1*lscm_weight[kk+1] + ww2*lscm_weight[kk]);
          hess_[row+1][col+1] += factor*( ww2*lscm_weight[kk+1]+ ww1*lscm_weight[kk]);
        }
        else if(row <= col+1){
          hess_[row][col+1] += factor*(-ww1*lscm_weight[kk+1] + ww2*lscm_weight[kk]);
        }
      }
    }
  }

  avg_area_ /= mesh_->num_faces_;

  return true;
}


// returns the size of the problem
bool lscm_ipopt_foldfree::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                       Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described has n variables, x[0] through x[n-1]
  n = hess_.size();

  // has m inequality constraints
  m = mesh_->num_faces_;
  //m = 0;

  // in this example the jacobian is dense and contains 8 nonzeros
  nnz_jac_g = m*6;

  // the hessian is also dense and has 16 total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag = 0;
  for(size_t i = 0;i < n;++i){
    nnz_h_lag += hess_[i].size();
  }

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool lscm_ipopt_foldfree::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                          Index m, Number* g_l, Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.

  // the variables have lower bounds of 1
  for (Index i=0; i<n; ++i) {
    x_l[i] = -2e19;
    x_u[i] =  2e19;
  }

  size_t num_feature = fix_ids_.size();
  for(size_t i=0;i < num_feature;++i){
    size_t index = 2 * back_map_[fix_ids_[i]];
    x_l[index]   = x_u[index]   = fix_coords_[2*i];
    x_l[index+1] = x_u[index+1] = fix_coords_[2*i+1];
  }

  // the first constraint g1 has NO upper bound, here we set it to 2e19.
  // Ipopt interprets any number greater than nlp_upper_bound_inf as
  // infinity. The default value of nlp_upper_bound_inf and nlp_lower_bound_inf
  // is 1e19 and can be changed through ipopt options.


  if(reverse_){
    for(size_t i=0;i < m;++i){
      g_l[i] = 1e-6; //
      //g_l[i] = 0;
      g_u[i] = 2e19;
    }
  }
  else{
    for(size_t i=0;i < m;++i){
      g_l[i] = -2e19;
      g_u[i] = -1e-6;
      //g_u[i] = 0;
    }
  }

  return true;
}

// returns the initial point for the problem
bool lscm_ipopt_foldfree::get_starting_point(Index n, bool init_x, Number* x,
                                             bool init_z, Number* z_L, Number* z_U,
                                             Index m, bool init_lambda,
                                             Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // initialize to the given starting point
  //    for (Index i=0; i<n; ++i) {
  //        x[i] = 1.0;
  //    }
  //texture_.resize(n);

  /*
    jy::lscm_knitro_problem0 lscm;
    if(!lscm.init(file_mesh.c_str(),file_fix_vert.c_str())){
        std::cerr<<"read error!"<<std::endl;
        return 0;
    }

    if(lscm.solve(texture_)){
        std::cout<<"solve successfully!"<<std::endl;
        return 1;
    }
*/
  //    if(!read_init_guess("../data/texture.uvw",x)){
  //        std::cout<<"read error!"<<std::endl;
  //        return false;
  //    }


  //  for (Index i=0; i<n; ++i) {
  //    x[i] = fix_coords_[0];
  //  }
  //  size_t num_feature = fix_ids_.size();
  //  for(size_t i=0;i < num_feature;++i){
  //    size_t index = 2 * back_map_[fix_ids_[i]];
  //    x[index]   = fix_coords_[2*i];
  //    x[index+1] = fix_coords_[2*i+1];
  //  }


  for (Index i=0; i < front_map_.size();++i) {
    Index j = 2 * front_map_[i];
    x[2*i]   = init_x_[j];
    x[2*i+1] = init_x_[j+1];
  }


  return true;
}

// returns the value of the objective function
bool lscm_ipopt_foldfree::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  //assert(n == 4);
  obj_value = 0.0;
  for(size_t i = 0;i < mesh_->num_faces_;++i){
    const double *w = &lscm_weight_[i*6];
    size_t *face = mesh_->faces_[i];
    size_t fv[3] = {2*back_map_[face[0]],2*back_map_[face[1]],2*back_map_[face[2]]};
    double v[6] = {x[fv[0]],x[fv[0]+1],x[fv[1]],x[fv[1]+1],x[fv[2]],x[fv[2]+1]};

    double fun1 = (w[0]*v[0] + w[2]*v[2] + w[4]*v[4] - w[1]*v[1] - w[3]*v[3]);
    double fun2 = (w[1]*v[0] + w[3]*v[2] + w[0]*v[1] + w[2]*v[3] + w[4]*v[5]);

    obj_value += int_weight_[i]*(fun1*fun1 + fun2*fun2);
  }

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool lscm_ipopt_foldfree::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  //assert(n == 4);
  for(size_t i = 0;i < n;++i){
    grad_f[i] = 0.0;
  }

  for(size_t i = 0;i < mesh_->num_faces_;++i){
    const double *w = &lscm_weight_[i*6];
    size_t *face = mesh_->faces_[i];
    size_t fv[3] = {2*back_map_[face[0]],2*back_map_[face[1]],2*back_map_[face[2]]};
    double v[6] = {x[fv[0]],x[fv[0]+1],x[fv[1]],x[fv[1]+1],x[fv[2]],x[fv[2]+1]};

    double f1 = 2*int_weight_[i]*(w[0]*v[0] - w[1]*v[1] + w[2]*v[2] - w[3]*v[3] + w[4]*v[4]);
    double f2 = 2*int_weight_[i]*(w[1]*v[0] + w[0]*v[1] + w[3]*v[2] + w[2]*v[3] + w[4]*v[5]);

    grad_f[fv[0]]   +=  f1*w[0] + f2*w[1];
    grad_f[fv[0]+1] += -f1*w[1] + f2*w[0];

    grad_f[fv[1]]   +=  f1*w[2] + f2*w[3];
    grad_f[fv[1]+1] += -f1*w[3] + f2*w[2];

    grad_f[fv[2]]   +=  f1*w[4];
    grad_f[fv[2]+1] +=  f2*w[4];
  }

  return true;
}

// return the value of the constraints: g(x)
bool lscm_ipopt_foldfree::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  for(Index i = 0;i < m;++i){//(w1*w2-w0*w3)*(x1*x2-x0*x3)+w1*w4*(x1*x4-x0*x5)+w3*w4*(x3*X4-x2*x5)
    //const double *w = lscm_weight_[i];
    size_t *face = mesh_->faces_[i];
    size_t fv[3] = {2*back_map_[face[0]],2*back_map_[face[1]],2*back_map_[face[2]]};
    //double v[6] = {x[fv[0]],x[fv[0]+1],x[fv[1]],x[fv[1]+1],x[fv[2]],x[fv[2]+1]};

    g[i] = (x[fv[1]]-x[fv[0]])*(x[fv[2]+1]-x[fv[0]+1])-
           (x[fv[2]]-x[fv[0]])*(x[fv[1]+1]-x[fv[0]+1]);
  }

  return true;
}

// return the structure or values of the jacobian
bool lscm_ipopt_foldfree::eval_jac_g(Index n, const Number* x, bool new_x,
                                     Index m, Index nele_jac, Index* iRow, Index *jCol,
                                     Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian
    Index idx = 0;
    for(Index i = 0;i < m;++i){
      size_t *f = mesh_->faces_[i];
      for(Index j = 0;j < 3;++j,idx += 2){
        iRow[idx] = i;
        jCol[idx] = 2*back_map_[f[j]];

        iRow[idx+1] = i;
        jCol[idx+1] = jCol[idx] + 1;
      }
    }
  }
  else {
    // return the values of the jacobian of the constraints
    for(Index i = 0,idx = 0;i < m;++i){
      size_t *face = mesh_->faces_[i];
      size_t fv[3] = {2*back_map_[face[0]],2*back_map_[face[1]],2*back_map_[face[2]]};
      double v[6] = {x[fv[0]],x[fv[0]+1],x[fv[1]],x[fv[1]+1],x[fv[2]],x[fv[2]+1]};

      values[idx++] = (v[3] - v[5]);
      values[idx++] = (v[4] - v[2]);
      values[idx++] = (v[5] - v[1]);
      values[idx++] = (v[0] - v[4]);
      values[idx++] = (v[1] - v[3]);
      values[idx++] = (v[2] - v[0]);
    }
  }

  return true;
}

//return the structure or values of the hessian
bool lscm_ipopt_foldfree::eval_h(Index n, const Number* x, bool new_x,
                                 Number obj_factor, Index m, const Number* lambda,
                                 bool new_lambda, Index nele_hess, Index* iRow,
                                 Index* jCol, Number* values)
{
  if(values == NULL){
    for(Index i = 0,idx = 0;i < n;++i){
      for(row_citerd it = hess_[i].begin();
          it != hess_[i].end();++it,++idx) {
        iRow[idx] = i;
        jCol[idx] = it->first;
      }
    }
    return true;
  }
  else{
    update_hessian(lambda,true);
    for(Index i = 0,idx = 0;i < n;++i){
      for(row_citerd it = hess_[i].begin();
          it != hess_[i].end();++it) {
        values[idx++] = it->second;
      }
    }
    update_hessian(lambda,false);
  }

  return true;
}

void lscm_ipopt_foldfree::finalize_solution(
    SolverReturn status,
    Index n, const Number* x, const Number* z_L, const Number* z_U,
    Index m, const Number* g, const Number* lambda,
    Number obj_value, const Ipopt::IpoptData* ip_data,
    Ipopt::IpoptCalculatedQuantities* ip_cq)
{
  //for(size_t i = 0;i < n;++i){
  check_fold_over(&x[0]);

  write_uv(output_uv_.c_str(),x);

  size_t sz = back_map_.size();
  vectord xx(2*sz);
  for(size_t i = 0,num_vert = n;i < sz;++i){
    size_t idx = 2*back_map_[i];
    if(idx < num_vert){
      xx[2*i]   = x[idx];
      xx[2*i+1] = x[idx+1];
    }
  }

  write_texture_mesh("texture.off",&xx[0]);

}

bool lscm_ipopt_foldfree::write_uv(const char *file_name, const double *uv)
{
  std::ofstream outfile(file_name);
  if (!outfile) {
    std::cerr << __FILE__ << " " << __LINE__ << " "
              << "load file " << file_name << " to write failed" << std::endl;
    return false;
  }

  for(size_t i = 0,num_vert = front_map_.size();i < num_vert;++i){
    outfile << front_map_[i]<<" ";
    outfile << uv[2*i] <<" " <<uv[2*i+1] <<" " <<0.0<<std::endl;
  }

  outfile.close();
}

void lscm_ipopt_foldfree::write_texture_mesh(const char *filename,const double *texture_uv)
{
  FILE *pf = fopen(filename,"w");
  fprintf(pf,"OFF\n%d %d 0\n",mesh_->num_verts_,mesh_->num_faces_);
  for(size_t i = 0;i < mesh_->num_verts_;i++){
    fprintf(pf,"%f %f 0\n",texture_uv[2*i],texture_uv[2*i+1]);
  }

  for(size_t i = 0;i < mesh_->num_faces_;i++){
    size_t *fv = mesh_->faces_[i];
    fprintf(pf,"3 %d %d %d\n",fv[0],fv[1],fv[2]);
  }
  fclose(pf);
}



size_t lscm_ipopt_foldfree::check_fold_over(const double *texture){
  size_t num_fold = 0;

  if(reverse_){
    for(size_t i = 0;i < mesh_->num_faces_;++i){
      size_t *face = mesh_->faces_[i];
      size_t fv[3] = {2*back_map_[face[0]],2*back_map_[face[1]],2*back_map_[face[2]]};

      if((texture[fv[1]]-texture[fv[0]])*(texture[fv[2]+1]-texture[fv[0]+1])
         - (texture[fv[2]]-texture[fv[0]])*(texture[fv[1]+1]-texture[fv[0]+1]) <= 0){
        ++num_fold;
      }
    }
  }
  else{
    for(size_t i = 0;i < mesh_->num_faces_;++i){
      size_t *face = mesh_->faces_[i];
      size_t fv[3] = {2*back_map_[face[0]],2*back_map_[face[1]],2*back_map_[face[2]]};

      if((texture[fv[1]]-texture[fv[0]])*(texture[fv[2]+1]-texture[fv[0]+1])
         - (texture[fv[2]]-texture[fv[0]])*(texture[fv[1]+1]-texture[fv[0]+1]) >= 0){
        ++num_fold;
      }
    }
  }


  std::cout<<"The number of fold-over triangles is:"<<num_fold<<std::endl;

  return num_fold;
}


void lscm_ipopt_foldfree::update_hessian(const double * lambda,bool is_add)
{
  for(int i = 0;i < mesh_->num_faces_;++i){
    size_t *face = mesh_->faces_[i];
    size_t fv[3] = {2*back_map_[face[0]],2*back_map_[face[1]],2*back_map_[face[2]]};

    size_t ff[6] = {fv[0],fv[0]+1,fv[1],fv[1]+1,fv[2],fv[2]+1};
    //(v[2]-v[0])*(v[5]-v[1])-(v[4]-v[0])*(v[3]-v[1])
    double lbd = (is_add) ? lambda[i] : -lambda[i];
    if(ff[0] <= ff[3]){
      hess_[ff[0]][ff[3]] += lbd;
    }
    if(ff[0] <= ff[5]){
      hess_[ff[0]][ff[5]] -= lbd;
    }

    if(ff[1] <= ff[2]){
      hess_[ff[1]][ff[2]] -= lbd;
    }
    if(ff[1] <= ff[4]){
      hess_[ff[1]][ff[4]] += lbd;
    }


    if(ff[2] <= ff[1]){
      hess_[ff[2]][ff[1]] -= lbd;
    }
    if(ff[2] <= ff[5]){
      hess_[ff[2]][ff[5]] += lbd;
    }

    if(ff[3] <= ff[0]){
      hess_[ff[3]][ff[0]] += lbd;
    }
    if(ff[3] <= ff[4]){
      hess_[ff[3]][ff[4]] -= lbd;
    }


    if(ff[4] <= ff[1]){
      hess_[ff[4]][ff[1]] += lbd;
    }
    if(ff[4] <= ff[3]){
      hess_[ff[4]][ff[3]] -= lbd;
    }


    if(ff[5] <= ff[0]){
      hess_[ff[5]][ff[0]] -= lbd;
    }
    if(ff[5] <= ff[2]){
      hess_[ff[5]][ff[2]] += lbd;
    }
  }
}


bool ipopt_foldfree_parameterization(const char *file_mesh,
                                     const char *file_fix_vert,
                                     const char *output_uv)
{
  SmartPtr<TNLP> lscm = new lscm_ipopt_foldfree(file_mesh,file_fix_vert,output_uv);


  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

  // Change some options
  // Note: The following choices are only examples, they might not be
  //       suitable for your optimization problem.
  app->Options()->SetNumericValue("tol", 1e-6);
  app->Options()->SetStringValue("mu_strategy", "adaptive");


  // Intialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
    return false;
  }

  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(lscm);

  if (status == Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;

    return true;
  }
  else {
    std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
    return false;
  }

  return true;
}

}
