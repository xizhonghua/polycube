#include "util.h"
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/linear_solver/linear_solver.h>
#include <zjucad/matrix/io.h>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseCholesky>

#include "../vol_param/descriptor/func_terms/arap.h"

int convert_rot_to_euler_zyx(const double *rot,
                             double *zyx)
{
  using namespace zjucad::matrix;
  const double PI = std::atan(1.0)*4;

  assert(rot && zyx);
  itr_matrix<const double*> R(3,3,rot);
  itr_matrix<double*> zyx_(3,2,zyx);

  double & phi1  = zyx_(0,0);
  double & phi2= zyx_(0,1);

  double & theta1 = zyx_(1,0);
  double & theta2 = zyx_(1,1);

  double & varphi1 = zyx_(2,0);
  double & varphi2 = zyx_(2,1);

  if(fabs(R(2,0) - 1) > 0 && fabs(R(2,0) + 1) > 0){
      theta1 = -1 * asin(std::min(1.0,std::max(-1.0,R(2,0))));
      theta2 = PI - theta1;
      varphi1 = atan2(R(2,1)/cos(theta1),R(2,2)/cos(theta1));
      varphi2 = atan2(R(2,1)/cos(theta2),R(2,2)/cos(theta2));

      phi1 = atan2(R(1,0)/cos(theta1),R(0,0)/cos(theta1));
      phi2 = atan2(R(1,0)/cos(theta2),R(0,0)/cos(theta2));
    }else{
      phi1 = phi2 = 0; // phi can be anything, so set them to be 0
      if(fabs(R(2,0) + 1) < 1e-8){
          theta1 = theta2 = PI/2;
          varphi2 = varphi1 = phi1 + atan2(R(0,1),R(0,2));
        }else{
          theta1 = theta2 = -1.0*PI/2;
          varphi2 = varphi1 = -1.0*phi1 + atan2(-1*R(0,1),-1*R(0,2));
        }
    }
  return 0;
}

int convert_euler_angle_to_rot(const double*zyx,double *rot)
{
  using namespace zjucad::matrix;
  assert(zyx && rot);
  const double & alpha = zyx[0];
  const double & belta = zyx[1];
  const double & gamma = zyx[2];

  itr_matrix<double*> rot_(3,3,rot);

  rot_(0,0) = cos(alpha) * cos(belta);
  rot_(1,0) = sin(alpha) * cos(belta);
  rot_(2,0) = -1 * sin(belta);

  rot_(0,1) = cos(alpha) * sin(belta) * sin(gamma) - sin(alpha) * cos(gamma);
  rot_(1,1) = sin(alpha) * sin(belta) * sin(gamma) + cos(alpha) * cos(gamma);
  rot_(2,1) = cos(belta) * sin(gamma);

  rot_(0,2) = cos(alpha) * sin(belta) * cos(gamma) + sin(alpha) * sin(gamma);
  rot_(1,2) = sin(alpha) * sin(belta) * cos(gamma) - cos(alpha) * sin(gamma);
  rot_(2,2) = cos(belta) * cos(gamma);
  return 0;
}

int convert_quat_to_euler_angle(const double *q, double *zyx)
{
  using namespace zjucad::matrix;
  assert(q && zyx);
  zyx[0] = atan2(2*(q[0] *q[3] + q[1] * q[2]),1-2*(q[2] * q[2] + q[3] * q[3]));
  zyx[1] = asin(std::min(1.0,std::max(-1.0,2*(q[0] * q[2] - q[3] * q[1]))));
  zyx[2] = atan2(2*(q[0] * q[1] + q[2] * q[3]),1-2*(q[1] * q[1] + q[2] * q[2]));
  return 0;
}

int convert_euler_angle_to_quat(const double *zyx, double *quat)
{
  using namespace zjucad::matrix;
  assert(quat && zyx);
  const double &z= zyx[0];
  const double &y= zyx[1];
  const double &x= zyx[2];

  quat[0] = cos(x/2)*cos(y/2)*cos(z/2)+sin(x/2)*sin(y/2)*sin(z/2);
  quat[1] = sin(x/2)*cos(y/2)*cos(z/2)-cos(x/2)*sin(y/2)*sin(z/2);
  quat[2] = cos(x/2)*sin(y/2)*cos(z/2)+sin(x/2)*cos(y/2)*sin(z/2);
  quat[3] = cos(x/2)*cos(y/2)*sin(z/2)-sin(x/2)*sin(y/2)*cos(z/2);

  return 0;
}

int from_angle_to_rotation_matrix(const double &angle,
                                  const matrixd &axis_in,
                                  matrixd &rot)
{
  using namespace zjucad::matrix;
  matrixd axis = axis_in / norm(axis_in);

  assert(rot.size() == 9); //3 * 3 matrix
  assert(axis.size() == 3);

  //angle *= -1.0;
  const double fHalfAngle = 0.5f*angle;
  const double fSin = sinf(fHalfAngle);
  matrixd quaternion(4,1); // x,y,z,w

  quaternion[0] = fSin*axis[0];
  quaternion[1] = fSin*axis[1];
  quaternion[2] = fSin*axis[2];
  quaternion[3] = cosf(fHalfAngle);

  const double &X = quaternion[0];
  const double &Y = quaternion[1];
  const double &Z = quaternion[2];
  const double &W = quaternion[3];

  rot(0,0) = 1.0f - 2.0f*Y*Y - 2.0f*Z*Z;
  rot(1,0) = 2.0f*X*Y + 2.0f*Z*W;
  rot(2,0) = 2.0f*X*Z - 2.0f*Y*W;

  rot(0,1) = 2.0f*X*Y - 2.0f*Z*W;;
  rot(1,1) = 1.0f - 2.0f*X*X - 2.0f*Z*Z;
  rot(2,1) = 2.0f*Z*Y + 2.0f*X*W;

  rot(0,2) = 2.0f*X*Z + 2.0f*Y*W;
  rot(1,2) = 2.0f*Z*Y - 2.0f*X*W;
  rot(2,2) = 1.0f - 2.0f*X*X - 2.0f*Y*Y;

  return 0;
}

void convert_quat_to_rotation_matrix(const double * quat,
                                     matrixd & rot)
{
  assert(quat);
  const double &qx = quat[0];
  const double &qy = quat[1];
  const double &qz = quat[2];
  const double &qw = quat[3];
  if(rot.size(1) != rot.size(2))
    rot.resize(3,3);
  rot(0,0) = 1 - 2*qy*qy - 2*qz*qz;
  rot(0,1) = 2*qx*qy - 2*qz*qw;
  rot(0,2) = 2*qx*qz + 2*qy*qw;
  rot(1,0) = 2*qx*qy + 2*qz*qw;
  rot(1,1) = 1 - 2*qx*qx - 2*qz*qz;
  rot(1,2) = 2*qy*qz - 2*qx*qw;
  rot(2,0) = 2*qx*qz - 2*qy*qw;
  rot(2,1) = 2*qy*qz + 2*qx*qw;
  rot(2,2) = 1 - 2*qx*qx - 2*qy*qy;
}

void convert_rotation_matrix_to_quat(const matrixd &rot, double *quat)
{
  assert(rot.size(1) == 3);
  assert(rot.size(2) == 3);
  assert(quat);

  double & x = quat[0];
  double & y = quat[1];
  double & z = quat[2];
  double & w = quat[3];

  w = sqrt(1.0+rot(0,0)+rot(1,1)+rot(2,2))/2.0;
  if(w < 1e-6){
      const double c = sqrt(rot(0,1)*rot(0,1)*rot(0,2)*rot(0,2)
                            + rot(0,1)*rot(0,1)*rot(1,2)*rot(1,2)
                            + rot(0,2)*rot(0,2)*rot(1,2)*rot(1,2));
      x = rot(0,2)*rot(0,1)/c;
      y = rot(0,1)*rot(1,2)/c;
      z = rot(0,2)*rot(1,2)/c;
    }else{
      x = (rot(2,1)-rot(1,2))/(4*w);
      y = (rot(0,2)-rot(2,0))/(4*w);
      z = (rot(1,0)-rot(0,1))/(4*w);
    }
}

int long_number_add(std::deque<size_t> & data, size_t  a)
{
  if(data.empty()){
      while(a > 0){
          data.push_back(a%10);
          a /= 10;
        }
      return 0;
    }else{
      data.front() += a;
      for(size_t i = 0; i < data.size()-1; ++i){
          if(data[i] > 9){
              data[i+1] += data[i]/10;
              data[i] = data[i]%10;
            }else
            break;
        }
      while(data.back() > 9){
          data.push_back(data.back()/10);
        }
      return 0;
    }
}

template <typename VAL_TYPE, typename INT_TYPE>
inline void add_to_diag(hj::sparse::csc<VAL_TYPE, INT_TYPE> &A, VAL_TYPE eps)
{
  for(INT_TYPE ci = 0; ci < A.size(2); ++ci) {
      INT_TYPE nzi = A.ptr()[ci];
      for(; nzi < A.ptr()[ci+1]; ++nzi) {
          if(A.idx()[nzi] == ci) {
              A.val()[nzi] += eps;
              break;
            }
        }
      if(nzi == A.ptr()[ci+1]) {
          std::cout << "bad add to diag." << std::endl;
        }
    }
}


void pow_method_for_smallest_eigenvector(
    const hj::sparse::csc<double, int32_t> & A,
    zjucad::matrix::matrix<double> & smallest_eig,
    const size_t max_iteration)
{
  hj::sparse::csc<double, int32_t> A_tmp(A);
  boost::property_tree::ptree pt;
//  pt.put("linear_solver/type.value","direct");
//  pt.put("linear_solver/name.value", "cholmod");

  std::shared_ptr<linear_solver> ls(linear_solver::create(&A_tmp.val()[0], &A_tmp.idx()[0],
      &A_tmp.ptr()[0], A_tmp.nnz(), A_tmp.size(1), A_tmp.size(2), pt));
  if(!ls.get()){
      size_t try_times = 0;
      while(!ls.get() && try_times <= 5){
          std::cerr << "# [info] add to dig of ATA, try: " << try_times << " " << std::endl;
          add_to_diag(A_tmp,1e-1*pow(2.0,try_times));
          ls.reset(linear_solver::create(&A_tmp.val()[0], &A_tmp.idx()[0],
              &A_tmp.ptr()[0], A_tmp.nnz(), A_tmp.size(1), A_tmp.size(2), pt));
          ++try_times;
        }
      if(try_times > 5){
          throw std::logic_error("[info] ATA condition number too bad, exceed try times");
        }
    }

  zjucad::matrix::matrix<double> uk = zjucad::matrix::ones<double>(A_tmp.size(1),1);
  zjucad::matrix::matrix<double> vk = uk;
  for(size_t k = 0; k < max_iteration; ++k){
      ls->solve(&uk[0], &vk[0], 1, pt);
      uk = vk/zjucad::matrix::max(zjucad::matrix::fabs(vk));
    }
  smallest_eig = uk/zjucad::matrix::norm(uk);
}


static double find_max_abs(const Eigen::VectorXcd &a){
  double max = -1*std::numeric_limits<double>::infinity();
  for(size_t i = 0; i < a.rows(); ++i){
      const std::complex<double> &v = a[i];
      const double norm_c = v.real()*v.real()+v.imag()*v.imag();
      if( norm_c > max){
          max = norm_c;
        }
    }
  return sqrt(max);
}

void pow_method_for_smallest_eigenvector(
    const hj::sparse::csc<std::complex<double>, int32_t> &A,
    zjucad::matrix::matrix<std::complex<double> > &eigen_vector,
    const size_t max_iteration)
{
  assert(A.size(1) == A.size(2));
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> A_(A.size(1), A.size(2));
  for(size_t ci = 0; ci < A.ptr().size()-1; ++ci){
      for(size_t off = A.ptr()[ci]; off != A.ptr()[ci+1]; ++off){
          A_.coeffRef(A.idx()[off], ci) = A.val()[off];
        }
    }

  for(size_t ci = 0; ci != A.size(1); ++ci){
      A_.coeffRef(ci, ci) += std::complex<double>(1e-8);
    }

  Eigen::SimplicialCholesky<Eigen::SparseMatrix<std::complex<double> > > chol(A_);
  Eigen::ComputationInfo info = chol.info();
  if(info == Eigen::NumericalIssue){
      throw std::logic_error("matrix is not SPD.");
    }

  Eigen::VectorXcd uk(A.size(1));
  for(size_t i = 0; i < A.size(1); ++i) uk[i] = std::complex<double>(1.0,0.0);
  Eigen::VectorXcd vk = uk;
  for(size_t k = 0; k < max_iteration; ++k){
      vk = chol.solve(uk);
      uk = vk/find_max_abs(vk);
    }
  eigen_vector.resize(A.size(1),1);
  for(size_t i = 0; i < A.size(1); ++i)
    eigen_vector[i] = uk[i];
}

void solve_smallest_eigenvalue_problem(const hj::sparse::csc<double, int32_t> &A,
                                       zjucad::matrix::matrix<double> &eigen_vector)
{
  assert(A.size(1) == A.size(2));
  pow_method_for_smallest_eigenvector(A,eigen_vector,20);

  {
    zjucad::matrix::matrix<double> Ax(eigen_vector.size(),1);
    hj::sparse::mv(false, A, eigen_vector, Ax);
    const double lambda = Ax[0]/eigen_vector[0];
    std::cerr << "# [info] lambda , min_energy " << lambda << " " << lambda*lambda*zjucad::matrix::norm(eigen_vector) << std::endl;
    std::cerr << "# [info] eigenvector differnce " << zjucad::matrix::norm(Ax-lambda*eigen_vector) << std::endl;
  }
}

void solve_smallest_eigenvalue_problem(const hj::sparse::csc<std::complex<double>, int32_t> &A,
                                       zjucad::matrix::matrix<std::complex<double> > &eigen_vector)
{
  assert(A.size(1) == A.size(2));
  pow_method_for_smallest_eigenvector(A,eigen_vector,20);
}


void cal_deformation_gradient(const zjucad::matrix::matrix<size_t> & mesh,
                              const zjucad::matrix::matrix<double> & node0,
                              const zjucad::matrix::matrix<double> & node1,
                              zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & deform_gradient)
{
  using namespace zjucad::matrix;
  matrix<double> tet_node0, tet_node1;
  matrix<double> grad_op(4,3);
  deform_gradient.resize(mesh.size(2),1);
  for(size_t ti = 0; ti < mesh.size(2); ++ti){
      tet_node0 = node0(colon(), mesh(colon(),ti));
      tet_node1 = node1(colon(), mesh(colon(),ti));
      calc_tet_def_grad_op(&tet_node0[0], &grad_op[0]);
      deform_gradient[ti] = tet_node1*grad_op;
    }
}
