#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <hjlib/math/polar.h>
#include "connection_generator.h"

namespace fxz {

#define CALL_FUNC(pn) std::cout << "# [ Run: " << #pn << " ]" << std::endl; \
  if(pn){ std::cerr << "# [ ERROR ] " << #pn << std::endl; return __LINE__;}

using zjucad::matrix::colon;

int connection_generator::run()
{
  // build fxz::tet_mesh
  mesh_ptr_.reset(new tet_mesh(tm_.tetmesh_.node_, tm_.tetmesh_.mesh_));

  if (size_field_.size() != mesh_ptr_->verts_num()) {
    std::cerr << "# size field should have the same size with verts num, "
              <<"but the size field size is " << size_field_.size() << std::endl;
    return 1;
  }
  CALL_FUNC(mesh_ptr_->build());
  CALL_FUNC(cal_tets_inv_grad());
  CALL_FUNC(cal_size_grad());
  CALL_FUNC(cal_faces_connection());

  return 0;
}

int connection_generator::cal_tets_inv_grad()
{
  size_t tets_num = mesh_ptr_->tets_num();
  grad_inv_.resize(tets_num);
  const matrixd& verts = mesh_ptr_->verts_coord();  
  for (size_t ti = 0; ti < tets_num; ++ti) {
    std::vector<size_t> tvs = mesh_ptr_->tet_verts(ti);
    matrixd p3(3,3);
    p3(colon(),0) = verts(colon(),tvs[1])-verts(colon(),tvs[0]);
    p3(colon(),1) = verts(colon(),tvs[2])-verts(colon(),tvs[0]);
    p3(colon(),2) = verts(colon(),tvs[3])-verts(colon(),tvs[0]);
    zjucad::matrix::inv(p3);
    grad_inv_[ti] = p3;
  }

  return 0;
}

int connection_generator::cal_size_grad()
{
  size_t tets_num = mesh_ptr_->tets_num();
  size_grad_.resize(tets_num);
  const matrixst& tets = mesh_ptr_->tets_verts();
  for (size_t ti = 0; ti < tets_num; ++ti) {
    matrixst vids = tets(colon(), ti);
    matrixd pss(1,3);
    for (size_t i=1; i<4; ++i) {
      pss[i-1] = size_field_[vids[i]]-size_field_[vids[0]];
    }
    size_grad_[ti] = pss*grad_inv_[ti];
  }

  return 0;
}

int connection_generator::cal_faces_connection()
{
  size_t faces_num = mesh_ptr_->direc_faces_num();
  const matrixd& verts = mesh_ptr_->verts_coord();  
  for (size_t fi = 0; fi < faces_num; fi+=2) {
    if (mesh_ptr_->is_surf_face(fi)) continue;
    std::vector<size_t> fvs = mesh_ptr_->face_verts(fi);
    
    std::pair<size_t,size_t> t2 = mesh_ptr_->face_adj_tets(fi);
    if (t2.first > t2.second) std::swap(t2.first, t2.second);
    
    matrixd tet_cen(3,1);
    matrixd face_cen = zjucad::matrix::zeros<double>(3,1);
    matrixd faces_p(3,3);
    for (size_t i=0; i<3; ++i) {
      face_cen += verts(colon(),fvs[i]);
    }
    face_cen /= 3.0;
    
    matrixd conn1;
    tet_cen = mesh_ptr_->tet_center(t2.first);    
    cal_W(t2.first, tet_cen, face_cen, conn1);
    matrixd conn2;
    tet_cen = mesh_ptr_->tet_center(t2.second);    
    cal_W(t2.second, tet_cen, face_cen, conn2);

    inv(conn1); // maybe need to modify
    matrixd R = conn1*conn2;
    hj::polar3d pol; pol(R,2);
    faces_conn_[t2] = R;
    faces_conn_[std::make_pair(t2.second,t2.first)] = trans(R);
  }

  return 0;
}

int connection_generator::cal_W(size_t tid, const matrixd& p0,
                                const matrixd& p1, matrixd& conn)
{
  double size_p0 = cal_point_size(tid, p0);
  double size_p1 = cal_point_size(tid, p1);

  matrixd e = p1 - p0;
  assert(zjucad::matrix::norm(e) > 1e-8);
  e /= zjucad::matrix::norm(e); // normalize
  e *= -1.0/(size_p0 + size_p1);
  e = zjucad::matrix::temp(trans(e));
  
  matrixd H0(3,3), H1(3,3), H2(3,3);
  const matrixd& sg = size_grad_[tid];
  
  H0(0,0) = 0.0; H0(0,1) =  sg[1]; H0(0,2) = sg[2];
  H0(1,0) = 0.0; H0(1,1) = -sg[0]; H0(1,2) =   0.0;
  H0(2,0) = 0.0; H0(2,1) =    0.0; H0(2,2) =-sg[0];

  H1(0,0) = -sg[1]; H1(0,1) = 0.0; H1(0,2) =   0.0;
  H1(1,0) =  sg[0]; H1(1,1) = 0.0; H1(1,2) = sg[2];  
  H1(2,0) =    0.0; H1(2,1) = 0.0; H1(2,2) =-sg[1];

  H2(0,0) = -sg[2]; H2(0,1) =    0.0;  H2(0,2) = 0.0;
  H2(1,0) =    0.0; H2(1,1) = -sg[2];  H2(1,2) = 0.0;
  H2(2,0) =  sg[0]; H2(2,1) =  sg[1];  H2(2,2) = 0.0;

  matrixd W(3,3);
  W(0,colon()) = e*H0;
  W(1,colon()) = e*H1;
  W(2,colon()) = e*H2;
  const double EPS = 1e-8;
  assert(fabs(W(0,0))<EPS && fabs(W(1,1))<EPS && fabs(W(2,2))<EPS);
  assert(fabs(W(0,1)+W(1,0))<EPS && fabs(W(0,2)+W(2,0))<EPS && fabs(W(1,2)+W(2,1))<EPS);

  // e^W => connection
  matrixd u(3,1);
  u[0] = -W(1,2); u[1] = W(0,2); u[2] = -W(0,1);
  double ang = zjucad::matrix::norm(u);
  if (fabs(ang) < EPS) {
    u[0] = u[1] = 0.0; u[2] = 1.0;
    //std::cerr << "-- Appear a zero rotation axis." << std::endl;
  } else u /= ang; // (x,y,z), x^2+y^2+z^2 = 1
  double cs = cos(ang*0.5), ss = sin(ang*0.5);
  conn.resize(3,3);
  conn(0,0) = 2.0*(u[0]*u[0]-1.0)*ss*ss+1.0;
  conn(0,1) = 2.0*u[0]*u[1]*ss*ss-2.0*u[2]*cs*ss;
  conn(0,2) = 2.0*u[0]*u[2]*ss*ss+2.0*u[1]*cs*ss;
  conn(1,0) = 2.0*u[0]*u[1]*ss*ss+2.0*u[2]*cs*ss;
  conn(1,1) = 2.0*(u[1]*u[1]-1.0)*ss*ss+1.0;
  conn(1,2) = 2.0*u[1]*u[2]*ss*ss-2.0*u[0]*cs*ss;
  conn(2,0) = 2.0*u[0]*u[2]*ss*ss-2.0*u[1]*cs*ss;
  conn(2,1) = 2.0*u[1]*u[2]*ss*ss+2.0*u[0]*cs*ss;
  conn(2,2) = 2.0*(u[2]*u[2]-1.0)*ss*ss+1.0;
  assert(norm(trans(conn)*conn-zjucad::matrix::eye<double>(3))<EPS);
  conn *= sqrt(size_p0/size_p1);
  
  return 0;
}

double connection_generator::cal_point_size(size_t tid, const matrixd& p)
{
  size_t first_vid = mesh_ptr_->tets_verts()(0,tid);
  matrixd P = p - mesh_ptr_->verts_coord()(colon(),first_vid);
  matrixd r = size_grad_[tid]*P;
  return (r[0]+size_field_[first_vid]);
}

} // namespace fxz

