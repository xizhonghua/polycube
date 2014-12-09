#include "Application2.h"
#include <fstream>

namespace DDG
{
  void Application2::compute_mean_curvature(DenseMatrix<Real> & rho, const DenseMatrix<Real> &x)
  {
    DenseMatrix<Real> HN = star0.inverse()*laplace * x;
    DenseMatrix<Real> H(HN.nRows(),1);
    for(int i = 0; i < H.nRows(); ++i)
      H(i,0) = 0.5*sqrt(HN(i,0)*HN(i,0)+HN(i,1)*HN(i,1)+HN(i,2)*HN(i,2));
    rho = H;

    N = HN;
    for(int i = 0; i < N.nRows(); ++i){
        const Real len = sqrt(N(i,0)*N(i,0)+N(i,1)*N(i,1)+N(i,2)*N(i,2));
        for(int j = 0; j < 3; ++j)
          N(i,j) /= len;
      }
  }

  void Application2::compute_flow_direction(Mesh &mesh, DenseMatrix<Real> &drho, const DenseMatrix<Real> &rho)
  {
    drho = rho;
    drho *= -1.0;
  }

  void Application2::build_constraint_basis(vector<DenseMatrix<Real> > &c, const Mesh &mesh)
  {
    vector<DenseMatrix<Real> > axis;
    DenseMatrix<Real> ones(mesh.vertices.size(),1);
    ones.zero(1.0);
    axis.push_back(ones);

    for(size_t i = 0; i < 3; ++i){
        DenseMatrix<Real> Ni = ones;
        for(int j = 0; j < N.nRows(); ++j)
          Ni(j,0) = N(j,i);
        axis.push_back(Ni);
      }

    gram_schmidt(axis, c);
  }

  void Application2::gram_schmidt(const vector<DenseMatrix<Real> > &axis, vector<DenseMatrix<Real> > &c)
  {
    c.resize(axis.size());
    c[0] = axis.front();
    c[0].normalize();

    for(size_t i = 1; i < axis.size(); ++i){
        c[i] = axis[i];
        for(size_t j = 0; j < i; ++j){
            DenseMatrix<Real> tmp = c[j];
            tmp *= inner(axis[i],c[j]);
            c[i] -= tmp;
          }
        c[i].normalize();
      }
  }

  void Application2::project_flow_onto_constraint(DenseMatrix<Real> &drho, const vector<DenseMatrix<Real> > &c)
  {
    DenseMatrix<Real> drho_new = drho;
    DenseMatrix<Real> new_c;
    for(size_t ci = 0; ci < c.size(); ++ci){
        new_c = c[ci];
        new_c *= inner(drho,c[ci]);
        drho_new -= new_c;
      }
    drho = drho_new;
  }

  void Application2::update_rho(const Real &step, DenseMatrix<Real> &rho, DenseMatrix<Real> &drho)
  {
    drho *= step;
    rho += drho;
  }

  void Application2::recover_tangents(const Mesh &mesh,
                                      const DenseMatrix<Real> &rho,
                                      DenseMatrix<Real> &lambda)
  {
    SparseMatrix<Real> D(4*rho.nRows(), 4*rho.nRows());

    for(FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f){
        Real area = f->area();
        HalfEdgeCIter he = f->he;
        do{
            int pj = he->vertex->index;
            int pi = he->flip->vertex->index;
            int pk = he->next->vertex->index;

            Vector ei = mesh.vertices[pk].position - mesh.vertices[pj].position;
            Vector ej = mesh.vertices[pi].position - mesh.vertices[pk].position;

            Quaternion qei(ei);
            Quaternion qej(ej);
            Quaternion xij = -0.25*qei * qej/area +
                (rho(pi,0)*qej-rho(pj,0)*qei)/6. +
                area*(rho(pi,0)*rho(pj,0))/9.;

            Real a = xij.re();
            const Vector &bcd = xij.im();
            const Real & b = bcd[0]; const Real & c = bcd[1]; const Real & d = bcd[2];

            D(4*pi+0,4*pj+0) = a;   D(4*pi+0,4*pj+1) = -b;  D(4*pi+0,4*pj+2) =-c;   D(4*pi+0,4*pj+3) = -d;
            D(4*pi+1,4*pj+0) = b;   D(4*pi+1,4*pj+1) = a;   D(4*pi+1,4*pj+2) =-d;   D(4*pi+1,4*pj+3) = c;
            D(4*pi+2,4*pj+0) = c;   D(4*pi+2,4*pj+1) = d;   D(4*pi+2,4*pj+2) = a;   D(4*pi+2,4*pj+3) = -b;
            D(4*pi+3,4*pj+0) = d;   D(4*pi+3,4*pj+1) = -c;  D(4*pi+3,4*pj+2) = b;   D(4*pi+3,4*pj+3) = a;

          }while(he != f->he);
      }

    smallestEig(D,lambda);
  }

  void Application2::recover_positions(const Mesh &mesh,
                                       DenseMatrix<Real> &x, DenseMatrix<Real> &lambda)
  {
    DenseMatrix<Real> New_T(mesh.edges.size(),3);
    for(EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e){
        int eidx = e->index;
        int pj = e->he->vertex->index;
        int pi = e->he->flip->vertex->index;
        Vector eij = mesh.vertices[pj].position-mesh.vertices[pi].position;
        Quaternion qeij(eij);
        Quaternion qi(lambda(4*pi,0), lambda(4*pi,1), lambda(4*pi,2), lambda(4*pi,3));
        Quaternion qj(lambda(4*pj,0), lambda(4*pj,1), lambda(4*pj,2), lambda(4*pj,3));
        Quaternion tmp = qi.conj()*qeij*qi/3.0+qi.conj()*qeij*qj/6.0
            +qj.conj()*qeij*qi/6.+qj.conj()*qeij*qj/3.;

        const Vector &new_eij = tmp.im();
        New_T(eidx,0) = new_eij[0];
        New_T(eidx,1) = new_eij[1];
        New_T(eidx,2) = new_eij[2];
      }
    DenseMatrix<Real> b = d0.transpose() * star1 * New_T;

    DenseMatrix<Real> b0(b.nRows(),1), b1(b.nRows(),1), b2(b.nRows(),1);
    DenseMatrix<Real> x0(b.nRows(),1), x1(b.nRows(),1), x2(b.nRows(),1);

    for(int i = 0; i < b.nRows(); ++i){
        b0(i,0) = b(i,0);
        b1(i,0) = b(i,1);
        b2(i,0) = b(i,2);
      }

    solvePositiveDefinite(laplace, x0, b0);
    solvePositiveDefinite(laplace, x1, b1);
    solvePositiveDefinite(laplace, x2, b2);

    for(int i = 0; i < b.nRows(); ++i){
        x(i,0) = x0(i,0);
        x(i,1) = x1(i,0);
        x(i,2) = x2(i,0);
      }
  }
}
