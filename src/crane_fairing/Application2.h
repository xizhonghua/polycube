#ifndef APPLICATION2_H
#define APPLICATION2_H

/*
 * This application is for robust fairing in
 * "Robust Fairing Via Conformal Curvature Flow"
*/


#include "Mesh.h"
#include "Real.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "DiscreteExteriorCalculus.h"

#include <fstream>
#include "../common/vtk.h"

namespace DDG
{
  class Application2
  {
  public:
    void run(const Real step, Mesh& mesh)
    {
      init(mesh);

      DenseMatrix<Real> x;
      getPositions(mesh, x);

      DenseMatrix<Real> rho, drho;
      vector<DenseMatrix<Real> > c;
      compute_mean_curvature(rho, x);

      {
        vector<double> H(N.nRows());
        vector<double> xp(3*x.nRows());
        for(int i = 0 ; i < N.nRows(); ++i)
          H[i] = rho(i,0);

        vector<size_t> face(mesh.faces.size() * 3);
        mesh.write(&face[0], &xp[0]);
        std::ofstream ofs("mean_curvature.vtk");
        tri2vtk(ofs, &xp[0], xp.size()/3, &face[0], face.size()/3);
        point_data(ofs, &H[0], H.size(), "mean-curvature");
      }


      compute_flow_direction(mesh, drho, rho);
      build_constraint_basis(c, mesh);
      project_flow_onto_constraint(drho, c);
      update_rho(step, rho, drho);

      DenseMatrix<Real> lambda(4*rho.nRows(),1);
      lambda.zero(1.0);
      recover_tangents(mesh, rho, lambda);
      recover_positions(mesh, x, lambda);
      setPositions(x, mesh);
    }

  protected:
    void gram_schmidt(const vector<DenseMatrix<Real> > &axis, vector<DenseMatrix<Real> > &c);

    void compute_mean_curvature(DenseMatrix<Real> & rho, const DenseMatrix<Real> &x);

    void compute_flow_direction(Mesh &mesh, DenseMatrix<Real> &drho, const DenseMatrix<Real> &rho);

    void build_constraint_basis(vector<DenseMatrix<Real> > &c, const Mesh &mesh);

    void project_flow_onto_constraint(DenseMatrix<Real> &drho, const vector<DenseMatrix<Real> > &c);

    void update_rho(const Real &step, DenseMatrix<Real> &rho, DenseMatrix<Real> &drho);

    void recover_tangents(const Mesh &mesh, const DenseMatrix<Real> &rho, DenseMatrix<Real> &lambda);

    void recover_positions(const Mesh &mesh, DenseMatrix<Real> &x, DenseMatrix<Real> &lambda);

    void init(Mesh &mesh)
    {
      HodgeStar0Form<Real>::build(mesh, star0);
      HodgeStar1Form<Real>::build(mesh, star1);
      ExteriorDerivative0Form<Real>::build(mesh, d0);
      //laplace = star0.inverse()*d0.transpose()*star1*d0;
      laplace = d0.transpose()*star1*d0;
    }

    void getPositions(const Mesh& mesh, DenseMatrix<Real>& x) const
    {
      x = DenseMatrix<Real>( mesh.vertices.size(), 3 );
      for ( VertexCIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v ++)
        {
          for( int i = 0; i < 3; ++i)
            x(v->index, i) = v->position[i];
        }
    }

    void setPositions(const DenseMatrix<Real>& x, Mesh& mesh)
    {
      for ( VertexIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v ++)
        {
          v->position = Vector(x(v->index, 0),
                               x(v->index, 1),
                               x(v->index, 2));
        }
    }
  private:
    SparseMatrix<Real> star0, star1, d0, laplace;
    DenseMatrix<Real> N;
  };
}

#endif // APPLICATION2_H
