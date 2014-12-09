#include <fstream>

#include <hjlib/math/polar.h>
#include <hjlib/math/blas_lapack.h>

#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/mesh.h>

#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/lapack.h>
#include <numeric>

#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"

using namespace std;
using namespace zjucad::matrix;

//////////////////////  interface  /////////////////////////////
////////////////////////////////////////////////////////////////
int defom_quality(const jtf::mesh::meshes &init_tm, const jtf::mesh::meshes &def_tm,
                  const char * output_vtk);

int defom_quality_raw(const matrix<double>& init_node,
                      const matrix<double>& def_node,
                      const matrix<size_t> & tet,
                      const char * output_vtk);
/////////////////////////////////////////////////////////////////

///////////////  calculate arap distortion //////////////////////
int calc_tet_def_grad_op(const double *tet, double *grad_op)
{
  const static double node2edge[] = {
    1, 0, 0, -1,
    0, 1, 0, -1,
    0, 0, 1, -1
  };
  const static itr_matrix<const double *> E(4, 3, node2edge);

  itr_matrix<const double *> T(3, 4, tet);
  itr_matrix<double *> G(4, 3, grad_op);
  matrix<double> A = T*E;
  if(inv(A))
    return __LINE__;
  G = E*A;
  return 0;
}

int calc_tet_arap_distortion(const double *orig_tet, const double * new_tet,
                             double *dis)
{
  matrix<double> grad_op = zeros(4,3) ;
  itr_matrix<const double *> orig_T(3, 4, orig_tet);
  itr_matrix<const double *> new_T(3, 4, new_tet);
  if(calc_tet_def_grad_op(&orig_T[0], &grad_op[0]))
    cerr << "# degenerated tet" << endl;
  itr_matrix<double*> f(3,3,dis);

  f = new_T * grad_op;
  matrix<double> R = f;
  hj::polar3d p;
  p(R,2);
  f -= R;
  return 0;
}

int calc_tri_arap_distortion(const double *orig_tri, const double * new_tri,
                             double *dis)
{
  itr_matrix<double *> dis_m(3,3,dis);
  itr_matrix<const double *> orig_tri_m(3,3, orig_tri);
  itr_matrix<const double *> new_tri_m(3,3, new_tri);
  dis_m = orig_tri_m;
  matrix<double > B = orig_tri_m;
  if(inv(B)) return __LINE__;
  //grad = (f0,f1,f2)*(p0,p1,p2)^-1
  matrix<double> A = new_tri_m;
  dis_m = temp(A * B);

  matrix<double> R = dis_m;
  hj::polar3d p;
  p(R,2);
  dis_m -= R;
  return 0;
}

/////////////////////////////////////////////////////////////////////


int deform_quality_raw(const matrix<double>& init_node,
                       const matrix<double>& def_node,
                       const matrix<size_t>& init_cell,
                       const matrix<size_t>& def_cell,
                       matrix<double> & arap_dis)
{
  assert(init_cell.size(1) == 3 || init_cell.size(1) == 4);
  assert(init_cell.size() == def_cell.size());

  arap_dis = zeros<double>(init_cell.size(2),1);
  matrix<double> dis = zeros<double>(3,3);
  matrix<double> init_cell_node, def_cell_node;
  for(size_t ti = 0; ti < arap_dis.size(); ++ti){

      init_cell_node = init_node(colon(), init_cell(colon(),ti));
      def_cell_node = def_node(colon(), def_cell(colon(),ti));

      if(init_cell.size(1) == 4)
        calc_tet_arap_distortion(&init_cell_node[0], &def_cell_node[0], &dis[0]);
      else if(init_cell.size(1) == 3)
        calc_tri_arap_distortion(&init_cell_node[0], &def_cell_node[0], &dis[0]);

      arap_dis[ti] = norm(dis) * norm(dis)/2.0;
    }

  cerr << "# [arap_dis] min/avg/max "
       << *min_element(arap_dis.begin(),arap_dis.end()) << "/"
       << std::accumulate(arap_dis.begin(),arap_dis.end(),0.0) /arap_dis.size()<< "/"
       << *max_element(arap_dis.begin(),arap_dis.end()) << endl;


  return 0;
}

int vis_tet_arap_distortion(int argc, char * argv[])
{
  if(argc != 5){
      cerr << "# [usage] vis_arap_distortion tet/tri init_cell def_cell output_vtk." << endl;
      return __LINE__;
    }

  const string cell_type = argv[1];
  if(cell_type != "tet") return __LINE__;

  jtf::mesh::meshes init_tm, def_tm;

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &init_tm.node_, &init_tm.mesh_))
    return __LINE__;

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[3], &def_tm.node_, &def_tm.mesh_))
    return __LINE__;

  matrix<double> arap_dis(init_tm.mesh_.size(2),1);
  deform_quality_raw(init_tm.node_, def_tm.node_, init_tm.mesh_, def_tm.mesh_, arap_dis);

  string output_name = argv[4];
  ofstream ofs_init((output_name + ".init_node.vtk").c_str());
  ofstream ofs_def((output_name + ".def_node.vtk").c_str());

  if(ofs_init.fail() || ofs_def.fail()){
      cerr << "# [error] can not open output vtk file." << endl;
      return __LINE__;
    }

  tet2vtk(ofs_init, &init_tm.node_[0], init_tm.node_.size(2),
      &init_tm.mesh_[0], init_tm.node_.size(2));
  cell_data(ofs_init, &arap_dis[0], arap_dis.size(), "arap_distortion");

  tet2vtk(ofs_def, &def_tm.node_[0], def_tm.node_.size(2),
      &def_tm.mesh_[0], def_tm.mesh_.size(2));
  cell_data(ofs_def, &arap_dis[0], arap_dis.size(), "arap_distortion");

  ofstream ofs_l("arap_dis_txt");
  for(size_t ti = 0; ti < arap_dis.size(); ++ti) ofs_l << arap_dis[ti] << endl;

  matrix<double> tet_vol(arap_dis.size(),1);

  for(size_t ti = 0; ti < init_tm.mesh_.size(2); ++ti){
      tet_vol[ti] = jtf::mesh::cal_tet_vol(init_tm.node_(colon(), init_tm.mesh_(colon(),ti)));
    }
  const double total_vol = std::accumulate(tet_vol.begin(), tet_vol.end(), 0.0);
  cerr << "# [info] vol_weighted_arap " << dot(tet_vol, arap_dis)/total_vol << endl;

  return 0;
}

int vis_tri_arap_distortion(int argc, char * argv[])
{
  if(argc != 5){
      cerr << "# [usage] vis_arap_distortion tet/tri init_cell def_cell output_vtk." << endl;
      return __LINE__;
    }

  const string cell_type = argv[1];
  if(cell_type != "tri") return __LINE__;

  jtf::mesh::meshes init_trimesh, def_trimesh;

  if(jtf::mesh::load_obj(argv[2], init_trimesh.mesh_, init_trimesh.node_))
    return __LINE__;

  if(jtf::mesh::load_obj(argv[3], def_trimesh.mesh_, def_trimesh.node_))
    return __LINE__;

  matrix<double> arap_dis;
  deform_quality_raw(init_trimesh.node_, def_trimesh.node_,
                     init_trimesh.mesh_, def_trimesh.mesh_, arap_dis);

  string output_name = argv[4];
  ofstream ofs_init((output_name + ".init_node.vtk").c_str());
  ofstream ofs_def((output_name + ".def_node.vtk").c_str());

  if(ofs_init.fail() || ofs_def.fail()){
      cerr << "# [error] can not open output vtk file." << endl;
      return __LINE__;
    }

  tri2vtk(ofs_init, &init_trimesh.node_[0], init_trimesh.node_.size(2),
      &init_trimesh.mesh_[0], init_trimesh.mesh_.size(2));
  cell_data(ofs_init, &arap_dis[0], arap_dis.size(), "arap_distortion");

  tri2vtk(ofs_def, &def_trimesh.node_[0], def_trimesh.node_.size(2),
      &def_trimesh.mesh_[0], def_trimesh.mesh_.size(2));
  cell_data(ofs_def, &arap_dis[0], arap_dis.size(), "arap_distortion");

  ofstream ofs_l("arap_dis_txt");
  for(size_t ti = 0; ti < arap_dis.size(); ++ti) ofs_l << arap_dis[ti] << endl;


  return 0;
}

//////////////////////////////////////////////////
// example:
//////////////////////////////////////////////////
int vis_arap_distortion(int argc, char * argv[])
{
  if(argc != 5){
      cerr << "# [usage] vis_arap_distortion tet/tri init_cell def_cell output_vtk." << endl;
      return __LINE__;
    }

  const string cell_type = argv[1];
  if(cell_type == "tet")
    return vis_tet_arap_distortion(argc, argv);
  else if(cell_type == "tri")
    return vis_tri_arap_distortion(argc, argv);
  else{
      cerr << "# [error] unknown cell type." << endl;
      return __LINE__;
    }
}
