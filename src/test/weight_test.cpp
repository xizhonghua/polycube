#include "weight_test.h"
#include "../tetmesh/tetmesh.h"
#include "../hex_ui/function_term.h"
#include <string>
using namespace std;

void calculate_init_zyz_scale(ptree& pt, const tetmesh& tm, matrixd & residual)
{
  auto_ptr<face2tet_adjacent> fa(face2tet_adjacent::create(tm.tet));

  matrixst out_face_idx;
  get_outside_face_idx(*fa, out_face_idx);

  std::vector<double> volume_of_tets(tm.tet.size(2));
  std::vector<double> area_of_surf_tris(out_face_idx.size());

  for(size_t ti = 0; ti < volume_of_tets.size(); ++ti)
  {
    volume_of_tets[ti] = fabs(tet_volume(tm.node(colon(), tm.tet(colon(), ti))));
  }
  for(size_t ti = 0; ti < area_of_surf_tris.size(); ++ti)
  {
    area_of_surf_tris[ti] = fabs(area_of_triangle(tm.node, &fa->faces_[out_face_idx[ti]][0]));
    CPPUNIT_ASSERT(area_of_surf_tris[ti] > 1e-8);
  }
  const double volume_of_object = accumulate(volume_of_tets.begin(), volume_of_tets.end(), 0.0);
  const double area_of_surface = accumulate(area_of_surf_tris.begin(), area_of_surf_tris.end(), 0.0);
  // smoothing term
  vector<boost::shared_ptr<function> > funcs;
  for(size_t i = 0; i < fa->faces_.size(); ++i) {// for each element, i.e. inner face
    if(fa->is_outside_face(fa->face2tet_[i])) continue;

    const size_t &vi = fa->face2tet_[i].first;
    const size_t &vj = fa->face2tet_[i].second;

    matrixd bary_i = zeros<double>(3,1);
    matrixd bary_j = zeros<double>(3,1);
    for(size_t v = 0; v < 4; ++v)
    {
      bary_i += tm.node(colon(), tm.tet(v,vi));
      bary_j += tm.node(colon(), tm.tet(v,vj));
    }
    bary_i /= 4;
    bary_j /= 4;
    const double dij = norm(bary_i - bary_j);
    const double weight = sqrt((volume_of_tets[vi] + volume_of_tets[vj] )/4)*pow(volume_of_object,-1.0/6.0) /dij  ;

    funcs.push_back(boost::shared_ptr<function>(
                      new inner_smooth(tm.tet,vi,vj,*fa,weight)));
    //sqrt((volume_of_tets[vi] + volume_of_tets[vj] )/(4*volume_of_object)))));
  }

  for(size_t i = 0; i < out_face_idx.size(); ++i) {
    const size_t &vi = fa->face2tet_[out_face_idx[i]].first;
    const size_t &vj = fa->face2tet_[out_face_idx[i]].second;
    const size_t v = (vi == -1)?vj:vi;

    const double weight = sqrt(area_of_surf_tris[i]/area_of_surface * pt.get<double>("align_w.value"));

    const vector<size_t> &vert_idx = fa->faces_[out_face_idx[i]];
    const matrixd edge[2] = {
      tm.node(colon(), vert_idx[1]) - tm.node(colon(), vert_idx[0]),
      tm.node(colon(), vert_idx[2]) - tm.node(colon(), vert_idx[0]),
    };

    matrixd normal = cross(edge[0], edge[1]);

    const double len = norm(normal);

    CPPUNIT_ASSERT(len > 1e-8);
    normal /= len;

    funcs.push_back(boost::shared_ptr<function>(
                      new align_sh_inner(tm.tet.size(2), v, &normal[0], weight)));
  }

  matrixd sh(9, tm.tet.size(2));
  std::auto_ptr<hj::function::function> all_f;
  all_f.reset(new_catenated_function<double, int32_t>(funcs));

  residual.resize(all_f->dim_of_f());

  cerr << "# begin to compute" << endl;
  sh = zeros<double>(9, tm.tet.size(2));
  zjucad::optimize(*all_f, sh, residual, pt);
}

void calculate_lbfgs_zyz_scale(ptree& pt, const tetmesh& tm, matrixd & residual,matrixd &zyz)
{
  sym_frame_opt  func;
  matrixd fixed_frame, aligned;
  matrixst fixed_frame_idx, aligned_idx;
  CPPUNIT_ASSERT(!load_from_tet(tm.node, tm.tet, fixed_frame, fixed_frame_idx,
                                aligned, aligned_idx,0,0,0));

  const double weight[2] = {
    pt.get<double>("align_w.value"),
    pt.get<double>("fix_w.value")};
  matrixd stiff = ones<double>(tm.tet.size(2), 1);
  const string smooth_strategy = pt.get<string>("inner_smooth.value","face");
  func.setup_equations_inner_zyz_barycenter_adjacent_face(
        tm.tet, tm.node, fixed_frame, fixed_frame_idx,
        aligned, aligned_idx, weight, stiff,smooth_strategy);
  zyz = zeros<double>(3, func.get()->dim_of_x()/3)+3.1415926/4.0;
  zjucad::optimize(*func.get(), zyz, residual, pt);
}
void weight_test::setUp()
{
  pt.put("tet.value","../../dat@jtf/fandisk/301k/tet/fandisk-301k.tet");
  pt.put("fix_w.value",1e3);
}

void weight_test::tearDown()
{}

void weight_test::weight_init_zyz_scale_test()
{
  pt.put("prog.value","init_zyz_inner");
  pt.put("align_w.value",100);
  pt.put("package.value","hj");
  pt.put("iter.value",10);
  tetmesh tm,tm_scale;
  CPPUNIT_ASSERT(!tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm.node, &tm.tet));
  CPPUNIT_ASSERT(!tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm_scale.node, &tm_scale.tet));
  tm_scale.node *= 10;
  matrixd residual,residual_scale;
  calculate_init_zyz_scale(pt,tm,residual);
  cerr << "# finish the scale = 1 " << endl;
  calculate_init_zyz_scale(pt,tm_scale,residual_scale);
  cerr << "# finish the scale = 10 " << endl;

  CPPUNIT_ASSERT_DOUBLES_EQUAL(norm(residual),norm(residual_scale),1e-8);

  for(size_t t = 0; t < residual.size(); ++t)
  {
    CPPUNIT_ASSERT(fabs(residual[t] - residual_scale[t]) < 1e-8);
  }
}

void weight_test::weight_lbfgs_zyz_scale_test()
{
  pt.put("prog.value","frame_inner");
  pt.put("align_w.value",0.1);
  pt.put("package.value","alglib");
  pt.put("alg.value","lbfgs");
  pt.put("lbfgs-len.value",7);
  pt.put("iter.value",100);
  tetmesh tm,tm_scale;
  CPPUNIT_ASSERT(!tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm.node, &tm.tet));
  CPPUNIT_ASSERT(!tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm_scale.node, &tm_scale.tet));
  tm_scale.node *= 5;
  matrixd residual,residual_scale;
  matrixd zyz,zyz_scale;

  calculate_lbfgs_zyz_scale(pt,tm,residual,zyz);
  cerr << "# finish the scale = 1 " << endl;

  calculate_lbfgs_zyz_scale(pt,tm_scale,residual_scale,zyz_scale);
  cerr << "# finish the scale = 5 " << endl;

  CPPUNIT_ASSERT_DOUBLES_EQUAL(norm(residual),norm(residual_scale),1e-8);

  for(size_t t = 0; t < residual.size(); ++t)
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(residual[t],residual_scale[t],1e-8);
  }

}
void weight_test::weight_lbfgs_zyz_subdivision_test()
{
  pt.put("prog.value","frame_inner");
  pt.put("align_w.value",100);
  pt.put("package.value","alglib");
  pt.put("alg.value","lbfgs");
  pt.put("lbfgs-len.value",7);
  pt.put("iter.value",100);
  tetmesh tm,tm_sub;
  matrixd residual,residual_sub;
  matrixd zyz,zyz_sub;

  pt.put("tet.value","../../dat@jtf/sphere/20k/tet/sphere-20k.tet");
  CPPUNIT_ASSERT(!tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm.node, &tm.tet));

  calculate_lbfgs_zyz_scale(pt,tm,residual,zyz);

  pt.put("tet.value","../../dat@jtf/sphere/194k/tet/sphere-194k.tet");
  CPPUNIT_ASSERT(!tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm_sub.node, &tm_sub.tet));

  calculate_lbfgs_zyz_scale(pt,tm_sub,residual_sub,zyz_sub);

  cerr << "# origin residual = " << norm(residual) << endl;
  cerr << "# subdivision residual = " << norm(residual_sub) << endl;
  //    matrixd zyz_from_sub = zeros<double>(3, tm.tet.size(2));
  //    matrixd barycenter = zeros<double>(3,tm.tet.size(2));
  //    matrixd barycenter_sub = zeros<double>(3,tm_sub.tet.size(2));
  //    for(size_t t = 0; t < tm.tet.size(2); ++t)
  //    {
  //        for(size_t  v = 0; v < 4; ++v)
  //            barycenter(colon(),t) += tm.node(colon(), tm.tet(v,t));
  //    }
  //    for(size_t t = 0; t < tm_sub.tet.size(2); ++t)
  //    {
  //        for(size_t  v = 0; v < 4; ++v)
  //            barycenter_sub(colon(),t) += tm_sub.node(colon(), tm_sub.tet(v,t));
  //    }
  //    for(size_t t = 0; t < tm.tet.size(2); ++t)
  //    {
  //        for(size_t t2 = 0; t2 < tm_sub.tet.size(2); ++t2)
  //        {
  //            if(norm(barycenter(colon(), t) - barycenter_sub(colon(), t2) < 1e-4))
  //                zyz_from_sub(colon(), t) += zyz_sub(colon(),t2);
  //        }
  //    }
  //    matrixd difference = zyz - zyz_from_sub;
  //    for(size_t t = 0; t < tm.tet.size(2); ++t)
  //    {
  //        CPPUNIT_ASSERT(norm(difference(colon(), t)) < 1e-4);
  //    }
}
