#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <map>
#include <omp.h>

#include <boost/property_tree/ptree.hpp>
#include <minpack.h>

#include <jtflib/mesh/io.h>
#include <hjlib/sparse/sparse.h>
#include <hjlib/function/function.h>
#include <zjucad/optimizer/optimizer.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/ptree/ptree.h>

#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/lapack.h>

#include "../hex_frame_opt/function_term.h"
#include "../common/zyz.h"
#include "../common/IO.h"
#include "../common/util.h"
#include "../common/vtk.h"
#include "../common/visualize_tool.h"
#include "../common/transition.h"
#include "../common/transition_type.h"
#include "../common/face_orient.h"



#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
#include "../tetmesh/util.h"
#include "../hex_param/hex_param.h"

#include "../hex_frame_opt/frame_function.h"
#include "../hex_frame_opt/util.h"
#include "../hex_frame_opt/hex_frame_opt.h"
#include "./common.h"
#include "../hex_param/io.h"
#include "../hex_param/cut_tet.h"
#include "../hex_param/hex_param.h"
#include "../hex_param/common.h"
#include "../hex_param/find_singularities.h"
#include "../hex_param/global_alignment.h"
#include "../hex_param/topology_analysis.h"
#include "../hex_param/topology_operation.h"
#include "../hex_param/wyz_format.h"
#include "../hex_param/remove_surface_wedge.h"

#include <jtflib/mesh/mesh.h>
#include <jtflib/math/math.h>
#include "../hexmesh/io.h"
#include "../hex_process/hex_process.h"

#include "../tetmesh_refine/tetmesh_refine.h"
#include "../tetmesh/util.h"

#include "../tet_mesh_sxx/tet_mesh_sxx.h"
#include "../hex_process/hex_process.h"
#include "../tetmesh/util.h"


using namespace std;
using namespace zjucad::matrix;
using namespace hj::sparse;
using namespace hj::function;

using boost::property_tree::ptree;

void surf2inner_dijkstra(const matrixst &prev_node,
                         const jtf::mesh::face2tet_adjacent &fa, matrixd &zyz)
{
  matrixst node2face(prev_node.size());
  for(size_t ni = 0; ni < zyz.size(2); ++ni) { // ni face
      node2face[fa.faces_[ni][0]] = ni;
    }
  for(size_t ni = 0; ni < zyz.size(2); ++ni) { // ni face
      size_t i;
      for(i = fa.faces_[ni][0]; i != prev_node[i]; i = prev_node[i]);
      zyz(colon(), ni) = zyz(colon(), node2face[i]);
    }
}

int surf2inner_dijkstra(ptree & pt)
{
  matrixd zyz;
  {
    ifstream ifs(pt.get<string>("surf-zyz.value").c_str(), ifstream::binary);
    if(!ifs.fail()) {
        if(jtf::mesh::read_matrix(ifs, zyz, 3)) {
            cerr << "# not a zyz file." << endl;
            return __LINE__;
          }
      }
    else {
        cerr << "open init fail." << endl;
        return __LINE__;
      }
  }
  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet").c_str(), &tm.node_, &tm.mesh_))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));

  ifstream ifs(pt.get<string>("surf2inner.value").c_str(), ifstream::binary);
  if(ifs.fail()) {
      cerr << "open surf2inner file fail." << endl;
      return __LINE__;
    }
  matrixd dist;
  matrixst prev_node;
  jtf::mesh::read_matrix(ifs, dist);
  jtf::mesh::read_matrix(ifs, prev_node);
  surf2inner_dijkstra(prev_node, *fa, zyz);
  cerr << "# interp to inner." << endl;

  ofstream ofs(pt.get<string>("output.value").c_str(), ofstream::binary);
  jtf::mesh::write_matrix(ofs, zyz);

  return 0;
}

int cut_tet_inner(ptree &pt)
{  
  jtf::tet_mesh tm(pt.get<string>("tet.value").c_str());

  matrixd tet_zyz(3, tm.tetmesh_.mesh_.size(2));
  pt.put("zyz.desc","tet zyz frame file");
  if(read_zyz(pt.get<string>("zyz.value").c_str(), tet_zyz))
    return __LINE__;

  tetmesh_cutter tmc(tm);
  tmc.cut(tet_zyz);

  pt.put("cut_tet.desc","output tet");

  if(jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("cut_tet.value").c_str(), &tmc.cut_tm_.node_, &tmc.cut_tm_.mesh_))
    return __LINE__;
  return 0;
}

int find_singularities2_inner(ptree & pt)
{
  jtf::tet_mesh tm(pt.get<string>("tet.value").c_str());

  matrixd tet_zyz;
  if(read_zyz(pt.get<string>("zyz.value").c_str(), tet_zyz))
    return __LINE__;

  if(tet_zyz.size(2) != tm.tetmesh_.mesh_.size(2)) {
      cerr << "# [error] incompatible tet zyz file." << endl;
      cerr << "# [error] tet num " << tm.tetmesh_.mesh_.size(2) << endl;
      cerr << "# [error] zyz num " << tet_zyz.size(2) << endl;
      return __LINE__;
    }

  pt.put("vtk_dump.desc","output singularity vtk file");
  const string vtk_file = pt.get<string>("vtk_dump.value");

  vector<deque<pair<size_t,size_t> > > singularities_chain;
  vector<deque<size_t> > singularities_type;

  singularity_extractor se(tm);
  se.extract(tet_zyz, singularities_chain, singularities_type);

  dump_singularity_to_vtk((vtk_file+"_edge_idx.vtk").c_str(),tm.tetmesh_.node_,singularities_chain );
  const string edge_vtk_file = vtk_file + "_edge.vtk";
  dump_singularity_chain_to_vtk_2(edge_vtk_file.c_str(), tm.tetmesh_.node_,singularities_chain, singularities_type);

  cerr << "# find " << singularities_chain.size() << " singularities chains." << endl;

  return 0;
}

int find_singularities_use_face_type(ptree & pt)
{
  pt.put("tet.desc","tet file");
  jtf::tet_mesh tm(pt.get<string>("tet.value").c_str());

  pt.put("vtk_dump.desc","output singularity vtk file");
  const string vtk_file = pt.get<string>("vtk_dump.value");

  pt.put("face_type.desc","face type file");
  boost::unordered_map<pair<size_t,size_t>,size_t>  inner_type;
  if(load_inner_face_jump_type( pt.get<string>("face_type.value").c_str(),
                                inner_type))
    return __LINE__;

  vector<deque<pair<size_t,size_t> > > singularities_chain;
  vector<deque<size_t> > singularities_type;

  singularity_extractor se(tm);
  se.extract(inner_type, singularities_chain, singularities_type);

  dump_singularity_to_vtk((vtk_file+"_edge_idx.vtk").c_str(),tm.tetmesh_.node_,singularities_chain );
  const string edge_vtk_file = vtk_file + "_edge.vtk";
  dump_singularity_chain_to_vtk_2(edge_vtk_file.c_str(), tm.tetmesh_.node_,singularities_chain, singularities_type);

  cerr << "# find " << singularities_chain.size() << " singularities chains." << endl;
  return 0;
}

int draw_sh_smooth(ptree &pt)
{
  matrixd zyz;
  jtf::tet_mesh tm(pt.get<string>("tet.value").c_str());

  pt.put("zyz.desc","need inner frame zyz file");
  if(read_zyz(pt.get<string>("zyz.value").c_str(), zyz)){
      cerr << "# [error] can not load zyz." << endl;
      return __LINE__;
    }
  if(zyz.size(2) != tm.tetmesh_.mesh_.size(2)) {
      cerr << "wrong zyz file." << endl;
      return __LINE__;
    }

  matrixd sh(9, zyz.size(2));
  for(size_t ti = 0; ti < sh.size(2); ++ti){
      calc_rot_cubic_f_sh_(&sh(0, ti), &zyz(0, ti));
    }

  std::string path = pt.get<string>("zyz.value");
  dump_sh_difference_to_vtk((path+".sh.residual.vtk").c_str(),
                            tm.tetmesh_.mesh_, tm.tetmesh_.node_, *tm.fa_, sh);
  cerr << "# sucessed" << endl;

  return 0;
}

int draw_inner_residual(ptree &pt)
{
  matrixd zyz;
  jtf::tet_mesh tm(pt.get<string>("tet.value").c_str());

  pt.put("zyz.desc","need inner frame zyz file");
  if(read_zyz(pt.get<string>("zyz.value").c_str(), zyz)){
      cerr << "# [error] can not load zyz." << endl;
      return __LINE__;
    }
  if(zyz.size(2) != tm.tetmesh_.mesh_.size(2)) {
      cerr << "wrong zyz file." << endl;
      return __LINE__;
    }

  matrixd tet_zyz_bkp = zyz;
  hj_frame_alignemt(tm.tetmesh_.mesh_,*tm.fa_,tet_zyz_bkp,zyz);

  if(zjucad::has("output/align_zyz.value",pt)){
      if(write_zyz(pt.get<string>("align_zyz.value").c_str(), zyz)){
          cerr << "# [error] can not open align_zyz file." << endl;
          return __LINE__;
        }
    }

  { // for vis
    string path = pt.get<string>("zyz.value");
    matrix<matrixd > frame(zyz.size(2));
    zyz2frame(zyz, frame);
    dump_frame_difference_to_vtk((path+".aligned.zyz.residual.vtk").c_str(),
                                 tm.tetmesh_.mesh_, tm.tetmesh_.node_, *tm.fa_, frame);
    cerr << "# sucessed" << endl;
  }
  return 0;
}

int label_polycube(ptree &pt)
{
  jtf::tet_mesh tm(pt.get<string>("tet.value").c_str());
  matrixd zyz;

  pt.put("zyz.desc","need inner frame zyz file");
  if(read_zyz(pt.get<string>("zyz.value").c_str(), zyz)){
      cerr << "# [error] can not load zyz." << endl;
      return __LINE__;
    }
  if(zyz.size(2) != tm.tetmesh_.mesh_.size(2)) {
      cerr << "wrong zyz file." << endl;
      return __LINE__;
    }

  matrixst outside_face_belong_to_tet_idx;
  matrixst outside_face_label;

  const size_t outside_face_num = tm.outside_face_idx_.size();
  outside_face_label.resize(outside_face_num);
  outside_face_belong_to_tet_idx.resize(outside_face_num);

  for(size_t fi = 0; fi < outside_face_num; ++fi) {
      const pair<size_t,size_t> &tet_adj = tm.fa_->face2tet_[tm.outside_face_idx_[fi]];
      outside_face_belong_to_tet_idx[fi] = (tet_adj.first == -1)? tet_adj.second:tet_adj.first;
    }

  matrix<matrixd > rot_frame_at_face(outside_face_num);
  matrix<matrixd > normal_rot_outside_face(outside_face_num);

  // change each outside tet zyz into rotation matrix format
  for(size_t i = 0; i < outside_face_num; ++i){
      rot_frame_at_face[i].resize(3,3);
      zyz_angle_2_rotation_matrix1(&zyz(0, outside_face_belong_to_tet_idx[i]), &rot_frame_at_face[i][0]);
    }

  // assign the stander axis
  matrixd stander_axis(3,6);
  {
    matrixd x = zeros<double>(3,1);
    matrixd y = zeros<double>(3,1);
    matrixd z = zeros<double>(3,1);
    x[0] = 1; y[1] = 1; z[2] = 1;

    stander_axis(colon(),0) = x;
    stander_axis(colon(),1) = -x;
    stander_axis(colon(),2) = y;
    stander_axis(colon(),3) = -y;
    stander_axis(colon(),4) = z;
    stander_axis(colon(),5) = -z;
  }

  matrix<float> scalar_normal_for_vtk(outside_face_num);
  // change each outside face normal into rotation matrix format
  vector<pair<double,size_t> > resi_label(6);
  matrixd tet_bary_center = zeros<double>(3,1);
  matrixd face_bary_center = zeros<double>(3,1);
  for(size_t i = 0; i < outside_face_num; ++i){
      const vector<size_t> &vert_idx = tm.fa_->faces_[tm.outside_face_idx_[i]];
      const matrixd edge[2] = {
        tm.tetmesh_.node_(colon(), vert_idx[1]) - tm.tetmesh_.node_(colon(), vert_idx[0]),
        tm.tetmesh_.node_(colon(), vert_idx[2]) - tm.tetmesh_.node_(colon(), vert_idx[0]),
      };
      matrixd normal = cross(edge[0], edge[1]);
      const double len = norm(normal);
      if(len > 1e-8) {
          normal /= len;
          {// test whether the nomral is right or oppoiste
            cal_average_node(tm.tetmesh_.node_(colon(),tm.tetmesh_.mesh_(colon(),outside_face_belong_to_tet_idx[i])),
                             tet_bary_center);
            cal_average_node(tm.tetmesh_.node_(colon(), tm.outside_face_(colon(),i)),
                             face_bary_center);
            const matrixd face_to_tet_bary = face_bary_center - tet_bary_center;

            if(dot(face_to_tet_bary,normal) < -1e-7)
              normal = -normal;
            scalar_normal_for_vtk[i] = dot(face_to_tet_bary,normal);
          }
          normal_rot_outside_face[i].resize(3,3);

          for(size_t li = 0; li < 6; ++li){
              resi_label[li]= make_pair(dot(normal,rot_frame_at_face[i] * stander_axis(colon(),li)),li);
            }
          sort(resi_label.begin(),resi_label.end());

          outside_face_label[i] = resi_label[5].second; // find the largest, which is the most aligned axis
        }else{
          cerr << "# degenerated surface triangle at face node:" << tm.outside_face_idx_[i] << endl;
        }
    }

  cerr << "# begin to label the surface." << endl;
  { // for vis
    string path = pt.get<string>("zyz.value");
    ofstream ofs((path + ".surface_label_for_polycube.vtk").c_str());
    tri2vtk(ofs, &tm.tetmesh_.node_[0], tm.tetmesh_.node_.size(2), &tm.outside_face_[0], tm.outside_face_.size(2));
    matrix<float> rgba,color_label;
    rgba.resize(4,outside_face_label.size());
    color_label = ones<float>(4,6);

    // x = whilt

    // -x = yellow
    color_label(2,1) = 0;

    // y = orange
    color_label(2,2) = 0;
    color_label(1,2) = 0.5;

    //-y = red
    color_label(1,3) = 0;
    color_label(2,3) = 0;

    // z = green
    color_label(0,4) = 0;
    color_label(2,4) = 0;

    // -z = blue
    color_label(0,5) = 0;
    color_label(1,5) = 0;

    for(size_t t = 0; t < outside_face_label.size(); ++t)
      {
        rgba(colon(),t) = color_label(colon(),outside_face_label[t]);
      }
    cell_data_rgba_and_scalar(ofs, &rgba[0], &scalar_normal_for_vtk[0], scalar_normal_for_vtk.size(),
        "polycube_label","normal_distance");
  }
  cerr << "# sucess labeling the surface." << endl;

  return 0;
}

static double h_scalar(const matrixd &coord)
{
  assert(coord.size() == 3);
  const double &x = coord[0];
  const double &y = coord[1];
  const double &z = coord[2];

  return (x * x) * ( y * y ) + (y * y) *( z * z) + (z * z) * ( x * x);
}

int refine_frame_field_after_aligned(ptree &pt)
{
  matrixd zyz;
  jtf::tet_mesh tm(pt.get<string>("tet.value").c_str());

  pt.put("zyz.desc","need inner frame zyz file");
  if(read_zyz(pt.get<string>("zyz.value").c_str(), zyz)){
      cerr << "# [error] read zyz fail." << endl;
      return __LINE__;
    }
  if(zyz.size(2) != tm.tetmesh_.mesh_.size(2)){
      cerr << "# [error] error zyz format." << endl;
      return __LINE__;
    }

  matrixd new_zyz;

  pt.put("frame_refine_strategy.desc","choose splitting strategy or original one to refine frame.[splitting/original]");
  string frame_refine_strategy = pt.get<string>("frame_refine_strategy.value","original");
  if(frame_refine_strategy == "splitting"){
      refine_frame_field_after_aligned_split(pt,tm,zyz,new_zyz);
    }else{
      refine_frame_field_after_aligned_new(pt,tm,zyz,new_zyz);
    }

  if(zjucad::has("output_zyz.value",pt)){
      write_zyz(pt.get<string>("output_zyz.value").c_str(), new_zyz);
    }
  return 0;
}

int split_tet_at_zigzag(ptree &pt)
{
  jtf::mesh::meshes tm;
  matrixd zyz;

  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),&tm.node_,&tm.mesh_)) {
      cerr << "# read tet fail." << endl;
      return __LINE__;
    }
  pt.put("zyz.desc","need inner frame zyz file");

  ifstream ifs(pt.get<string>("zyz.value").c_str(),ifstream::binary);
  if(!ifs.fail())
    {
      jtf::mesh::read_matrix(ifs,zyz);
      if(!zyz.size()) {
          cerr << "#read zyz fail." << endl;
          return __LINE__;
        }else if(zyz.size(2) != tm.mesh_.size(2)){
          cerr << "# error zyz format." << endl;
          return __LINE__;
        }
    }else{
      cerr << "# open zyz file fail." << endl;
      return __LINE__;
    }

  matrix<matrix<double> > frame_inner;
  zyz2frame(zyz, frame_inner);

  while(1){
      unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));

      jtf::mesh::one_ring_tet_at_edge ortae;
      ortae.add_tets(tm.mesh_, *fa);

      if(ortae.sort_into_loop(tm.mesh_,tm.node_)) {
          cerr << "# sort error." << endl;
          return __LINE__;
        }

      vector<vector<size_t> > singularities_tet_loop;
      matrixst outside_face;
      get_outside_face(*fa,outside_face);

      boost::unordered_set<pair<size_t,size_t> > singularities_segment;
      set<size_t> singularity_points;
      // find out all the singulatiy segments
      for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator eti = ortae.e2t_.begin();
          eti != ortae.e2t_.end(); ++eti ) {
          const vector<size_t> &loop = eti->second;
          const pair<size_t,size_t> & one_edge = eti->first;

          if(loop.front() != loop.back()) // error format or surface edge or surface edge
            continue;
          if(loop.front() == -1) // open loop
            continue;

          matrixd transition = eye<double>(3);
          matrixd shuffle_rot(3,3); // what's the use of shuffle_rot

          for(size_t i = 0; i < loop.size()-1; ++i) {
              get_best_alignment(&frame_inner[loop[i]][0], &frame_inner[loop[i+1]][0],
                  &shuffle_rot[0]);

              transition = temp(transition * shuffle_rot);
            }

          // label the singularity edge
          if(norm(transition - eye<double>(3)) > 1e-8) {
              if(one_edge.first < one_edge.second)
                singularities_segment.insert(one_edge);
              else
                singularities_segment.insert(make_pair(one_edge.second, one_edge.first));
              singularity_points.insert(one_edge.first);
              singularity_points.insert(one_edge.second);
            }
        }


      set<pair<size_t,size_t> > split_edges;
      //      vector<int> singularity_edge_flag(3);
      //      for(size_t fi = 0; fi < fa->faces_.size(); ++fi){
      //          const vector<size_t> & one_face = fa->faces_[fi];

      //          for(size_t i = 0; i < one_face.size(); ++i){
      //              pair<size_t,size_t> one_edge(one_face[i], one_face[(i+1)%one_face.size()]);
      //              if(one_edge.first > one_edge.second)
      //                swap(one_edge.first, one_edge.second);
      //              if(singularities_segment.find(one_edge) != singularities_segment.end())
      //                singularity_edge_flag[i] = 1;
      //              else
      //                singularity_edge_flag[i] = 0;
      //            }
      //          if(std::accumulate(singularity_edge_flag.begin(),
      //                             singularity_edge_flag.end(),
      //                             static_cast<int>(0)) == 2){
      //              const size_t ei =
      //                  min_element(singularity_edge_flag.begin(),
      //                              singularity_edge_flag.end())
      //                  - singularity_edge_flag.begin();
      //              pair<size_t,size_t> one_split_edge(one_face[ei], one_face[(ei+1)%one_face.size()]);
      //              if(one_split_edge.first > one_split_edge.second)
      //                swap(one_split_edge.first, one_split_edge.second);
      //              split_edges.insert(one_split_edge);
      //            }
      //        }
      for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator eti = ortae.e2t_.begin();
          eti != ortae.e2t_.end(); ++eti  ){
          const pair<size_t,size_t> &one_edge = eti->first;

          const vector<size_t> &loop = eti->second;

          if(loop.front() != loop.back()) // error format or surface edge or surface edge
            continue;
          if(loop.front() == -1) // open loop
            continue;

          matrixd transition = eye<double>(3);
          matrixd shuffle_rot(3,3); // what's the use of shuffle_rot

          for(size_t i = 0; i < loop.size()-1; ++i) {
              get_best_alignment(&frame_inner[loop[i]][0], &frame_inner[loop[i+1]][0],
                  &shuffle_rot[0]);

              transition = temp(transition * shuffle_rot);
            }

          // label the singularity edge
          if(norm(transition - eye<double>(3)) < 1e-8) {
              if(singularity_points.find(one_edge.first) != singularity_points.end() &&
                 singularity_points.find(one_edge.second) != singularity_points.end())
                {
                  split_edges.insert(one_edge);
                }
            }

        }

      cerr << "# [info] find zigzag edges " << split_edges.size() << endl;
      if(split_edges.size() == 0) break;

      sxx::tet_mesh stm;
      stm.create_tetmesh(tm.node_, tm.mesh_);
      for(const auto & one_edge: split_edges){
          stm.split_edge(one_edge);
        }

      matrix<size_t> tet_mapping;
      stm.write_tetmesh_to_matrix(tm.node_, tm.mesh_);
      stm.get_tet2orginal_index(tm.mesh_, tet_mapping);

      matrix<matrix<double> > new_frame(tm.mesh_.size(2),1);
      for(size_t ti = 0; ti < tm.mesh_.size(2); ++ti){
          new_frame[ti] = frame_inner[tet_mapping[ti]];
        }
      swap(new_frame, frame_inner);

    }

  {// dump out the new tet and new zyz
    pt.put("output_tet.desc","split output tet mesh.");
    pt.put("output_zyz.desc","split output zyz.");

    string output_tet = pt.get<string>("output_tet.value");

    if(jtf::mesh::tet_mesh_write_to_zjumat(output_tet.c_str(),&tm.node_, &tm.mesh_))
      return __LINE__;

    matrix<double> new_zyz(3,frame_inner.size());
    frame2zyz(frame_inner, new_zyz);

    write_zyz(pt.get<string>("output_zyz.value").c_str(), new_zyz);

    cerr << "# finish dump out the split tet mesh and zyz ." << endl;
  }

  return 0;
}


int check_jump_type_and_folding(ptree &pt)
{
  matrixd zyz;
  jtf::tet_mesh tm(pt.get<string>("tet.value").c_str());

  pt.put("zyz.desc","need inner frame zyz file");
  if(read_zyz(pt.get<string>("zyz.value").c_str(), zyz)){
      cerr << "# [error] wrong zyz file." << endl;
      return __LINE__;
    }
  if(zyz.size(2) != tm.tetmesh_.mesh_.size(2)){
      cerr << "# error zyz format." << endl;
      return __LINE__;
    }

  matrix<matrixd > frame_inner;
  zyz2frame(zyz, frame_inner);

  matrixd rot = zeros<double>(3,3);
  matrixd diff_angle = zeros<double>(3,1);

  vector<double> error_angle;
  string rotation_jump_str = pt.get<string>("zyz.value");
  rotation_jump_str += ".jump_type_no_identity";
  ofstream ofs(rotation_jump_str.c_str());
  error_angle.reserve(tm.fa_->face2tet_.size());
  for(size_t t = 0; t < tm.fa_->face2tet_.size(); ++t){
      const pair<size_t,size_t> &tet_pair = tm.fa_->face2tet_[t];
      if(tet_pair.first == -1 || tet_pair.second == -1) continue;
      get_best_alignment(&frame_inner[tet_pair.first](0,0),
          &frame_inner[tet_pair.second](0,0),
          &rot[0]);

      if(norm(rot - eye<double>(3)) > 1e-8){
          ofs << tet_pair.first << " " << tet_pair.second << " " << type_transition1(rot) << endl;
          ofs << tet_pair.second << " " << tet_pair.first << " " << type_transition1(trans(rot)) << endl;
        }
      const matrixd aligned_left = frame_inner[tet_pair.first] * rot;
      for(size_t i = 0; i < 3; ++i){
          diff_angle[i] = fabs(
                jtf::math::safe_acos(
                  dot(aligned_left(colon(),i),
                      frame_inner[tet_pair.second](colon(),i))
                /(norm(aligned_left(colon(),i))*norm(frame_inner[tet_pair.second](colon(),i)))))
              * 180.0 / My_PI();
        }
      error_angle.push_back(*max_element(diff_angle.begin(),diff_angle.end()));
    }
  assert(error_angle.size());
  const double max_err_angle = *max_element(error_angle.begin(),error_angle.end());
  const double average_err_anle = std::accumulate(error_angle.begin(),error_angle.end(),0.0) / error_angle.size();

  cerr << "# max_err_angle: " << max_err_angle << endl;
  cerr << "# average_err_angle: " << average_err_anle << endl;

  return 0;
}

int pare_hexmesh_from_surface(ptree &pt)
{
  using namespace jtf::hexmesh;
  pt.put("hex.desc","file name to the input hex mesh");
  pt.put("hex_format.desc","hex format, default is 1");
  jtf::mesh::meshes hm;
  if(hex_mesh_read_from_wyz(pt.get<string>("hex.value").c_str(),
                            hm.mesh_, hm.node_,
                            pt.get<size_t>("hex_format.value"))){
      cerr << "# [error]  can not load hex mesh." << endl;
      return __LINE__;
    }

  matrixst hex0 = hm.mesh_;
  matrixst hex1 = hex0;

  pt.put("lay_n.desc","layer num need to be removed");
  const size_t layer_num = pt.get<size_t>("lay_n.value");
  string str = "_pare_hex_mesh.vtk";
  string frame_idx;
  stringstream ss;
  for(size_t t = 0; t < layer_num; ++t){
      frame_idx.clear();
      ss.str("");
      ss << t;
      //frame_idx = ss.str();
      string file_name = ss.str() + str;
      cerr << "# [info] vtk name: " << file_name << endl;
      ofstream ofs(file_name.c_str());

      hex2vtk(ofs,&hm.node_[0],hm.node_.size(2),&hex1[0],hex1.size(2));
      swap(hex0,hex1);
      pare_hexmesh_from_surface(hex0,hm.node_,hex1);

    }
  return 0;
}


int extend_tetmesh(ptree &pt)
{
  jtf::mesh::meshes tm;

  pt.put("tet.desc","tet mesh filename");
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),&tm.node_,&tm.mesh_)) {
      cerr << "# [error] can not read tet." << endl;
      return __LINE__;
    }

  pt.put("zyz.desc", "input zyz");
  matrix<double> zyz;
  if(zjucad::has("zyz.value", pt)){
      if(jtf::mesh::read_matrix(pt.get<string>("zyz.value").c_str(),zyz))
        return __LINE__;
    }

  pt.put("lay_n.desc","the num of extended tet mesh layer");
  const size_t layer_num = pt.get<size_t>("lay_n.value");
  jtf::mesh::meshes tm0 = tm,tm1 = tm;

  map<size_t,size_t> surface_point_map;
  for(size_t t = 0; t < layer_num; ++t){
      swap(tm1,tm0);
      jtf::tetmesh::extend_tetmesh(tm0.mesh_,tm0.node_,tm1.mesh_,tm1.node_,surface_point_map, zyz.size()==0?0:&zyz);
    }

  pt.put("output_tet.desc","output tet mesh file name");

  if(zjucad::has("zyz.value",pt)){
      const string output_zyz_str = pt.get<string>("output_zyz.value");
      if(jtf::mesh::write_matrix(output_zyz_str.c_str(),zyz))
        return __LINE__;
    }
  string output_tet_str = pt.get<string>("output_tet.value");

  {// deform new surface points to original ones
    jtf::tetmesh::deform_tet_accordint_to_surface_constraints(tm1.mesh_, tm1.node_, tm0.mesh_, tm0.node_, surface_point_map,pt);
  }

  jtf::mesh::tet_mesh_write_to_zjumat(output_tet_str.c_str(),&tm1.node_,&tm1.mesh_);

  return 0;
}

int pare_hex_outside_tet_surface(ptree &pt)
{
  using namespace jtf::hexmesh;
  pt.put("tet.desc","original tet filename");
  pt.put("hex.desc","input hex filename");

  jtf::mesh::meshes hm,new_hm;
  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),&tm.node_,&tm.mesh_))
    return __LINE__;

  if(hex_mesh_read_from_wyz(pt.get<string>("hex.value").c_str(),
                            hm.mesh_,hm.node_,
                            pt.get<size_t>("hex_format.value")))
    return __LINE__;

  trivial_hit_func thf;
  pare_hex_outside_tet_surface(tm.mesh_,tm.node_,hm.mesh_,hm.node_,thf,new_hm.mesh_,new_hm.node_);

  pt.put("new_hex.desc","new hex mesh after pared");
  hex_mesh_write_to_wyz(pt.get<string>("new_hex.value").c_str(),new_hm.mesh_,new_hm.node_);

#if 1 // visual
  if(zjucad::has("vtk.value",pt)){
      ofstream ofs(pt.get<string>("vtk.value").c_str());
      hex2vtk(ofs,&new_hm.node_[0],new_hm.node_.size(2),&new_hm.mesh_[0],new_hm.mesh_.size(2));
    }
#endif

  return 0;
}

int relax_singularities(ptree &pt)
{
  jtf::mesh::meshes tm;
  pt.put("tet.desc","tet mesh filename");
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),&tm.node_,&tm.mesh_)) {
      cerr << "# [error] can not read tet." << endl;
      return __LINE__;
    }

  pt.put("singularity_chain.desc","file name of singularity chain");
  std::vector<std::deque<std::pair<size_t,size_t> > > singularities_chain;
  std::vector<size_t> singularities_type;
  if(load_singularity_chain(pt.get<string>("singularity_chain.value").c_str(),
                            singularities_chain,
                            singularities_type))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  matrixst outside_face;
  get_outside_face(*fa,outside_face);
  pt.put("iter_relax.desc","relaxing iteration num");
  pt.put("new_tet.desc","file name of new_tet after relax");
  const string str_new_tet = pt.get<string>("new_tet.value");
  const size_t iter_relax = pt.get<size_t>("iter_relax.value");
  for(size_t i = 0; i < iter_relax; ++i)
    relax_singularities(tm.mesh_,tm.node_,singularities_chain,*fa,outside_face,pt);

  { // dump out new tet
    jtf::mesh::tet_mesh_write_to_zjumat(str_new_tet.c_str(),
                                        &tm.node_,&tm.mesh_);
    dump_singularity_to_vtk("singularity_after_relax.vtk", tm.node_,
                            singularities_chain);
  }
  return 0;
}

int calc_tet_rotation(ptree &pt)
{  
  jtf::tet_mesh tm(pt.get<string>("tet.value").c_str());

  pt.put("inner_face_jump_type.desc","file name of inner face jump type");
  boost::unordered_map<pair<size_t,size_t>,size_t> inner_face_jump_type;

  if(load_inner_face_jump_type(pt.get<string>("inner_face_jump_type.value").c_str(),
                               inner_face_jump_type))
    throw std::invalid_argument("can not load inner type.");

  pt.put("dump_tet_rot.desc","dump out tet rotation filename");
  const string dump_tet_rot_str = pt.get<string>("dump_tet_rot.value");

  vector<size_t> tet_rot_to_root;
  calc_tet_rotation(tm,inner_face_jump_type,tet_rot_to_root);

  ofstream ofs(dump_tet_rot_str.c_str());
  if(ofs.fail()){
      cerr << "# [error] can not open tet rot file." << endl;
      return __LINE__;
    }
  ofs << "tet " << tet_rot_to_root.size() << " root " << 0 << endl;
  for(size_t t = 0; t < tet_rot_to_root.size(); ++t){
      ofs << tet_rot_to_root[t] << " ";
    }

  return 0;
}

int arap_to_concentrate_points(ptree &pt)
{
  jtf::mesh::meshes tm;
  pt.put("tet.desc","tet mesh filename");
  pt.put("new_tet.desc","tet mesh after arap deormation");
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),
                                          &tm.node_,&tm.mesh_)) {
      cerr << "# [error] can not read tet." << endl;
      return __LINE__;
    }

  arap_to_concentrate_points(tm.mesh_,tm.node_,pt);

  if(jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("new_tet.value").c_str(),
                                         &tm.node_,&tm.mesh_))
    return __LINE__;
  return 0;
}

int remove_degenerated_edges(ptree &pt)
{
  jtf::mesh::meshes tm;
  matrixd zyz;

  pt.put("tet.desc","tet mesh filename");
  pt.put("zyz.desc","zyz field filename");
  pt.put("new_tet.desc","new tet mesh after modifying");
  pt.put("new_zyz.desc","new zyz after modifying");

  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),
                                          &tm.node_,&tm.mesh_)) {
      cerr << "# [error] can not read tet." << endl;
      return __LINE__;
    }

  if(read_zyz(pt.get<string>("zyz.value").c_str(),
              zyz))
    return __LINE__;

  vector<size_t> new_tet_vec(tm.mesh_.size());
  vector<double> new_node_vec(tm.node_.size());
  vector<matrixd > new_zyz_vec(tm.mesh_.size(2));

  copy(tm.mesh_.begin(),tm.mesh_.end(),new_tet_vec.begin());
  copy(tm.node_.begin(),tm.node_.end(),new_node_vec.begin());
  for(size_t t = 0; t < tm.mesh_.size(2); ++t){
      new_zyz_vec[t] = zyz(colon(),t);
    }

  jtf::mesh::one_ring_tet_at_edge ortae;
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_,"topology"));
  ortae.add_tets(tm.mesh_,*fa);
  ortae.sort_into_loop(tm.mesh_,tm.node_);

  vector<deque<pair<size_t,size_t> > > degenerated_chains;
  pt.put("strategy.desc",
         "strategy to choose handle original black edge or integer conflicts"
         "[black]");
  const string strategy = pt.get<string>("strategy.value");

  if(strategy == "black"){
      pt.put("singularity_chain.desc","singularity chain file");
      vector<deque<pair<size_t,size_t> > > singularity_chain;
      vector<deque<size_t> > singularity_type;
      if(load_singularity_chain_new(
           pt.get<string>("singularity_chain.value").c_str(),
           singularity_chain,singularity_type))
        return __LINE__;

      for(size_t t = 0; t < singularity_chain.size(); ++t){
          const deque<pair<size_t,size_t> > & one_chain = singularity_chain[t];
          vector<size_t> degenerated_points_vec;
          if(is_black_line_new(singularity_type[t].front()))
            {
              degenerated_chains.push_back(one_chain);
              continue;
            }
        }
    }else{
      cerr << "# [error] can not handle such strategy." << endl;
      return __LINE__;
    }

  for(size_t t = 0; t < degenerated_chains.size(); ++t){
      const deque<pair<size_t,size_t> > & one_chain = degenerated_chains[t];
      remove_one_chain_and_update_info(one_chain,new_tet_vec,new_node_vec,
                                       new_zyz_vec,ortae);
    }

  itr_matrix<size_t*> new_tet_mat(4,new_tet_vec.size()/4,&new_tet_vec[0]);
  itr_matrix<double*> new_node_mat(3,new_node_vec.size()/3,&new_node_vec[0]);

  matrixst new_tet = new_tet_mat;
  matrixd new_node = new_node_mat;

  jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("new_tet.value").c_str(),
                                      &new_node,&new_tet);

  assert(new_tet_mat.size(2) == new_zyz_vec.size());

  matrixd new_zyz_mat(3,new_zyz_vec.size());

  for(size_t t = 0; t < new_zyz_vec.size(); ++t)
    new_zyz_mat(colon(),t) = new_zyz_vec[t];

  write_zyz(pt.get<string>("new_zyz.value").c_str(),
            new_zyz_mat);
  return 0;
}

int minimal_cut_tet(ptree &pt)
{
  jtf::tet_mesh tm(pt.get<string>("tet.value").c_str());

//  vector<deque<pair<size_t,size_t> > > singularity_edges;
//  vector<deque<size_t> > singularity_type;

//  pt.put("singularity_chain.desc", "input singularity chian");
//  if(load_singularity_chain_new(
//       pt.get<string>("singularity_chain.value").c_str(),
//       singularity_edges, singularity_type))
//    return __LINE__;


  pt.put("zyz.desc", "input zyz file");
  matrixd zyz(3, tm.tetmesh_.mesh_.size(2));
  if(read_zyz(pt.get<string>("zyz.value").c_str(), zyz)){
      cerr << "# [error] read zyz error" << endl;
      return  __LINE__;
    }

  if(zyz.size(2) != tm.tetmesh_.mesh_.size(2)) {
      cerr << "# [error] wrong zyz." << endl;
      return __LINE__;
    }

  boost::unordered_map<pair<size_t,size_t>, size_t> inner_face_jump_type;
  pt.put("inner_face_jump_type.desc", " input inner_face_jump_type");

  if(zjucad::has("inner_face_jump_type.value",pt)){
      if(load_inner_face_jump_type(
           pt.get<string>("inner_face_jump_type.value").c_str(),
           inner_face_jump_type))
        return __LINE__;
    }else{
      extract_inner_face_type_from_zyz(*tm.fa_, zyz, inner_face_jump_type);
    }

  matrixst cut_tet;
  matrixst face_pair;
  matrixst face_type;

  // test param node
  pt.put("param_strategy.desc", "param node strategy [wyz/jtf]");
  if(zjucad::has("param_strategy.value", pt)){
      const string strategy_ = pt.get<string>("param_strategy.value");
      matrixd param(3, tm.tetmesh_.mesh_.size(2) * 4);

      if(strategy_ == "wyz"){
          pt.put("param.desc", "input param nodes");
          if(zjucad::has("param.value",pt)){
              if(load_wyz_param_file(pt.get<string>("param.value").c_str(),
                                     tm.tetmesh_.mesh_.size(2),
                                     param))
                return __LINE__;
            }
        }else if(strategy_ == "jtf"){
          pt.put("param.desc", "input deformated tet mesh");
          matrixst cut_tet;
          if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("param.value").c_str(),
                                                  &param, &cut_tet))
            return __LINE__;
          if(cut_tet.size() != tm.tetmesh_.mesh_.size())
            return __LINE__;
        }
      minimal_cut_tet(tm, zyz, cut_tet, face_pair, face_type, & param);
    }else
    minimal_cut_tet(tm, zyz, cut_tet, face_pair, face_type);

  pt.put("cut_tet_mat.desc", "only cut_tet mesh mtrix");
  pt.put("face_pair_with_type.desc", "output face pair with type");

  matrix<int> cut_tet_int = cut_tet;
  jtf::mesh::write_matrix(pt.get<string>("cut_tet_mat.value").c_str(),
                          cut_tet_int);

  if(cut_tet.size() != tm.tetmesh_.mesh_.size())
    cerr << "# [error] wrong cut_tet size" << endl;

  assert(face_pair.size() == face_type.size());
  matrixst face_pair_with_type(3, face_type.size());
  for(size_t t = 0; t < face_pair.size(); ++t){
      face_pair_with_type(0,t) = t;
      face_pair_with_type(1,t) = face_pair[t];
      face_pair_with_type(2,t) = face_type[t];
    }

  jtf::mesh::write_matrix(pt.get<string>("face_pair_with_type.value").c_str(),
                          face_pair_with_type);

  matrixst cut_tet2tet(max(cut_tet) + 1);
  cut_tet2tet(cut_tet) = tm.tetmesh_.mesh_(colon());
  matrixd cut_node(3, cut_tet2tet.size());
  for(size_t t = 0; t < cut_tet2tet.size(); ++t){
      cut_node(colon(), t) = tm.tetmesh_.node_(colon(), cut_tet2tet[t]);
    }

  pt.put("cut_tet.desc", "output cut tet");
  if(jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("cut_tet.value").c_str(),
                                         &cut_node, &cut_tet));
  return 0;
}

int tet_mesh_inflation_for_ball(ptree &pt)
{
  jtf::mesh::meshes tm;
  pt.put("tet.desc", "input tet mesh");
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),
                                          &tm.node_, &tm.mesh_))
    return __LINE__;

  tet_mesh_inflation_for_ball(tm.mesh_, tm.node_, pt);

  pt.put("output_tet.desc", "output tet mesh");
  jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("output_tet.value").c_str(),
                                      &tm.node_, &tm.mesh_);
  return 0;
}

int check_bad_rounding(ptree & pt)
{
  jtf::mesh::meshes tm;
  pt.put("tet.desc", "input tetmesh");
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),
                                          &tm.node_, &tm.mesh_))
    return __LINE__;

  vector<deque<pair<size_t,size_t> > > singularity;
  vector<deque<size_t> > singularity_type;
  pt.put("singulairty_chain.desc", "singularity_chains filename");

  if(load_singularity_chain_new(
       pt.get<string>("singularity_chain.value").c_str(),
       singularity,singularity_type))
    return __LINE__;

  set<size_t> singularity_points;
  for(size_t t = 0; t < singularity.size(); ++t){
      const deque<pair<size_t,size_t> > & one_chain = singularity[t];
      for(size_t i = 0; i < one_chain.size(); ++i){
          singularity_points.insert(one_chain[i].first);
          singularity_points.insert(one_chain[i].second);
        }
    }

#if 0
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_,"topology"));
  pt.put("wyz_face_mapping.desc", "input face index from wyz");
  matrixst wyz_face_idx_to_jtf;
  map<size_t,size_t> jtf_face_idx_to_wyz;
  if(convert_wyz_face_index_mapping(
       pt.get<string>("wyz_face_mapping.value").c_str(),
       *fa,wyz_face_idx_to_jtf,jtf_face_idx_to_wyz))
    return __LINE__;
#endif
  pt.put("wyz_param.desc","wyz's parameterization file");

  matrixd param;

  if(load_wyz_param_file(pt.get<string>("wyz_param.value").c_str(),
                         tm.mesh_.size(2),param))
    return __LINE__;

#if 1 // check
  {
    // TODO: there may be a bug in parameterization construction
    // some data becomes very huge
    for(size_t t = 0; t < param.size(); ++t){
        if(fabs(param[t]) > 2e2)
          param[t] = param[t]>0?100:-100;
      }
    cerr << "# [info] max element " << *max_element(param.begin(),param.end()) << endl;
    cerr << "# [info] min element " << *min_element(param.begin(),param.end()) << endl;

    matrixst param_tet(4,tm.mesh_.size(2));
    param_tet(colon()) = colon(0,4 * tm.mesh_.size(2)-1);

    set<size_t> singularity_point_in_param_tet;
    for(set<size_t>::const_iterator scit = singularity_points.begin();
        scit != singularity_points.end(); ++scit){
        for(size_t i = 0; i < tm.mesh_.size(2); ++i){
            for(size_t t = 0; t < tm.mesh_.size(1); ++t){
                if(tm.mesh_(t,i) == *scit)
                  singularity_point_in_param_tet.insert(param_tet(t,i));
              }
          }
      }

    vector<size_t> singularity_points_vec(singularity_point_in_param_tet.size());
    copy(singularity_point_in_param_tet.begin(),singularity_point_in_param_tet.end(),
         singularity_points_vec.begin());
    ofstream ofs_s("singularity_point_in_param.vtk");
    point2vtk(ofs_s,&param[0],param.size(2),&singularity_points_vec[0],singularity_points_vec.size());

    ofstream ofs("param.vtk");
    tet2vtk(ofs,&param[0],param.size(2),&param_tet[0],param_tet.size(2));
  }
#endif
  return 0;
}

int check_polycube_validation_with_equation_graph(ptree &pt)
{
  jtf::mesh::meshes tm, cut_tm;
  pt.put("tet.desc", "input tetmesh");
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),
                                          &tm.node_, &tm.mesh_))
    return __LINE__;

  pt.put("cut_tet.desc", "input cut_tet tetmesh");
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("cut_tet.value").c_str(),
                                          &cut_tm.node_, &cut_tm.mesh_))
    return __LINE__;

  if(tm.mesh_.size() != cut_tm.mesh_.size()){
      cerr << "# [error] input cut_tet has different size with original tet." << endl;
      return __LINE__;
    }

  boost::unordered_map<pair<size_t,size_t>,size_t> inner_face_jump_type_cut;
  boost::unordered_map<size_t,size_t> surface_type_cut;

  if(zjucad::has("inner_face_jump_type.value",pt)){
      if(load_inner_face_jump_type(
           pt.get<string>("inner_face_jump_type.value").c_str(),
           inner_face_jump_type_cut))
        return __LINE__;
    }else{
      cerr << "# [info] pure polycube structure." << endl;
    }

  pt.put("no_surface.desc", "[y/n] whether need surface type constraints.");
  bool no_surface = false;

  const string yn = pt.get<string>("no_surface.value");
  if(yn == "Yes" || yn == "y" || yn == "Y" || yn == "YES" || yn == "yes")
    no_surface = true;

  pt.put("surface_type.desc", "restricted surface type on orig mesh");

  if(!no_surface){
      if(load_surface_type(
           pt.get<string>("surface_type.value").c_str(), surface_type_cut))
        return __LINE__;
    }

  aggressive_assign_transition_with_type(
        tm.mesh_, cut_tm.mesh_, tm.node_, inner_face_jump_type_cut, surface_type_cut,no_surface);

  return 0;
}

int aggressive_extract_jump_type(ptree &pt)
{
  jtf::mesh::meshes tm;
  pt.put("tet.desc", "input tetmesh");
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),
                                          &tm.node_, &tm.mesh_))
    return __LINE__;

  boost::unordered_map<pair<size_t,size_t>,size_t> face_type;
  //load_inner_face_jump_type(pt.get<string>("inner_face_jump_type.value").c_str(),
  //                           face_type);
  // aggressive_extract_jump_type(tm.mesh_, tm.node_, zyz, face_type);

  bool no_surface = false;
  pt.put("no_surface.desc", "[y/n] whether need surface constriants.");
  if(zjucad::has("no_surface.value", pt)){
      const string no_surface_str = pt.get<string>("no_surface.value");
      if(no_surface_str == "yes" || no_surface_str == "Yes" ||
         no_surface_str == "YES" || no_surface_str == "y" ||
         no_surface_str == "Y")
        no_surface = true;
    }

  cerr << "# [info] no surface constraints: " << (no_surface?"yes":"no") << endl;

  if(zjucad::has("inner_face_jump_type.value",pt) ||
     zjucad::has("surface_type.value",pt)){
      load_inner_face_jump_type(
            pt.get<string>("inner_face_jump_type.value").c_str(), face_type);

      pt.put("zyz.desc", "input zyz file");
      matrixd zyz(3, tm.mesh_.size(2));
      if(zjucad::has("zyz.value", pt)){
          ifstream ifs(pt.get<string>("zyz.value").c_str(), ifstream::binary);
          if(jtf::mesh::read_matrix(ifs, zyz)){
              cerr << "# [error] read zyz error" << endl;
              return  __LINE__;
            }

          if(zyz.size(2) != tm.mesh_.size(2)) {
              cerr << "# [error] wrong zyz." << endl;
              return __LINE__;
            }
        }

      aggressive_assign_transition_with_origin_type(
            tm.mesh_, tm.node_, face_type,pt, no_surface, &zyz);

    }else{
      pt.put("zyz.desc", "input zyz file");
      matrixd zyz(3, tm.mesh_.size(2));
      ifstream ifs(pt.get<string>("zyz.value").c_str(), ifstream::binary);
      if(jtf::mesh::read_matrix(ifs, zyz)){
          cerr << "# [error] read zyz error" << endl;
          return  __LINE__;
        }

      if(zyz.size(2) != tm.mesh_.size(2)) {
          cerr << "# [error] wrong zyz." << endl;
          return __LINE__;
        }
      aggressive_assign_transition(tm.mesh_,tm.node_,zyz,no_surface,face_type);
    }
  return 0;
}

int deform_hex_to_original_tet(ptree &pt)
{
  using namespace jtf::hexmesh;
  jtf::mesh::meshes original_tet;
  pt.put("tet.desc", "input tetmesh");
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),
                                          &original_tet.node_, &original_tet.mesh_))
    return __LINE__;

  jtf::mesh::meshes input_hex;
  pt.put("hex_format.desc", "hex format of wyz [1/2/3], 1: original format;"
         "2: hex point with parameter; 3: invalid point idx is -1 ");

  const size_t hex_format = pt.get<size_t>("hex_format.value");

  pt.put("hex.desc", "input hexmesh");
  if (hex_mesh_read_from_wyz(pt.get<string>("hex.value").c_str(),
                             input_hex.mesh_, input_hex.node_, hex_format, true))
    return __LINE__;

  cerr << "# [info] hex element num " << input_hex.mesh_.size(2) << endl;
  deform_hex_to_original_tet(input_hex.mesh_, input_hex.node_, original_tet.mesh_,
                             original_tet.node_, pt);

  ofstream ofs("hex_deformed.vtk");
  hex2vtk(ofs, &input_hex.node_[0], input_hex.node_.size(2),
      &input_hex.mesh_[0], input_hex.mesh_.size(2));
  cout << " hex_for_wyz2 will generate" <<endl;
  hex_mesh_write_to_wyz("hex_for_wyz2.hex", input_hex.mesh_, input_hex.node_);
  return 0;
}

int split_tet_around_singularities(ptree & pt)
{
  matrixst split_tet;
  matrixd split_node, zyz;

  pt.put("tet.desc", "input tet mesh");
  pt.put("zyz.desc", "input zyz frame field");

  jtf::tet_mesh tm(pt.get<string>("tet.value").c_str());

  if(jtf::mesh::read_matrix(pt.get<string>("zyz.value").c_str(), zyz))
    return __LINE__;

  matrix<matrixd> frame_inner(zyz.size(2));
  zyz2frame(zyz, frame_inner);

  vector<deque<pair<size_t,size_t> > > chain_list;
  vector<deque<size_t> > singularities_type;

  singularity_extractor se(tm);
  se.extract(frame_inner, chain_list, singularities_type);

  boost::unordered_set<size_t> singularity_points;
  for(size_t ci = 0; ci < chain_list.size(); ++ci){
      const deque<pair<size_t,size_t> > & one_chain = chain_list[ci];
      for(size_t ei = 0; ei < one_chain.size(); ++ei){
          singularity_points.insert(one_chain[ei].first);
          singularity_points.insert(one_chain[ei].second);
        }
    }

  boost::unordered_set<pair<size_t,size_t> > edges_need_to_split;

  jtf::tetmesh::one_ring_edge_at_point oreap;
  oreap.add_tets(tm.tetmesh_.mesh_);

  for(boost::unordered_set<size_t>::const_iterator bucit
      = singularity_points.begin(); bucit != singularity_points.end(); ++bucit){
      const boost::unordered_set<size_t> & linked_points = oreap.p2p_[*bucit];
      for(boost::unordered_set<size_t>::const_iterator bucitj =
          linked_points.begin(); bucitj != linked_points.end(); ++bucitj){
          pair<size_t,size_t> edge(*bucit, *bucitj);
          if(edge.first > edge.second)
            swap(edge.first, edge.second);
          edges_need_to_split.insert(edge);
        }
    }

  // split tet
  {
    split_tet = tm.tetmesh_.mesh_;
    split_node = tm.tetmesh_.node_;
    sxx::tet_mesh stm;
    stm.create_tetmesh(tm.tetmesh_.node_, tm.tetmesh_.mesh_);
    //    std::tuple<size_t,size_t,size_t> trash;
    size_t mid_index;
    matrix<size_t> tet_idxmap;
    for(boost::unordered_set<pair<size_t,size_t> >::const_iterator eit =
        edges_need_to_split.begin(); eit != edges_need_to_split.end(); ++eit){
        mid_index = stm.split_edge(*eit);
      }
    stm.write_tetmesh_to_matrix(split_node, split_tet);
    stm.get_tet2orginal_index(split_tet, tet_idxmap);
    matrix<double> old_zyz = zyz;
    zyz.resize(3, split_tet.size(2));
    for(size_t ti = 0; ti < zyz.size(2); ++ti){
        zyz(colon(), ti) = old_zyz(colon(), tet_idxmap[ti]);
      }
  }


  {
    // visual tet
    ofstream ofs("tet_after_split.vtk");
    tet2vtk(ofs, &split_node[0], split_node.size(2), &split_tet[0], split_tet.size(2));
    vector<double> rev_volume_tet(split_tet.size(2));
    for(size_t ti = 0; ti < split_tet.size(2); ++ti){
        rev_volume_tet[ti] =
            jtf::mesh::cal_tet_vol(split_node(colon(), split_tet(colon(),ti)));
        rev_volume_tet[ti] = 1.0/rev_volume_tet[ti];
      }
    cell_data(ofs, &rev_volume_tet[0], rev_volume_tet.size(), "rev_volume");
  }

  jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("output/tet.value").c_str(), &split_node, &split_tet);
  jtf::mesh::write_matrix(pt.get<string>("output/zyz.value").c_str(), zyz);
  return 0;
}

int remove_near_miss_by_minimal_cut(ptree & pt)
{
  matrixst cut_tet; // the input tet mush be cutted
  matrixd cut_node;

  matrixst tet;
  matrixd node;

  pt.put("cut_tet.desc", "input cut tet");
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("cut_tet.value").c_str(),
                                          &cut_node, &cut_tet)){
      cerr << "# [error] can not read cut_tet." << endl;
      return __LINE__;
    }

  pt.put("tet.desc", "input tet");
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),
                                          &node, &tet)){
      cerr << "# [error] can not read cut_tet." << endl;
      return __LINE__;
    }

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  unique_ptr<jtf::mesh::face2tet_adjacent> fa_cut(jtf::mesh::face2tet_adjacent::create(cut_tet));

  if(!fa.get() || !fa_cut.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }

  boost::unordered_map<pair<size_t,size_t>,size_t> inner_face_jump_type;
  boost::unordered_map<size_t,size_t> surface_type;

  pt.put("inner_face_jump_type.desc", "input inner face jump type");
  pt.put("surface_type.desc", "input surface type which is extracted from graph,0/1/2");
  if(load_inner_face_jump_type(
       pt.get<string>("inner_face_jump_type.value").c_str(),
       inner_face_jump_type))
    return __LINE__;

  int rtn = load_surface_type(pt.get<string>("surface_type.value").c_str(),
                              surface_type);
  if(rtn != 0 && rtn != 1)
    return __LINE__;

  pt.put("g_unknown_face.desc", "input g_unknown_face file.");
  vector<pair<size_t,size_t> > g_unknown_faces;
  if(load_g_unknown_face(pt.get<string>("g_unknown_face.value").c_str(),
                         g_unknown_faces)) {
      return __LINE__;
    }

  vector<size_t> fnode;
  pt.put("fnode_group.desc", "fnode_group info");
  if(load_fnode_group(pt.get<string>("fnode_group.value").c_str(),
                      cut_tet, fnode))
    return __LINE__;


  const string strategy = pt.get<string>("strategy.value","new");
  if(strategy == "new")
    remove_near_miss_by_minimal_cut_new(
          cut_tet, tet, cut_node, node, inner_face_jump_type, surface_type,
          g_unknown_faces,*fa, *fa_cut, fnode);
  else
    remove_near_miss_by_minimal_cut_old(
          cut_tet, tet, cut_node, node, inner_face_jump_type, surface_type,
          g_unknown_faces,*fa, *fa_cut, fnode);

  dump_inner_face_jump_type("inner_face_jump_type_after_removing_near_miss",
                            inner_face_jump_type);

  return 0;
}

int remove_surface_wedge(ptree & pt)
{
  jtf::mesh::meshes tm, cut_tet;

  pt.put("tet.desc", "input tet mesh.");
  pt.put("cut_tet.desc", "input cut_tet mesh");
  pt.put("surface_type_orient.desc", "surface type with orientation [0-5]: +u,-u,+v,-v,+w,-w");
  pt.put("surface_type.desc", "input surface restricted type on cut_tetmesh");
  pt.put("inner_face_jump_type.desc", "inner face jump type");
  pt.put("loop_points.desc", "loop points file.");
  pt.put("g_unknown_face.desc", "input g_unknown_faces");
  pt.put("degenerated_edges.desc", "input degenerated edges");

  if(jtf::mesh::tet_mesh_read_from_zjumat(
       pt.get<string>("tet.value").c_str(), &tm.node_, &tm.mesh_)){
      cerr << "# [error] can not read tet." << endl;
      return __LINE__;
    }

  if(jtf::mesh::tet_mesh_read_from_zjumat(
       pt.get<string>("cut_tet.value").c_str(), &cut_tet.node_, &cut_tet.mesh_)){
      cerr << "# [error] can not read cut_tet." << endl;
      return __LINE__;
    }

  boost::unordered_map<std::pair<size_t,size_t>,size_t> inner_face_jump_type;
  if(load_inner_face_jump_type(
       pt.get<string>("inner_face_jump_type.value").c_str(), inner_face_jump_type))
    {
      return __LINE__;
    }

  boost::unordered_map<size_t,size_t> surface_type_cut;
  bool use_uvw_surface_type = true; // only uvw, no orientation
  if(zjucad::has("surface_type_orient.value", pt)){
      int rtn =load_surface_type(
            pt.get<string>("surface_type_orient.value").c_str(), surface_type_cut);
      if(rtn != 1)
        return __LINE__;
      use_uvw_surface_type = false;
    }else{
      if(!zjucad::has("surface_type.value",pt)){
          cerr << "# [error] You Must choose surface_type_orient [with +/- direction] "
               << " or surface_type [no +/- direction]." << endl;
          return __LINE__;
        }
      if(load_surface_type(
           pt.get<string>("surface_type.value").c_str(), surface_type_cut)){
          return __LINE__;
        }
    }

  vector<pair<size_t,size_t> > g_unknown_face_pair;
  if(zjucad::has("g_unknown_face.value", pt)){
      if(load_g_unknown_face(pt.get<string>("g_unknown_face.value").c_str(),
                             g_unknown_face_pair)){
          return __LINE__;
        }
      cerr << "# [info] use g_unknown_face." << endl;
    }

  boost::unordered_map<std::pair<size_t,size_t>,size_t> restricted_edges;
  if(zjucad::has("restricted_edge.value",pt)){
      if(load_restricted_edge(pt.get<string>("restricted_edge.value").c_str(),
                              restricted_edges))
        return __LINE__;
      cerr << "# [info] use restricted_edges." << endl;
    }

  matrixd zyz;
  if(zjucad::has("zyz.value", pt)){
      if(jtf::mesh::read_matrix(pt.get<string>("zyz.value").c_str(), zyz))
        return __LINE__;
      if(zyz.size(2) != tm.mesh_.size(2)){
          cerr << "# [error] imcompatible zyz size." << endl;
          return __LINE__;
        }
    }

  jtf::mesh::meshes new_tm = tm;
  if(zjucad::has("zyz.value", pt))
    remove_surface_wedge(
          tm.mesh_, tm.node_, cut_tet.mesh_, cut_tet.node_, g_unknown_face_pair,
          restricted_edges,use_uvw_surface_type, inner_face_jump_type,
          new_tm.mesh_, new_tm.node_, surface_type_cut, &zyz);
  else
    remove_surface_wedge(
          tm.mesh_, tm.node_, cut_tet.mesh_, cut_tet.node_, g_unknown_face_pair,
          restricted_edges, use_uvw_surface_type, inner_face_jump_type,
          new_tm.mesh_, new_tm.node_, surface_type_cut, 0);


  { // dump out data.
    pt.put("new_tet.desc", "output new_tet");
    if(jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("new_tet.value").c_str(),
                                           &new_tm.node_, &new_tm.mesh_)){
        return __LINE__;
      }

    if(zjucad::has("zyz.value",pt)){
        if(jtf::mesh::write_matrix(pt.get<string>("zyz.value").c_str(),zyz))
          return __LINE__;
        cerr << "# [info] write zyz." << endl;
      }
  }

  return 0;
}


int glue_cut_tet_according_to_cut_face(ptree & pt)
{
  pt.put("tet.desc", "input tet mesh");
  pt.put("cut_tet.desc", "input cut_tet mesh");
  pt.put("new_cut_tet.desc", "output new cut tet");

  jtf::tet_mesh tm(pt.get<string>("tet.value").c_str()), cut_tm(pt.get<string>("cut_tet.value").c_str());

  if(tm.tetmesh_.mesh_.size(2) != cut_tm.tetmesh_.mesh_.size(2)){
      cerr << "# [error] strange, tet size is not compatiable with cut_tet." << endl;
      return __LINE__;
    }
  vector<pair<size_t,size_t> > g_unknown_face_pair;

  pt.put("g_unknown_face.desc", "input g_unknown_face");
  if(load_g_unknown_face(pt.get<string>("g_unknown_face.value").c_str(),
                         g_unknown_face_pair)){
      return __LINE__;
    }

  vector<pair<size_t,size_t> > cut_tet_pair(g_unknown_face_pair.size());
  { // convert to tet pair
    for(size_t i = 0; i < g_unknown_face_pair.size(); ++i){
        const pair<size_t,size_t> & face_pair = g_unknown_face_pair[i];
        const pair<size_t,size_t> & tet_pair_0 = cut_tm.fa_->face2tet_.at(face_pair.first);
        const pair<size_t,size_t> & tet_pair_1 = cut_tm.fa_->face2tet_.at(face_pair.second);
        if(!cut_tm.fa_->is_outside_face(tet_pair_0) ||
           !cut_tm.fa_->is_outside_face(tet_pair_1)){
            cerr << "# [error] strange, face " << face_pair.first << " "
                 << face_pair.second << " should be outside face." << endl;
            return __LINE__;
          }

        cut_tet_pair[i] = make_pair(
              (tet_pair_0.first == -1?tet_pair_0.second:tet_pair_0.first),
              (tet_pair_1.first == -1?tet_pair_1.second:tet_pair_1.first));
      }
  }

  matrixst new_cut_tet;
  glue_tet_with_given_cut_faces(tm.tetmesh_.mesh_, tm.tetmesh_.node_, cut_tet_pair,new_cut_tet);

  matrixst new_cut2tet(max(new_cut_tet)+1);
  new_cut2tet(new_cut_tet) = tm.tetmesh_.mesh_(colon());

  matrixd new_cut_node = zeros<double>(3, new_cut2tet.size());
  for(size_t pi = 0; pi < new_cut2tet.size(); ++pi){
      new_cut_node(colon(), pi) = tm.tetmesh_.node_(colon(), new_cut2tet[pi]);
    }

  if(jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("new_cut_tet.value").c_str(),
                                         &new_cut_node, &new_cut_tet)){
      return __LINE__;
    }

  return 0;
}

int extract_relax_surface(ptree & pt)
{
  jtf::mesh::meshes tm;

  pt.put("tet.desc", "input tetmesh");
  pt.put("surface_type.desc", "input surface type");
  pt.put("loop_points.desc", "input loop points");
  pt.put("restricted_surface_type.desc", "output restricted surface type");

  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(),
                                          &tm.node_, &tm.mesh_))
    return __LINE__;

  boost::unordered_map<size_t,size_t> surface_type;
  int rtn = load_surface_type(pt.get<string>("surface_type.value").c_str(),
                              surface_type);
  if(rtn != 0 && rtn != 1)
    return __LINE__;

  boost::unordered_set<size_t> loop_points;
  if(load_loop_points(pt.get<string>("loop_points.value").c_str(), loop_points))
    return __LINE__;

  // this surface_type will relax parts of surface, those should be ignored in
  // polycube L1 optimization
  boost::unordered_map<size_t,size_t> restricted_surface_type;
  extract_relax_surface(tm.mesh_, tm.node_, restricted_surface_type, loop_points,
                        restricted_surface_type);

  dump_surface_restricted_type("after_relax_surface_type", restricted_surface_type);
  return 0;
}

int remove_surface_degenerated_patch_by_collapsing(ptree & pt)
{
  jtf::mesh::meshes tm;
  pt.put("tet.desc", "input tetmesh");
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm.node_, &tm.mesh_ ))
    return __LINE__;

  boost::unordered_map<size_t,size_t> surface_type;
  pt.put("surface_type.desc", "input surface_type");
  int rtn = load_surface_type(
        pt.get<string>("surface_type.value").c_str(), surface_type);
  if(rtn != 0)
    return __LINE__;

  boost::unordered_map<pair<size_t,size_t>, size_t> inner_face_jump_type;
  pt.put("inner_face_jump_type.desc", "input inner_face_jump_type");
  if(zjucad::has("inner_face_jump_type.value",pt)){
      if(load_inner_face_jump_type(
           pt.get<string>("inner_face_jump_type.value").c_str(),
           inner_face_jump_type))
        return __LINE__;
    }else{
      cerr << "# [info] pure polycube structure." << endl;
    }

  //  vector<pair<size_t,size_t> > g_unknown_face_pair;

  //  pt.put("g_unknown_face.desc", "input g_unknown_face");
  //  if(zjucad::has("g_unknown_face.value",pt)){
  //    if(load_g_unknown_face(pt.get<string>("g_unknown_face.value").c_str(),
  //                           g_unknown_face_pair)){
  //      return __LINE__;
  //    }
  //  }

  iterately_remove_surface_degenerated_patch(
        tm.mesh_, tm.node_, inner_face_jump_type, surface_type);

  {
    unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
    if(!fa.get()){
        cerr << "# [error] can not build facetet_adjacent." << endl;
        return __LINE__;
      }
    dump_surface_restricted_type_to_vtk(
          "after_iterately_remove_surface_type.vtk", "patch", tm.node_,
          *fa, surface_type);
  }

  {
    if(jtf::mesh::tet_mesh_write_to_zjumat("new_tet.tet", &tm.node_, &tm.mesh_))
      return __LINE__;
    dump_surface_restricted_type("surface_type_after_automatically_modification",
                                 surface_type);
    dump_inner_face_jump_type("inner_face_jump_type_after_automatically_modification",
                              inner_face_jump_type);
  }
  return 0;
}

int visualize_normal(ptree &pt) {
  string input_tet_file = pt.get<string>("input_tet_vtk.value");
  string output_surface_color_file = pt.get<string>("output_surface_color_vtk.value");
  string output_edge_color_file = pt.get<string>("output_edge_color_vtk.value");
  cout << output_surface_color_file << endl;
  cout << output_edge_color_file << endl;
  jtf::mesh::meshes tm;
  jtf::mesh::tet_mesh_read_from_vtk(input_tet_file.c_str(), &tm.node_, &tm.mesh_);
  if(jtf::mesh::tet_mesh_write_to_zjumat("vtk2tet.tet", &tm.node_, &tm.mesh_))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  matrixst outside_face, outside_face_idx;
  matrixd face_normal;
  get_outside_face(*fa,outside_face);
  get_outside_face_idx(*fa, outside_face_idx);
  jtf::mesh::cal_face_normal(outside_face, tm.node_, face_normal);
  jtf::tetmesh::orient_face_normal_outside_tetmesh(tm.mesh_,tm.node_,outside_face, outside_face_idx,*fa,face_normal);
  ofstream ofs(output_surface_color_file.c_str());
  tri2vtk(ofs, &tm.node_[0], tm.node_.size(2), &outside_face[0], outside_face.size(2));
  cout << "save surface triangle face " << tm.node_.size(2) << endl;
  matrixd rgba = ones<double>(4, outside_face.size(2));
  for(size_t i = 0; i < rgba.size(2); ++i) {
      for(size_t j = 0; j < 3; ++j) {
          rgba(j, i) = fabs(face_normal(j, i));//(face_normal(j, i) + 1.0) / 2.0;
        }
    }
  ofs << "CELL_DATA " << rgba.size(2) << endl;
  vtk_data_rgba(ofs, &rgba[0],rgba.size(2),"normal_color");
  ofs.close();
  cout << "save surface over!" << endl;
  unique_ptr<jtf::mesh::edge2cell_adjacent> tet_surface(jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!tet_surface.get()) {
      cerr << "# [error] build edge2cell_adjacent fail." << endl;
      return __LINE__;
    }
  vector<deque<pair<size_t, size_t> > > visualize_edge;
  const vector<pair<size_t, size_t> >& edge2tri = tet_surface->edge2cell_;
  const vector<pair<size_t, size_t> >& edge2vertex = tet_surface->edges_;
  std::vector<pair<size_t, size_t> > edge_info;
  for(size_t i = 0; i < edge2tri.size(); ++i) {
      if(edge2tri[i].first == -1 || edge2tri[i].second == -1) {
          edge_info.push_back(edge2vertex[i]);
          continue;
        }
      if(!dzw::is_same_orient_two_faces(face_normal(colon(), edge2tri[i].first), face_normal(colon(), edge2tri[i].second))) {
          edge_info.push_back(edge2vertex[i]);
        }
    }
  jtf::util::extract_chain_from_edges(edge_info, visualize_edge);
  cout << "save surface edge" << endl;

  dump_singularity_to_cylinder(output_edge_color_file.c_str(), tm.node_, visualize_edge, 0.002);
  cout << "save surface edge over!" << endl;
  return 0;
}
