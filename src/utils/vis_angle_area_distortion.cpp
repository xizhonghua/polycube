#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <stack>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <jtflib/util/util.h>
#include <jtflib/util/vertex_connection.h>

#include "../numeric/util.h"
#include "../common/vtk.h"

using namespace std;
using namespace zjucad::matrix;

int cal_angle_distortion(const jtf::mesh::meshes & init_tri,
                         const jtf::mesh::meshes & def_tri,
                         matrix<double> & angle_distortion)
{
  const double PI = 2 * acos(0);
  assert(init_tri.mesh_.size(2) == def_tri.mesh_.size(2));
  angle_distortion.resize(def_tri.mesh_.size(2),1);
  vector<double> angle,abc(3);
  double init_area = 0, def_area = 0, total_def_area = 0;
  matrix<double> edges(3,3);

  for(size_t fi = 0; fi < init_tri.mesh_.size(2); ++fi){
      jtf::mesh::cal_face_angle(init_tri.mesh_(colon(),fi), init_tri.node_, angle);
      def_area = jtf::mesh::cal_face_area(def_tri.mesh_(colon(),fi), def_tri.node_);

      for(size_t i = 0; i < init_tri.mesh_.size(1); ++i)
        edges(colon(),i) =
            init_tri.node_(colon(), init_tri.mesh_((i+1)%init_tri.mesh_.size(1),fi))
            - init_tri.node_(colon(), init_tri.mesh_((i+2)%init_tri.mesh_.size(1),fi));

      angle_distortion[fi] = 1/(4.0* def_area);
      double a = 0;
      for(size_t i = 0; i < 3; ++i){
          const double len = dot(edges(colon(),i),edges(colon(),i));
          const double cot_angle = 1.0/tan(angle[i]/180.0*PI);
          a += cot_angle * len;
        }
      angle_distortion[fi] *= a;
      if(angle_distortion[fi] < 1){
          //cerr << "# [error] strange angle dis < 1." << endl;
          angle_distortion[fi] = 1.0; // may be numeric problem
        }
    }

  return 0;
}

int cal_area_distortion(const jtf::mesh::meshes & init_tri,
                        const jtf::mesh::meshes & def_tri,
                        matrix<double> & area_distortion)
{
  assert(init_tri.mesh_.size(2) == def_tri.mesh_.size(2));
  area_distortion.resize(def_tri.mesh_.size(2),1);
  double init_area = 0, def_area = 0;

  for(size_t fi = 0; fi < init_tri.mesh_.size(2); ++fi){
      init_area = jtf::mesh::cal_face_area(init_tri.mesh_(colon(),fi), init_tri.node_);
      def_area = jtf::mesh::cal_face_area(def_tri.mesh_(colon(),fi), def_tri.node_);
      area_distortion[fi] = 0.5*(init_area/def_area + def_area/init_area);
      if(area_distortion[fi] == std::numeric_limits<double>::infinity())
        area_distortion[fi] = 0;
    }

  return 0;
}


int vis_angle_area_distortion(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] vis_angle_area_distortion init_tri def_tri vtk_name" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes init_tri, def_tri;
  if(jtf::mesh::load_obj(argv[1], init_tri.mesh_, init_tri.node_))
    return __LINE__;

  if(jtf::mesh::load_obj(argv[2], def_tri.mesh_, def_tri.node_))
    return __LINE__;

  matrix<double> angle_distortion = zeros<double>(init_tri.mesh_.size(2),1);
  matrix<double> area_distortion = zeros<double>(init_tri.mesh_.size(2),1);
  matrix<double> area_init = zeros<double>(init_tri.mesh_.size(2),1);

  cal_angle_distortion(init_tri, def_tri, angle_distortion);
  cal_area_distortion(init_tri, def_tri, area_distortion);

  for(size_t fi = 0; fi < init_tri.mesh_.size(2); ++fi){
      area_init[fi] = jtf::mesh::cal_face_area(init_tri.mesh_(colon(),fi), init_tri.node_);
      if(area_init[fi] == std::numeric_limits<double>::infinity())
        area_init[fi] = 0;
    }

  string init_dis = argv[3];
  init_dis +=  "_init_dis.vtk";
  string def_dis = argv[3];
  def_dis += "_def_dis.vtk";
  ofstream ofs_init(init_dis.c_str());

  tri2vtk(ofs_init, &init_tri.node_[0], init_tri.node_.size(2), &init_tri.mesh_[0], init_tri.mesh_.size(2));
  cell_data(ofs_init, &angle_distortion[0], angle_distortion.size(), "angle_distion", "angle_table");
  vtk_data(ofs_init, &area_distortion[0], area_distortion.size(), "area_distion", "area_table");

  ofstream ofs_def(def_dis.c_str());
  tri2vtk(ofs_def, &def_tri.node_[0], def_tri.node_.size(2), &def_tri.mesh_[0], def_tri.mesh_.size(2));
  cell_data(ofs_def, &angle_distortion[0], angle_distortion.size(), "angle_distion", "angle_table");
  vtk_data(ofs_def, &area_distortion[0], area_distortion.size(), "area_distion", "area_table");

  const double total_area = std::accumulate(area_init.begin(), area_init.end(), 0.0);
  const double angle_dis_all = dot(angle_distortion, area_init) / total_area;
  const double avg_angle_dis =
      std::accumulate(angle_distortion.begin(), angle_distortion.end(), 0.0)
      / angle_distortion.size();
  const double area_dis_all = dot(area_distortion, area_init) / total_area;
  const double avg_area_dis =
      std::accumulate(area_distortion.begin(), area_distortion.end(), 0.0)
      / area_distortion.size();
  cerr << "# [info] angle_distortion " << angle_dis_all << endl;
  cerr << "# [info] area_distion " << area_dis_all << endl;
  return 0;

}

int find_shortest_path(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [usage] find_shortest_path input_obj input_fl" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes trimesh;
  vector<vector<size_t> > feature_line;
  if(jtf::mesh::load_obj(argv[1], trimesh.mesh_, trimesh.node_))
    return __LINE__;
  if(jtf::mesh::load_feature_line(argv[2], feature_line))
    return __LINE__;

  std::shared_ptr<jtf::mesh::edge2cell_adjacent> ec(
        jtf::mesh::edge2cell_adjacent::create(trimesh.mesh_));
  if(!ec.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }

  map<pair<size_t,size_t>,double> edge_weight;
  for(size_t ei = 0; ei < ec->edges_.size(); ++ei){
      const pair<size_t,size_t> & one_edge = ec->edges_[ei];
      edge_weight[ec->edges_[ei]] =
          sqrt(dot(trimesh.node_(colon(), one_edge.first) - trimesh.node_(colon(), one_edge.second),
                   trimesh.node_(colon(), one_edge.first) - trimesh.node_(colon(), one_edge.second)));
    }

  std::shared_ptr<vertex_connection<UNDIRECT> > vc(
        vertex_connection<UNDIRECT>::create(edge_weight));
  if(!vc.get()){
      cerr << "# [error] can not build vertex_connection." << endl;
      return __LINE__;
    }

  vector<vector<size_t> > new_feature_lines(feature_line.size());
  for(size_t fi = 0; fi < feature_line.size(); ++fi){
      const vector<size_t> & one_feature_line = feature_line[fi];
      vc->get_shortest_path(one_feature_line.front(),one_feature_line.back(),
                            new_feature_lines[fi]);
    }

  ofstream ofs("new_feature_line.fl");
  ofs << new_feature_lines.size() << endl;
  for(size_t fi = 0; fi < new_feature_lines.size(); ++fi){
      ofs << new_feature_lines[fi].size() << endl;
      for(size_t i = 0; i < new_feature_lines[fi].size(); ++i){
          ofs << new_feature_lines[fi][i] << " ";
        }
      ofs << endl;
    }
  return 0;
}

static int load_surface_patch(
    const char * surface_patch_file,
    const matrix<size_t> & trimesh,
    std::vector<matrix<size_t> > &patch_faces)
{
  ifstream ifs_patch(surface_patch_file);
  if(ifs_patch.fail()){
      cerr << "# [error] can not open patch surface file." << endl;
      return __LINE__;
    }

  size_t patch_num = 0;
  ifs_patch >> patch_num;

  cerr << "# [info] patch number " << patch_num << endl;

  if(patch_num == 0) return 0;

  patch_faces.resize(patch_num);
  size_t one_patch_face_num, face_idx, trash,trash_b;
  for(size_t pi = 0; pi < patch_num; ++pi){
      ifs_patch >> one_patch_face_num >> trash;
      patch_faces[pi].resize(trimesh.size(1), one_patch_face_num);
      for(size_t i = 0; i < one_patch_face_num; ++i){
          ifs_patch >> face_idx;
          patch_faces[pi](colon(),i) = trimesh(colon(),face_idx);
        }
      for(size_t i = 0; i < trash; ++i) ifs_patch >> trash_b;
    }

  return 0;
}

int snap_each_patch(const matrix<size_t> &patch_faces,
                    const map<size_t, matrix<double> > &handle_target_coord,
                    const jtf::mesh::meshes &init_tri,
                    jtf::mesh::meshes &output_tri)
{
  using namespace jtf::mesh;
  std::shared_ptr<edge2cell_adjacent> ea(edge2cell_adjacent::create(patch_faces));
  if(!ea.get()) {
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }
  matrix<size_t> all_edges;
  get_boundary_edge(*ea, all_edges);
  vector<pair<size_t,size_t> > edge2pair;
  for(size_t ei = 0; ei < all_edges.size(2); ++ei)
    edge2pair.push_back(make_pair(all_edges(0,ei), all_edges(1,ei)));

  vector<deque<pair<size_t,size_t> > > chain;
  jtf::util::extract_chain_from_edges(edge2pair, chain);
  for(map<size_t,matrix<double> >::const_iterator cit = handle_target_coord.begin();
      cit != handle_target_coord.end(); ++cit){
      output_tri.node_(colon(), cit->first) = cit->second;
    }

  set<size_t> fixed_points;

  for(size_t ei = 0; ei < chain.front().size(); ++ei){
      const pair<size_t,size_t> & one_edge = chain.front()[ei];
      fixed_points.insert(one_edge.first);
      fixed_points.insert(one_edge.second);
    }

  deque<pair<size_t,size_t> > edges_chain_copy = chain.front();
  deque<size_t> one_queue;

  while(!edges_chain_copy.empty()){
      pair<size_t,size_t> & one_edge = edges_chain_copy.front();
      map<size_t,matrix<double> >::const_iterator cit =
          handle_target_coord.find(one_edge.first);
      if(cit == handle_target_coord.end()){
          if(one_queue.empty()){
              edges_chain_copy.push_back(one_edge);
              edges_chain_copy.pop_front();
              continue;
            }
          one_queue.push_back(one_edge.first);
          edges_chain_copy.pop_front();

          if(edges_chain_copy.empty()
             && handle_target_coord.find(one_edge.second) != handle_target_coord.end())
            edges_chain_copy.push_back(make_pair(one_edge.second, -1));
          continue;
        }else{
          if(one_queue.empty()){
              one_queue.push_back(one_edge.first);
              edges_chain_copy.pop_front();
              continue;
            }
          {// one queue is finished
            const size_t N = one_queue.size();
            const double unit_len = norm(output_tri.node_(colon(),one_queue.front())
                                         - output_tri.node_(colon(), one_edge.first))
                / N;
            matrix<double> direction = output_tri.node_(colon(),one_edge.first)
                - output_tri.node_(colon(), one_queue.front());
            direction /= norm(direction);
            for(size_t i = 1; i < one_queue.size(); ++i)
              output_tri.node_(colon(), one_queue[i]) = output_tri.node_(colon(),one_queue.front())
                  + direction * unit_len * i;
            one_queue.clear();
            one_queue.push_back(one_edge.first);
            edges_chain_copy.pop_front();
          }
        }

    }

  std::shared_ptr<one_ring_point_at_point> orpap(
        one_ring_point_at_point::create(patch_faces));
  if(!orpap.get()){
      cerr << "# [error] can not build one_ring_point_at_point." << endl;
      return __LINE__;
    }

  vector<double> weight;
  matrix<double> edges(3,3);
  for(size_t it = 0; it < 1000; ++it)  {
      for(one_ring_point_at_point::p2p_type::const_iterator cit =
          orpap->p2p_.begin(); cit != orpap->p2p_.end(); ++cit){
          if(fixed_points.find(cit->first) != fixed_points.end()) continue;
          const vector<size_t> & one_ring_points = cit->second;
          weight.resize(one_ring_points.size());

          for(size_t i = 0; i < one_ring_points.size(); ++i){
              size_t edge_idx = ea->get_edge_idx(cit->first, one_ring_points[i]);
              assert(edge_idx != -1);
              const pair<size_t,size_t> &adj_faces = ea->edge2cell_[edge_idx];
              assert(adj_faces.first != -1 && adj_faces.second != -1);
              size_t other_point_idx_0, other_point_idx_1;
              other_point_idx_0 = std::accumulate(
                    patch_faces(colon(),adj_faces.first).begin(),
                    patch_faces(colon(),adj_faces.first).end(),0)
                  - cit->first - one_ring_points[i];

              other_point_idx_1 = std::accumulate(
                    patch_faces(colon(),adj_faces.second).begin(),
                    patch_faces(colon(),adj_faces.second).end(),0)
                  - cit->first - one_ring_points[i];

              edges(colon(),0) = init_tri.node_(colon(),one_ring_points[i])
                  - init_tri.node_(colon(),cit->first);
              edges(colon(),1) = init_tri.node_(colon(),other_point_idx_0)
                  - init_tri.node_(colon(),cit->first);
              edges(colon(),2) = init_tri.node_(colon(),other_point_idx_1)
                  - init_tri.node_(colon(),cit->first);

              // tan(\theta/2) = sin(\theta)/(cos(\theta)+1)
              const double cos_01 = dot(edges(colon(),0),edges(colon(),1))
                  /(norm(edges(colon(),0)) * norm(edges(colon(),1)));
              const double cos_02 = dot(edges(colon(),0),edges(colon(),2))
                  /(norm(edges(colon(),0)) * norm(edges(colon(),2)));
              const double sin_01 = sqrt(1.0-cos_01*cos_01);
              const double sin_02 = sqrt(1.0-cos_02*cos_02);

              weight[i] = (sin_01/(cos_01+1.0) + sin_02/(cos_02+1.0))/(norm(edges(colon(),0)));
            }

          output_tri.node_(colon(),cit->first) *= 0;
          const double total_w = std::accumulate(weight.begin(), weight.end(),0.0);
          for(size_t i = 0; i < weight.size(); ++i)
            output_tri.node_(colon(),cit->first) +=
                weight[i] / total_w
                * output_tri.node_(colon(),one_ring_points[i]);
        }
    }

  return 0;
}

int lin_para(int argc, char * argv[])
{
  using namespace jtf::mesh;
  if(argc != 5){
      cerr << "# [usage] lin_para input_obj input_polycube handle_map patch" << endl;
      return __LINE__;
    }

  meshes init_tri, init_polycube;
  map<size_t,size_t> handle_map;
  if(load_obj(argv[1], init_tri.mesh_, init_tri.node_))
    return __LINE__;

  if(load_obj(argv[2], init_polycube.mesh_, init_polycube.node_))
    return __LINE__;

  //read handle map
  {
    ifstream ifs(argv[3]);
    if(ifs.fail()) {
        cerr << "# [error] can not open file " << argv[3] << endl;
        return __LINE__;
      }
    size_t handle_number = 0;
    ifs >> handle_number;
    size_t a,b;
    for(size_t i = 0; i < handle_number; ++i){
        ifs >> a >> b;
        handle_map[a] = b;
      }
  }
  map<size_t, matrix<double> > handle_target_coord;
  for(map<size_t,size_t>::const_iterator cit = handle_map.begin();
      cit != handle_map.end(); ++cit){
      handle_target_coord[cit->second] = init_polycube.node_(colon(), cit->first);
    }

  std::vector<matrix<size_t> > patch_faces;

  if(load_surface_patch(argv[4], init_tri.mesh_, patch_faces))
    return __LINE__;

  {
    ofstream ofs("patch.vtk");
    vector<size_t> faces;
    vector<size_t> faces_type;
    for(size_t pi = 0; pi < patch_faces.size(); ++pi){
        faces.insert(faces.end(), patch_faces[pi].begin(), patch_faces[pi].end());
        for(size_t fi = 0; fi < patch_faces[pi].size(2); ++fi)
          faces_type.push_back(pi);
      }
    tri2vtk(ofs, &init_tri.node_[0], init_tri.node_.size(), &faces[0], faces.size()/3);
    cell_data(ofs, &faces_type[0], faces_type.size(), " face_type");
  }

  // for each patch
  meshes output_tri = init_tri;
  for(size_t pi = 0; pi < patch_faces.size(); ++pi){
      snap_each_patch(patch_faces[pi], handle_target_coord, init_tri, output_tri);
    }

  if(save_obj("output.obj", output_tri.mesh_, output_tri.node_))
    return  __LINE__;

  return 0;
}
