#include <numeric>
#include <fstream>

#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/io.h>

#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"
#include "../common/vtk.h"
#include "../hex_param/io.h"

using namespace std;

int check_polycube_surface(const matrixst & tet,
                           const matrixd & original_node,
                           const matrixd & polycube_node,
                           const matrixst * orig_tet_ptr = 0)
{
  using namespace zjucad::matrix;
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
    cerr << "# [error] can not build face2tet_adjacnet." << endl;
    return __LINE__;
  }

  matrixst outside_face, outside_face_idx;
  get_outside_face(*fa, outside_face);
  get_outside_face_idx(*fa, outside_face_idx);
  matrixd face_normal_in_polycube;

  jtf::mesh::cal_face_normal(outside_face, polycube_node, face_normal_in_polycube);
  jtf::tetmesh::orient_face_normal_outside_tetmesh(
        tet, polycube_node, outside_face, outside_face_idx, *fa,
        face_normal_in_polycube);

  // if fa is inner face in orig_tet, then is_orig_outside_face will be -1
  // or will be the face idx in orig tet
  matrix<size_t> is_orig_outside_face = outside_face_idx;

  if(orig_tet_ptr){
    if((*orig_tet_ptr).size(2) != tet.size(2)){
      cerr << "# [error] orig_tet is not compatible with given tet." << endl;
      return __LINE__;
    }
    matrixst cut_tet2tet(max(tet)+1);
    cut_tet2tet(tet) = (*orig_tet_ptr)(colon());
    unique_ptr<jtf::mesh::face2tet_adjacent> fa_orig(jtf::mesh::face2tet_adjacent::create(*orig_tet_ptr));
    if(!fa_orig.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }

    for(size_t fi = 0; fi < outside_face.size(2); ++fi){
      const size_t face_idx =
          fa_orig->get_face_idx(
            cut_tet2tet[outside_face(0,fi)],  cut_tet2tet[outside_face(1,fi)],
            cut_tet2tet[outside_face(2,fi)]);
      if(face_idx == -1){
        cerr << "# [error] strange face " << outside_face(0,fi) << " "
             << outside_face(1,fi) << " " << outside_face(2,fi)
             << " with orig_idx " << cut_tet2tet[outside_face(0,fi)] << " "
             << cut_tet2tet[outside_face(1,fi)] << " "
             << cut_tet2tet[outside_face(2,fi)]  << "is not in orig tet." << endl;
        return __LINE__;
      }
      const pair<size_t,size_t> & tet_pair = fa_orig->face2tet_[face_idx];
      if(!fa_orig->is_outside_face(tet_pair))
        is_orig_outside_face[fi] = -1;
      else
        is_orig_outside_face[fi] = face_idx;
    }
  }

  const matrixd axes = eye<double>(3);
  vector<pair<double,int> > axis_choice(6);
  //vector<pair<double>>
  matrixst outside_face_type(outside_face.size(2));
  matrixst outside_face_type_with_orientation(outside_face.size(2));

  for(size_t fi = 0; fi < outside_face.size(2); ++fi){
    if(is_orig_outside_face[fi] == -1) continue;
    for(size_t ai = 0; ai < 3; ++ai){
      axis_choice[2*ai] =
          make_pair(
            dot(face_normal_in_polycube(colon(),fi), axes(colon(),ai)),
            2*ai);
      axis_choice[2*ai+1] =
          make_pair(
            -1*dot(face_normal_in_polycube(colon(),fi), axes(colon(),ai)),
            2*ai+1);
    }
    sort(axis_choice.begin(), axis_choice.end());
    outside_face_type[fi] = axis_choice.back().second/2;
    outside_face_type_with_orientation[fi] = axis_choice.back().second;
  }

  matrixd outside_face_volume = zeros<double>(outside_face.size(2), 1);
  for(size_t fi = 0; fi < outside_face.size(2); ++fi){
    if(fi == 7)
      cerr << endl;
    const size_t & face_idx = outside_face_idx[fi];
    const pair<size_t,size_t> tet_pair = fa->query(&outside_face(0,fi));
    assert(fa->is_outside_face(tet_pair));
    const size_t & tet_idx = (tet_pair.first == -1? tet_pair.second:tet_pair.first);
    outside_face_volume[fi] = jtf::mesh::cal_tet_vol(polycube_node(colon(), tet(colon(),tet_idx)));
  }

  const double average_vol =
      std::accumulate(outside_face_volume.begin(), outside_face_volume.end(), 0.0)
      / outside_face.size(2);
  for(size_t fi = 0; fi < outside_face.size(2); ++fi)
    outside_face_volume[fi] /= average_vol ;

  {
    vector<size_t> packaged_outside_face, packaged_outside_face_type,
        packaged_outside_face_type_with_orientation;
    packaged_outside_face.reserve(outside_face.size());
    packaged_outside_face_type.reserve(outside_face.size(2));
    packaged_outside_face_type_with_orientation.reserve(outside_face.size(2));

    for(size_t fi = 0; fi < is_orig_outside_face.size(); ++fi){
      if(is_orig_outside_face[fi] == -1) continue;
      packaged_outside_face.insert(packaged_outside_face.end(),
                                   outside_face(colon(),fi).begin(),
                                   outside_face(colon(),fi).end());
      //if(outside_face_volume[fi] > 0){
        packaged_outside_face_type.push_back(outside_face_type[fi]);
        packaged_outside_face_type_with_orientation.push_back(outside_face_type_with_orientation[fi]);
//      }else if(outside_face_volume[fi] < 0){
//        packaged_outside_face_type.push_back(3);
//        packaged_outside_face_type_with_orientation.push_back(3);
//      }else if(fabs(outside_face_volume[fi]) < 1e-3){
//        packaged_outside_face_type.push_back(4);
//        packaged_outside_face_type_with_orientation.push_back(4);
//      }
    }
    ofstream ofs("surface_type_in_original_tet.vtk");
    tri2vtk(ofs, &original_node[0], original_node.size(2),
            &packaged_outside_face[0], packaged_outside_face.size()/3);
    cell_data(ofs, &packaged_outside_face_type[0],
              packaged_outside_face_type.size(), "surface_type");

    ofstream ofs_p("surface_type_in_polycube_tet.vtk");
    tri2vtk(ofs_p, &polycube_node[0], polycube_node.size(2),
            &packaged_outside_face[0], packaged_outside_face.size()/3);
    cell_data(ofs_p, &packaged_outside_face_type[0],
              packaged_outside_face_type.size(), "surface_type");

    ofstream ofs_wo("surface_type_with_orientation_in_original_tet.vtk");
    tri2vtk(ofs_wo, &original_node[0], original_node.size(2),
            &packaged_outside_face[0], packaged_outside_face.size()/3);
    cell_data(ofs_wo, &packaged_outside_face_type_with_orientation[0],
              packaged_outside_face_type_with_orientation.size(), "surface_type");

    ofstream ofs_p_wo("surface_type_with_orientation_in_polycube_tet.vtk");
    tri2vtk(ofs_p_wo, &polycube_node[0], polycube_node.size(2),
            &packaged_outside_face[0], packaged_outside_face.size()/3);
    cell_data(ofs_p_wo, &packaged_outside_face_type_with_orientation[0],
              packaged_outside_face_type_with_orientation.size(), "surface_type");

  }

  { // dump to file
    ofstream ofs("surface_type_from_polycube");
    for(size_t fi = 0; fi  < outside_face_type.size(); ++fi){
      if(is_orig_outside_face[fi] == -1) continue;
      ofs << outside_face_idx[fi] << " " << outside_face_type[fi] << endl;
    }
    if(orig_tet_ptr){
      ofstream ofs_orig("surface_type_from_polycube_in_orig_tet");
      for(size_t fi = 0; fi  < outside_face_type.size(); ++fi){
        if(is_orig_outside_face[fi] == -1) continue;
        ofs_orig << is_orig_outside_face[fi] << " " << outside_face_type[fi] << endl;
      }
    }
  }

  { // dump to file
    ofstream ofs("surface_type_with_orientation_from_polycube");
    for(size_t fi = 0; fi  < outside_face_type_with_orientation.size(); ++fi){
      if(is_orig_outside_face[fi] == -1) continue;
      ofs << outside_face_idx[fi] << " " << outside_face_type_with_orientation[fi] << endl;
    }
    if(orig_tet_ptr){
      ofstream ofs_orig("surface_type_with_orientation_from_polycube_in_orig_tet");
      for(size_t fi = 0; fi  < outside_face_type_with_orientation.size(); ++fi){
        if(is_orig_outside_face[fi] == -1) continue;
        ofs_orig << is_orig_outside_face[fi] << " " << outside_face_type_with_orientation[fi] << endl;
      }
    }
  }

  return 0;
}

int check_polycube_surface(int argc, char *argv[])
{
  if(argc != 3 && argc != 4){
    cerr << "# [usage] check_polycube_surface cut_tet polycube_tet [orig_tet]." << endl;
    return __LINE__;
  }

  jtf::mesh::meshes cut_tm, polycube_tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &cut_tm.node_, &cut_tm.mesh_)){
    return __LINE__;
  }

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &polycube_tm.node_, &polycube_tm.mesh_)){
    return __LINE__;
  }

  if(cut_tm.mesh_.size() != polycube_tm.mesh_.size()){
    cerr << "# [error] input tet is not the same as polycube_tet in topology. "
         << fabs(norm(cut_tm.mesh_ - polycube_tm.mesh_)) << endl;
    return __LINE__;
  }

  if(argc == 3)
    check_polycube_surface(cut_tm.mesh_, cut_tm.node_, polycube_tm.node_);
  else{
    jtf::mesh::meshes orig_tet;
    if(jtf::mesh::tet_mesh_read_from_zjumat(argv[3], &orig_tet.node_, &orig_tet.mesh_))
      return __LINE__;
    check_polycube_surface(cut_tm.mesh_, cut_tm.node_, polycube_tm.node_, &orig_tet.mesh_);
  }

  return 0;
}
