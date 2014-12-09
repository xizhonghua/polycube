#include <fstream>
#include <iostream>

#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>

#include "../common/vtk.h"
#include "../common/def.h"

using namespace std;
using namespace zjucad::matrix;

int dump_fv_to_cross(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] dump_fv_to_cross obj fv vtk" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes trim;
  if(jtf::mesh::load_obj(argv[1], trim.mesh_, trim.node_)){
      return __LINE__;
    }

  vector<matrixd> face_vector;
  ifstream ifs(argv[2]);
  if(ifs.fail()){
      cerr << "# [error] can not open fv file." << endl;
      return __LINE__;
    }

  size_t face_num = 0;
  ifs >> face_num;
  if(face_num != trim.mesh_.size(2)){
      cerr << "# [error] wrong fv." << endl;
      cerr << "# tri size " << trim.mesh_.size(2) << " fv size " << face_num << endl;
      return __LINE__;
    }

  face_vector.resize(face_num);
  for(size_t fi = 0; fi < face_num; ++fi){
      face_vector[fi].resize(3,2);
      for(size_t i = 0; i < 2; ++i){
          ifs >> face_vector[fi](0,i) >> face_vector[fi](1,i) >> face_vector[fi](2,i);
        }
    }

  const double av_len = jtf::mesh::cal_average_edge(trim.mesh_,trim.node_);
  vector<size_t> edges;
  edges.reserve(trim.mesh_.size(2) * 4);
  matrixd new_node(3, trim.mesh_.size(2) * 4);
  matrixd center;
  for(size_t tri = 0; tri < trim.mesh_.size(2); ++tri){
      center = zeros<double>(3,1);
      for(size_t pi = 0; pi < trim.mesh_.size(1); ++pi){
          center += trim.node_(colon(), trim.mesh_(pi, tri));
        }
      center /= trim.mesh_.size(1);
      new_node(colon(), 4 * tri + 0) = center + face_vector[tri](colon(),0) * av_len/3;
      new_node(colon(), 4 * tri + 1) = center - face_vector[tri](colon(),0) * av_len/3;
      new_node(colon(), 4 * tri + 2) = center + face_vector[tri](colon(),1) * av_len/3;
      new_node(colon(), 4 * tri + 3) = center - face_vector[tri](colon(),1) * av_len/3;
      edges.push_back(4 * tri + 0);
      edges.push_back(4 * tri + 1);
      edges.push_back(4 * tri + 2);
      edges.push_back(4 * tri + 3);
    }

  {
    ofstream ofs(argv[3]);
    if(ofs.fail()){
        cerr << "# [error] can not open vtk file." << endl;
        return __LINE__;
      }
    line2vtk(ofs , &new_node[0], new_node.size(2),&edges[0], edges.size()/2);
  }

  return 0;
}

int dump_pv_to_cross(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] dump_pv_to_cross obj pv vtk" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes trim;
  if(jtf::mesh::load_obj(argv[1], trim.mesh_, trim.node_)){
      return __LINE__;
    }

  vector<matrixd> point_vector;
  ifstream ifs(argv[2]);
  if(ifs.fail()){
      cerr << "# [error] can not open fv file." << endl;
      return __LINE__;
    }

  size_t point_num = 0;
  size_t N = 0;
  ifs >> point_num >> N;
  if(point_num != trim.node_.size(2)){
      cerr << "# [error] wrong pv." << endl;
      cerr << "# point size " << trim.node_.size(2) << " pv size " << point_num << endl;
      return __LINE__;
    }

  matrix<double> point_normal;
  jtf::mesh::cal_point_normal(trim.mesh_, trim.node_, point_normal);

  point_vector.resize(point_num);
  for(size_t pi = 0; pi < point_num; ++pi){
      point_vector[pi].resize(3,2);
      ifs >> point_vector[pi][0] >> point_vector[pi][1] >> point_vector[pi][2];
      point_vector[pi](colon(),1) = cross(point_normal(colon(),pi), point_vector[pi](colon(),0));
    }

  const double av_len = jtf::mesh::cal_average_edge(trim.mesh_,trim.node_);
  vector<size_t> edges;
  edges.reserve(trim.mesh_.size(2) * 4);
  matrixd new_node(3, trim.mesh_.size(2) * 4);
  matrixd center;
  for(size_t pi = 0; pi < trim.node_.size(2); ++pi){
      center = trim.node_(colon(),pi);
      new_node(colon(), 4 * pi + 0) = center + point_vector[pi](colon(),0) * av_len/3;
      new_node(colon(), 4 * pi + 1) = center - point_vector[pi](colon(),0) * av_len/3;
      new_node(colon(), 4 * pi + 2) = center + point_vector[pi](colon(),1) * av_len/3;
      new_node(colon(), 4 * pi + 3) = center - point_vector[pi](colon(),1) * av_len/3;
      edges.push_back(4 * pi + 0);
      edges.push_back(4 * pi + 1);
      edges.push_back(4 * pi + 2);
      edges.push_back(4 * pi + 3);
    }

  {
    ofstream ofs(argv[3]);
    if(ofs.fail()){
        cerr << "# [error] can not open vtk file." << endl;
        return __LINE__;
      }
    line2vtk(ofs , &new_node[0], new_node.size(2),&edges[0], edges.size()/2);
  }

  return 0;
}

int dump_fv_to_N_field(int argc, char *argv[])
{
  if(argc != 3){
      cerr << "# [usage] dump_fv_to_N_field fv N_field" << endl;
      return __LINE__;
    }
  vector<matrixd> face_vector;
  ifstream ifs(argv[1]);
  if(ifs.fail()){
      cerr << "# [error] can not open fv file." << endl;
      return __LINE__;
    }

  size_t face_num = 0;
  ifs >> face_num;
  face_vector.resize(face_num);
  for(size_t fi = 0; fi < face_num; ++fi){
      face_vector[fi].resize(3,2);
      for(size_t i = 0; i < 2; ++i){
          ifs >> face_vector[fi](0,i) >> face_vector[fi](1,i) >> face_vector[fi](2,i);
        }
    }

  ofstream ofs(argv[2]);
  ofs << face_num << " " << 4 << endl;
  for(size_t fi = 0; fi < face_vector.size(); ++fi)
      ofs << face_vector[fi][0] << " " << face_vector[fi][1] << " " << face_vector[fi][2] << endl;
  return 0;
}
