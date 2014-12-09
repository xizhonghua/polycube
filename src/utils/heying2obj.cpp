#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>

using namespace std;
using namespace jtf::mesh;
using namespace zjucad::matrix;

int load_m_file(const char * file,
                jtf::mesh::meshes &polycube_m,
                jtf::mesh::meshes &orig_m)
{
  ifstream ifs(file);
  if(ifs.fail()){
      cerr << "# [error] can not open m file." << endl;
      return __LINE__;
    }

  string trash;
  size_t vertex_num = 0, face_num = 0;
  while(!ifs.eof()){
      ifs >> trash;
      if(ifs.eof()) break;
      if(trash == "Vertex") ++vertex_num;
      if(trash == "Face") ++face_num;
    }

  polycube_m.node_.resize(3, vertex_num);
  polycube_m.mesh_.resize(3, face_num);

  orig_m.node_.resize(3, vertex_num);
  orig_m.mesh_.resize(3, face_num);

  ifs.clear();
  ifs.seekg(0,std::ios::beg);

  size_t vid = 0;
  size_t vid_o = 0;
  size_t fid = 0;
  while(!ifs.eof()){
      ifs >> trash;
      if(ifs.eof()) break;
      if(trash == "Vertex"){
          ifs >> vid;
          ifs >> polycube_m.node_(0,vid-1)
              >> polycube_m.node_(1,vid-1)
              >> polycube_m.node_(2,vid-1);
          continue;
        }
      if(trash == "Face"){
          ifs >> fid;
          ifs >> polycube_m.mesh_(0,fid-1)
              >> polycube_m.mesh_(1,fid-1)
              >> polycube_m.mesh_(2,fid-1);
          continue;
        }
      // find "Opos=("
      size_t pos = trash.find("Opos=(");
      if(pos == std::string::npos) continue;
      {
        std::stringstream os_value;
        os_value << trash.substr(pos+6, trash.size()-6);
        os_value >> orig_m.node_(0,vid_o) ;

        ifs >> orig_m.node_(1, vid_o);
        ifs >> trash;

        //assert(trash.back() == ')');
        os_value.clear();
        os_value << trash.substr(0, trash.size()-1);
        os_value >> orig_m.node_(2, vid_o);
        ++vid_o;
        continue;
      }
    }

  polycube_m.mesh_ -= 1;
  orig_m.mesh_ = polycube_m.mesh_;

  return 0;
}

int heying2obj(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [usage] heying2obj heying.m output_obj" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes polycube_m, orig_m;
  if(load_m_file(argv[1], polycube_m, orig_m))
    return __LINE__;

  string obj_file = argv[2];
  string polycube_obj = "polycube_" + obj_file;
  string orig_obj = "orig_" + obj_file;

  matrix<size_t> quad_mesh = zeros<size_t>(4, polycube_m.mesh_.size(2)/2);
  if(save_obj(polycube_obj.c_str(), polycube_m.mesh_, polycube_m.node_))
    return __LINE__;

  if(save_obj(orig_obj.c_str(), orig_m.mesh_, orig_m.node_))
    return __LINE__;

  return 0;
}
