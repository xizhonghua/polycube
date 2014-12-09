#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
#include "../common/zyz.h"
#include "../common/IO.h"
#include <fstream>
#include <istream>

using namespace std;

int tet_tf2yf(const string & input_tet,
              const string & output_tet)
{
  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(input_tet.c_str(), &tm.node_, &tm.mesh_))
    return __LINE__;

  ofstream ofs(output_tet.c_str());
  ofs << tm.node_.size(2) << " vertices" << endl;
  ofs << tm.mesh_.size(2) << " cells" << endl;

  for(size_t ni = 0; ni < tm.node_.size(2); ++ni){
    for(size_t di = 0; di < 3; ++di){
      ofs << tm.node_(di, ni) << " ";
    }
    ofs << endl;
  }

  for(size_t ti = 0; ti < tm.mesh_.size(2); ++ti){
    ofs << tm.mesh_.size(1) << " ";
    for(size_t di = 0; di < tm.mesh_.size(1); ++di){
      ofs << tm.mesh_(di,ti) << " ";
    }
    ofs << endl;
  }

  return 0;
}

int fld_tf2yf(const string & input_fld,
              const string & output_fld)
{
  matrixd zyz;
  if(jtf::mesh::read_matrix(input_fld.c_str(), zyz)) {
    cerr << "# [error] can not read zyz file." << endl;
    return __LINE__;
  }

  ofstream ofs(output_fld.c_str());
  ofs << zyz.size(2) << endl;

  matrixd frame(3,3);
  for(size_t ti = 0; ti < zyz.size(2); ++ti){
    zyz_angle_2_rotation_matrix1(&zyz(0,ti), &frame[0]);
    for(size_t i = 0; i < frame.size(); ++i){
      ofs << frame[i] << " " ;
    }
    ofs << endl;
  }

  return 0;
}

int tf2yf(int argc, char *argv[])
{
  if(argc != 4){
    cerr << "# [usage] tf2yf tet/fld input_tet[input_fld] output_tet[output_fld]."  << endl;
    return __LINE__;
  }

  const std::string file = argv[1];

  if(file == "tet"){
    tet_tf2yf(argv[2], argv[3]);
  }else if(file == "fld"){
    fld_tf2yf(argv[2], argv[3]);
  }else {
    cerr << "# [error] unsupported file." << endl;
    return __LINE__;
  }

  return 0;
}
