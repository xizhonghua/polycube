#include "../common/def.h"
#include "../common/transition_type.h"
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
#include "../common/vtk.h"
#include <vector>
#include <fstream>
#include <set>

#include <jtflib/mesh/io.h>

using namespace std;

int validate_equation_file(int argc, char *argv[])
{
  if(argc != 3){
    cerr << "# [usage] validate_equation_file file cut_tet" << endl;
    return __LINE__;
  }

  matrixst cut_tet;
  matrixd cut_node;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &cut_node, &cut_tet)){
    cerr << "# [error] can not read tet file." << endl;
    return __LINE__;
  }

  ifstream ifs(argv[1]);
  if(ifs.fail()){
    cerr << "# [error] can not read equation file." << endl;
    return __LINE__;
  }

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(cut_tet));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  typedef vector<pair<size_t,double> > equation_type;
  vector<equation_type> equations;

  vector<size_t> face_to_vtk;
  vector<size_t> face_type_to_vtk;
  string trash;
  size_t item_num = 0;
  while(!ifs.eof()){
    ifs >> trash >> trash >> item_num;
    if(trash.empty())
      break;
    if(item_num != 4){
      cerr << "# [error] unsupport this equation." << endl;
      return __LINE__;
    }
    equation_type  temp(item_num);
    for(size_t i = 0; i < item_num; ++i){
      ifs >> temp[i].first;
    }
    for(size_t i = 0; i < item_num; ++i){
      ifs >> temp[i].second;
    }
    equations.push_back(temp);
    trash.clear();
  }

  assert(equations.size() % 6 == 0);
  matrixd from(3,2), to(3,2);

  for(size_t t = 0; t < equations.size() / 6; ++t){
    set<size_t> face_set;
    for(size_t i = 0; i < 3; ++i){
      from(i,0) = equations[t * 6 + 2 * i + 0][0].second;
      from(i,1) = equations[t * 6 + 2 * i + 0][1].second;
      to(i,0) = equations[t * 6 + 2 * i + 0][2].second;
      to(i,1) = equations[t * 6 + 2 * i + 0][3].second;
    }
    face_set.insert(equations[t * 6 + 0][0].first/3);
    face_set.insert(equations[t * 6 + 0][1].first/3);
    face_set.insert(equations[t * 6 + 1][0].first/3);
    face_set.insert(equations[t * 6 + 1][1].first/3);
    size_t type = -1;
    for(size_t i = 0; i < 24; ++i){
      if(norm(from - trans(type_transition2(i)) * to) < 1e-6){
        type = i;
        break;
      }
    }

    face_to_vtk.insert(face_to_vtk.end(), face_set.begin(), face_set.end());
    face_type_to_vtk.push_back(type);
  }

  ofstream ofs("equation_type.vtk");
  tri2vtk(ofs, &cut_node[0], cut_node.size(2), &face_to_vtk[0], face_to_vtk.size()/3);
  cell_data(ofs, &face_type_to_vtk[0], face_type_to_vtk.size(), "face_type");
  return 0;
}
