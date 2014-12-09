#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"
using namespace std;
using namespace zjucad::matrix;

int load_range_file(const char * filename,
                    zjucad::matrix::matrix<double> &bb)
{
  bb.resize(3,2);
  bb(colon(),0) = -1*std::numeric_limits<double>::min();
  bb(colon(),1) = std::numeric_limits<double>::max();

  ifstream ifs(filename);
  string xyz, flag;
  double r;
  size_t xyz_idx = -1, range_idx = -1;
  while(!ifs.eof()){
      ifs >> xyz >> flag >> r;
      if(ifs.eof()) return 0;
      if(xyz == "x" || xyz == "X") xyz_idx = 0;
      if(xyz == "y" || xyz == "Y") xyz_idx = 1;
      if(xyz == "z" || xyz == "Z") xyz_idx = 2;
      if(flag == ">"){ range_idx = 0;}
      else if(flag == "<") range_idx = 1;
      else {
          cerr << "# [error] can not handle " << flag << "symbol." << endl;
          return __LINE__;
        }
      bb(xyz_idx,range_idx) = r;
    }
  return 0;
}

int select_edges(int argc, char *argv[])
{
  if(argc != 4){
      cerr << "# [usage] select_edges input_tet range_file select_edges" << endl;
      return __LINE__;
    }

  jtf::tet_mesh tm(argv[1]);
  matrix<double> bb(3,2);
  if(load_range_file(argv[2], bb))
    return __LINE__;

  cerr << "range " << bb << endl;
  vector<size_t> select_edges;
  matrix<double> mid_point;
  for(const auto & one_edge : tm.ortae_.e2t_){
      const pair<size_t,size_t> &one_edge_idx = one_edge.first;
      mid_point = (tm.tetmesh_.node_(colon(), one_edge_idx.first)
                   + tm.tetmesh_.node_(colon(), one_edge_idx.second))/2.0;
      bool selected = true;
      for(size_t j = 0; j < tm.tetmesh_.node_.size(1); ++j){
          if(mid_point[j] < bb(j,0) || mid_point[j] > bb(j,1)){
              selected = false; break;
            }
        }
      if(selected){
          select_edges.push_back(one_edge_idx.first);
          select_edges.push_back(one_edge_idx.second);
        }
    }

  {
    ofstream ofs("select_edge.vtk");
    line2vtk(ofs, &tm.tetmesh_.node_[0], tm.tetmesh_.node_.size(2),
        &select_edges[0], select_edges.size()/2);
  }
  ofstream ofs("select_edges");
  ofs << select_edges.size()/2 << endl;
  for(size_t i = 0 ;i < select_edges.size()/2; ++i){
      ofs << select_edges[i*2+0] << " " << select_edges[i*2+1] << endl;
    }
  return 0;
}
