#include "../tetmesh/tetmesh.h"

using namespace std;
using namespace zjucad::matrix;


int load_singularity_vtk(const char * file,
                         vector<pair<size_t,size_t> > & singularities,
                         matrix<double> & node)
{
  ifstream ifs(file);
  if(ifs.fail()) {
      cerr << "[info] " << "can not open file" << file << endl;
      return __LINE__;
    }


  string str;
  int point_num = 0,cell_num = 0;

  while(!ifs.eof()){
      ifs >> str;
      if(str == "POINTS"){
          ifs >> point_num >> str;
          node.resize(3, point_num);
          for(size_t i = 0;i < point_num; ++i){
              for(size_t j = 0;j < 3; ++j)
                ifs >> node(j, i);
            }
          continue;
        }
      if(str == "CELLS"){
          ifs >> cell_num >> str;
          singularities.reserve(cell_num);
          int point_number_of_cell = 0;
          pair<size_t,size_t> edge_temp;
          for(size_t ci = 0; ci < cell_num; ++ci){
              ifs >> point_number_of_cell;
              if(point_number_of_cell != 2){
                  for(size_t i = 0; i < point_number_of_cell; ++i)
                    ifs >> str;
                }else{
                  ifs >> edge_temp.first >> edge_temp.second;
                }
              singularities.push_back(edge_temp);
            }
        }
    }
  return 0;

}

int save_size_file(const char * file, const zjucad::matrix::matrix<double> & node,
                   const vector<pair<size_t,size_t> > & singularities,
                   const double size)
{
  ofstream ofs(file);
  if(ofs.fail()) {
      cerr << "# [error] can not open size file." << endl;
      return __LINE__;
    }
  ofs << 0 << endl;
  ofs << singularities.size() << endl;
  for(size_t ei = 0; ei < singularities.size(); ++ei){
      const pair<size_t,size_t> & one_edge = singularities[ei];
      for(size_t di = 0; di < node.size(1); ++di)
        ofs << node(di, one_edge.first) << " ";
      for(size_t di = 0; di < node.size(1); ++di)
        ofs << node(di, one_edge.second) << " ";
      ofs << size << endl;
    }
  ofs << endl;
  return 0;
}

int dump_singularity_vtk_to_netgen_size(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] dump_singularity_vtk_to_netgen_size singularity_vtk size netgen_size_file"  << endl;
      return __LINE__;
    }

  vector<pair<size_t,size_t> > singularities;
  matrix<double> node;
  if(load_singularity_vtk(argv[1], singularities, node))
    return __LINE__;

  const double size = atof(argv[2]);

  if(save_size_file(argv[3], node, singularities, size))
    return __LINE__;
  return 0;
}
