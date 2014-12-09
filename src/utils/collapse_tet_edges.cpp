#include "../tetmesh/tetmesh.h"
#include "../tet_mesh_sxx/tet_mesh_sxx.h"
#include "../common/util.h"
#include <set>

using namespace std;
using namespace zjucad::matrix;

int load_edge_file(const char * file, set<pair<size_t,size_t> >  & edges)
{
  ifstream ifs(file);
  if(ifs.fail()) {
      cerr << "# [error] can not open file." << endl;
      return __LINE__;
    }
  size_t edge_num;
  ifs >> edge_num;
  pair<size_t,size_t> one_edge;
  for(size_t ei = 0; ei < edge_num; ++ei){
      ifs >> one_edge.first >> one_edge.second;
      if(one_edge.first > one_edge.second) swap(one_edge.second, one_edge.first);
      edges.insert(one_edge);
    }
  return 0;
}

int collapse_tet_edge_with_zyz(jtf::tet_mesh & tm, const set<pair<size_t,size_t> > & edges,
                               zjucad::matrix::matrix<double> &zyz)
{
  sxx::tet_mesh stm;
  stm.create_tetmesh(tm.tetmesh_.node_, tm.tetmesh_.mesh_);
  for(const auto & one_edge : edges){
      stm.collapse_edge(one_edge);
    }
  stm.write_tetmesh_to_matrix(tm.tetmesh_.node_, tm.tetmesh_.mesh_,false);
  if(zyz.size() != 0){
      matrix<size_t> tet_map2original;
      stm.get_tet2orginal_index(tm.tetmesh_.mesh_, tet_map2original);

      matrix<double>  old_zyz = zyz;
      zyz.resize(3, tm.tetmesh_.mesh_.size(2));
      for(size_t ti = 0; ti < tet_map2original.size(); ++ti){
          zyz(colon(),ti) = old_zyz(colon(), tet_map2original[ti]);
        }
    }
  remove_extra_node(tm.tetmesh_.mesh_, tm.tetmesh_.node_);
  return 0;
}

int collapse_tet_edges(int argc, char * argv[])
{
  if(argc != 4 && argc != 6){
      cerr << "# [usage] collapse_tet_edges tet edge_file output_tet [input_zyz] [output_zyz]" << endl;
      return __LINE__;
    }

  jtf::tet_mesh tm(argv[1]);
  set<pair<size_t,size_t> > collapse_edges;
  if(load_edge_file(argv[2], collapse_edges))
    return __LINE__;

  matrix<double>  zyz;
  if(argc == 6){
      if(jtf::mesh::read_matrix(argv[4], zyz))
        return __LINE__;
      if(zyz.size(2) != tm.tetmesh_.mesh_.size(2)){
          cerr << "# [error] zyz is not compatible with tet." << endl;
          return __LINE__;
        }
    }

  collapse_tet_edge_with_zyz(tm, collapse_edges, zyz);

  jtf::mesh::tet_mesh_write_to_zjumat(argv[3], &tm.tetmesh_.node_, &tm.tetmesh_.mesh_);
  if(argc == 6){
      jtf::mesh::write_matrix(argv[5], zyz);
    }
  cerr << "# [info] success." << endl;

  return 0;
}

