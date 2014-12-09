#include "../tetmesh/tetmesh.h"
#include "../tet_mesh_sxx/tet_mesh_sxx.h"

int load_edge_file(const char * file,
                   vector<pair<size_t,size_t> > & split_edges)
{
  ifstream ifs(file);
  if(ifs.fail()){
      cerr << "# [error] can not open edge file." << endl;
      return __LINE__;
    }
  size_t edge_num = 0;
  ifs >> edge_num;
  split_edges.clear();
  pair<size_t,size_t> one_edge;
  for(size_t ei = 0; ei < edge_num; ++ei){
      ifs >> one_edge.first >> one_edge.second;
      split_edges.push_back(one_edge);
    }

  return 0;
}

int split_tets_on_edges(int argc, char * argv[])
{
  using  namespace zjucad::matrix;
  if(argc != 4 && argc != 6){
      cerr << "# [usage] split_tets_on_edges input_tet input_edges_file output_tet [input_zyz] [output_zyz]" << endl;
      return __LINE__;
    }

  vector<pair<size_t,size_t> > split_edges;
  if(load_edge_file(argv[2], split_edges)){
      return __LINE__;
    }

  sxx::tet_mesh stm;
  stm.create_tetmesh(argv[1]) ;

  for(const auto & one_edge : split_edges){
      stm.split_edge(one_edge);
    }

  stm.write_tetmesh_to_file(argv[3]);

  if(argc == 6){
      zjucad::matrix::matrix<size_t> new_tet;
      zjucad::matrix::matrix<double> new_node;
      stm.write_tetmesh_to_matrix(new_node, new_tet);
      zjucad::matrix::matrix<double> zyz;
      jtf::mesh::read_matrix(argv[4], zyz);
      zjucad::matrix::matrix<size_t> tet_map;
      stm.get_tet2orginal_index(new_tet, tet_map);
      zjucad::matrix::matrix<double> old_zyz = zyz;
      zyz.resize(3, new_tet.size(2));
      for(size_t ti = 0; ti < new_tet.size(2); ++ti){
          zyz(colon(),ti) = old_zyz(colon(), tet_map[ti]);
        }
      jtf::mesh::write_matrix(argv[5], zyz);
    }
  return 0;
}

