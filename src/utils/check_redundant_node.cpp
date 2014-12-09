#include <iostream>
#include <map>
#include <jtflib/mesh/io.h>
#include <zjucad/matrix/matrix.h>
#include "../hex_param/io.h"

using namespace std;
using namespace zjucad::matrix;

int check_redundant_node(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] check_redundant_node cut_tet uncut_tet node_group " << endl;
      return __LINE__;
    }

  matrix<size_t> cut_tet;
  matrix<double> cut_node;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &cut_node, &cut_tet))
    return __LINE__;

  matrix<size_t> uncut_tet;
  matrix<double> uncut_node;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &uncut_node, &uncut_tet))
    return __LINE__;

  matrix<size_t> cut_tet2tet(max(cut_tet) + 1);
  cut_tet2tet(cut_tet) = uncut_tet(colon());


  vector<vector<size_t> > node_groups;
  if(load_group_file(argv[3], node_groups))
    return __LINE__;

  map<size_t,size_t> variant2group;
  for(size_t gi = 0; gi < node_groups.size(); ++gi){
      const vector<size_t> & one_group = node_groups[gi];
      for(size_t ni = 0; ni < one_group.size(); ++ni){
          variant2group[one_group[ni]] = gi;
        }
    }

  matrix<size_t>  cut_node_group_varint(3, cut_node.size(2));
  for(size_t i = 0; i < cut_node_group_varint.size(); ++i){
      cut_node_group_varint[i] = variant2group[i];
    }

  map<vector<size_t>, vector<size_t> > coord_points;
  vector<size_t> one_node(3);
  for(size_t i = 0; i < cut_node_group_varint.size(2); ++i){
      std::copy(cut_node_group_varint(colon(),i).begin(),
                cut_node_group_varint(colon(),i).end(), one_node.begin());
      coord_points[one_node].push_back(i);
    }

  size_t point_number = 0;
  for(const auto & one_coord: coord_points){
      if(one_coord.second.size() < 2) continue;
      cerr << "points with same coord: ";
      for(size_t i = 0; i < one_coord.second.size(); ++i)
        cerr << one_coord.second[i] << " " ;
      cerr << norm(cut_node(colon(),one_coord.second.front()) - cut_node(colon(),one_coord.second.back())) << endl;
      cerr << endl;
      ++point_number;
    }
  cerr << "same points " << point_number << endl;
  return 0;
}

//int remove_redundant_node(int argc, char * argv[])
//{
//  if(argc != 5){
//      cerr << "# [usage] check_redundant_node cut_tet uncut_tet node_group node_equation " << endl;
//      return __LINE__;
//    }

//  matrix<size_t> cut_tet;
//  matrix<double> cut_node;
//  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &cut_node, &cut_tet))
//    return __LINE__;

//  matrix<size_t> uncut_tet;
//  matrix<double> uncut_node;
//  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &uncut_node, &uncut_tet))
//    return __LINE__;

//  matrix<size_t> cut_tet2tet(max(cut_tet) + 1);
//  cut_tet2tet(cut_tet) = uncut_tet(colon());


//  vector<vector<size_t> > node_groups;
//  if(load_group_file(argv[3], node_groups))
//    return __LINE__;

//  map<size_t,size_t> variant2group;
//  for(size_t gi = 0; gi < node_groups.size(); ++gi){
//      const vector<size_t> & one_group = node_groups[gi];
//      for(size_t ni = 0; ni < one_group.size(); ++ni){
//          variant2group[one_group[ni]] = gi;
//        }
//    }

//  matrix<size_t>  cut_node_group_varint(3, cut_node.size(2));
//  for(size_t i = 0; i < cut_node_group_varint.size(); ++i){
//      cut_node_group_varint[i] = variant2group[i];
//    }

//  map<vector<size_t>, vector<size_t> > coord_points;
//  vector<size_t> one_node(3);
//  for(size_t i = 0; i < cut_node_group_varint.size(2); ++i){
//      std::copy(cut_node_group_varint(colon(),i).begin(),
//                cut_node_group_varint(colon(),i).end(), one_node.begin());
//      coord_points[one_node].push_back(i);
//    }

//  //////////////////////////////////////////////////////////////////////////////
//  /// find points with the same coordinates, update them

//  matrix<size_t> update_cut_tet = cut_tet;
//  for(const auto & one_coord : coord_points){
//      if(one_coord.second.size() < 2) continue;
//      const vector<size_t> & same_points = one_coord.second;

//    }
//  return 0;
//}
