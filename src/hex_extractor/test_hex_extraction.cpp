#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>


#include <zjucad/ptree/ptree.h>
#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/io.h>
#include <ctime>
#include <boost/filesystem.hpp>

#include "../common/vtk.h"
#include "extractor.h"
#include <omp.h>
#include "../hexmesh/io.h"

using namespace std;
using namespace zjucad::matrix;
using boost::property_tree::ptree;
using namespace jtf::mesh;
using namespace extractor;
using namespace boost;

#ifdef USE_PTREE
int main(int argc, char *argv[])
{
  
  ptree pt;

  try {

    zjucad::read_cmdline(argc, argv, pt);
    string input_file = pt.get<string>("input_file.value");

    matrixd node;
    matrixst tet;
    Imatrix int_point;

    tet_mesh_read_from_txt(input_file.c_str(), &node, &tet, int_point);

    generate_all_integer_grid(&node, &tet, int_point);

    print("../text.txt", int_point);


  }catch(std::exception &e){
    cerr << endl;
    cerr << "[error] "<<e.what() << endl;
    cerr << "Usage:" << endl;
    zjucad::show_usage_info(std::cerr, pt);
  }

  return 0;
}
#else

void point2vtk(ofstream &out, vector<boost::tuple<double, double, double> > &param_nodes, vector<point_info> &param_mesh)
{
  out << "# vtk DataFile Version 2.0\nTRI\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";
  out << "POINTS " << param_nodes.size() << " float\n";
  for (size_t i = 0; i < param_nodes.size(); i++)
    out << param_nodes[i].get<0>() << " " << param_nodes[i].get<1>() << " " << param_nodes[i].get<2>() << endl;

  out << "CELLS " << param_mesh.size() << " " << param_mesh.size()*2 << "\n";
  for(size_t i = 0; i < param_mesh.size(); ++i)
    out << 1 << " " << param_mesh[i].coordinate_idx << "\n";

  out << "CELL_TYPES " << param_mesh.size() << "\n";
  for(size_t i = 0; i < param_mesh.size(); ++i)
    out << 1 << "\n";
}

int read_inner_type(const char* path, std::map<pair<size_t, size_t>, size_t> &inner_type)
{
  ifstream os(path);
  if ( !os ) {
      cerr << "[INFO]no inner trasition\n";
      return 1;
    }
  size_t s, t, r;
  while ( os >> s >> t >> r )
    inner_type.insert(std::make_pair(std::make_pair(s, t), r));
  return 0;
}

void see_final_nodes(ofstream &out, const zjucad::matrix::matrix<double> &nodes)
{
  out << "# vtk DataFile Version 2.0\nTRI\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";
  out << "POINTS " << nodes.size(2) << " float\n";
  for (size_t i = 0; i < nodes.size(2); i++)
    out << nodes(0, i) << " " << nodes(1, i) << " " << nodes(2, i) << endl;

  out << "CELLS " << nodes.size(2) << " " << nodes.size(2)*2 << "\n";
  for(size_t i = 0; i < nodes.size(2); ++i)
    out << 1 << " " << i << "\n";

  out << "CELL_TYPES " << nodes.size() << "\n";
  for(size_t i = 0; i < nodes.size(); ++i)
    out << 1 << "\n";
}



int main(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [usage] test_hex_extraction input_path scale " << endl;
      return __LINE__;
    }
  zjucad::matrix::matrix<double> orig_nodes;
  zjucad::matrix::matrix<size_t> orig_mesh;
  zjucad::matrix::matrix<double> cut_nodes;
  zjucad::matrix::matrix<size_t> cut_mesh;
  std::vector<boost::tuple<double, double, double> > i_nodes;
  std::vector<point_info> i_mesh;
  std::vector<point_info> c_mesh;
  std::vector<boost::tuple<double, double, double> > c_nodes;
  std::map<pair<size_t, size_t>, size_t> inner_type;
  zjucad::matrix::matrix<double> hex_nodes;
  zjucad::matrix::matrix<size_t> hex_mesh;
  zjucad::matrix::matrix<double> final_nodes;
  zjucad::matrix::matrix<size_t> final_mesh;
  zjucad::matrix::matrix<double> inv_base;
  zjucad::matrix::matrix<bool>   invertible;

  hex_extractor handle;

  string path = argv[1];
  double scale = atof(argv[2]);

  string path0 = path + "/orig.tet.vtk";
  string path1 = path + "/param.tet.vtk";
  string path2 = path + "/orig.hex.vtk";
  string path3 = path + "/param.hex.vtk";

  std::ofstream OUT0(path2);
  std::ofstream OUT1(path3);

  string inner_type_path = path + "/inner_type";

  tet_mesh_read_from_vtk(path0.c_str(), &orig_nodes, &orig_mesh);
  tet_mesh_read_from_vtk(path1.c_str(), &cut_nodes, &cut_mesh);

 read_inner_type(inner_type_path.c_str(), inner_type);

  clock_t start = clock();

  handle.SCALE = scale;

  handle.mesh_scaling(cut_nodes, scale);

  handle.calculate_tet_base_inverse(cut_mesh, cut_nodes, inv_base, invertible);

  //extract integer points
  handle.extract_points(cut_mesh, cut_nodes, 0, inv_base, invertible, i_mesh, i_nodes);

  handle.merge_points_with_same_orig_coord(cut_mesh, i_nodes, i_mesh);

  //extract half integer points
  handle.extract_points(cut_mesh, cut_nodes, 0.5, inv_base, invertible, c_mesh, c_nodes);


  handle.extract_mesh(orig_mesh, orig_nodes, cut_mesh, cut_nodes, inv_base, invertible, i_mesh, i_nodes,
                      c_mesh, c_nodes, inner_type, hex_mesh, hex_nodes);

  handle.mapping_to_orig_space(orig_mesh, orig_nodes, i_mesh, hex_mesh, hex_nodes, final_mesh, final_nodes);

  clock_t end = clock();

  cout << "[INFO]THE TIME COST : " << (double) ( end - start ) / CLOCKS_PER_SEC << endl;

  //    zjucad::matrix::matrix<double> orig_c_nodes;
  //    orig_c_nodes.resize(3, c_nodes.size());
  //    for (size_t i = 0; i < c_mesh.size(); ++i)
  //    {
  //        zjucad::matrix::matrix<double> tet_nodes;
  //        tet_nodes = orig_nodes(colon(), orig_mesh(colon(), c_mesh[i].tet_idx));
  //        orig_c_nodes(colon(), c_mesh[i].coordinate_idx) = tet_nodes * c_mesh[i].b_coeff;
  //    }
  //    see_final_nodes(out6, final_nodes);
  //    see_final_nodes(out10, orig_c_nodes);

  hex2vtk(OUT1, &hex_nodes[0], hex_nodes.size(2), &hex_mesh[0], hex_mesh.size(2));
  hex2vtk(OUT0, &final_nodes[0], final_nodes.size(2), &final_mesh[0], final_mesh.size(2));

  cout << "ALL DONE!" << endl;

  return 0;

}

#endif
