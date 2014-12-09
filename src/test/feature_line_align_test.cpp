#include "feature_line_align_test.h"
#include "../hexmesh/io.h"
#include "../hexmesh/hexmesh.h"
#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"

#include <iostream>

using namespace std;

void feature_line_align_test::setUp()
{

}

void feature_line_align_test::tearDown()
{

}

void feature_line_align_test::feature_line_align()
{
  /* size_t triface_mat[] = { 5, 1, 4,
                                                        5, 4, 8,                                                                                      3, 7, 8,
                                                        3, 8, 4,
                                                        2, 6, 3,
                                                        6, 7, 3,
                                                        1, 5, 2,
                                                        5, 6, 2,
                                                        5, 8, 6,
                                                        8, 7, 6,
                                                        1, 2, 3,
                                                        1, 3, 4 };
    double tri_node_mat[] = { 1.000000, -1.000000, -1.000000,
                                                          1.000000, -1.000000, 1.000000,
							  -1.000000, -1.000000, 1.000000,
                                                          -1.000000, -1.000000, -1.000000,
                                                          1.000000, 1.000000, -0.999999,
                                                          0.999999, 1.000000, 1.000001,
                                                          -1.000000, 1.000000, 1.000000,
                                                          -1.000000, 1.000000, -1.000000};
    size_t quadface_mat[] = { 1, 2, 3, 4,
                                                    5, 8, 7, 6,
                                                    1, 5, 6, 2,
                                                    2, 6, 7, 3,
                                                    3, 7, 8, 4,
                                                    5, 1, 4, 8 };
    double quad_node_mat[] = { 1.000000, -1.000000, -1.000000,
                                                          1.000000, -1.000000, 1.000000,
							  -1.000000, -1.000000, 1.000000,
                                                          -1.000000, -1.000000, -1.000000,
                                                          1.000000, 1.000000, -0.999999,
                                                          0.999999, 1.000000, 1.000001,
                                                          -1.000000, 1.000000, 1.000000,
                                                          -1.000000, 1.000000, -1.000000}   zjucad::matrix::itr_matrix< size_t*> triface(3, 12, &triface_mat[0]);
   zjucad::matrix::itr_matrix< double*> tri_node(3, 8, &tri_node_mat[0]);
   zjucad::matrix::itr_matrix< size_t*> quadface(4, 6, &quadface_mat[0]);
   zjucad::matrix::itr_matrix< double*> quad_node(3, 8, &quad_node_mat[0]);
   triface = triface -1;
   quadface = quadface - 1;
   cout << triface << endl;
   cout << quadface << endl;
   cout << tri_node << endl;*/
   matrixst hex;
   matrixd quad_node;
   matrixst tet;
   matrixd tri_node;
   matrixst triface;
   matrixst quadface;
   matrixst tri_feature;
   matrixst quad_feature;
   boost::unordered_map<size_t, matrixd >  q2t;

   tet_mesh_read_from_zjumat("fandisk-301k.tet", &tri_node, &tet);
   cout << tet.size(1) <<endl;
   static face2tet_adjacent* tet_face = tet_face -> create(tet);
   get_outside_face(*tet_face, triface);
   cout << triface.size(2) <<endl;
   jtf::trimesh::extract_trimesh_feature_line(triface, tri_node, tri_feature);
   cout << tri_feature << endl;
   cout << "the num of trimesh feature line is:" << tri_feature.size(2) << endl;
   ofstream ofs("tri_feature.vtk");
   line2vtk(ofs, &tri_node[0],tri_node.size(2), &tri_feature[0], tri_feature.size(2));

   jtf::hexmesh::hex_mesh_read_from_wyz("fandisk2.hex", hex, quad_node,1);
   cout << hex.size(2) << endl;
   static jtf::hexmesh::face2hex_adjacent* hex_face = hex_face -> create(hex);
   jtf::hexmesh:: get_outside_face(*hex_face, quadface);
   cout << quadface.size(2) <<endl;

   jtf::quadmesh::extract_quadmesh_feature_line(quadface, quad_node, quad_feature);
   cout << quad_feature << endl;
   cout << "the num of quadmesh feature line is:" << quad_feature.size(2) << endl;
   ofstream ofs2("quad_feature.vtk");
   line2vtk(ofs2, &quad_node[0],quad_node.size(2), &quad_feature[0], quad_feature.size(2));
   project_quad_feature_to_tri_feature(quad_feature, quad_node,
                                       tri_feature, tri_node, q2t);
   boost::unordered_map<size_t, matrixd >::const_iterator
     map_it = q2t.begin();
   while(map_it != q2t.end())
     {
       cout << map_it -> first << endl;
       cout << map_it -> second << endl;
       ++ map_it;
     }
}
