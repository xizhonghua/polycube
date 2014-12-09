/*=============================================================================*\
 *                                                                             *
 *                               Volume Frame                                  *
 *         Copyright (C) 2011-2021 by State Key Lab of CAD&CG 410              *
 *			Zhejiang University, shizeyun                          *
 *                        http://www.cad.zju.edu.cn/                           *
 *                                                                             *
 *-----------------------------------------------------------------------------*
 *  This file is part of Volume Frame                                          *
 *  Created by shizeyun                                                        *
 *                                                                             *
\*=============================================================================*/

//== INCLUDES ===================================================================

//STANDARD
#include<algorithm>
#include<cassert>
#include<iostream>
#include<fstream>
#include<vector>
#include<zjucad/matrix/io.h>

//LOCAL
#include"hex_process.h"
#include"../common/vtk.h"

//== NAMESPACES =================================================================

using zjucad::matrix::colon;
using zjucad::matrix::cross;
using zjucad::matrix::dot;
using zjucad::matrix::matrix;
using zjucad::matrix::norm;
using zjucad::matrix::zeros;

//===============================================================================
int create_tet_from_hex_cut_raw(
    const matrixst & hex,
    const matrixd & node,
    matrixst & tet)
{
  tet.resize(4, 6 * hex.size(2));

  const size_t tet_indice [24] = {
    7,5,0,1,
    7,5,4,0,
    7,4,6,0,
    7,1,0,3,
    7,0,6,2,
    7,0,2,3
  };
  zjucad::matrix::itr_matrix<const size_t*> tet_indice_mat(4, 6, &tet_indice[0]);
  for(size_t hi = 0; hi < hex.size(2); ++hi){
    tet(colon(), colon(6 * hi, 6 * hi + 5))(colon()) = hex(tet_indice_mat, hi);
  }

  orient_tet(node, tet);
  
#if 1 // test
  {
    matrixst tet_faces;
    jtf::mesh::get_outside_face(*jtf::mesh::face2tet_adjacent::create(tet), tet_faces);

    //check valid by write to vtk
    std::ofstream outfile("tet.vtk");
    tet2vtk(outfile, &node[0], node.size(2), &tet[0], tet.size(2));
  }
#endif
  return 0;
}

//int create_tet_from_hex_cut(const hexmesh &my_hex, jtf::mesh::meshes &my_tet) {
//  //tet has same nodes with hex
//  my_tet.node = my_hex.node;
//  const size_t hex_num = my_hex.mesh_.size(2);
//  const size_t tet_num = 6*hex_num;
//  my_tet.mesh_.resize(4, tet_num);

//  std::cout << "# [info] node number of tet: " << my_tet.node_.size(2) << std::endl;
//  std::cout << "# [info] tet number of tet: " << my_tet.mesh_.size(2) << std::endl;
//  std::cout << "# [info] node number of hex: " << my_hex.node_.size(2) << std::endl;
//  std::cout << "# [info] hex number of hex: " << my_hex.mesh_.size(2) << std::endl;

//  //construct tet
//  const size_t tet_indice[6][4] = {
//    {7,5,0,1},
//    {7,5,4,0},
//    {7,4,6,0},
//    {7,1,0,3},
//    {7,0,6,2},
//    {7,0,2,3}
//  };
//  for (size_t hid = 0; hid < hex_num; ++hid) {
//    for (size_t i = 0; i < 6; ++i) {
//      for (size_t j = 0; j < 4; ++j) {
//        my_tet.mesh_(j, hid*6+i) = my_hex.hex(tet_indice[i][j], hid);
//      }
//    }
//  }

//  orient_tet(my_tet.node_, my_tet.mesh_);
//  matrixst tet_faces;
//  get_outside_face(*face2tet_adjacent::create(my_tet.mesh_), tet_faces);
  
//  //check valid by write to vtk
//  std::ofstream outfile("tet.vtk");
//  tet2vtk(outfile,
//          &my_tet.node[0], my_tet.node_.size(2),
//          &my_tet.mesh[0], tet_num);

//  return 0;
//}

//===============================================================================

//===============================================================================
