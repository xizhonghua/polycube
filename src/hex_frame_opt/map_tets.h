#ifndef MAP_TETS_H
#define MAP_TETS_H

#include "../tetmesh/tetmesh.h"
void map_tets(jtf::tet_mesh &tm0, const jtf::tet_mesh &tm1,
              boost::property_tree::ptree &pt);
void map_tris(const zjucad::matrix::matrix<size_t> &tri,
             const zjucad::matrix::matrix<double> &node0,
              const zjucad::matrix::matrix<double> &node1,
              boost::property_tree::ptree &pt);
#endif // MAP_TETS_H
