#include "hex_process.h"
#include <jtflib/util/container_operation.h>

#include <iostream>
#include <numeric>
#include <zjucad/matrix/matrix.h>

using namespace std;
using namespace zjucad::matrix;
#define NUMERIC_ERROR 1e-6

int calc_point_normal(const size_t & tri_face_idx,
                      const matrixd & uvw,
                      const jtf::mesh::one_ring_face_at_point & orfap,
                      const matrixst & tri_faces,
                      const matrixd & tri_nodes,
                      const matrixd & tri_face_normal,
                      matrixd & normal)
{
  assert(tri_faces.size(1) == 3 && uvw.size() == 3);
  assert(tri_face_idx < tri_faces.size(2));
  assert(fabs(uvw[0] + uvw[1] + uvw[2] - 1) < NUMERIC_ERROR); // u+v+w = 1

  typedef jtf::mesh::one_ring_face_at_point::p2f_type::const_iterator p2fit;

  vector<size_t> zero_idx;
  for(size_t t = 0; t < uvw.size(); ++t){
    if(fabs(uvw[t]) < NUMERIC_ERROR)
      zero_idx << t;
  }

  if(zero_idx.empty()){ // point locates inside the triangle
    normal = tri_face_normal(colon(),tri_face_idx);
  }else if(zero_idx.size() == 1){// point locates on one edge
    vector<vector<size_t>  > adj_faces;
    for(size_t t = 0; t < tri_faces.size(1); ++t){
      if(t != zero_idx.front()){
        p2fit it = orfap.p2f_.find(tri_faces(t, tri_face_idx));
        if(it == orfap.p2f_.end()){
          cerr << "# [error] strange can not find point "
               << tri_faces(t, tri_face_idx) << " in one_ring_face_at_point."
               << endl;
          return __LINE__;
        }
        adj_faces.push_back(it->second);
      }
    }
    vector<size_t> shared_faces(adj_faces[0].size());
    vector<size_t>::iterator shared_end =
        find_intersection_set(adj_faces.front().begin(),
                              adj_faces.front().end(),
                              adj_faces.back().begin(),
                              adj_faces.back().end(),
                              shared_faces.begin());
    if(shared_end == shared_faces.begin()){
      cerr << "# [error] strange can not find common faces for points: "
           << tri_faces((zero_idx.front() + 1)%tri_faces.size(1), tri_face_idx)
           <<  ","
            << tri_faces((zero_idx.front() + 2)%tri_faces.size(1), tri_face_idx)
            << endl;
      return __LINE__;
    }

    if(find(shared_faces.begin(),shared_end, -1) != shared_end){
      cerr << "# [error] strane, found -1 as shared faces." << endl;
      return __LINE__;
    }

    normal = zeros<double>(3,1);
    for(vector<size_t>::iterator it = shared_faces.begin(); it != shared_end;
        ++it){
      const size_t &face_idx = *it;

      const double area =
          norm(cross(tri_nodes(colon(), tri_faces(1, face_idx))
                     - tri_nodes(colon(), tri_faces(0, face_idx)),
                     tri_nodes(colon(), tri_faces(2, face_idx))
                     - tri_nodes(colon(), tri_faces(1, face_idx))));
      normal += area * tri_face_normal(colon(),face_idx);
    }
    const double length = norm(normal);
    if(length > NUMERIC_ERROR)
      normal /= length;

  }else { // locate at points
    assert(zero_idx.size() == 2);
    const size_t point_idx =
        std::accumulate(tri_faces(colon(), tri_face_idx).begin(),
                        tri_faces(colon(), tri_face_idx).end(),  0)
        - tri_faces(zero_idx.front(), tri_face_idx)
        - tri_faces(zero_idx.back(), tri_face_idx);
    p2fit it = orfap.p2f_.find(point_idx);
    if(it == orfap.p2f_.end()){
      cerr << "# [error] can not find point " << point_idx
           << " in one_ring_face_at_point" << endl;
      return __LINE__;
    }

    normal = zeros<double>(3,1);
    for(size_t fi = 0; fi < it->second.size(); ++fi){
      if(it->second[fi] == -1) continue; // outside face
      const size_t &face_idx = it->second[fi];
      const double area =
          norm(cross(tri_nodes(colon(), tri_faces(1, face_idx))
                     - tri_nodes(colon(), tri_faces(0, face_idx)),
                     tri_nodes(colon(), tri_faces(2, face_idx))
                     - tri_nodes(colon(), tri_faces(1, face_idx))));
      normal += area * tri_face_normal(colon(), face_idx);
    }
    const double length = norm(normal);
    if(length > NUMERIC_ERROR)
    normal /= length;
  }
  return 0;
}
