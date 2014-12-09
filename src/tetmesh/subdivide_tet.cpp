#include "subdivide_tet.h"

#include <iostream>
#include <set>

#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>

#include "tetmesh.h"

using namespace std;
using namespace zjucad::matrix;

//! edge lenth = 2, center at zero
class regular_tetrahedron
{
public:
  typedef pair<matrix<double>, matrix<size_t> > sym_t;
  regular_tetrahedron()
    : nodes_(3, 4) {
    const double coor[12] = {
      -1, 0, -1/sqrt(2.0),
       1, 0, -1/sqrt(2.0),
       0, 1,  1/sqrt(2.0),
       0,-1,  1/sqrt(2.0)
    };
    copy(coor, coor+12, nodes_.begin());

    // get symmetric transfrom
    matrix<double> A = nodes_(colon(), colon(0, 2));
    inv(A);
    sym_rot_.resize(4*3);
    matrix<double> rot_nodes(3, 3), rot;
    for(int i = 0; i < 4; ++i) {
      for(int j = 0; j < 3; ++j) {
        matrix<size_t> &perm = sym_rot_[i*3+j].second;
        perm.resize(4, 1);
        perm[0] = i;
        perm[1] = (i+j+1)%4;
        matrix<double> &R = sym_rot_[i*3+j].first;
        int k;
        for(k = 0; k < 4; ++k) {
          if(k == perm[0] || k == perm[1]) continue;
          perm[2] = k;
          rot_nodes = nodes_(colon(), perm(colon(0, 2)));
          R = trans(rot_nodes*A);
          //          R = floor(R+0.5);
          matrix<double> tmp = R;
          if(det(tmp) < 0 || max(fabs(R*trans(R)-eye<double>(3))) > 1e-8) continue;
          break; // find the rot
        }
        if(k == 4) {
          cerr << nodes_ << perm << endl;
        }
        perm[3] = 6-sum(perm(colon(0, 2)));
      }
    }
#if 0 // validate
    for(size_t i = 0; i < 12; ++i) {
      if(norm(nodes_ - sym_rot_[i].first*nodes_(colon(), sym_rot_[i].second)) > 1e-8) {
        cerr << "error in regular tet: " << i << endl;
        cerr << sym_rot_[i].first << sym_rot_[i].second;
        exit(2);
      }
    }
#endif
  }
  const matrix<double> &nodes(void) const { return nodes_; }
  const vector<pair<matrix<double>, matrix<size_t> > > &sym_rot() const { return sym_rot_; }
private:
  matrix<double> nodes_;
  vector<pair<matrix<double>, matrix<size_t> > > sym_rot_;
};

// NOTICE: if a face has and only has two edge to be subdivded, it's
// not robust because it may lead to inconsistent tessellation. We
// should subdivide face first to ensure the consistency for
// neighboring tets, and then go here with hint points of new node and
// new face center.
class subdivide_regular_tetrahedron
{
public:
  typedef pair<matrix<double>, matrix<size_t> > pattern_t;
  subdivide_regular_tetrahedron() {
    pts_.resize(4, 4);
    for(int i = 0; i < 4; ++i) {
      for(int j = 0; j < 4; ++j) {
        pts_(i, j) = (tet_.nodes()(colon(), i) + tet_.nodes()(colon(), j))/2;
      }
    }

    { // one edge
      const size_t idx[] = {
        0, 0, 0, 1, 2, 2, 3, 3,
        0, 1, 1, 1, 2, 2, 3, 3
      };
      insert_table2pattern(idx, sizeof(idx)/sizeof(size_t));
    }
    { // two opposite edges
      const size_t idx[] = {
        0, 0, 0, 1, 2, 3, 3, 3,
        0, 0, 0, 1, 2, 2, 2, 3,
        0, 1, 1, 1, 2, 3, 3, 3,
        0, 1, 1, 1, 2, 2, 2, 3
      };
      insert_table2pattern(idx, sizeof(idx)/sizeof(size_t));
    }

    { // three edges in a face
      const size_t idx[] = {
        0, 0, 1, 1, 1, 2, 1, 3,
        0, 0, 1, 2, 2, 2, 2, 3,
        0, 0, 2, 3, 3, 3, 1, 3,
        0, 0, 1, 2, 2, 3, 1, 3
      };
      insert_table2pattern(idx, sizeof(idx)/sizeof(size_t));
    }

    { // six edges
      const size_t idx[] = {
        0, 0, 0, 1, 0, 2, 0, 3,
        1, 1, 1, 2, 0, 1, 1, 3,
        2, 2, 1, 2, 2, 3, 0, 2,
        3, 3, 1, 3, 0, 3, 2, 3,
        0, 1, 1, 2, 0, 2, 1, 3,
        0, 1, 1, 3, 0, 2, 0, 3,
        2, 3, 1, 3, 0, 3, 0, 2,
        2, 3, 1, 3, 0, 2, 1, 2        
      };
      insert_table2pattern(idx, sizeof(idx)/sizeof(size_t));
    }
#if 0 // validate
    const double regular_tet_vol = tet_volume(tet_.nodes());
    cerr << "regular_tet_vol: " << regular_tet_vol << endl;
    matrix<double> tet(3, 4);
    for(size_t i = 0; i < patterns_.size(); ++i) {
      double total_vol = 0;
      for(size_t j = 0; j < patterns_[i].second.size(2); j += 4) {
        for(int k = 0; k < 4; ++k) {
          tet(colon(), k) = pts()(patterns_[i].second(0, j+k), patterns_[i].second(1, j+k));
        }
        const double v = tet_volume(tet);
        if(v < 0) {
          cerr << "# negative tet vol in subdivide_regular_tetrahedron."
               << i << " " << j/4 << endl;
          exit(3);
        }
        //        cerr << "sub vol: " << j/4 << " " << v << endl;
        total_vol += v;
      }
      if(fabs(total_vol - regular_tet_vol) > 1e-8) {
          cerr << "# missing tet vol in subdivide_regular_tetrahedron: "
               << i << " " << total_vol << endl;
          exit(3);
      }
    }
    cerr << "# validate finished." << endl;
#endif
  }

  const regular_tetrahedron tet_;
  const matrix<matrix<double> > &pts() const { return pts_; }
  const vector<pattern_t> &patterns() const { return patterns_; }
private:
  void insert_table2pattern(const size_t *idx, size_t len) {
      vector<pair<size_t, size_t> > nodes_idx;
      for(size_t j = 0; j < len/2; ++j) { // for each node
        if(idx[j*2] == idx[j*2+1]) continue;
        nodes_idx.push_back(make_pair(idx[j*2], idx[j*2+1]));
      }
      sort(nodes_idx.begin(), nodes_idx.end());
      size_t num = unique(nodes_idx.begin(), nodes_idx.end())-nodes_idx.begin();
      matrix<double> nodes(3, num);
      for(size_t j = 0; j < num; ++j) {
        nodes(colon(), j) = pts_(nodes_idx[j].first, nodes_idx[j].second);
      }
      patterns_.push_back
        (make_pair(nodes,
                   itr_matrix<const size_t *>(2, len/2, idx)));
  }
  // the first is used to match configuration
  vector<pattern_t> patterns_;
  matrix<matrix<double> > pts_;
};

static const subdivide_regular_tetrahedron g_srt;

void add_edge_node4subdivide(const matrix<size_t> &edge_node, matrix<size_t> &new_node)
{
  const static size_t all_edges[] = {0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3};
  const size_t N = edge_node.size(2);
  int op_type = -1; // 0: no change, 1, already done, 6: all edges
  int visited[4] = {0, 0, 0, 0};
  if(N < 2 || N == 6)
    op_type = 0;
  else if(N > 3)
    op_type = 6;
  else { // N is in 2, 3
    for(int i = 0; i < edge_node.size(); ++i)
      ++visited[edge_node[i]];
    int missing_node;
    for(missing_node = 0; missing_node < 4; ++missing_node)
      if(visited[missing_node] == 0)
        break;
    if(N == 2) {
      if(missing_node == 4)
        op_type = 0;
      else {
        new_node.resize(2, 3);
        new_node(colon(), colon(0, 1)) = edge_node;
        for(int i = 0, j = 0; i < 4; ++i) {
          if(visited[i] == 1) new_node(j++, 2) = i;
        }
//        /assert(j == 2);
        op_type = 1;
      }
    }
    else { // N == 3
      op_type = (missing_node == 4)?6:0;
    }
  }
  switch(op_type) {
  case 0:
    new_node = edge_node;
    break;
  case 1:
    break;
  case 6:
    new_node.resize(2, 6);
    copy(all_edges, all_edges+12, new_node.begin());
    break;
  default:
    cerr << "# cannot be here in add_edge_node4subdivide." << endl;
    exit(5);
  }
}

//! @return: children is a 2xn matrix point to local node_idx
int subdivide_tet(const matrix<size_t> &edge_node,
                   matrix<size_t> &children)
{
  matrix<double> P(3, edge_node.size(2));
  for(size_t i = 0; i < P.size(2); ++i)
    P(colon(), i) = g_srt.pts()(edge_node(0, i), edge_node(1, i));
  size_t i;
  for(i = 0; i < g_srt.patterns().size(); ++i) {
    const subdivide_regular_tetrahedron::pattern_t &pat = g_srt.patterns()[i];
    if(P.size(2) != pat.first.size(2)) continue;
    size_t r;
    for(r = 0; r < g_srt.tet_.sym_rot().size(); ++r) { // for each possibility
      const regular_tetrahedron::sym_t & s = g_srt.tet_.sym_rot()[r];
      matrix<double> rot_P = s.first * P;
      int j;
      for(j = 0; j < rot_P.size(2); ++j) {
        int k;
        for(k = 0; k < rot_P.size(2); ++k) {
          if(max(fabs(rot_P(colon(), j) - pat.first(colon(), k))) < 1e-1) {
            rot_P(colon(), j) = 0; // find
            break;
          }
        }
        if(k == rot_P.size(2)) // cannot find
          break;
      }
      if(max(fabs(rot_P)) > 1e-1) // cannot find match configuration
        continue;
      children.resize(2, pat.second.size(2));
      children(colon()) = s.second(pat.second);
      break;
    }
    if(r != g_srt.tet_.sym_rot().size())
      break;
  }
  if(i == g_srt.patterns().size()) {
    cerr << "cannot find matched pattern." << edge_node << endl;
    exit(4);
  }
}

int subdivide_tet(const matrix<size_t> &tet, const matrix<size_t> &edge_node,
                   matrix<size_t> &children)
{
  matrix<size_t> local_edge_node = edge_node;
  for(size_t i = 0; i < local_edge_node.size(); ++i) {
    int j;
    for(j = 0; j < 4; ++j) {
      if(tet[j] == edge_node[i]) {
        local_edge_node[i] = j;
        break;
      }
    }
    assert(j != 4 && "not in this tet");
  }
  subdivide_tet(local_edge_node, children);
  children(colon()) = tet(children);
}

typedef vector<size_t> id_type;

//! @ param: eles is in num x dim, but cells is in dim x num, TODO: fix it
void create_ordered_table(const zjucad::matrix::matrix<size_t> &cells,
                     const zjucad::matrix::matrix<size_t> &pattern,
                     matrix<size_t> &eles)
{
  matrix<size_t> id(pattern.size(1));
  ordered_table_builder b(id.size());
  for(size_t ti = 0; ti < cells.size(2); ++ti) {
    for(int ei = 0; ei < pattern.size(2); ++ei) {
      id = cells(pattern(colon(), ei), ti);
      b << &id[0];
    }
  }
  b >> eles;
}

void create_edge2face(const matrix<size_t> &edge_set,
                      const matrix<size_t> &face_set,
                      vector<vector<size_t> > &edge2face)
{
  static const size_t tri_edge_idx[] = {0, 1, 0, 2, 1, 2};
  id_type edge_id(2);
  edge2face.resize(edge_set.size(1));
  for(size_t ti = 0; ti < face_set.size(1); ++ti) {
    for(int ei = 0; ei < 3; ++ei) {
      edge_id[0] = face_set(ti, tri_edge_idx[ei*2]);
      edge_id[1] = face_set(ti, tri_edge_idx[ei*2+1]);
      assert(edge_id[0] < edge_id[1]);
      const size_t edge_idx = find(edge_set, &edge_id[0]);
      assert(edge_idx != -1);
      edge2face[edge_idx].push_back(ti);
    }
  }
}

bool is_ordered_table(const zjucad::matrix::matrix<size_t> &ids)
{
  if(ids.size(1) < 2)
    return true;
  vector<size_t> rows[2];
  for(int i = 0; i < 2; ++i)
    rows[i].resize(ids.size(2));
  for(int i = 0; i < ids.size(2); ++i)
    rows[0][i] = ids(0, i);
  for(size_t r = 1; r < ids.size(1); ++r) {
    for(size_t c = 0; c < ids.size(2); ++c)
      rows[1][c] = ids(r, c);
    if(rows[0] > rows[1])
      return false;
    swap(rows[0], rows[1]);
  }
  return true;
}


void validate_edge_nodes(const zjucad::matrix::matrix<size_t> &tetmesh,
                         const zjucad::matrix::matrix<size_t> &edge_node,
                         zjucad::matrix::matrix<size_t> &new_nodes)
{
  assert(edge_node.size(2) == 2);
  static const size_t edge_idx[] = {0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3};
  itr_matrix<const size_t *> tet2edge(2, 6, edge_idx);
  matrix<size_t> edge_set;
  create_ordered_table(tetmesh, tet2edge, edge_set);

  static const size_t face_idx[] = {
    0, 1, 2,
    0, 1, 3,
    0, 2, 3,
    1, 2, 3
  };
  itr_matrix<const size_t *> tet2face(3, 4, face_idx);
  matrix<size_t> face_set;
  create_ordered_table(tetmesh, tet2face, face_set);

  vector<vector<size_t> > edge2face;
  create_edge2face(edge_set, face_set, edge2face);

  id_type edge_id(2);
  vector<set<size_t> > edge_node_in_face(face_set.size(1));
  for(size_t ei = 0; ei < edge_node.size(1); ++ei) {
    edge_id[0] = edge_node(ei, 0);
    edge_id[1] = edge_node(ei, 1);
    make_id(&edge_id[0], 2);
    const size_t edge_idx = find(edge_set, &edge_id[0]);
    assert(edge_idx != -1);
    const vector<size_t> &face_idx = edge2face[edge_idx];
    for(size_t fi = 0; fi < face_idx.size(); ++fi)
      edge_node_in_face[face_idx[fi]].insert(edge_idx);
  }

  set<size_t> invalid_face;
  for(size_t fi = 0; fi < edge_node_in_face.size(); ++fi) {
    if(edge_node_in_face[fi].size() == 2)
      invalid_face.insert(fi);
  }
  // it must converge because in the worest case, all edges are
  // subdivided.
  while(!invalid_face.empty()) { 
    const size_t x = *invalid_face.begin();
    set<size_t> &edges = edge_node_in_face[x];
    // find the duplicated node
    size_t node_idx_sum = 0;
    for(set<size_t>::const_iterator i = edges.begin(); i != edges.end(); ++i) {
      node_idx_sum += (edge_set(*i, 0)+edge_set(*i, 1));
    }
    node_idx_sum -= sum(face_set(x, colon()));
    for(int i = 0, j = 0; i < 3; ++i) {
      if(face_set(x, i) == node_idx_sum) continue;
      edge_id[j++] = face_set(x, i);
    }
    //assert(j == 2);
    make_id(&edge_id[0], 2);
    const size_t edge_idx = find(edge_set, &edge_id[0]); // the edge to be added
    const vector<size_t> &nb_faces = edge2face[edge_idx];
    for(size_t fi = 0; fi < nb_faces.size(); ++fi) {
      set<size_t> &edges = edge_node_in_face[nb_faces[fi]];
      edges.insert(edge_idx);
      if(edges.size() == 2)
        invalid_face.insert(nb_faces[fi]);
      else
        invalid_face.erase(nb_faces[fi]);
    }
  }
  cerr << "validate over." << endl;
  set<size_t> all_edges;
  for(size_t fi = 0; fi < edge_node_in_face.size(); ++fi) {
    all_edges.insert(edge_node_in_face[fi].begin(), edge_node_in_face[fi].end());
  }
  new_nodes.resize(all_edges.size(), 2);
  int k = 0;
  for(set<size_t>::const_iterator i = all_edges.begin();
      i != all_edges.end(); ++i, ++k)
    new_nodes(k, colon()) = edge_set(*i, colon());

  assert(is_ordered_table(new_nodes));
}

void create_edge2tet(const zjucad::matrix::matrix<size_t> &tets,
                     const matrix<size_t> &edge_set,
                     vector<vector<size_t> > &edge2tet)
{
  static const size_t local_edge_idx[] = {0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3};
  id_type edge_id(2);
  edge2tet.resize(edge_set.size(1));
  for(size_t ti = 0; ti < tets.size(2); ++ti) {
    for(int ei = 0; ei < 6; ++ei) {
      edge_id[0] = tets(local_edge_idx[ei*2], ti);
      edge_id[1] = tets(local_edge_idx[ei*2+1], ti);
      make_id(&edge_id[0], 2);
      const size_t edge_idx = find(edge_set, &edge_id[0]);
      assert(edge_idx != -1);
      edge2tet[edge_idx].push_back(ti);
    }
  }
}

int subdivide_tetmesh(const zjucad::matrix::matrix<size_t> &tetmesh,
                      size_t parent_node_num,
                      const zjucad::matrix::matrix<size_t> &edge_node,
                      zjucad::matrix::matrix<size_t> &children)
{
  assert(edge_node.size(2) == 2 && is_ordered_table(edge_node));
  static const size_t edge_idx[] = {0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3};
  itr_matrix<const size_t *> tet2edge(2, 6, edge_idx);
  matrix<size_t> edge_set;
  create_ordered_table(tetmesh, tet2edge, edge_set);
  vector<vector<size_t> > edge2tet;
  create_edge2tet(tetmesh, edge_set, edge2tet);
  vector<vector<size_t> > edge_in_tet(tetmesh.size(2));
  id_type edge_id(2);
  for(size_t ei = 0; ei < edge_node.size(1); ++ei) {
    edge_id[0] = edge_node(ei, 0);
    edge_id[1] = edge_node(ei, 1);
    const size_t edge_idx = find(edge_set, &edge_id[0]);
    assert(edge_idx != -1);
    const vector<size_t> &nb_tets = edge2tet[edge_idx];
    for(size_t ti = 0; ti < nb_tets.size(); ++ti) {
      edge_in_tet[nb_tets[ti]].push_back(edge_idx);
    }
  }
  vector<pair<matrix<size_t>, size_t> > new_tet_vec;
  matrix<size_t> tet_edge_node;
  for(size_t ti = 0; ti < tetmesh.size(2); ++ti) {
    const vector<size_t> &edges = edge_in_tet[ti];
    if(edges.empty()) {
      new_tet_vec.push_back(make_pair(tetmesh(colon(), ti), ti));
      continue;
    }
    tet_edge_node.resize(2, edges.size());
    for(size_t ei = 0; ei < edges.size(); ++ei) {
      tet_edge_node(colon(), ei) = trans(edge_set(edges[ei], colon()));
    }
    matrix<size_t> tet = tetmesh(colon(), ti);
    subdivide_tet(tet, tet_edge_node, children);
    matrix<size_t> t(4);
    for(size_t i = 0; i < children.size(2)/4; ++i) {
      for(int j = 0; j < 4; ++j) {
        if(children(0, i*4+j) == children(1, i*4+j)) { // parent node
          t[j] = children(0, i*4+j);
        }
        else {
          make_id(&children(0, i*4+j), 2);
          t[j] = find(edge_node, &children(0, i*4+j));
          assert(t[j] != -1);
          t[j] += parent_node_num;
        }
      }
      new_tet_vec.push_back(make_pair(t, ti));
    }
  }
  children.resize(4, new_tet_vec.size());
  for(size_t i = 0; i < new_tet_vec.size(); ++i) {
    children(colon(), i) = new_tet_vec[i].first;
  }
}

ordered_table_builder::ordered_table_builder(size_t dim)
  :dim_(dim) {
}

ordered_table_builder &ordered_table_builder::operator << (const size_t *id)
{
  vector<size_t> v(dim_);
  copy(id, id+dim_, v.begin());
  sort(v.begin(), v.end());
  buf_.push_back(v);
  return *this;
}

void ordered_table_builder::operator >> (zjucad::matrix::matrix<size_t> &tab)
{
  sort(buf_.begin(), buf_.end());
  size_t num = unique(buf_.begin(), buf_.end())-buf_.begin();
  buf_.resize(num);
  tab.resize(num, dim_);
  for(size_t i = 0; i < num; ++i)
    for(size_t j = 0; j < dim_; ++j)
      tab(i, j) = buf_[i][j];
}

void subdivide_top2geo(const zjucad::matrix::matrix<double> &parent_node,
                       const zjucad::matrix::matrix<size_t> &new_node_parent_id, //nx2 to ori node
                       zjucad::matrix::matrix<double> &child_node
                       )
{
  child_node.resize(3, parent_node.size(2)+new_node_parent_id.size(1));
  child_node(colon(), colon(0, parent_node.size(2)-1)) = parent_node;
  for(size_t i = 0; i < new_node_parent_id.size(1); ++i) {
    child_node(colon(), parent_node.size(2)+i)
      = parent_node(colon(), new_node_parent_id(i, colon()))*ones<double>(2, 1)/2.0;
  }
}
