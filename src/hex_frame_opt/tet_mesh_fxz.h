#ifndef TET_MESH_FXZ_H
#define TET_MESH_FXZ_H

#include <jtflib/mesh/mesh.h>

namespace fxz {

  const size_t INVALID_NUM = (size_t)(-1); 
  const size_t UNSIGN_NEG_ONE = (size_t)(-1);

  typedef zjucad::matrix::matrix<double> matrixd;
  typedef zjucad::matrix::matrix<size_t> matrixst;

  class tet_mesh
  {
 public:
    enum status{ SURF_TET=1, SURF_VERT=(1<<1), SURF_FACE=(1<<2),
                 SURF_EDGE=(1<<3),
                 INVALID_ID=INVALID_NUM};


 tet_mesh(matrixd& vs, matrixst& ts): verts_(vs), tets_verts_(ts) {}

    ~tet_mesh() { }

    int build();
    
    jtf::mesh::face2tet_adjacent& face2tet() { return *fa_;}
    
    inline const matrixd& verts_coord() const { return verts_; }
    inline matrixd& verts_coord() { return verts_; }
    inline const matrixst& tets_verts() const { return tets_verts_; }
    inline matrixst& tets_verts() { return tets_verts_; }
    inline const jtf::mesh::face2tet_adjacent& face2tet() const { return *fa_;}

    inline size_t tets_num() const { return tets_verts_.size(2); }
    inline size_t verts_num() const { return verts_coord().size(2); }
    inline size_t direc_faces_num() const { return faces_verts().size()*2; }
    inline size_t undir_faces_num() const { return faces_verts().size(); }
    inline size_t direc_edges_num() const { return edges_verts_.size()*2; }
    inline size_t undir_edges_num() const { return edges_verts_.size(); }

    inline size_t pair_edge(size_t index) const { return (index^1); }
    inline size_t pair_face(size_t index) const { return (index^1); }

    inline std::vector<double> vert_coord(size_t index) const
    {
      assert(index < verts_num());
      std::vector<double> r(3);
      r[0] = verts_(0, index);
      r[1] = verts_(1, index);
      r[2] = verts_(2, index);
      return r;
    }

    // get the vert index of element: face, tet, edge
    //! TODO  the vertices should be right hand sort, but not be done.
    inline std::vector<size_t> face_verts(size_t index) const
    {
      std::vector<size_t> r(3);
      assert((index>>1) < faces_verts().size());
      r[0] = faces_verts()[index>>1][0];
      r[1] = faces_verts()[index>>1][1];
      r[2] = faces_verts()[index>>1][2];
      return r;
    }
    
    inline std::vector<size_t> tet_verts(size_t index) const
    {
      assert(index < tets_num());
      std::vector<size_t> r(4);
      r[0] = tets_verts()(0, index);
      r[1] = tets_verts()(1, index);
      r[2] = tets_verts()(2, index);
      r[3] = tets_verts()(3, index);
      return r;
    }

    inline std::pair<size_t,size_t> edge_verts(size_t index) const
    {
      assert((index>>1)<edges_verts_.size());
      std::pair<size_t, size_t> ev = edges_verts_[index>>1];
      assert(ev.first < ev.second);
      if (index%2==1) std::swap(ev.first, ev.second);
      return ev;
    }

    // get index : face, edge
    inline size_t face_index(size_t a, size_t b) const // for two tet indices
    {
      assert(a < t2t_.size());
      for (size_t i=0; i<t2t_[a].size(); ++i)
        if (t2t_[a][i] == b) return t2f_[a][i];
      return INVALID_ID;
    }
    
    inline size_t edge_index(size_t a, size_t b) const // for two vert indices
    {
      assert(a < v2v_.size());
      for (size_t i=0; i<v2v_[a].size(); ++i)
        if (v2v_[a][i] == b) return v2e_[a][i];
      return INVALID_ID;
    }
    
    const std::vector<size_t>& tet_adj_tets(size_t index) const
    {
      assert(index < t2t_.size());
      return t2t_[index];
    }

    const std::vector<size_t>& tet_adj_faces(size_t index) const
    {
      assert(index < t2f_.size());
      return t2f_[index];
    }
    
    std::vector<size_t> edge_adj_loop_tets(size_t index) const
    {
      assert(index < direc_edges_num());
      const std::vector<size_t>& loop = e2t_[index>>1];
      std::vector<size_t> r(loop.size());
      if (index%2 == 1) {
        // assert that the loop is right hand when (a<b) (a,b) is the edge verts
        std::copy(loop.rbegin(), loop.rend(), r.begin());
        return r;
      }
      return loop;
    }

    std::vector<size_t> edge_adj_loop_faces(size_t index) const
    {
      assert(index < direc_edges_num());
      if (is_surf_edge(index)) return surf_edge_adj_loop_faces(index);
      const std::vector<size_t>& loop = e2t_[index>>1];
      std::vector<size_t> fs;
      if (index%2 == 1) {
        for (size_t i=loop.size()-1; i>0; --i)
          fs.push_back(face_index(loop[i], loop[i-1])); // assert first==last in loop
      } else {
        for (size_t i=0; i<loop.size()-1; ++i)
          fs.push_back(face_index(loop[i], loop[i+1]));
      }
      return fs;
    }

    inline std::pair<size_t,size_t> face_adj_tets(size_t index) const {
      assert((index>>1)<faces_tets().size());
      std::pair<size_t,size_t> t2 = faces_tets()[index>>1];
      if (t2.first==UNSIGN_NEG_ONE) t2.first = INVALID_ID;
      if (t2.second==UNSIGN_NEG_ONE) t2.second = INVALID_ID;
      if (index%2==1) {
        if (t2.first < t2.second) std::swap(t2.first, t2.second);
      } else {
        if (t2.first > t2.second) std::swap(t2.first, t2.second);
      }
      return t2;
    }

    // set status

    inline void set_face_status(size_t index, status s)
    {
      faces_status_[index>>1] |= s; // 2*i and (2*i+1) -> i
    }

    inline void set_vert_status(size_t index, status s)
    {
      verts_status_[index] |= s;
    }

    inline void set_tet_status(size_t index, status s)
    {
      tets_status_[index] |= s;
    }

    inline void set_edge_status(size_t index, status s)
    {
      edges_status_[index>>1] |= s; // 2*i and (2*i+1) -> i
    }

    inline bool is_surf_tet(size_t index) const
    {
      return ((tets_status_[index]&SURF_TET)>0);
    }

    inline bool is_surf_face(size_t index) const
    {
      return ((faces_status_[index>>1]&SURF_FACE)>0);
    }

    inline bool is_surf_vert(size_t index) const
    {
      return ((verts_status_[index]&SURF_VERT)>0);
    }

    inline bool is_surf_edge(size_t index) const
    {
      return ((edges_status_[index>>1]&SURF_EDGE)>0);
    }

    inline matrixd tet_center(size_t index) const
    {
      assert(index < tets_num());
      matrixd c = zjucad::matrix::zeros<double>(3,1);
      std::vector<size_t> tet_vs = tet_verts(index);
      assert(tet_vs.size() == 4);
      for (size_t vi = 0; vi < 4; ++vi)
        c += verts_coord()(zjucad::matrix::colon(), tet_vs[vi]);
      c *= 0.25;
      return c;
    }

    inline matrixd face_center(size_t index) const
    {
      assert(index < direc_faces_num());
      matrixd c = zjucad::matrix::zeros<double>(3,1);
      std::vector<size_t> vid = face_verts(index);
      assert(vid.size() == 3);
      for (size_t pi=0; pi<vid.size(); ++pi)
        c += verts_(zjucad::matrix::colon(), vid[pi]);
      return c/3.0;
    }

    inline double cal_tet_vol(size_t index) const
    {
      //double r = 0.0;
      using zjucad::matrix::colon;
      matrixst vid = tets_verts()(colon(), index);
      matrixd p1 = verts_(colon(), vid[0]);
      matrixd p2 = verts_(colon(), vid[1]);
      matrixd p3 = verts_(colon(), vid[2]);
      matrixd p4 = verts_(colon(), vid[3]);

      using zjucad::matrix::cross;
      using zjucad::matrix::dot;
      return (dot(cross(p2-p1,p3-p1),(p4-p1)))/6.0;
    }

    inline int cal_face_norm(size_t index, matrixd& norm) const
    {
      using zjucad::matrix::colon;
      assert(index < direc_faces_num());
      const std::vector<size_t>& vid = faces_verts()[index>>1];
      matrixd p1 = verts_(colon(), vid[0]);
      matrixd p2 = verts_(colon(), vid[1]);
      matrixd p3 = verts_(colon(), vid[2]);
      p2 -= p1; p3 -= p1;
      p1 = zjucad::matrix::cross(p2, p3);
      double len = zjucad::matrix::norm(p1);
      assert(len>1e-10);
      if (len < 1e-10) norm = zjucad::matrix::zeros<double>(3,1);
      else norm = p1/len;
      return 0;
    }

    const std::vector<std::pair<size_t,size_t> >& faces_tets() const {return fa_->face2tet_;}    
    std::vector<std::pair<size_t,size_t> >& faces_tets() {return fa_->face2tet_;}
    
 private:
    
    int build_face_tet_relation();
    int build_vert_edge_relation();
    std::vector<size_t> surf_edge_adj_loop_faces(size_t index) const;

    const std::vector<std::vector<size_t> >& faces_verts() const { return fa_->faces_; }
    std::vector<std::vector<size_t> >& faces_verts() { return fa_->faces_; }
    /// s->t or t->s, not all


    
 private:

    std::shared_ptr<jtf::mesh::face2tet_adjacent> fa_;
    matrixd& verts_;
    matrixst& tets_verts_;
    std::vector<std::pair<size_t,size_t> > edges_verts_;

    std::vector<std::vector<size_t> > t2t_;
    std::vector<std::vector<size_t> > t2f_;
    std::vector<std::vector<size_t> > v2v_;
    std::vector<std::vector<size_t> > v2e_;
    
    //std::vector<std::vector<size_t> > f2t_; // store the half data ( a<b, not b>a)
    std::vector<std::vector<size_t> > e2t_; // store the half data ( a<b, not b>a)
    //std::vector<std::vector<size_t> > e2f_; // store the half data ( a<b, not b>a)

    std::vector<size_t> tets_status_;
    std::vector<size_t> faces_status_;
    std::vector<size_t> verts_status_;
    std::vector<size_t> edges_status_;
  };
  
} // namespace 

#endif
