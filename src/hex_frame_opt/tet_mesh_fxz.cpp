#include "tet_mesh_fxz.h"

namespace fxz {

inline void print(const std::string& s) { std::cout << s << std::endl; }
inline void print_err(const std::string& s) { std::cerr << s <<std::endl;}


int tet_mesh::build()
{
  print("# [ TET MESH ] Start to build...");
  if (build_face_tet_relation()) return __LINE__;
  if (build_vert_edge_relation()) return __LINE__;
  print("# [ TET MESH ] End of building.");
  return 0;
}

int tet_mesh::build_face_tet_relation()
{
  print("# [ FACE2TET ] begin...");
  fa_.reset(jtf::mesh::face2tet_adjacent::create(tets_verts_));
  print("# [ FACE2TET ] end.");


  print("# Build face tet relation...");
  tets_status_.resize(tets_num(), 0);
  verts_status_.resize(verts_num(), 0);
  faces_status_.resize(undir_faces_num(), 0);

  t2t_.resize(tets_num());
  t2f_.resize(tets_num());

  for (size_t i=0; i<tets_num(); ++i) {
    t2t_[i].reserve(4);
    t2f_[i].reserve(4);
  }

  for (size_t fi=0; fi<faces_tets().size(); ++fi) {
    size_t a = faces_tets()[fi].first;
    size_t b = faces_tets()[fi].second;
    if (a==UNSIGN_NEG_ONE || b==UNSIGN_NEG_ONE) {
      set_face_status(2*fi, SURF_FACE); // directed face
      set_face_status(2*fi+1, SURF_FACE);
      if (a != UNSIGN_NEG_ONE) set_tet_status(a, SURF_TET);
      if (b != UNSIGN_NEG_ONE) set_tet_status(b, SURF_TET);
      for (size_t vi=0; vi<faces_verts()[fi].size(); ++vi) {
        set_vert_status(faces_verts()[fi][vi], SURF_VERT);
      }
    }
    
    if (a==b) continue;

    if (a == UNSIGN_NEG_ONE) {
      a = INVALID_ID;
      std::swap(a, b);
   } else if (b == UNSIGN_NEG_ONE) {
      b = INVALID_ID;
    } else {
      if (a > b) std::swap(a, b);
    }
    
    std::vector<size_t>::iterator it;
    if (b!=INVALID_ID) it = std::find(t2t_[a].begin(), t2t_[a].end(), b);
    else it = t2t_[a].end();
    if (it == t2t_[a].end()) {
      t2t_[a].push_back(b);
      t2f_[a].push_back(2*fi);
    }

    if (b != INVALID_ID) {
      it = std::find(t2t_[b].begin(), t2t_[b].end(), a);
      if (it == t2t_[b].end()) {
        t2t_[b].push_back(a);
        t2f_[b].push_back(2*fi+1);
      }
    }
  }
  
  return 0;
}

int tet_mesh::build_vert_edge_relation()
{
  jtf::mesh::one_ring_tet_at_edge ortae;
  ortae.add_tets(tets_verts_, *fa_);

  if(ortae.sort_into_loop(tets_verts_, verts_)) {
    std::cerr << "# sort error." << std::endl;
    return __LINE__;
  }

  print("# Build vert edge relation...");
  edges_verts_.resize(ortae.e2t_.size());
  e2t_.resize(edges_verts_.size());

  jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator eit;
  size_t cnt = 0;
  for (eit=ortae.e2t_.begin(); eit != ortae.e2t_.end(); ++eit, ++cnt) {
    const std::vector<size_t>& loop = eit->second;
    const std::pair<size_t, size_t>& e = eit->first;
    assert(e.first < e.second);
    edges_verts_[cnt] = e;
    e2t_[cnt] = loop;
  }

  v2v_.resize(verts_num());
  v2e_.resize(verts_num());
  for (size_t i=0; i<verts_num(); ++i) {
    v2v_[i].reserve(6); v2e_[i].reserve(6);
  }

  for (size_t ei=0; ei<undir_edges_num(); ++ei) {
    size_t a = edges_verts_[ei].first;
    size_t b = edges_verts_[ei].second;
    std::vector<size_t>::iterator it;
    it = std::find(v2v_[a].begin(), v2v_[a].end(), b);
    if (it==v2v_[a].end()) {
      v2v_[a].push_back(b); v2e_[a].push_back(2*ei);
    }
    it = std::find(v2v_[b].begin(), v2v_[b].end(), a);
    if (it==v2v_[b].end()) {
      v2v_[b].push_back(a); v2e_[b].push_back(2*ei+1);
    }
  }

  edges_status_.resize(ortae.e2t_.size());
  assert(edges_status_.size() == undir_edges_num());

  for (size_t i=0; i<direc_faces_num(); i+=2) {
    if (!is_surf_face(i)) continue;
    const std::vector<size_t>& fid = face_verts(i);
    assert(fid.size() == 3);
    for (size_t j=0; j<3; ++j) {
      size_t eid = edge_index(fid[j],fid[(j+1)%3]);
      set_edge_status(eid, SURF_EDGE);
    }
  }

  print("# CHECK ...");
  for (eit=ortae.e2t_.begin(); eit != ortae.e2t_.end(); ++eit, ++cnt) {
    const std::vector<size_t>& loop = eit->second;
    const std::pair<size_t, size_t>& e = eit->first;

    if (!is_surf_edge(edge_index(e.first,e.second))) {
      if (loop[0]==UNSIGN_NEG_ONE) {
        print("# [ ERROR ]");
        return __LINE__;
      }
    }
  }
  print("# CHECK end.");
  
  
  return 0;
}


//!TODO this implementation isn't sorted
std::vector<size_t> tet_mesh::surf_edge_adj_loop_faces(size_t index) const
{
  std::vector<size_t> fs;
  const std::vector<size_t>& loop = e2t_[index>>1];

  for (size_t i=0; i<loop.size()-1; ++i) {
    if (loop[i]!=UNSIGN_NEG_ONE && loop[i+1]!=UNSIGN_NEG_ONE) {
      fs.push_back(face_index(loop[i],loop[i+1]));
    }
  }
  
  std::pair<size_t,size_t> evs = edge_verts(index);
  size_t x=INVALID_ID;
  if (loop[0]==UNSIGN_NEG_ONE && loop[1]!=UNSIGN_NEG_ONE) x = 1;
  else if (loop[0] != UNSIGN_NEG_ONE) x = 0;

  if (x!=INVALID_ID) {
    const std::vector<size_t>& tfs = tet_adj_faces(loop[x]);
    const std::vector<size_t>& tts = tet_adj_tets(loop[x]);
    assert(tts.size() == 4);
    for (size_t i=0; i<tts.size(); ++i)
      if (tts[i]==INVALID_ID) {
        const std::vector<size_t> fvs = face_verts(tfs[i]);
        size_t cnt = 0;
        for (size_t j=0; j<3; ++j) {
          if (fvs[j]==evs.first || fvs[j]==evs.second) ++cnt;
        }
        if (cnt == 2) fs.push_back(tfs[i]);
      }
  }
  x = INVALID_ID;
  if (loop[loop.size()-1]==UNSIGN_NEG_ONE && loop[loop.size()-2]!=UNSIGN_NEG_ONE) x = loop.size()-2;
  else if (loop[loop.size()-1]!=UNSIGN_NEG_ONE) x = loop.size()-1;

  if (x!=INVALID_ID) {
    size_t tid = loop[x];
    const std::vector<size_t>& tfs = tet_adj_faces(tid);
    const std::vector<size_t>& tts = tet_adj_tets(tid);
    for (size_t i=0; i<tts.size(); ++i)
      if (tts[i]==INVALID_ID) {
        const std::vector<size_t> fvs = face_verts(tfs[i]);
        size_t cnt = 0;
        for (size_t j=0; j<3; ++j) {
          if (fvs[j]==evs.first || fvs[j]==evs.second) ++cnt;
        }
        if (cnt == 2) fs.push_back(tfs[i]);
      }
  }

  return fs;
}

} // namespace fxz
