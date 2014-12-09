#include "type_patch_graph.h"
#include "remove_surface_wedge.h"

using namespace std;

int type_patch_graph::build_graph(
    const matrixst &tet,
    const matrixst &cut_tet,
    const matrixd &node,
    const matrixst &cut_tet2tet,
    const matrixst &outside_face,
    const matrixst &outside_face_idx,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const jtf::mesh::edge2cell_adjacent &ea,
    const vector<pair<size_t,size_t> > &g_unknown_face_pair)
{
  return extract_surface_patch_graph(
        tet, cut_tet, node, cut_tet2tet, outside_face, outside_face_idx, fa,
        fa_cut, ea, g_unknown_face_pair, surface_type_, patches_,
        patch_linking_info_);
}

int type_patch_graph::merge_group(const size_t &g0, const size_t &g1)
{
  if(g0 >= patches_.size() || g1 >= patches_.size()){
    cerr << "# [error] wrong group idx." << endl;
    return __LINE__;
  }

  patches_[g1].insert(patches_[g0].begin(), patches_[g0].end());

  boost::unordered_map<size_t, boost::unordered_set<size_t> >::iterator it_g0
      = patch_linking_info_.find(g0);
  boost::unordered_map<size_t, boost::unordered_set<size_t> >::iterator it_g1
      = patch_linking_info_.find(g1);

  if(it_g0 == patch_linking_info_.end() || it_g1 == patch_linking_info_.end()){
    cerr << "# [info] already merged." << endl;
    return 0;
  }

  it_g1->second.insert(it_g0->second.begin(), it_g0->second.end());

  return del_group(g0);
}

int type_patch_graph::del_group(const size_t &g)
{
  if(g >= patches_.size()){
    cerr << "# [error] wrong patch number." << endl;
    return  __LINE__;
  }

  boost::unordered_map<size_t, boost::unordered_set<size_t> >::iterator it_g
      = patch_linking_info_.find(g);

  if(it_g == patch_linking_info_.end() && patches_[g].empty()){
    cerr << "# [info] already deleated." << endl;
    return 0;
  }else if(it_g == patch_linking_info_.end() && !patches_[g].empty()){
    cerr << "# [error] strange: patch linking info broken." << endl;
    return __LINE__;
  }else if(it_g != patch_linking_info_.end() && patches_[g].empty()){
    cerr << "# [error] strange: patch linking info broken." << endl;
    return __LINE__;
  }else {
    for(boost::unordered_set<size_t>::const_iterator cit = it_g->second.begin();
        cit != it_g->second.end(); ++cit){
      const size_t & other_patch_idx = *cit;
      boost::unordered_map<size_t, boost::unordered_set<size_t> >::iterator other_it
          = patch_linking_info_.find(other_patch_idx);
      if(other_it == patch_linking_info_.end()){
        cerr << "# [error] can not find linking info of " << other_patch_idx << endl;
        return __LINE__;
      }
      boost::unordered_set<size_t> & other_linking = other_it->second;
      boost::unordered_set<size_t>::iterator it = other_linking.find(g);
      if(it == other_linking.end()){
        cerr << "# [error] stragne can not find linking from " << other_patch_idx
             << " to " << g << endl;
        return __LINE__;
      }
      other_linking.erase(it);
    }
    patch_linking_info_.erase(it_g);
  }

  return 0;
}

size_t type_patch_graph::get_group_type(
    const size_t &g,
    const bool check) const
{
  if(g >= patches_.size() || patches_[g].empty()){
    return -1;
  }

  const size_t & first_face_idx = *patches_[g].begin();
  boost::unordered_map<size_t,size_t>::const_iterator cit =
      surface_type_.find(first_face_idx);
  if(cit == surface_type_.end()){
    cerr << "# [error] strange can not find face type of " << first_face_idx << endl;
    return __LINE__;
  }
  return cit->second;

  if(check){
    const boost::unordered_set<size_t> & group = patches_[g];
    for(boost::unordered_set<size_t>::const_iterator ccit = group.begin();
        ccit != group.end(); ++ccit){
      boost::unordered_map<size_t,size_t>::const_iterator tit =
          surface_type_.find(*ccit);
      if(tit == surface_type_.end()){
        cerr << "# [error] can not find surface type of " << *ccit << endl;
        return __LINE__;
      }
      if(tit->second != cit->second){
        cerr << "# [error] this group contains different type faces." << endl;
        return __LINE__;
      }
    }
  }
}

const boost::unordered_set<size_t> & type_patch_graph::get_linking_of_group(
    const size_t &g) const
{
  boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
      cit = patch_linking_info_.find(g);
  if(cit == patch_linking_info_.end()){
    throw "# [error] can not find group." ;
  }

  return cit->second;
}

size_t type_patch_graph::get_patch_degree(const size_t & g) const
{
  linking_type::const_iterator cit = patch_linking_info_.find(g);
  if(cit == patch_linking_info_.end())
    return -1;
  else
    return cit->second.size();
}
