#include <boost/tuple/tuple_comparison.hpp>
#include <memory>
#include <numeric>
#include <stack>
#include <fstream>

#include <jtflib/util/container_operation.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>


#include "remove_surface_wedge.h"
#include "type_patch_graph.h"
#include "cut_tet.h"
#include "io.h"
#include "../tetmesh/tetmesh.h"
#include "../common/util.h"
#include "../common/vtk.h"
#include "../common/visualize_tool.h"
#include "../common/transition_type.h"
#include "../tet_mesh_sxx/tet_mesh_sxx.h"
#include "../tet_mesh_sxx/edge_split_in_mid.h"

using namespace std;
using namespace zjucad::matrix;

//////////////// these two items are used to count type issue /////////////
int boundary_issue_num = 0;
int degree_issue_num = 0;

int remove_isolated_patch(
    const matrixst & tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent & fa,
    const vector<boost::unordered_set<size_t> > patches,
    const boost::unordered_map<size_t, boost::unordered_set<size_t> > &patch_links_info,
    boost::unordered_map<size_t,size_t> & surface_type_original)
{
  for(size_t gi = 0; gi < patches.size(); ++gi){
    boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
        cit = patch_links_info.find(gi);
    if(cit == patch_links_info.end()){
      cerr << "# [error] strange can not find patch " << gi << endl;
      return __LINE__;
    }
    if(cit->second.size() == 1){ // isolated patch
      const size_t &other_group_idx = *cit->second.begin();
      const size_t &other_group_face_idx = *patches[other_group_idx].begin();
      boost::unordered_map<size_t,size_t>::const_iterator mcit =
          surface_type_original.find(other_group_face_idx);
      if(mcit == surface_type_original.end()){
        cerr << "# [error] strange can not find surface type of face "
             << other_group_face_idx << endl;
        return __LINE__;
      }
      const size_t & modified_type = mcit->second;
      ++degree_issue_num;
      for(boost::unordered_set<size_t>::const_iterator bumcit = patches[gi].begin();
          bumcit != patches[gi].end(); ++bumcit){
        boost::unordered_map<size_t,size_t>::iterator fit =
            surface_type_original.find(*bumcit);
        if(fit == surface_type_original.end()){
          cerr << "# [error] strange can not find face " << *bumcit << endl;
          return __LINE__;
        }
        fit->second = modified_type;
      }
    }
  }

  dump_surface_restricted_type_to_vtk(
        "surface_type_original_after_remove_isolated_patch.vtk",
        "patch", node, fa, surface_type_original);

  return 0;
}

static int use_surface_normal_to_modify_surface_type(
    const boost::unordered_set<size_t> & patches_along_chains,
    const boost::unordered_map<size_t, matrixd> & surface_normal,
    const size_t &current_patch_type,
    boost::unordered_map<size_t,size_t> &surface_type_original)
{
  matrixd combined_face_normal = zjucad::matrix::zeros<double>(3);
  for(boost::unordered_set<size_t>::const_iterator cit
      = patches_along_chains.begin(); cit != patches_along_chains.end(); ++cit){
    boost::unordered_map<size_t, matrixd>::const_iterator mcit =
        surface_normal.find(*cit);
    if(mcit == surface_normal.end()){
      cerr << "# [error] strange can not find surface normal of face " << *cit << endl;
      return __LINE__;
    }
    combined_face_normal += mcit->second;
  }

  const double len = norm(combined_face_normal);
  if(len > 1e-6)
    combined_face_normal /= len;
  if(current_patch_type > 2){
    cerr << "# [error] strange patch type should only be [0/1/2] not "
         << current_patch_type << endl;
    return __LINE__;
  }

  const matrixd I = eye<double>(3);
  vector<pair<double,int> > res;
  for(size_t i = 1; i < 3; ++i){
    res.push_back(make_pair(fabs(dot(combined_face_normal,
                                     I(colon(), (i + current_patch_type)%3))),
                            (i + current_patch_type)%3));
  }

  sort(res.begin(), res.end());

  const size_t &modified_type = res.back().second;
  assert(modified_type != current_patch_type);

  for(boost::unordered_set<size_t>::const_iterator cit = patches_along_chains.begin();
      cit != patches_along_chains.end(); ++cit){
    boost::unordered_map<size_t,size_t>::iterator mit =
        surface_type_original.find(*cit);
    if(mit == surface_type_original.end()){
      cerr << "# [error] strange can not find surface type of face " << *cit << endl;
      return __LINE__;
    }
    mit->second = modified_type;
  }

  return 0;
}


int modify_surface_type_to_remove_separated_edges(
    const std::deque<std::pair<size_t,size_t> > & one_chain,
    const boost::unordered_set<size_t> &one_patch,
    const matrixst & outside_face_idx,
    const boost::unordered_map<size_t, matrixd> & surface_normal,
    boost::unordered_map<size_t,size_t> & surface_type_original,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::N_ring_face_at_point &nrfap,
    const size_t N_ring_to_modify)
{
  boost::unordered_set<size_t> patches_along_chains;
  boost::unordered_set<size_t> chain_points;

  for(size_t ei = 0; ei < one_chain.size(); ++ei){
    chain_points.insert(one_chain[ei].first);
    chain_points.insert(one_chain[ei].second);
  }

  // gather N ring faces along the separaed chain
  for(boost::unordered_set<size_t>::const_iterator cit = chain_points.begin();
      cit != chain_points.end(); ++cit){
    jtf::mesh::N_ring_face_at_point::p2f_type::const_iterator pcit
        = nrfap.p2f_.find(*cit);
    if(pcit == nrfap.p2f_.end()){
      cerr << "# [error] strange can not find surface of point " << *cit << endl;
      return __LINE__;
    }
    const vector<vector<size_t> > & N_ring_faces = pcit->second;
    for(size_t ni = 0; ni < N_ring_faces.size(); ++ni)
      patches_along_chains.insert(N_ring_faces[ni].begin(), N_ring_faces[ni].end());
  }

  // collect those faces which are inside this patch
  {
    boost::unordered_set<size_t> temp;
    for(boost::unordered_set<size_t>::const_iterator cit = patches_along_chains.begin();
        cit != patches_along_chains.end(); ++cit){
      if(one_patch.find(outside_face_idx[*cit]) != one_patch.end()){
        temp.insert(outside_face_idx[*cit]);}
    }
    patches_along_chains = temp;
  }

  boost::unordered_map<size_t,size_t>::const_iterator cit =
      surface_type_original.find(*one_patch.begin());
  if(cit == surface_type_original.end()){
    cerr << "# [error] strange can not find surface type of face "
         << *one_patch.begin() << endl;
    return __LINE__;
  }

  const size_t current_patch_type = cit->second;

  if(one_chain.front().first == one_chain.back().second){

    use_surface_normal_to_modify_surface_type(
          patches_along_chains, surface_normal, cit->second, surface_type_original);
    return 0;

  }else{
    jtf::mesh::N_ring_face_at_point::p2f_type::const_iterator cit_first =
        nrfap.p2f_.find(one_chain.front().first);
    jtf::mesh::N_ring_face_at_point::p2f_type::const_iterator cit_second =
        nrfap.p2f_.find(one_chain.back().second);
    if(cit_first == nrfap.p2f_.end() || cit_second == nrfap.p2f_.end()){
      cerr << "# [error] can not find point in nrfap." << endl;
      return __LINE__;
    }

    boost::unordered_set<size_t> uvw_type_around_chain_ends;
    for(size_t fi = 0; fi < cit_first->second[0].size(); ++fi){
      boost::unordered_map<size_t,size_t>::const_iterator cit =
          surface_type_original.find(outside_face_idx[cit_first->second[0][fi]]);
      if(cit == surface_type_original.end()){
        cerr << "# [error] strange can not find surface type of face "
             << outside_face_idx[cit_first->second[0][fi]] << endl;
        return __LINE__;
      }
      uvw_type_around_chain_ends.insert(cit->second);
    }

    for(size_t fi = 0; fi < cit_second->second[0].size(); ++fi){
      boost::unordered_map<size_t,size_t>::const_iterator cit =
          surface_type_original.find(outside_face_idx[cit_second->second[0][fi]]);
      if(cit == surface_type_original.end()){
        cerr << "# [error] strange can not find surface type of face "
             << outside_face_idx[cit_second->second[0][fi]] << endl;
        return __LINE__;
      }
      uvw_type_around_chain_ends.insert(cit->second);
    }

    if(uvw_type_around_chain_ends.size() == 3){
      use_surface_normal_to_modify_surface_type(
            patches_along_chains, surface_normal, current_patch_type, surface_type_original);
    }else if(uvw_type_around_chain_ends.size() == 2){
      const size_t modified_type =
          (0+1+2) - std::accumulate(uvw_type_around_chain_ends.begin(),
                                    uvw_type_around_chain_ends.end(),
                                    static_cast<size_t>(0));
      for(boost::unordered_set<size_t>::const_iterator buscit
          = patches_along_chains.begin(); buscit != patches_along_chains.end();
          ++buscit){
        boost::unordered_map<size_t,size_t>::iterator mit
            = surface_type_original.find(*buscit);
        if(mit == surface_type_original.end()){
          cerr << "# [error] strange can not find surface type of face "
               << *buscit << endl;
          return __LINE__;
        }
        assert(mit->second != modified_type);
        mit->second = modified_type;
      }
    }
  }

  return 0;
}

int remove_multi_orientation_group(
    const matrixd & node,
    const jtf::mesh::edge2cell_adjacent & ea,
    const jtf::mesh::face2tet_adjacent & fa,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const boost::unordered_map<size_t, matrixd> & surface_normal,
    const boost::unordered_map<size_t,size_t> &surface_type_orient_orig,
    const std::vector<boost::unordered_set<size_t> > &patches,
    const boost::unordered_map<size_t, boost::unordered_set<size_t> > &patch_links_info,
    const size_t N_ring_to_modify,
    boost::unordered_map<size_t,size_t> &surface_type_original)
{
  jtf::mesh::N_ring_face_at_point nrfap;
  nrfap.add_all_faces(outside_face, ea, N_ring_to_modify);

#define debug
#ifdef debug
  vector<deque<pair<size_t,size_t> > > chains_vec;
#endif

  for(size_t pi = 0; pi < patches.size(); ++pi){
    const boost::unordered_set<size_t> & one_patch = patches[pi];
    size_t orientation = -1;
    bool is_multi_orientation_group = false;

    // check wether is a multi_orientation_patch
    for(boost::unordered_set<size_t>::const_iterator cit = one_patch.begin();
        cit != one_patch.end(); ++cit){
      boost::unordered_map<size_t,size_t>::const_iterator mcit =
          surface_type_orient_orig.find(*cit);
      if(mcit == surface_type_orient_orig.end()){
        cerr << "# [error] can not find " << *cit << " in surface_type_orient_orig." << endl;
        return __LINE__;
      }

      if(orientation == -1)
        orientation = mcit->second;
      if(mcit->second / 2 != orientation / 2){
        cerr << "# [error] strange this patch contains multi uvw axis restricted type." << endl;
        return __LINE__;
      }

      if(mcit->second != orientation){
        is_multi_orientation_group = true;
        break;
      }
    }

    if(is_multi_orientation_group){
      // find separated edges which adjacent to different orientations
      //boost::unordered_set<size_t> gathered_orients;
      boost::unordered_map<pair<size_t,size_t>, boost::unordered_set<size_t> > edge_with_orient;
      for(boost::unordered_set<size_t>::const_iterator buscit = one_patch.begin();
          buscit != one_patch.end(); ++buscit){
        boost::unordered_map<size_t,size_t>::const_iterator itt =
            surface_type_orient_orig.find(*buscit);
        if(itt == surface_type_orient_orig.end()){
          cerr << "# [error] strange can not find surface type of face "
               << *buscit << endl;
          return __LINE__;
        }
        //gathered_orients.insert(itt->second);
        const vector<size_t> & face_vec = fa.faces_.at(*buscit);
        for(size_t p = 0; p < face_vec.size(); ++p){
          pair<size_t,size_t>  edge(face_vec[p], face_vec[(p+1)%face_vec.size()]);
          if(edge.first > edge.second)
            swap(edge.first, edge.second);
          edge_with_orient[edge].insert(itt->second);
        }
      }

      vector<pair<size_t,size_t> > separaet_edges;
      for(boost::unordered_map<std::pair<size_t,size_t>,
          boost::unordered_set<size_t> >::const_iterator
          eit = edge_with_orient.begin(); eit != edge_with_orient.end(); ++eit){
        if(eit->second.size() == 2) separaet_edges.push_back(eit->first);
      }

      vector<deque<pair<size_t,size_t> > > chains;
      jtf::util::extract_chain_from_edges(separaet_edges,chains);

#ifdef debug
      chains_vec.insert(chains_vec.end(), chains.begin(), chains.end());
#endif

      for(size_t ci = 0 ; ci < chains.size(); ++ci)
        modify_surface_type_to_remove_separated_edges(
              chains[ci], one_patch, outside_face_idx, surface_normal,
              surface_type_original, ea, nrfap, N_ring_to_modify);
    }
  }

#ifdef debug
  dump_singularity_to_vtk("separated_edges.vtk",node, chains_vec);
#endif
  return 0;
}

int remove_fractional_degree_two_patches(
    const matrixst & tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent & fa,
    const vector<boost::unordered_set<size_t> > patches,
    const boost::unordered_map<size_t, boost::unordered_set<size_t> > &patch_links_info,
    boost::unordered_map<size_t,size_t> & surface_type_original,
    const boost::unordered_map<pair<size_t,size_t>,double> &patch2patch_len,
    const size_t fractional_cell_num )
{
  //const size_t fractional_cell_num = 10;
  for(size_t gi = 0; gi < patches.size(); ++gi){
    boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
        cit = patch_links_info.find(gi);
    if(cit == patch_links_info.end()){
      cerr << "# [error] strange can not find patch " << gi << endl;
      return __LINE__;
    }
    if((cit->second.size() == 2) && (patches[gi].size() < fractional_cell_num)){
      size_t other_group_idx = -1;
      double max_len = 0;
      const boost::unordered_set<size_t> & other_linkes  = cit->second;
      for(boost::unordered_set<size_t>::const_iterator cit = other_linkes.begin();
          cit != other_linkes.end(); ++cit){
        boost::unordered_map<pair<size_t,size_t>,double>::const_iterator p2pit
            = patch2patch_len.find(make_pair(gi, *cit));
        if(p2pit == patch2patch_len.end()){
          cerr << "# [error] can not find patch " << gi << " to " << *cit
               << " length" << endl;
          return __LINE__;
        }
        if(p2pit->second > max_len)
          other_group_idx = *cit;
      }

      const size_t &other_group_face_idx = *patches[other_group_idx].begin();
      boost::unordered_map<size_t,size_t>::const_iterator mcit =
          surface_type_original.find(other_group_face_idx);
      if(mcit == surface_type_original.end()){
        cerr << "# [error] strange can not find surface type of face "
             << other_group_face_idx << endl;
        return __LINE__;
      }
      ++degree_issue_num;
      const size_t & modified_type = mcit->second;
      for(boost::unordered_set<size_t>::const_iterator bumcit = patches[gi].begin();
          bumcit != patches[gi].end(); ++bumcit){
        boost::unordered_map<size_t,size_t>::iterator fit =
            surface_type_original.find(*bumcit);
        if(fit == surface_type_original.end()){
          cerr << "# [error] strange can not find face " << *bumcit << endl;
          return __LINE__;
        }
        fit->second = modified_type;
      }
    }
  }

  return 0;
}


int remove_fractional_patches_by_degenerated_edges(
    const matrixst & tet,
    const matrixd & node,
    const matrixst & outside_face_idx,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::edge2cell_adjacent & ea,
    const boost::unordered_set<std::pair<size_t,size_t> > & degenerated_edges,
    boost::unordered_map<size_t,size_t> &surface_type_original)
{
  boost::unordered_map<size_t, boost::unordered_set<size_t> > face2edges;
  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit =
      degenerated_edges.begin(); cit != degenerated_edges.end(); ++cit){
    const pair<size_t,size_t> & edge = *cit;
    const size_t &edge_idx = ea.get_edge_idx(edge.first,edge.second);
    if(edge_idx == -1){
      continue;
    }

    const pair<size_t,size_t> & face_pair = ea.edge2cell_[edge_idx];
    face2edges[outside_face_idx[face_pair.first]].insert(edge_idx);
    face2edges[outside_face_idx[face_pair.second]].insert(edge_idx);
  }

  while(1){
    bool do_fix_a_degeneration = false;
    for(boost::unordered_map<size_t,boost::unordered_set<size_t> >::const_iterator
        cit = face2edges.begin(); cit != face2edges.end(); ++cit){
      if(cit->second.size() == 3){
        boost::unordered_map<size_t,size_t>::iterator bumit =
            surface_type_original.find(cit->first);
        if(bumit == surface_type_original.end()){
          cerr << "# [error] strange can not find surface type of "
               << bumit->first << endl;
          return __LINE__;
        }
        const size_t & current_type = bumit->second;
        vector<size_t> other_types(3,0);
        const boost::unordered_set<size_t> & adjacent_edges = cit->second;
        for(boost::unordered_set<size_t>::const_iterator scit =
            adjacent_edges.begin(); scit != adjacent_edges.end(); ++scit){
          const pair<size_t,size_t> & tri_pair = ea.edge2cell_[*scit];
          const size_t other_face_idx =
              outside_face_idx[tri_pair.first] +
              outside_face_idx[tri_pair.second] - cit->first;
          boost::unordered_map<size_t,size_t>::const_iterator other_face_it =
              surface_type_original.find(other_face_idx);
          if(other_face_it == surface_type_original.end()){
            cerr << "# [error] can not find surface type of " << other_face_idx
                 << endl;
            return __LINE__;
          }
          if(other_face_it->second != current_type){
            ++other_types[other_face_it->second];
          }
        }
        vector<size_t>::const_iterator max_e =
            max_element(other_types.begin(), other_types.end());
        if(*max_e < 2) continue;
        bumit->second = max_e - other_types.begin(); // modify the surface type
        do_fix_a_degeneration = true;
      }
    }
    if(!do_fix_a_degeneration) break;
  }

  return 0;
}

static int surface_type_cut2orig(
    const boost::unordered_map<size_t,size_t> & surface_type_cut,
    boost::unordered_map<size_t,size_t> & surface_type_orig,
    const matrixst & cut_tet2tet,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    boost::unordered_map<size_t,size_t> * orig_surf2cut_surf_map_ptr = 0)
{
  for(boost::unordered_map<size_t,size_t>::const_iterator cit =
      surface_type_cut.begin(); cit != surface_type_cut.end(); ++cit){
    const vector<size_t> & face_vec = fa_cut.faces_[cit->first];
    const size_t & face_idx =
        fa.get_face_idx(cut_tet2tet[face_vec[0]],
                        cut_tet2tet[face_vec[1]],
                        cut_tet2tet[face_vec[2]]);
    if(face_idx == -1){
      cerr << "# [error] strange can not find orig face of " << cit->first << endl;
      return __LINE__;
    }
    surface_type_orig[face_idx] = cit->second;
    if(orig_surf2cut_surf_map_ptr)
      (*orig_surf2cut_surf_map_ptr)[face_idx] = cit->first;
  }
  return 0;
}

int remove_surface_wedge(
    const matrixst & tet,
    const matrixd &node,
    const matrixst &cut_tet,
    const matrixd & cut_node,
    const vector<pair<size_t,size_t> > &g_unknown_face_pair,
    const boost::unordered_map<pair<size_t,size_t>,size_t > & restricted_edge_from_graph,
    const bool use_uvw_restricted_surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    matrixst & new_tet,
    matrixd & new_node,
    boost::unordered_map<size_t,size_t> & surface_type_cut,
    matrixd * zyz_ptr)
{
  boost::unordered_map<size_t,size_t> surface_type_orient_cut;
  if(!use_uvw_restricted_surface_type){
    surface_type_orient_cut = surface_type_cut;
    for(boost::unordered_map<size_t,size_t>::iterator it = surface_type_cut.begin();
        it != surface_type_cut.end(); ++it){
      it->second /= 2;
    }
  }


  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  unique_ptr<jtf::mesh::face2tet_adjacent> fa_cut(jtf::mesh::face2tet_adjacent::create(cut_tet));
  if(!fa.get() || !fa_cut.get()){
    cerr << "# [error] can not build face2tet adjacent." << endl;
    return __LINE__;
  }

  matrixst cut_tet2tet(max(cut_tet)+1);
  boost::unordered_map<size_t,size_t> original_surface2cut_surface;
  boost::unordered_map<size_t,size_t> surface_type_original;
  matrixst outside_face, outside_face_idx;
  {
    cut_tet2tet(cut_tet) = tet(colon());

    // the input surface type is defined on cut_tet, since the cut tet mesh is
    // very fractional, and original surface triangle will not be overlapped, so
    // we can transfor these types to original tet, thus we can separate these faces
    // into groups

    if(surface_type_cut2orig(surface_type_cut, surface_type_original, cut_tet2tet,
                             *fa, *fa_cut, &original_surface2cut_surface))
      return __LINE__;

    get_outside_face(*fa, outside_face);
    get_outside_face_idx(*fa, outside_face_idx);

    if(surface_type_original.size() != outside_face.size(2)){
      cerr << "# [error] uncompatabile surface type, surface_type size "
           << surface_type_original.size() << " outside_face_num "
           << outside_face.size(2) << endl;
      return __LINE__;
    }
  }

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  straighten_patch_boundary(outside_face, outside_face_idx, *ea,*fa,
                            surface_type_original, inner_face_jump_type);

  dump_surface_restricted_type_to_vtk(
        "surface_type_original_after_straighten_boundary.vtk",
        "patch", node, *fa, surface_type_original);

  vector<boost::unordered_set<size_t> > patches;
  boost::unordered_map<size_t, boost::unordered_set<size_t> > patch_links_info;
  boost::unordered_map<pair<size_t,size_t>, double> patch_to_patch_boundary_len;
  extract_surface_patch_graph(
        tet,cut_tet,node,cut_tet2tet, outside_face, outside_face_idx,
        *fa,*fa_cut,*ea, g_unknown_face_pair, surface_type_original,
        patches, patch_links_info, &patch_to_patch_boundary_len);

  cerr << "# [info] patch number " << patches.size() << endl;

  if(!use_uvw_restricted_surface_type){
    // use orientation to help remove multi_orientation group
    boost::unordered_map<size_t,size_t> surface_type_orient_orig;
    if(surface_type_cut2orig(surface_type_orient_cut,surface_type_orient_orig,
                             cut_tet2tet, *fa, *fa_cut))
      return __LINE__;

    matrixd face_normal;

    jtf::mesh::cal_face_normal(outside_face, node, face_normal);
    jtf::tetmesh::orient_face_normal_outside_tetmesh(
          tet, node, outside_face, outside_face_idx, *fa, face_normal);

    boost::unordered_map<size_t, matrixd > surface_normal;
    for(size_t fi = 0; fi < outside_face_idx.size(); ++fi)
      surface_normal[outside_face_idx[fi]] = face_normal(colon(), fi);

    remove_multi_orientation_group(
          node,
          *ea,*fa,outside_face, outside_face_idx, surface_normal,
          surface_type_orient_orig,patches, patch_links_info, 1,
          surface_type_original);

    extract_surface_patch_graph(
          tet,cut_tet,node,cut_tet2tet, outside_face, outside_face_idx,
          *fa,*fa_cut,*ea, g_unknown_face_pair, surface_type_original,
          patches, patch_links_info);
  }

  remove_isolated_patch(tet, node, *fa, patches,patch_links_info,
                        surface_type_original);

  cerr << "# [info] finish remove isolated patch." << endl;

  //  extract_surface_patch_graph(
  //        tet,cut_tet,node,cut_node,cut_tet2tet, outside_face, outside_face_idx,
  //        *fa,*fa_cut,*ea, g_unknown_face_pair, surface_type_original,
  //        patches, patch_links_info);

  remove_fractional_degree_two_patches(
        tet, node, *fa,  patches,patch_links_info,surface_type_original,patch_to_patch_boundary_len);

  cerr << "# [info] finish remove fractional_degree_two_patches." << endl;


  //  {
  //    surface_type_cut.clear();
  //    for(boost::unordered_map<size_t,size_t>::const_iterator cit =
  //        surface_type_original.begin(); cit != surface_type_original.end();
  //        ++cit){
  //      surface_type_cut[original_surface2cut_surface[cit->first]] = cit->second;
  //    }
  //    dump_surface_restricted_type(
  //          "surface_type_cut_after_remove_degeneration",
  //          surface_type_cut);
  //  }

  boost::unordered_map<pair<size_t,size_t>,size_t > restricted_edges;
  gather_restricted_edges_on_polycube(
        tet, *ea, outside_face_idx, surface_type_original,restricted_edges);

#if 0
  {
    ofstream ofs("surface_restricted_edge.vtk");
    vector<size_t> edges;
    vector<size_t> edge_type;
    for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit
        = restricted_edges.begin(); cit != restricted_edges.end(); ++cit){
      edges.push_back(cit->first.first);
      edges.push_back(cit->first.second);
      edge_type.push_back(cit->second);
    }
    line2vtk(ofs, &node[0], node.size(2), &edges[0], edges.size()/2);
    cell_data(ofs, &edge_type[0], edge_type.size(), "type");
  }
#endif

  cerr << "# [info] restricted_edges number " << restricted_edges.size() << endl;
  remove_one_face_degeneration(tet, outside_face, outside_face_idx,
                               restricted_edges, surface_type_original);

  // since the one_face_removing will change the restricted edges
  gather_restricted_edges_on_polycube(
        tet, *ea, outside_face_idx, surface_type_original,restricted_edges);

  cerr << "# [info] finish remove one_face_degeneration." << endl;

  // the straighten step can handel surface zigzag.
  //  if(!restricted_edge_from_graph.empty()){
  //    boost::unordered_map<pair<size_t,size_t>,size_t> restricted_edge_from_graph_orig;
  //    for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator
  //        cit = restricted_edge_from_graph.begin();
  //        cit != restricted_edge_from_graph.end(); ++cit){
  //      restricted_edge_from_graph_orig[make_pair(cut_tet2tet[cit->first.first],
  //                                                cut_tet2tet[cit->first.second])]
  //          = cit->second;
  //    }

  //    remove_surface_zigzag_by_restricted_edges_from_graph(
  //          tet,node, outside_face_idx,*fa, *ea, restricted_edge_from_graph_orig,
  //          surface_type_original, new_tet, new_node, inner_face_jump_type, zyz_ptr);

  //    cerr << "# [info] finish remove zigzag edges by restricted edges from graph." << endl;
  //  }else{
  //    remove_surface_zigzag_by_restricted_edges(
  //          tet, node, outside_face_idx, *ea, restricted_edges, *fa,new_tet,
  //          new_node, surface_type_original, inner_face_jump_type, zyz_ptr);
  //    cerr << "# [info] finish remove zigzag edges by restricted edges from surface patch." << endl;
  //  }

  unique_ptr<jtf::mesh::face2tet_adjacent> fa_new(jtf::mesh::face2tet_adjacent::create(new_tet));
  if(!fa_new.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent for new_tet." << endl;
    return __LINE__;
  }

  //  // if there is one edge does not exist on any patch boundary, but both of
  //  // it's ending points are on patch boundary, this edge should be splitted
  //  // or it will result in degeneration
  //  remove_face_degeneration_by_splitting(*fa_new, surface_type_original, new_node,
  //                                        new_tet, inner_face_jump_type);

  dump_surface_restricted_type_to_vtk(
        "surface_type_original_after_modified.vtk",
        "patch", new_node, *fa_new, surface_type_original);

  dump_surface_restricted_type(
        "surface_type_original_after_remove_degeneration",
        surface_type_original);

  if(!inner_face_jump_type.empty())
    dump_inner_face_jump_type("inner_face_jump_type_after_remove_degeneration",
                              inner_face_jump_type);

  cerr << "# [info] +++++++++++++type issue statistics+++++++++++++" << endl;
  cerr << "# [info] boundary issue number " << boundary_issue_num << endl;
  cerr << "# [info] degree issue number " << degree_issue_num << endl;
  return 0;
}

//!!!! arounding_faces must be ordered, or the result is not right!!!!
bool is_critical_point(
    const vector<size_t> & arounding_faces,
    const matrixst & outside_face_idx,
    const boost::unordered_map<size_t,size_t> & surface_type)
{
  vector<size_t> arounding_type;
  for(size_t i = 0; i < arounding_faces.size(); ++i){
    const size_t &face_idx  = outside_face_idx[arounding_faces[i]];
    boost::unordered_map<size_t,size_t>::const_iterator cit =
        surface_type.find(face_idx);
    if(cit == surface_type.end()){
      cerr << "# [error] can not find surface idx " << face_idx << endl;
      return false;
    }
    const size_t &face_type = cit->second;
    if(arounding_type.empty() || arounding_type.back() != face_type)
      arounding_type.push_back(face_type);
  }
  assert(!arounding_type.empty());
  if(arounding_type.back() == arounding_type.front())
    arounding_type.pop_back();

  assert(arounding_type.size() != 1);

  if(arounding_type.empty()) return false;
  if(arounding_type.size() == 2) return false;

  if(arounding_type.size() > 3) return true; // since type can only be 0 ,1 ,2

  if(arounding_type.size() == 3){
    boost::unordered_set<size_t> type_set(
          arounding_type.begin(), arounding_type.end());
    if(type_set.size() == arounding_type.size())
      return false;
    return true;
  }

  return false;
}

int remove_face_degeneration_by_splitting(
   jtf::mesh::face2tet_adjacent &fa,
    boost::unordered_map<size_t,size_t> &surface_type_original,
    matrix<double> &new_node,
    matrix<size_t> &new_tet,
    boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type)
{
  matrix<size_t> outside_face, outside_face_idx;
  get_outside_face(fa, outside_face);
  get_outside_face_idx(fa, outside_face_idx);

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  matrix<size_t> surface_type(outside_face_idx.size(),1);
  for(size_t fi = 0; fi < surface_type.size(); ++fi){
    boost::unordered_map<size_t,size_t>::const_iterator cit =
        surface_type_original.find(outside_face_idx[fi]);
    if(cit == surface_type_original.end()){
      cerr << "# [error] can not find surface type of " << outside_face_idx[fi]
           << endl;
      return __LINE__;
    }
    surface_type[fi] = cit->second;
  }

  vector<pair<size_t,size_t> > patch_boundary_edges;
  boost::unordered_set<size_t> boundary_point;

  vector<pair<size_t,size_t> > not_boundary_edges;
  for(size_t ei = 0; ei < ea->edge2cell_.size(); ++ei){
    const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
    if(surface_type[tri_pair.first] == surface_type[tri_pair.second]) {
      not_boundary_edges.push_back(ea->edges_[ei]);
      continue;
    }
    pair<size_t,size_t> one_edge = ea->edges_[ei];
    if(one_edge.first > one_edge.second)
      swap(one_edge.first, one_edge.second);
    patch_boundary_edges.push_back(one_edge);
    boundary_point.insert(one_edge.first);
    boundary_point.insert(one_edge.second);
  }

  vector<deque<pair<size_t,size_t> > > chains;
  jtf::util::extract_chain_from_edges(patch_boundary_edges, chains);

  boost::unordered_map<size_t, boost::unordered_set<size_t> > bounary_point2chains;
  for(size_t ci = 0; ci < chains.size(); ++ci){
    const deque<pair<size_t,size_t> > & one_chain = chains[ci];

    for(size_t ei = 0; ei < one_chain.size(); ++ei){
      const pair<size_t,size_t> & one_edge = one_chain[ei];
      bounary_point2chains[one_edge.first].insert(ci);
      bounary_point2chains[one_edge.second].insert(ci);
    }
  }

  boost::unordered_set<pair<size_t,size_t> > edges_need_to_split;

  for(size_t ei = 0; ei < not_boundary_edges.size(); ++ei){
    const pair<size_t,size_t> &one_edge = not_boundary_edges[ei];
    if(boundary_point.find(one_edge.first) == boundary_point.end()) continue;
    if(boundary_point.find(one_edge.second) == boundary_point.end()) continue;
    boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
        cit_first = bounary_point2chains.find(one_edge.first);
    boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
        cit_second = bounary_point2chains.find(one_edge.second);
    if(cit_first == bounary_point2chains.end() ||
       cit_second == bounary_point2chains.end()){
      cerr << "# [error] boundary point can not find linking chain." << endl;
      return __LINE__;
    }
    vector<size_t> first_link(cit_first->second.begin(), cit_first->second.end());
    vector<size_t> second_link(cit_second->second.begin(), cit_second->second.end());
    vector<size_t> intersection_set(first_link.size() + second_link.size());
    vector<size_t>::iterator inter_end =
        find_intersection_set(first_link.begin(), first_link.end(),
                              second_link.begin(), second_link.end(),
                              intersection_set.begin());
    //if(inter_end == intersection_set.begin()) continue;
    edges_need_to_split.insert(one_edge);
  }

  ofstream ofs("edges_need_to_split");
  ofs << edges_need_to_split.size() << endl;
  cerr << "# [info] remove face degeneration, number " << edges_need_to_split.size() << endl;
  sxx::tet_mesh stm;
  stm.create_tetmesh(new_node, new_tet);
  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit
      = edges_need_to_split.begin(); cit != edges_need_to_split.end(); ++cit){
    ofs << cit->first << " " << cit->second << endl;
    stm.split_edge(*cit);
  }

  stm.write_tetmesh_to_matrix(new_node, new_tet);

  // update
  unique_ptr<jtf::mesh::face2tet_adjacent> fa_new(jtf::mesh::face2tet_adjacent::create(new_tet));
  if(!fa_new.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  //fa_new = *fa;
  matrix<size_t> tet_idx_map;
  boost::unordered_map<vector<size_t>, vector<size_t> > f2omap;
  sxx::get_face2orginal(stm, f2omap);
  sxx::get_tet2orginal_index(stm, new_tet, tet_idx_map);
  update_surface_type(fa, *fa_new, f2omap, surface_type_original);
  update_inner_face_jump_type(*fa_new, tet_idx_map, inner_face_jump_type);
  fa = *fa_new;
  return 0;
}

// this function is used to build a graph only with inner face jump type
//int remove_surface_degenerated_patch_prev_process(
//    const matrixst & tet,
//    const matrixst & cut_tet,
//    const matrixd & node,
//    const jtf::mesh::face2tet_adjacent & fa,
//    const jtf::mesh::face2tet_adjacent & fa_cut,
//    const jtf::mesh::one_ring_tet_at_edge & ortae,
//    const matrixst & outside_face,
//    const matrixst & outside_face_cut,
//    const matrixst & outside_face_idx,
//    const matrixst & cut_tet2tet,
//    const boost::unordered_set<size_t> & loop_points,
//    const boost::unordered_map<size_t,size_t> & inner_face_jump_type,
//    singularity_graph & sg_without_surface)
//{
//  matrixst face_pair_cut;
//  analysis_transition_raw(tet, node, fa, fa_cut, cut_tet, face_pair_cut);

//  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
//        jtf::mesh::edge2cell_adjacent::create(outside_face));
//  if(!ea.get()){
//    cerr << "# [error] can not build edge2tri adjacent" << endl;
//    return __LINE__;
//  }

//  jtf::mesh::N_ring_face_at_point nrfap;
//  const size_t ring_num = 2;
//  nrfap.add_all_faces(outside_face, ea, ring_num); // find 2 ring faces for each point
//  vector<size_t> possible_degenerated_faces;

//  for(boost::unordered_set<size_t>::const_iterator cit = loop_points.begin();
//      cit != loop_points.end(); ++cit){
//    jtf::mesh::N_ring_face_at_point::p2f_type::const_iterator ccit =
//        nrfap.p2f_.find(*cit);
//    const vector<vector<size_t> > & N_ring_faces = ccit->second;
//    assert(ring_num == N_ring_faces.size());
//    for(size_t ri = 0; ri < ring_num; ++ri){
//      for(size_t fi = 0; fi < N_ring_faces[ri].size(); ++fi){
//        if(find(possible_degenerated_faces.begin(),
//                possible_degenerated_faces.end(),
//                outside_face_idx[N_ring_faces[ri][fi]])
//           == possible_degenerated_faces.end()){
//          possible_degenerated_faces.push_back(
//                outside_face_idx[N_ring_faces[ri][fi]]);
//        }
//      }
//    }
//  }
//  vector<size_t> other_faces;
//  find


//  return 0;
//}

int extract_surface_patch_graph(
    const matrixst & tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const matrixst & cut_tet2tet,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const jtf::mesh::edge2cell_adjacent & ea,
    const std::vector<std::pair<size_t,size_t> > & g_unknown_face_pair,
    const boost::unordered_map<size_t,size_t> & surface_type_original,
    std::vector<boost::unordered_set<size_t> > & groups,
    boost::unordered_map<size_t,boost::unordered_set<size_t> > & group_linking_info,
    boost::unordered_map<std::pair<size_t,size_t>,double> *patch_to_patch_boundary_len)
{
  // gather edges which is offered by g_unknown_faces, these edges will
  // be the obstacle between patches
  boost::unordered_set<pair<size_t,size_t> > surface_cut_edges;
  {
    for(size_t fi = 0; fi < g_unknown_face_pair.size(); ++fi){
      const vector<size_t> & face_vec =
          fa_cut.faces_.at(g_unknown_face_pair[fi].first);
      for(size_t pi = 0; pi < face_vec.size(); ++pi){
        const size_t edge_idx =
            ea.get_edge_idx(cut_tet2tet[face_vec[pi]],
                            cut_tet2tet[face_vec[(pi+1)%face_vec.size()]]);
        if(edge_idx != -1){
          pair<size_t,size_t> cut_edge(
                cut_tet2tet[face_vec[pi]],
                cut_tet2tet[face_vec[(pi+1)%face_vec.size()]]);
          if(cut_edge.first > cut_edge.second)
            swap(cut_edge.first, cut_edge.second);
          surface_cut_edges.insert(cut_edge);
        }
      }
    }
  }

  vector<bool> is_visited_face(outside_face_idx.size(), false);

  boost::unordered_map<size_t,size_t> face_idx2vec_idx;
  {
    for(size_t fi = 0; fi < outside_face_idx.size(); ++fi){
      face_idx2vec_idx[outside_face_idx[fi]] = fi;
    }
  }

  stack<size_t> face_stack;
  groups.clear();
  vector<bool>::const_iterator it =
      find(is_visited_face.begin(), is_visited_face.end(), false);

  while(it != is_visited_face.end()){
    if(face_stack.empty()){
      face_stack.push(outside_face_idx[it-is_visited_face.begin()]);
      is_visited_face.at(it-is_visited_face.begin()) = true;
    }

    boost::unordered_map<size_t,size_t>::const_iterator bumcit =
        surface_type_original.find(face_stack.top());
    if(bumcit == surface_type_original.end()){
      cerr << "# [error] can not find surface in surface type." << endl;
      return __LINE__;
    }

    const size_t group_type = bumcit->second;
    //boost::unordered_set<size_t> boundary_points;
    boost::unordered_set<size_t> one_group;
    one_group.insert(face_stack.top());
    while(!face_stack.empty()){
      const size_t current_face = face_stack.top();
      face_stack.pop();
      const vector<size_t> & face_vec = fa.faces_[current_face];
      for(size_t i = 0; i < face_vec.size(); ++i){
        pair<size_t,size_t> edge(face_vec[i], face_vec[(i+1)%face_vec.size()]);
        if(edge.first > edge.second)
          swap(edge.first, edge.second);

        if(surface_cut_edges.find(edge) != surface_cut_edges.end())
          continue;
        const pair<size_t,size_t> face_pair = ea.query(edge.first,edge.second);
        assert(outside_face_idx[face_pair.first] == current_face ||
               outside_face_idx[face_pair.second] == current_face);

        const size_t other_face_idx =
            (outside_face_idx[face_pair.first] + outside_face_idx[face_pair.second])
            - current_face;

        boost::unordered_map<size_t,size_t>::const_iterator cit =
            face_idx2vec_idx.find(other_face_idx);
        if(cit == face_idx2vec_idx.end()){
          cerr << "# [error] strange can not find face_idx " << other_face_idx
               << " to idx_vec." << endl;
          return __LINE__;
        }
        if(is_visited_face.at(cit->second) == true) {
          continue; // boundary edge
        }

        boost::unordered_map<size_t,size_t>::const_iterator type_cit =
            surface_type_original.find(other_face_idx);
        if(type_cit == surface_type_original.end()){
          cerr << "# [error] strange can not find surface type of "
               << other_face_idx << endl;
          return __LINE__;
        }
        if(type_cit->second != group_type){
          continue; // boundary edge
        }else{
          is_visited_face.at(cit->second) = true;
          face_stack.push(other_face_idx);
          one_group.insert(other_face_idx);
        }
      }
    }

    groups.push_back(one_group);
    it = find(is_visited_face.begin(), is_visited_face.end(), false);
  }


  // after grouping, next step is to build linking graph
  boost::unordered_map<size_t,size_t> face_idx2group;
  for(size_t gi = 0; gi < groups.size(); ++gi){
    const boost::unordered_set<size_t> & one_group = groups[gi];
    for(boost::unordered_set<size_t>::const_iterator cit = one_group.begin();
        cit != one_group.end(); ++cit){
      face_idx2group[*cit] = gi;
    }
  }

  group_linking_info.clear();
  for(size_t ei = 0; ei < ea.edge2cell_.size(); ++ei){
    pair<size_t,size_t> one_edge = ea.edges_[ei];
    if(one_edge.first > one_edge.second)
      swap(one_edge.first, one_edge.second);
    if(surface_cut_edges.find(one_edge) != surface_cut_edges.end()) continue;
    pair<size_t,size_t> tri_pair(
          outside_face_idx[ea.edge2cell_[ei].first],
          outside_face_idx[ea.edge2cell_[ei].second]);
    //pair<size_t,size_t> group_pair;
    boost::unordered_map<size_t,size_t>::const_iterator f2g_it_0
        = face_idx2group.find(tri_pair.first);
    boost::unordered_map<size_t,size_t>::const_iterator f2g_it_1
        = face_idx2group.find(tri_pair.second);
    if(f2g_it_0  == face_idx2vec_idx.end() ||
       f2g_it_1  == face_idx2vec_idx.end()){
      cerr << "# [error] can not find face2group." << endl;
      return __LINE__;
    }
    if(f2g_it_0->second == f2g_it_1->second) continue;
    group_linking_info[f2g_it_0->second].insert(f2g_it_1->second);
    group_linking_info[f2g_it_1->second].insert(f2g_it_0->second);

    if(patch_to_patch_boundary_len->find(
         make_pair(f2g_it_0->second, f2g_it_1->second))
       == patch_to_patch_boundary_len->end()){
      (*patch_to_patch_boundary_len)[make_pair(f2g_it_0->second, f2g_it_1->second)]
          = norm(node(colon(), one_edge.first) - node(colon(), one_edge.second));
      (*patch_to_patch_boundary_len)[make_pair(f2g_it_1->second, f2g_it_0->second)]
          = norm(node(colon(), one_edge.first) - node(colon(), one_edge.second));
    }else{
      (*patch_to_patch_boundary_len)[make_pair(f2g_it_0->second, f2g_it_1->second)]
          += norm(node(colon(), one_edge.first) - node(colon(), one_edge.second));
      (*patch_to_patch_boundary_len)[make_pair(f2g_it_1->second, f2g_it_0->second)]
          += norm(node(colon(), one_edge.first) - node(colon(), one_edge.second));
    }
  }

  { // visual
    vector<size_t> patches_face;
    vector<size_t> patches_face_type;

    for(size_t gi = 0; gi < groups.size(); ++gi){
      const boost::unordered_set<size_t> & one_group = groups[gi];
      const size_t link_num = group_linking_info[gi].size();
      for(boost::unordered_set<size_t>::const_iterator scit = one_group.begin();
          scit != one_group.end(); ++scit){
        const vector<size_t> & face_vec = fa.faces_[*scit];
        patches_face.insert(patches_face.end(), face_vec.begin(),
                            face_vec.end());
        patches_face_type.push_back(gi);
      }
    }

    ofstream ofs("pacthes_with_idx.vtk");
    tri2vtk(ofs, &node[0], node.size(2), &patches_face[0], patches_face.size()/3);
    cell_data(ofs, &patches_face_type[0], patches_face_type.size(), "patch");
  }

  return 0;
}


int gather_restricted_edges_on_polycube(
    const matrixst & tet,
    const jtf::mesh::edge2cell_adjacent & ea,
    const matrixst & outside_face_idx,
    const boost::unordered_map<size_t,size_t> & surface_type_original,
    boost::unordered_map<std::pair<size_t,size_t>,size_t > & restricted_edges)
{
  for(size_t ei = 0; ei < ea.edge2cell_.size(); ++ei){
    const pair<size_t,size_t> & tri_pair = ea.edge2cell_[ei];
    pair<size_t,size_t> two_type(-1,-1);
    boost::unordered_map<size_t,size_t>::const_iterator cit_0 =
        surface_type_original.find(outside_face_idx[tri_pair.first]);
    boost::unordered_map<size_t,size_t>::const_iterator cit_1 =
        surface_type_original.find(outside_face_idx[tri_pair.second]);
    if(cit_0 != surface_type_original.end())
      two_type.first = cit_0->second;
    if(cit_1 != surface_type_original.end())
      two_type.second = cit_1->second;
    if(two_type.first != -1 && two_type.second != -1 &&
       two_type.first != two_type.second){
      if(ea.edges_[ei].first < ea.edges_[ei].second)
        restricted_edges[ea.edges_[ei]] =
            (0+1+2) - (two_type.first + two_type.second);
      else
        restricted_edges[make_pair(ea.edges_[ei].second,
                                   ea.edges_[ei].first)] =
            (0+1+2) - (two_type.first + two_type.second);
    }
  }
  return 0;
}

int remove_one_face_degeneration(
    const matrixst &tet,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t > & restricted_edges,
    boost::unordered_map<size_t,size_t> &surface_type_original)
{
  vector<size_t> face_edge_type(3);
  for(size_t fi = 0; fi < outside_face.size(2); ++fi){
    face_edge_type.clear();
    boost::unordered_map<size_t,size_t>::iterator fcit =
        surface_type_original.find(outside_face_idx[fi]);
    if(fcit == surface_type_original.end()){
      cerr << "# [error] strange can not find surface type of "
           << outside_face_idx[fi] << endl;
      return __LINE__;
    }

    const size_t & current_type = fcit->second;

    for(size_t pi = 0; pi < outside_face.size(1) ;++pi){
      pair<size_t,size_t> edge(outside_face(pi,fi),
                               outside_face((pi+1)%outside_face.size(1),fi));
      if(edge.first > edge.second)
        swap(edge.first, edge.second);

      boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit
          = restricted_edges.find(edge);
      if(cit == restricted_edges.end())
        cit = restricted_edges.find(make_pair(edge.second, edge.first));
      if(cit == restricted_edges.end())
        continue;
      face_edge_type.push_back(cit->second);
    }
    if(face_edge_type.size() == 3){ // all edges are restricted
      // since face type can only be 0/1/2, three edges of this face can only
      // be with two/one type
      boost::unordered_set<size_t> edge_type(face_edge_type.begin(), face_edge_type.end());
      if(edge_type.size() == 1){
        assert(*edge_type.begin() != current_type); // all edge type is 0, this face must be 1 or 2
        fcit->second = (0+1+2) - *edge_type.begin() - current_type;
      }else if(edge_type.size() == 2){
        // denote the three edge type as 0,0,1,
        // then the arounding face type of such face must be 1,1,0,
        // while the current type must be 2, so this face should be modified to 1
        if(count(face_edge_type.begin(),face_edge_type.end(),
                 *edge_type.begin()) == 1)
          fcit->second = *edge_type.begin();
        else{
          boost::unordered_set<size_t>::const_iterator it = edge_type.begin();
          ++it;
          fcit->second = *it;
        }
      }
      ++degree_issue_num;
    }
  }

  return 0;
}

static int remove_degree_two_patch_by_erasing_or_reconnecting(
    const size_t & g,
    type_patch_graph &tpg)
{
  if(tpg.get_patch_degree(g) != 2){
    cerr << "# [error] not degree two patch." << endl;
    return __LINE__;
  }

  const boost::unordered_set<size_t> & linking = tpg.get_linking_of_group(g);
  boost::unordered_set<size_t> two_ring_patches;

  double erase_cost = 0.0, reconnect_cost = 0.0;
  for(boost::unordered_set<size_t>::const_iterator cit = linking.begin();
      cit != linking.end(); ++cit){
    const size_t degree_other = tpg.get_patch_degree(*cit);
    if(degree_other < 5)
      erase_cost = std::numeric_limits<double>::infinity();
    if(degree_other == 2){
      cerr << "# [error] degree 2 patch linking degree 2 patch, need "
           << "special addressing." << endl;
      return __LINE__;
    }
    const boost::unordered_set<size_t> & one_ring_to_cit =
        tpg.get_linking_of_group(*cit);
    two_ring_patches.insert(one_ring_to_cit.begin(), one_ring_to_cit.end());
  }

  // if erase_cost == infinity, this patch must not been erased, or it will
  // introduce new degeneration.
  // if erase_cost == 0, this patch can be erased or reconnect by extend current patch
  // or inserting a new one.

  {
    // remove one_ring_patches from two_ring_patches, this step is to ensure
    // two_ring_patches contains pure two_rings
    for(boost::unordered_set<size_t>::const_iterator cit = linking.begin();
        cit != linking.end(); ++cit){
      boost::unordered_set<size_t>::iterator it_two_ring = two_ring_patches.find(*cit);
      if(it_two_ring != two_ring_patches.end())
        two_ring_patches.erase(it_two_ring);
    }
  }

  boost::unordered_set<size_t> degree_two_of_two_ring_patchs;
  for(boost::unordered_set<size_t>::const_iterator cit = two_ring_patches.begin();
      cit != two_ring_patches.end(); ++cit){
    const size_t degree = tpg.get_patch_degree(*cit);
    if(degree == 2)
      degree_two_of_two_ring_patchs.insert(*cit);
  }



  return 0;
}

int remove_degenerated_patch_and_break_through(
    const matrixst &tet,
    const matrixd &node,
    const matrixst &cut_tet,
    const matrixst &cut_tet2tet,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const jtf::mesh::edge2cell_adjacent & ea,
    const vector<pair<size_t,size_t> > & g_unknown_face_pair,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const matrixd & polycube_node,
    boost::unordered_map<size_t,size_t> &surface_type)
{
  cerr << "# [error] not finished." << endl;
  return __LINE__;

  if(node.size() != polycube_node.size()){
    cerr << "# [error] uncompatiable node size: " << node.size() << " "
         << polycube_node.size() << endl;
    return __LINE__;
  }

  type_patch_graph tpg(surface_type);
  tpg.build_graph(
        tet, cut_tet, node, cut_tet2tet, outside_face, outside_face_idx, fa,
        fa_cut, ea, g_unknown_face_pair);

  boost::unordered_set<size_t> degree_two_patches;

  while(1){
    degree_two_patches.clear();
    const type_patch_graph::linking_type & lt = tpg.get_linking_info();
    for(type_patch_graph::linking_type::const_iterator
        cit = lt.begin(); cit != lt.end(); ++cit){
      if(cit->second.size() == 2)
        degree_two_patches.insert(cit->first);
    }
    const size_t degree_two_number = degree_two_patches.size();
    if(degree_two_number == 0) return 0; // I can not handle degree 3 patch now

    for(boost::unordered_set<size_t>::const_iterator cit
        = degree_two_patches.begin(); cit != degree_two_patches.end(); ++cit){
      remove_degree_two_patch_by_erasing_or_reconnecting(*cit, tpg);
    }
  }

  return 0;
}


int remove_surface_zigzag_by_restricted_edges_from_graph(
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face_idx,
   jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::edge2cell_adjacent &ea,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &restricted_edge_from_graph_orig,
    boost::unordered_map<size_t,size_t> &surface_type_original,
    matrixst &new_tet,
    matrixd &new_node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    matrixd * zyz_ptr)
{
  sxx::tet_mesh stm;
  stm.create_tetmesh(node, tet);

  for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit
      = restricted_edge_from_graph_orig.begin();
      cit != restricted_edge_from_graph_orig.end(); ++cit){
    const pair<size_t,size_t> & edge = cit->first;
    const size_t edge_idx = ea.get_edge_idx(edge.first, edge.second);
    if(edge_idx == -1){
      continue; // since it's an inner edge
      //      cerr << "# [error] can not find edge " << edge.first << " " << edge.second << endl;
      //      return __LINE__;
    }
    const pair<size_t,size_t> & tri_pair = ea.edge2cell_[edge_idx];
    boost::unordered_map<size_t,size_t>::const_iterator f0 =
        surface_type_original.find(outside_face_idx[tri_pair.first]);
    boost::unordered_map<size_t,size_t>::const_iterator f1 =
        surface_type_original.find(outside_face_idx[tri_pair.second]);
    if(f0 == surface_type_original.end() || f1 == surface_type_original.end())
    {
      cerr << "# [error] strange can not find face type." << endl;
      return __LINE__;
    }

    if(f0->second == f1->second){
      stm.split_edge(edge);
    }
  }

  stm.write_tetmesh_to_matrix(new_node, new_tet);

  unique_ptr<jtf::mesh::face2tet_adjacent> fa_new(jtf::mesh::face2tet_adjacent::create(new_tet));
  if(!fa_new.get()){
    cerr << "# [error] can not build face2tet_adjacnet." << endl;
    return __LINE__;
  }

  boost::unordered_map<vector<size_t>, vector<size_t> > new_face2orig_face_mapping;

  sxx::get_face2orginal(stm, new_face2orig_face_mapping);

  matrixst outside_face_new, outside_face_idx_new;
  get_outside_face(*fa_new, outside_face_new);
  get_outside_face_idx(*fa_new, outside_face_idx_new);

  update_surface_type(fa, *fa_new, new_face2orig_face_mapping, surface_type_original);

  if(!inner_face_jump_type.empty()){
    matrixst tet_idx_map;
    sxx::get_tet2orginal_index(stm, new_tet, tet_idx_map);
    update_inner_face_jump_type(*fa_new, tet_idx_map, inner_face_jump_type);
  }

  if(zyz_ptr){
    zjucad::matrix::matrix<size_t> tet_index_map;
    sxx::get_tet2orginal_index(stm,new_tet, tet_index_map);
    matrixd new_zyz(3, new_tet.size(2)) ;
    for(size_t ti = 0; ti < tet_index_map.size(); ++ti){
      new_zyz(colon(), ti) = (*zyz_ptr)(colon(), tet_index_map[ti]);
    }
    *zyz_ptr = new_zyz;
  }

  fa = *fa_new;
  return 0;
}

int remove_surface_critcial_points(
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    const matrixst &outside_face_idx,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::edge2cell_adjacent &ea,
    boost::unordered_map<size_t,size_t> & surface_type)
{
  return __LINE__;
  //  jtf::mesh::one_ring_face_at_point orfap;
  //  orfap.add_all_faces(outside_face);
  //  orfap.sort_int_loop(outside_face, node);

  //  for(jtf::mesh::one_ring_face_at_point::p2f_type::const_iterator pcit =
  //      orfap.p2f_.begin(); pcit != orfap.p2f_.end(); ++pcit){
  //    const size_t & point_idx = pcit->first;
  //    const vector<size_t> & around_faces = pcit->second;
  //    if(is_critical_point(around_faces,outside_face_idx,surface_type)){

  //    }
  //  }

  return 0;
}

int remove_surface_zigzag_by_restricted_edges(
    const matrixst &tet,
    const matrixd & node,
    const matrixst & outside_face_idx,
    const jtf::mesh::edge2cell_adjacent & ea,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &restricted_edges,
    const jtf::mesh::face2tet_adjacent &fa,
    matrixst &new_tet,
    matrixd &new_node,
    boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inne_face_jump_type,
    matrixd * zyz_ptr)
{
  boost::unordered_set<pair<size_t,size_t> > edges_need_to_split;

  for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator
      cit = restricted_edges.begin(); cit != restricted_edges.end(); ++cit){
    boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator citt = cit;
    ++citt;
    const pair<size_t,size_t> & current_edge = cit->first;
    for(; citt != restricted_edges.end(); ++citt){
      const pair<size_t,size_t> & next_edge = citt->first;

      if(next_edge.first == current_edge.first ||
         next_edge.first == current_edge.second ||
         next_edge.second == current_edge.first ||
         next_edge.second == current_edge.second){ // at least they have one shared point
        deque<size_t> three_points;
        three_points.push_back(current_edge.first);
        three_points.push_back(current_edge.second);
        if(next_edge.first == three_points.back())
          three_points.push_back(next_edge.second);
        else if(next_edge.second == three_points.back())
          three_points.push_back(next_edge.first);
        else if(next_edge.first == three_points.front())
          three_points.push_front(next_edge.second);
        else if(next_edge.second == three_points.front())
          three_points.push_front(next_edge.first);

        const size_t face_idx =
            fa.get_face_idx(three_points[0], three_points[1], three_points[2]);
        if(face_idx != -1){
          edges_need_to_split.insert(make_pair(three_points.front(),
                                               three_points.back()));
        }
      }
    }
  }

  sxx::tet_mesh stm;
  stm.create_tetmesh(node, tet);
  vector<std::tuple<size_t,size_t,size_t> > split_edge_info;
  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit
      = edges_need_to_split.begin(); cit != edges_need_to_split.end(); ++cit){
    //cerr << cit->first << " " << cit->second << endl;
    const size_t new_vertex_ind = stm.split_edge(*cit);
    if(new_vertex_ind == -1) continue;
    split_edge_info.push_back(
          std::make_tuple( new_vertex_ind, cit->first, cit->second));
  }

  {
    vector<size_t> edges_need_to_split_vec;
    for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit
        = edges_need_to_split.begin(); cit != edges_need_to_split.end(); ++cit){
      edges_need_to_split_vec.push_back(cit->first);
      edges_need_to_split_vec.push_back(cit->second);
    }

    ofstream ofs("edges_need_to_split.vtk");
    line2vtk(ofs, &node[0], node.size(2), &edges_need_to_split_vec[0],
             edges_need_to_split_vec.size()/2);
  }
  stm.write_tetmesh_to_matrix(new_node, new_tet);

  unique_ptr<jtf::mesh::face2tet_adjacent> fa_new(jtf::mesh::face2tet_adjacent::create(new_tet));
  if(!fa_new.get()){
    cerr << "# [error] can not build face2tet_adjacnet." << endl;
    return __LINE__;
  }

  boost::unordered_map<vector<size_t>, vector<size_t> > new_face2orig_face_mapping;

  sxx::get_face2orginal(stm, new_face2orig_face_mapping);

  matrixst outside_face_new, outside_face_idx_new;
  get_outside_face(*fa_new, outside_face_new);
  get_outside_face_idx(*fa_new, outside_face_idx_new);

  update_surface_type(fa, *fa_new, new_face2orig_face_mapping, surface_type);

  if(!inne_face_jump_type.empty()){
    matrixst tet_idx_map;
    sxx::get_tet2orginal_index(stm, new_tet, tet_idx_map);
    update_inner_face_jump_type(*fa_new, tet_idx_map, inne_face_jump_type);
  }
  //  boost::unordered_map<size_t,size_t> new_surface_type;
  //  vector<size_t> one_face(outside_face_new.size(1));
  //  for(size_t fi = 0; fi < outside_face_new.size(2); ++fi){
  //    copy(outside_face_new(colon(),fi).begin(),
  //         outside_face_new(colon(), fi).end(),  one_face.begin());
  //    sort(one_face.begin(), one_face.end());
  //    boost::unordered_map<vector<size_t>,vector<size_t> >::const_iterator cit
  //        = new_face2orig_face_mapping.find(one_face);
  //    if(cit == new_face2orig_face_mapping.end()){
  //      cerr << "# [error] strange can not find face mapping. " << endl;
  //      return __LINE__;
  //    }
  //    const vector<size_t> & orig_face = cit->second;
  //    const size_t & orig_face_idx = fa.get_face_idx(&orig_face[0]);
  //    if(orig_face_idx == -1){
  //      cerr << "# [error] strange can not find original face idx."  << endl;
  //      return __LINE__;
  //    }
  //    boost::unordered_map<size_t,size_t>::const_iterator bumcit =
  //        surface_type.find(orig_face_idx);
  //    if(bumcit == surface_type.end()){
  //      cerr << "# [error] can not find surface type of orig face "
  //           << orig_face_idx << endl;
  //      return __LINE__;
  //    }
  //    new_surface_type[outside_face_idx_new[fi]] = bumcit->second;
  //  }

  //  surface_type = new_surface_type;

  if(zyz_ptr){
    zjucad::matrix::matrix<size_t> tet_index_map;
    sxx::get_tet2orginal_index(stm,new_tet, tet_index_map);
    matrixd new_zyz(3, new_tet.size(2)) ;
    for(size_t ti = 0; ti < tet_index_map.size(); ++ti){
      new_zyz(colon(), ti) = (*zyz_ptr)(colon(), tet_index_map[ti]);
    }
    *zyz_ptr = new_zyz;
  }

  return 0;
}

int extract_relax_surface(
    const matrixst & tet,
    const matrixd & node,
    const boost::unordered_map<size_t,size_t> &surface_type,
    const boost::unordered_set<size_t> &loop_points,
    boost::unordered_map<size_t,size_t> &restricted_surface_type)
{
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  matrixst outside_face, outside_face_idx;
  get_outside_face(*fa, outside_face);
  get_outside_face_idx(*fa, outside_face_idx);

  if(surface_type.size() != outside_face.size(2)){
    cerr << "# [error] imcompatiable surface type setting." << endl;
    return __LINE__;
  }
  // first, find the the faces which will is 2-ring near the degenerated points,
  // relax their surface restricted types

  boost::unordered_set<size_t> relaxed_face_idx;
  jtf::mesh::N_ring_face_at_point nrfap;
  const size_t N_ring = 2;
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }


  nrfap.add_all_faces(outside_face, *ea, N_ring);
  for(boost::unordered_set<size_t>::const_iterator cit = loop_points.begin();
      cit != loop_points.end(); ++cit){
    const size_t & point_idx = *cit;
    jtf::mesh::N_ring_face_at_point::p2f_type::const_iterator p2cit =
        nrfap.p2f_.find(point_idx);
    if(p2cit == nrfap.p2f_.end()){
      cerr << "# [error] strange can not find point " << point_idx
           << " on surface." << endl;
      return __LINE__;
    }
    const vector<vector<size_t> > & N_ring_faces = p2cit->second;
    for(size_t ri = 0; ri < N_ring_faces.size(); ++ri){
      for(size_t fi = 0; fi < N_ring_faces[ri].size(); ++fi){
        relaxed_face_idx.insert( outside_face_idx[N_ring_faces[ri][fi]]);
      }
    }
  }

  { //
    ofstream ofs("relaxed_surface.vtk");
    vector<size_t> face_vec;
    for(boost::unordered_set<size_t>::const_iterator cit = relaxed_face_idx.begin();
        cit != relaxed_face_idx.end(); ++cit){
      face_vec.insert(face_vec.end(), fa->faces_[*cit].begin(),
                      fa->faces_[*cit].end());
    }
    tri2vtk(ofs, &node[0], node.size(2), &face_vec[0], face_vec.size()/3);
  }

  restricted_surface_type.clear();
  for(boost::unordered_map<size_t,size_t>::const_iterator cit = surface_type.begin();
      cit != surface_type.end(); ++cit){
    // is not inside the relaxed face set
    if(relaxed_face_idx.find(cit->first) == relaxed_face_idx.end()){
      restricted_surface_type[cit->first] = cit->second;
    }
  }

  return 0;
}

int update_surface_type(
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_new,
    const boost::unordered_map<vector<size_t>,vector<size_t> > & f2omap,
    boost::unordered_map<size_t,size_t> & surface_type)
{
  boost::unordered_map<size_t,size_t> new_surface_type;

  for(boost::unordered_map<vector<size_t>,vector<size_t> >::const_iterator
      cit = f2omap.begin(); cit != f2omap.end(); ++cit){
    const vector<size_t> & new_face = cit->first;
    const vector<size_t> & old_face = cit->second;
    const size_t & new_face_idx = fa_new.get_face_idx(&new_face[0]);
    const size_t & old_face_idx = fa.get_face_idx(&old_face[0]);

    const pair<size_t,size_t> & tet_pair = fa.face2tet_[old_face_idx];

    if(new_face_idx == -1 || old_face_idx == -1){
      cerr << "# [error] can not find face idx." << endl;
      return __LINE__;
    }

    if(!fa.is_outside_face(tet_pair)) continue;

    boost::unordered_map<size_t,size_t>::const_iterator ccit =
        surface_type.find(old_face_idx);
    if(ccit == surface_type.end()){
      cerr << "# [error] strange can not find face type of "
           << old_face_idx << " in original surface type." << endl;
      return __LINE__;
    }
    new_surface_type[new_face_idx] = ccit->second;
  }

  surface_type = new_surface_type;
  return 0;
}

int update_inner_face_jump_type(
    const jtf::mesh::face2tet_adjacent &fa_new,
    const matrixst & tet_idx_map,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type)
{
  boost::unordered_map<pair<size_t,size_t>,size_t> new_tet_pair_type;

  for(size_t fi = 0; fi < fa_new.face2tet_.size(); ++fi){
    const pair<size_t,size_t> & tet_pair = fa_new.face2tet_[fi];
    if(fa_new.is_outside_face(tet_pair)) continue;
    const pair<size_t,size_t> orig_tet_pair(
          tet_idx_map[tet_pair.first], tet_idx_map[tet_pair.second]);
    boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit
        = inner_face_jump_type.find(orig_tet_pair);
    if(cit == inner_face_jump_type.end()) continue;
    new_tet_pair_type[tet_pair] = cit->second;
    new_tet_pair_type[
        make_pair(tet_pair.second, tet_pair.first)] =
        get_trans_type(cit->second);
  }

  inner_face_jump_type = new_tet_pair_type;

  return 0;
}

int collapse_degenerated_edges(
    matrixst & tet,
    matrixd & node,
   jtf::mesh::face2tet_adjacent & fa,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & restricted_edges)
{
  for(boost::unordered_map<pair<size_t,size_t>,size_t>::iterator it =
      inner_face_jump_type.begin(); it != inner_face_jump_type.end(); ){
    if(is_trivial_type(it->second)) inner_face_jump_type.erase(it++);
    else
      ++it;
  }

  boost::unordered_set<pair<size_t,size_t> > edges_need_to_collapse;
  for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit =
      restricted_edges.begin(); cit != restricted_edges.end(); ++cit){
    if(is_black_line_new(cit->second))
      edges_need_to_collapse.insert(cit->first);
  }

  if(edges_need_to_collapse.empty())
    return 1;

  sxx::tet_mesh stm;
  stm.create_tetmesh(node, tet);

  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit =
      edges_need_to_collapse.begin(); cit != edges_need_to_collapse.end(); ++cit){
    //cerr << cit->first << " " << cit->second << endl;
    stm.collapse_edge(*cit);
    //cerr << "# [error] edge delete is not finished." << endl;
  }

  stm.write_tetmesh_to_matrix(node, tet);

  boost::unordered_map<vector<size_t>,vector<size_t> > f2omap;
  matrixst tet_idx_map;
  sxx::get_face2orginal(stm, f2omap);
  sxx::get_tet2orginal_index(stm, tet, tet_idx_map);

  unique_ptr<jtf::mesh::face2tet_adjacent> fa_new(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa_new.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  boost::unordered_map<pair<size_t,size_t>,size_t> new_tet_pair_type;

  if(!inner_face_jump_type.empty()){
    update_inner_face_jump_type(*fa_new, tet_idx_map, inner_face_jump_type)  ;
  }

  update_surface_type(fa,*fa_new, f2omap, surface_type);

  fa = *fa_new;
  return 0;
}

int straighten_patch_boundary(
    const matrixst &outside_face,
    const matrixst &outside_face_idx,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::face2tet_adjacent &fa,
    boost::unordered_map<size_t,size_t> &surface_type_original,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & rot_type)
{
  vector<vector<boost::unordered_map<size_t,size_t>::iterator> > face2adj(outside_face_idx.size());

  for(size_t fi = 0; fi < outside_face_idx.size(); ++fi){
    boost::unordered_map<size_t,size_t>::iterator it =
        surface_type_original.find(outside_face_idx[fi]);
    const pair<size_t,size_t> & this_tet_pair = fa.face2tet_[outside_face_idx[fi]];
    assert(this_tet_pair.first == -1 || this_tet_pair.second == -1);
    const size_t this_tet_idx =
        (this_tet_pair.first == -1?this_tet_pair.second:this_tet_pair.first);

    if(it == surface_type_original.end()){
      cerr << "# [error] strange can not find surface type of face "
           << outside_face_idx[fi] << endl;
      return __LINE__;
    }
    face2adj[fi].push_back(it);

    for(size_t pi = 0; pi < outside_face.size(1); ++pi){
      const size_t edge_idx =
          ea.get_edge_idx(outside_face(pi, fi),
                          outside_face((pi+1)%outside_face.size(1),fi));
      if(edge_idx == -1) {
        cerr << "# [error] strange can not find edge "
             << outside_face(pi, fi)
             << outside_face((pi+1)%outside_face.size(1),fi) << endl;
        return __LINE__;
      }

      const pair<size_t,size_t> & tri_pair = ea.edge2cell_[edge_idx];
      assert(tri_pair.first == fi || tri_pair.second == fi);

      const size_t & other_tri = tri_pair.first + tri_pair.second - fi;
      const size_t & other_tri_idx = outside_face_idx[other_tri];

      const pair<size_t,size_t> & other_tet_pair = fa.face2tet_[other_tri_idx];
      assert(other_tet_pair.first == -1 || other_tet_pair.second == -1);

      const size_t other_tet_idx =
          (other_tet_pair.first == -1?other_tet_pair.second:other_tet_pair.first);
      boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit =
          rot_type.find(make_pair(this_tet_idx, other_tet_idx));

      if(cit != rot_type.end() && !is_trivial_type(cit->second)) continue;
      boost::unordered_map<size_t,size_t>::iterator it_other =
          surface_type_original.find(other_tri_idx);
      if(it_other == surface_type_original.end()){
        cerr << "# [error] strange can not find surface type of face "
             << other_tri_idx << endl;
        return __LINE__;
      }
      face2adj[fi].push_back(it_other);
    }
    //assert(face2adj[fi].size() == 4); // on closed surface , each triangle connect four faces
  }

  vector<size_t> uvw_type(3,0);

  while(1){
    size_t number = 0;
    for(size_t fi = 0; fi < face2adj.size(); ++fi){
      fill(uvw_type.begin(), uvw_type.end(),0);
      size_t &current_face_type = face2adj[fi].front()->second;
      for(size_t ci = 1; ci < face2adj[fi].size(); ++ci){
        ++uvw_type[face2adj[fi][ci]->second];
      }
      if(uvw_type[0] == uvw_type[1] && uvw_type[1] == uvw_type[2])
        continue;
      if(std::accumulate(uvw_type.begin(), uvw_type.end(),0.0) < 3)
        //means a face beside seam, do not have enough information to determin
        continue;

      const size_t max_type = max_element(uvw_type.begin(), uvw_type.end()) - uvw_type.begin();
      if(max_type == current_face_type) continue;
      current_face_type = max_type;
      ++number;
    }
    if(number == 0) break;
  }
  return 0;
}

int straighten_patch_boundary(
    const matrixst &outside_face,
    const matrixst &outside_face_idx,
    const jtf::mesh::edge2cell_adjacent &ea,
    boost::unordered_map<size_t,size_t> &surface_type_original)
{
  vector<vector<boost::unordered_map<size_t,size_t>::iterator> > face2adj(outside_face_idx.size());

  for(size_t fi = 0; fi < outside_face_idx.size(); ++fi){
    boost::unordered_map<size_t,size_t>::iterator it =
        surface_type_original.find(outside_face_idx[fi]);
    if(it == surface_type_original.end()){
      cerr << "# [error] strange can not find surface type of face "
           << outside_face_idx[fi] << endl;
      return __LINE__;
    }
    face2adj[fi].push_back(it);

    for(size_t pi = 0; pi < outside_face.size(1); ++pi){
      const size_t edge_idx =
          ea.get_edge_idx(outside_face(pi, fi),
                          outside_face((pi+1)%outside_face.size(1),fi));
      if(edge_idx == -1) {
        cerr << "# [error] strange can not find edge "
             << outside_face(pi, fi)
             << outside_face((pi+1)%outside_face.size(1),fi) << endl;
        return __LINE__;
      }
      const pair<size_t,size_t> & tri_pair = ea.edge2cell_[edge_idx];
      assert(tri_pair.first == fi || tri_pair.second == fi);
      const size_t & other_tri = tri_pair.first + tri_pair.second - fi;
      const size_t & other_tri_idx = outside_face_idx[other_tri];
      boost::unordered_map<size_t,size_t>::iterator it_other =
          surface_type_original.find(other_tri_idx);
      if(it_other == surface_type_original.end()){
        cerr << "# [error] strange can not find surface type of face "
             << other_tri_idx << endl;
        return __LINE__;
      }
      face2adj[fi].push_back(it_other);
    }
    assert(face2adj[fi].size() == 4); // on closed surface , each triangle connect four faces
  }

  vector<size_t> uvw_type(3,0);

  while(1){
    size_t number = 0;
    for(size_t fi = 0; fi < face2adj.size(); ++fi){
      fill(uvw_type.begin(), uvw_type.end(),0);
      size_t &current_face_type = face2adj[fi].front()->second;
      for(size_t ci = 1; ci < face2adj[fi].size(); ++ci){
        ++uvw_type[face2adj[fi][ci]->second];
      }
      if(uvw_type[0] == uvw_type[1] && uvw_type[1] == uvw_type[2])
        continue;
      const size_t max_type = max_element(uvw_type.begin(), uvw_type.end()) - uvw_type.begin();
      if(max_type == current_face_type) continue;
      current_face_type = max_type;
      ++number;
    }
    boundary_issue_num += number;
    if(number == 0) break;
  }
  return 0;
}

int remove_dgenerated_patch_by_modify_type(
    const matrixst &tet,
    const matrixst & cut_tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent & fa,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const jtf::mesh::edge2cell_adjacent & ea,
    const std::vector<std::pair<size_t,size_t> > & g_unknown_face_pair,
    boost::unordered_map<size_t,size_t> &surface_type)
{
  // first step is to detect the patch whose boundary is loop points, and
  // then we modify their types to merge them to their neighboring patches.

  //  matrixst outside_face, outside_face_idx;
  //  get_outside_face(fa, outside_face);
  //  get_outside_face_idx(fa, outside_face_idx);

  //  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
  //        jtf::mesh::edge2cell_adjacent::create(outside_face));
  //  if(!ea.get()){
  //    cerr << "# [error] strange can not build edge2cell_adjacent." << endl;
  //    return __LINE__;
  //  }

  vector<boost::unordered_set<size_t> > groups;
  boost::unordered_map<size_t, boost::unordered_set<size_t> > group_linking_info;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa_cut(jtf::mesh::face2tet_adjacent::create(cut_tet));
  if(!fa_cut.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  matrixst cut_tet2tet_mapping(max(cut_tet)+1);
  cut_tet2tet_mapping(cut_tet) = tet(colon());
  boost::unordered_map<pair<size_t,size_t>,double> p2p_len;
  extract_surface_patch_graph(tet, cut_tet, node, cut_tet2tet_mapping,
                              outside_face, outside_face_idx, fa, *fa_cut, ea,
                              g_unknown_face_pair, surface_type, groups,
                              group_linking_info, &p2p_len);

  cerr << "# [info] finish extract surface_patch_graph." << endl;

  remove_isolated_patch(tet, node, fa, groups, group_linking_info,
                        surface_type);

  cerr << "# [info] finish remove ioslated patch." << endl;

  remove_fractional_degree_two_patches(
        tet, node, fa, groups, group_linking_info, surface_type, p2p_len);

  cerr << "# [info] finish remove fractional degree two patches." << endl;

  return 0;
}

