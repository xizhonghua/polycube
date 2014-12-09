#include <stack>
#include <fstream>
#include <sstream>
#include <numeric>

#include <verdict/verdict.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <jtflib/util/util.h>
#include <hjlib/math/polar.h>
#include <hjlib/math/blas_lapack.h>

#include "../mesh_func/tri-normal.h"
#include "util.h"
#include "tet_func.h"
#include "../common/vtk.h"
#include "../tetmesh/util.h"
#include "../common/util.h"
#include "../tetmesh/tetmesh.h"
#include "../tet_mesh_sxx/tet_mesh_sxx.h"
#include "adaptive.h"
#include "../tetmesh/subdivide_tet.h"

using namespace std;
using namespace zjucad::matrix;

bool is_degenerate(const zjucad::matrix::matrix<double> &tri)
{
    matrix<double> e1 = tri(colon(), 0)-tri(colon(), 2),
            e2 = tri(colon(), 1)-tri(colon(), 2);
    if(norm(cross(e1, e2)) < 1e-8) {
        return true;
    }
    return false;
}

int extra_surface_patch_according_to_type(
        const boost::unordered_map<size_t,size_t> &surface_type,
        const matrix<size_t> & outside_face_idx,
        const jtf::mesh::edge2cell_adjacent & ea,
        const jtf::mesh::face2tet_adjacent & fa,
        vector<boost::unordered_set<size_t> >  &surface_patches,
        boost::unordered_set<pair<size_t,size_t> > &patch_boundary)
{
    vector<bool> is_visited_face(outside_face_idx.size(), false);

    boost::unordered_map<size_t,size_t> face_idx2vec_idx;
    for(size_t fi = 0; fi < outside_face_idx.size(); ++fi){
        face_idx2vec_idx[outside_face_idx[fi]] = fi;
    }

    stack<size_t> face_stack;
    surface_patches.clear();
    vector<bool>::const_iterator it =
            find(is_visited_face.begin(), is_visited_face.end(), false);

    while(it != is_visited_face.end()){
        if(face_stack.empty()){
            face_stack.push(outside_face_idx[it-is_visited_face.begin()]);
            is_visited_face.at(it-is_visited_face.begin()) = true;
        }

        boost::unordered_map<size_t,size_t>::const_iterator bumcit =
                surface_type.find(face_stack.top());
        if(bumcit == surface_type.end()){
            cerr << "# [error] can not find surface in surface type." << endl;
            return __LINE__;
        }

        const size_t group_type = bumcit->second;

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
                        surface_type.find(other_face_idx);
                if(type_cit == surface_type.end()){
                    cerr << "# [error] strange can not find surface type of "
                         << other_face_idx << endl;
                    return __LINE__;
                }
                if(type_cit->second != group_type){
                    patch_boundary.insert(edge);
                    continue; // boundary edge
                }else{
                    is_visited_face.at(cit->second) = true;
                    face_stack.push(other_face_idx);
                    one_group.insert(other_face_idx);
                }
            }
        }

        surface_patches.push_back(one_group);
        it = find(is_visited_face.begin(), is_visited_face.end(), false);
    }
    return 0;
}

int sort_edge(pair<size_t,size_t> & edge)
{
    if(edge.first > edge.second)
        swap(edge.first, edge.second);
    return 0;
}

int collapse_tet_edge(matrix<size_t> & tet,
                      matrix<double> & orig_node,
                      matrix<double> & deform_node,
                      const boost::unordered_set<pair<size_t,size_t> > & edges_need_to_collapse)
{
    sxx::tet_mesh stm, stm_defor;
    stm.create_tetmesh(orig_node, tet);
    stm_defor.create_tetmesh(deform_node, tet);

    //for(size_t ei = 0; ei < edges_need_to_collapse.size(); ++ei){
    for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit =
        edges_need_to_collapse.begin(); cit != edges_need_to_collapse.end(); ++cit){
        stm.collapse_edge(*cit);
        stm_defor.collapse_edge(*cit);
    }

    stm.write_tetmesh_to_matrix(orig_node, tet);
    stm_defor.write_tetmesh_to_matrix(deform_node, tet);

    return 0;
}

// return dot(p1-p0,p2-p0)/(norm(p1-p0, p2-p0))
double get_cos_angle(const matrix<double> & node,
                     const size_t p0, const size_t p1, const size_t p2)
{
    matrix<double> e1 = node(colon(), p1) - node(colon(),p0);
    matrix<double> e2 = node(colon(), p2) - node(colon(),p0);

    double len1 = norm(e1);
    if(len1 < 1e-6) len1 = 1.0;

    double len2 = norm(e2);
    if(len2 < 1e-6) len2 = 1.0;

    return dot(e1,e2)/(len1 * len2);
}

int split_tet_at_large_arap_distortion(
        matrix<size_t> &orig_tet,
        matrix<double> &orig_node,
        matrix<double> &deform_node,
        matrix<size_t> &new_node_parent_id,
        const double split_percent)
{
    vector<pair<double, size_t> > arap_distortion_with_index;
    arap_distortion_with_index.reserve(orig_tet.size(2));
    matrix<double> orig_tet_node, deform_tet_node;
    matrix<double> arap_dis = zeros<double>(3,3);
    for(size_t ti = 0; ti < orig_tet.size(2); ++ti){
        orig_tet_node = orig_node(colon(), orig_tet(colon(),ti));
        deform_tet_node = deform_node(colon(), orig_tet(colon(),ti));
        calc_tet_arap_distortion(&orig_tet_node[0], &deform_tet_node[0], &arap_dis[0]);
        arap_distortion_with_index.push_back(make_pair(norm(arap_dis), ti));
    }
    sort(arap_distortion_with_index.begin(), arap_distortion_with_index.end());

    // rules: find the largest 5%100 tets, and for each tet, find the largest deformed edge, split.

    size_t begin_idx = arap_distortion_with_index.size() -
            arap_distortion_with_index.size() * (split_percent / 100.0);
    //  size_t begin_idx = -1;

    //  for(size_t ti = 0; ti < arap_distortion_with_index.size(); ++ti)
    //    if(arap_distortion_with_index[ti].first > 1.0) {
    //      begin_idx = ti;
    //      break;
    //    }
    //  if(begin_idx == -1) return 0;

    const size_t edge2idx_[] = {0,1,0,2,0,3,1,2,1,3,2,3};
    vector<pair<double, size_t> > edge_deform_idx(6);
    const itr_matrix<const size_t*> edge2idx_mat(2,6, &edge2idx_[0]);
    vector<size_t> edges_need_to_split;
    for(size_t i = begin_idx; i < arap_distortion_with_index.size(); ++i){
        const size_t &tet_idx = arap_distortion_with_index[i].second;

        for(size_t ei = 0; ei < edge2idx_mat.size(2); ++ei){
            const double deform_len =
                    norm(deform_node(colon(), orig_tet(edge2idx_mat(0,ei), tet_idx))
                         - deform_node(colon(), orig_tet(edge2idx_mat(1,ei), tet_idx)));
            const double orig_len =
                    norm(orig_node(colon(), orig_tet(edge2idx_mat(0,ei), tet_idx))
                         - orig_node(colon(), orig_tet(edge2idx_mat(1,ei), tet_idx)));
            edge_deform_idx[ei] = make_pair(deform_len/orig_len, ei);
        }

        sort(edge_deform_idx.begin(), edge_deform_idx.end());
        //if(edge_deform_idx.back().first < 5) continue; // if the scaled_len/original_len < 3, it seems not such terrible
        const size_t &large_deform_edge_idx = edge_deform_idx.back().second;
        // to determine whether this edge should be split. we should calculate the opposite angle
        // if two opposite angle assoicates with this edge > 120, we split it or, ignore it
        const pair<size_t,size_t> edge(orig_tet(edge2idx_mat(0,large_deform_edge_idx), i),
                                       orig_tet(edge2idx_mat(1,large_deform_edge_idx), i));
        vector<size_t> other_edge;
        for(size_t p = 0; p < orig_tet.size(1); ++p){
            if(orig_tet(p,i) != edge.first && orig_tet(p,i) != edge.second)
                other_edge.push_back(orig_tet(p,i));
        }
        assert(other_edge.size() == 2);

        const double angle0 = get_cos_angle(deform_node, other_edge.front(), edge.first, edge.second);
        const double angle1 = get_cos_angle(deform_node, other_edge.back(), edge.first, edge.second);

        if(angle0 < -0.5 && angle1 < -0.5){ // two angels are larger than 120
            edges_need_to_split.push_back(orig_tet(edge2idx_mat(0,large_deform_edge_idx),tet_idx));
            edges_need_to_split.push_back(orig_tet(edge2idx_mat(1,large_deform_edge_idx),tet_idx));
        }
    }

    cerr << "# [info] edges need to split: " << edges_need_to_split.size() << endl;
    if(edges_need_to_split.empty()) return 0;

    ordered_table_builder b(2); // split edges

    for(size_t ei = 0; ei < edges_need_to_split.size()/2; ++ei)
        b << &edges_need_to_split[2 * ei];
    b >> new_node_parent_id;

    matrix<size_t> before_validate = new_node_parent_id;
    validate_edge_nodes(orig_tet, before_validate, new_node_parent_id);
    cerr << "# validate subdivide: " << before_validate.size(1)
         << " " << new_node_parent_id.size(1) << endl;

    cerr << "# [info] tet size before sub " << orig_tet.size(2) << endl;

    matrix<size_t> children_tet;
    matrix<double> children_orig_node, children_defom_node;
    subdivide_tetmesh(orig_tet, orig_node.size(2), new_node_parent_id, children_tet);

    subdivide_top2geo(orig_node, new_node_parent_id, children_orig_node);
    subdivide_top2geo(deform_node, new_node_parent_id, children_defom_node);

    cerr << " [info] sub tet size: " << children_tet.size(2) << endl;

    orig_tet = children_tet;
    orig_node = children_orig_node;
    deform_node = children_defom_node;

#if  0 // debug
    {
        // check if there are extra nodes;
        boost::unordered_set<size_t> used_node(orig_tet.begin(), orig_tet.end());
        if(used_node.size() != orig_node.size(2)){
            cerr << "# tet node does not match node number: " << used_node.size()
                 << " node used, total node number " << orig_node.size(2) << endl;
        }
    }
#endif
    return 0;
}

int local_remesh_tet(matrix<size_t> & tet,
                     matrix<double> & orig_node,
                     matrix<double> & deform_node)
{
    matrix<double> scaled_jacobian = zeros<double>(tet.size(2),1);
#if 1 // collapse
    {
        // phase 1: remove |scaled_jacobian| < 0.05 tet.
        for(size_t ti = 0; ti < tet.size(2); ++ti)
            scaled_jacobian[ti] = tet_scaled_jacobian(deform_node(colon(), tet(colon(),ti)));

        cerr << "# [scaled_jac] before remesh, min/avg/max "
             << *min_element(scaled_jacobian.begin(), scaled_jacobian.end()) << "/"
             << std::accumulate(scaled_jacobian.begin(), scaled_jacobian.end(),
                                0.0)/scaled_jacobian.size() << "/"
             << *max_element(scaled_jacobian.begin(), scaled_jacobian.end()) << endl;

        boost::unordered_set<pair<size_t,size_t> > edges_need_to_remove;
        const size_t edge2idx_[] = {0,1,0,2,0,3,1,2,1,3,2,3};

        vector<pair<double, size_t> > edge_def2idx(6);
        for(size_t ti = 0; ti < tet.size(2); ++ti){
            if(fabs(scaled_jacobian[ti]) < 1e-3){
                for(size_t ei = 0; ei < 6; ++ei){
                    const double def_len =
                            norm(deform_node(colon(), tet(edge2idx_[2*ei+1] ,ti))
                                 - deform_node(colon(), tet(edge2idx_[2*ei] ,ti)));
                    const double orig_len =
                            norm(orig_node(colon(), tet(edge2idx_[2*ei+1] ,ti))
                                 - orig_node(colon(), tet(edge2idx_[2*ei] ,ti)));
                    edge_def2idx[ei] = make_pair(def_len/orig_len, ei);
                }
                sort(edge_def2idx.begin(), edge_def2idx.end());
                pair<size_t,size_t> select_edge(
                            tet(edge2idx_[2 * edge_def2idx.front().second],ti),
                            tet(edge2idx_[2 * edge_def2idx.front().second+1],ti));
                sort_edge(select_edge);
                edges_need_to_remove.insert(select_edge);
            }
        }

        cerr << "# [info] edges_need_to_collapse: " << edges_need_to_remove.size() << endl;
        if(!edges_need_to_remove.empty()){

            collapse_tet_edge(tet, orig_node, deform_node,  edges_need_to_remove);

            scaled_jacobian = zeros<double>(tet.size(2),1);
            for(size_t ti = 0; ti < tet.size(2); ++ti)
                scaled_jacobian[ti] = tet_scaled_jacobian(deform_node(colon(), tet(colon(),ti)));

            cerr << "# [scaled_jac] after collapse, min/avg/max "
                 << *min_element(scaled_jacobian.begin(), scaled_jacobian.end()) << "/"
                 << std::accumulate(scaled_jacobian.begin(), scaled_jacobian.end(),
                                    0.0)/scaled_jacobian.size() << "/"
                 << *max_element(scaled_jacobian.begin(), scaled_jacobian.end()) << endl;
        }
    }
#endif

#if 1 // split
    {
        matrix<size_t> new_node_parent_idx;
        const double split_percent = 1.0;
        cerr << "# [befor split] tet " << tet.size() << " node " << orig_node.size(2) << endl;
        split_tet_at_large_arap_distortion(tet, orig_node, deform_node,
                                           new_node_parent_idx, split_percent);
        cerr << "# [after split] tet " << tet.size() << " node " << orig_node.size(2) << endl;


        scaled_jacobian = zeros<double>(tet.size(2),1);
        for(size_t ti = 0; ti < tet.size(2); ++ti)
            scaled_jacobian[ti] = tet_scaled_jacobian(deform_node(colon(), tet(colon(),ti)));

        cerr << "# [scaled_jac] after split, min/avg/max "
             << *min_element(scaled_jacobian.begin(), scaled_jacobian.end()) << "/"
             << std::accumulate(scaled_jacobian.begin(), scaled_jacobian.end(),
                                0.0)/scaled_jacobian.size() << "/"
             << *max_element(scaled_jacobian.begin(), scaled_jacobian.end()) << endl;
        //  {
        //    ofstream ofs_tet("tet_scaled_jacobian.vtk");
        //    tet2vtk(ofs_tet, &orig_node[0], orig_node.size(2), &tet[0], tet.size(2));
        //    cell_data(ofs_tet, &scaled_jacobian[0], scaled_jacobian.size(), "scaled_jacobian");
        //  }

        cerr << "# [validate begin] " << endl;
        set<size_t> tet_nodes(tet.begin(), tet.end());
        if(tet_nodes.size() != orig_node.size(2) ||
                tet_nodes.size() != deform_node.size(2))
            cerr << "# [error] not compressed tet " << endl;
        cerr << "# [validate end]" << endl;
    }
#endif
    return 0;
}

int convert_surface_type_to_surface_patches(
        const matrix<size_t> &tet,
        const boost::unordered_map<size_t,size_t> &surface_type,
        vector<matrix<size_t> > &surface_patches,
        boost::unordered_set<pair<size_t,size_t> > &patch_boundary)
{
    unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
    if(!fa.get()){
        cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
        return __LINE__;
    }

    matrix<size_t> outside_face, outside_face_idx;
    get_outside_face(*fa, outside_face,true);
    get_outside_face_idx(*fa, outside_face_idx);

    unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
                jtf::mesh::edge2cell_adjacent::create(outside_face));
    if(!ea.get()){
        cerr << "# [error] can not build edge2cell_adjacent." << endl;
        return __LINE__;
    }

    vector<boost::unordered_set<size_t> > surface_patches_set;

    extra_surface_patch_according_to_type(surface_type, outside_face_idx, *ea, *fa,
                                          surface_patches_set, patch_boundary);

    surface_patches.resize(surface_patches_set.size());
    for(size_t pi = 0; pi < surface_patches.size(); ++pi){
        surface_patches[pi].resize(3,surface_patches_set[pi].size());
        const boost::unordered_set<size_t> & one_patch = surface_patches_set[pi];
        size_t fi = 0;
        for(boost::unordered_set<size_t>::const_iterator cit = one_patch.begin();
            cit != one_patch.end(); ++cit, ++fi){
            const vector<size_t> & one_face = fa->faces_[*cit];
            copy(one_face.begin(), one_face.end(), surface_patches[pi](colon(),fi).begin());
        }
    }

    return 0;
}

int smooth_boundary(
        const boost::unordered_set<pair<size_t,size_t> > & boundary,
        const matrix<double> & orig_node,
        matrix<double> & polycube_node,
        const size_t iter)
{
    if(iter == 0) return 0;
    vector<deque<pair<size_t,size_t> > > chains;
    vector<pair<size_t,size_t> > boundary_edges(boundary.begin(), boundary.end());
    jtf::util::extract_chain_from_edges(boundary_edges, chains);

    // smooth each boundary
    //for(size_t it = 0; it < iter; ++it){
    for(size_t li = 0; li < chains.size(); ++li){
        const deque<pair<size_t,size_t> > & one_chain = chains[li];
        vector<double> length_seg(one_chain.size()+1,0);
        for(size_t ei = 0; ei < one_chain.size(); ++ei) {
            length_seg[ei+1] = norm(orig_node(colon(), one_chain[ei].first) -
                                    orig_node(colon(), one_chain[ei].second)) +
                    length_seg[ei];
        }

        matrix<double> dis = polycube_node(colon(), one_chain.back().second)
                - polycube_node(colon(), one_chain.front().first);
        for(size_t ei = 0; ei < one_chain.size(); ++ei){
            const pair<size_t,size_t> & edge = one_chain[ei];
            //const pair<size_t,size_t> & edge_next = one_chain[ei+1];

            polycube_node(colon(), edge.second) =
                    polycube_node(colon(), one_chain.front().first) +
                    dis*length_seg[ei+1]/length_seg.back();
        }
    }
    //}
    return 0;
}

int smooth_patch(const vector<matrixst> & surface_patch,
                 const boost::unordered_set<pair<size_t,size_t> > &boundary,
                 matrixd & node,
                 const size_t iter)
{
    if(iter == 0) return 0;
    boost::unordered_set<size_t> fix_node;
    for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit =
        boundary.begin();  cit != boundary.end(); ++cit){
        fix_node.insert(cit->first);
        fix_node.insert(cit->second);
    }

    typedef boost::unordered_map<size_t, vector<size_t> > p2p_type;
    vector<p2p_type> p2p_group(surface_patch.size());

    for(size_t pi = 0; pi < surface_patch.size(); ++pi){
        unique_ptr<jtf::mesh::one_ring_point_at_point> orpap(
                    jtf::mesh::one_ring_point_at_point::create(surface_patch[pi]));
        p2p_group[pi] = orpap->p2p_;
    }

    // remove fix_node
    for(size_t pi = 0; pi < p2p_group.size(); ++pi){
        p2p_type & p2p_ = p2p_group[pi];
        for(boost::unordered_set<size_t>::const_iterator cit = fix_node.begin();
            cit != fix_node.end(); ++cit){
            const size_t & node_idx = *cit;
            p2p_type::iterator pit = p2p_.find(node_idx);
            if(pit != p2p_.end()) p2p_.erase(pit);
        }
    }


    matrixd center_node = zeros<double>(3,1);
    for(size_t it = 0; it < iter; ++it){
        for(size_t pi = 0; pi < p2p_group.size(); ++pi){
            const p2p_type & p2p_ = p2p_group[pi];
            for(p2p_type::const_iterator pcit = p2p_.begin();
                pcit != p2p_.end(); ++pcit) {
                center_node = zeros<double>(3,1);
                const vector<size_t> & linked_node = pcit->second;
                for(size_t scidx = 0; scidx < linked_node.size(); ++scidx){
                    center_node += node(colon(), linked_node[scidx]);
                }
                center_node /= linked_node.size();
                node(colon(), pcit->first) = center_node;
            }
        }
    }
    return 0;
}

int smooth_volume(
        const matrix<size_t> & tet,
        const vector<matrix<size_t> > &surface_patches,
        matrix<double> &node,
        const size_t iter_v )
{
    boost::unordered_set<size_t> fix_node;
    for(size_t pi = 0; pi < surface_patches.size(); ++pi){
        fix_node.insert(surface_patches[pi].begin(), surface_patches[pi].end());
    }

    jtf::tetmesh::one_ring_point_at_point orpap;
    orpap.build(tet);

    boost::unordered_map<size_t, boost::unordered_set<size_t> > p2p ;

    //p2p = orpap.p2p_;
    // only smooth flipped tet
    {
        boost::unordered_set<size_t> associated_nodes;
        double vol;
        for(size_t ti = 0; ti < tet.size(2); ++ti){
            vol = jtf::mesh::cal_tet_vol(node(colon(), tet(colon(),ti)));
            if(vol < 0){
                // flipped_tet.push_back(ti);
                associated_nodes.insert(tet(colon(),ti).begin(), tet(colon(),ti).end());
            }
        }
        cerr << "# [info] left " << associated_nodes.size() << " flipped nodes need"
             << " to be smoothed." << endl;
        for(boost::unordered_set<size_t>::const_iterator cit = associated_nodes.begin();
            cit != associated_nodes.end(); ++cit){
            boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
                    mscit = orpap.p2p_.find(*cit);
            if(mscit == orpap.p2p_.end()){
                cerr << "# [error] can not find p2p mapping." << endl;
                return __LINE__;
            }
            p2p[*cit] = mscit->second;
        }

    }


    for(boost::unordered_set<size_t>::const_iterator cit = fix_node.begin();
        cit != fix_node.end(); ++cit){
        boost::unordered_map<size_t, boost::unordered_set<size_t> >::iterator it
                = p2p.find(*cit);
        if(it  != p2p.end()) p2p.erase(it);
    }

    matrix<double> center_node = zeros<double>(3,1);
    for(size_t it = 0; it < iter_v; ++it){
        for(boost::unordered_map<size_t,boost::unordered_set<size_t> >::const_iterator
            cit = p2p.begin(); cit != p2p.end(); ++cit){
            center_node = zeros<double>(3,1);
            const boost::unordered_set<size_t> & linked_nodes = cit->second;
            for(boost::unordered_set<size_t>::const_iterator lit = linked_nodes.begin();
                lit != linked_nodes.end(); ++lit){
                center_node += node(colon(), *lit);
            }
            center_node /= linked_nodes.size();
            node(colon(), cit->first) = center_node;
        }
    }
    return 0;
}

int get_face_patches_according_to_boundary(
        const zjucad::matrix::matrix<size_t> &outside_face,
        const jtf::mesh::edge2cell_adjacent & ea,
        const boost::unordered_set<std::pair<size_t,size_t> > &boundary_edges,
        std::vector<std::vector<size_t> > &patches)
{
    // this function assume each edge in boundary_edges is increasing order
    vector<bool> is_face_visited(outside_face.size(2), false);

    stack<size_t> face_stack;

    patches.clear();

    vector<bool>::const_iterator cit = find(is_face_visited.begin(),
                                            is_face_visited.end(), false);
    while(cit != is_face_visited.end()){
        vector<size_t> one_patch;
        face_stack.push(cit-is_face_visited.begin());
        while(!face_stack.empty()){
            const size_t f_idx = face_stack.top();
            face_stack.pop();
            is_face_visited[f_idx] = true;
            one_patch.push_back(f_idx);

            for(size_t pi = 0; pi < outside_face.size(1); ++pi){
                pair<size_t,size_t> one_edge(
                            outside_face(pi,f_idx),
                            outside_face((pi+1)%outside_face.size(1),f_idx));
                if(one_edge.first > one_edge.second)
                    swap(one_edge.first, one_edge.second);
                if(boundary_edges.find(one_edge) != boundary_edges.end()) continue;

                const size_t edge_idx = ea.get_edge_idx(one_edge.first, one_edge.second);
                if(edge_idx == -1){
                    cerr << "# [error] strange can not find edge idx of "
                         << one_edge.first << " " << one_edge.second << endl;
                    return __LINE__;
                }
                const pair<size_t,size_t> & tri_pair = ea.edge2cell_[edge_idx];
                if(ea.is_boundary_edge(tri_pair)) continue;
                if(f_idx != tri_pair.first && f_idx != tri_pair.second){
                    cerr << "# [error] strange face pair of edge "
                         << one_edge.first << " " << one_edge.second << " is "
                         << tri_pair.first << " " << tri_pair.second << " without "
                         << f_idx << endl;
                    return __LINE__;
                }
                const size_t other_face_idx = tri_pair.first + tri_pair.second - f_idx;
                if(is_face_visited[other_face_idx]) continue;
                face_stack.push(other_face_idx);
            }
        }
        patches.push_back(one_patch);
        cit = find(is_face_visited.begin(), is_face_visited.end(), false);
    }

    return 0;
}

int remove_one_flipped_tet(
        matrixd & node,
        matrixst & tet,
        const matrixst & one_tet)
{
    sxx::tet_mesh stm;
    stm.create_tetmesh(node, tet);
    int rtn = stm.split_tet(&one_tet[0]);
    if(rtn == -1) return __LINE__;
    const size_t new_idx = rtn;

    for(size_t pi = 0; pi < one_tet.size(); ++pi){
        int rtn = stm.collapse_edge(make_pair(one_tet[pi], new_idx),true);
        if(rtn == -1) return __LINE__;
    }

    stm.write_tetmesh_to_matrix(node, tet);

    return 0;
}

int remove_flipped_tet(
        matrixst &orig_tet,
        matrixd &node,
        matrixst &polycube_tet,
        matrixd &polycube_node)
{
    vector<size_t> flipped_tet;
    for(size_t ti = 0; ti < polycube_tet.size(2); ++ti){
        double vol = jtf::mesh::cal_tet_vol(polycube_node(colon(), polycube_tet(colon(),ti)));
        if(vol < 0)
            flipped_tet.push_back(ti);
    }

    cerr << "# [info] flipped_tet number " << flipped_tet.size() << endl;
    for(size_t ti = 0; ti < flipped_tet.size(); ++ti){
        cerr << "# " << ti << endl;
        remove_one_flipped_tet(node, orig_tet, orig_tet(colon(), flipped_tet[ti]));
        remove_one_flipped_tet(polycube_node, polycube_tet, polycube_tet(colon(), flipped_tet[ti]));
    }


    {
        // recheck
        flipped_tet.clear();
        for(size_t ti = 0; ti < polycube_tet.size(2); ++ti){
            double vol = jtf::mesh::cal_tet_vol(polycube_node(colon(), polycube_tet(colon(),ti)));
            if(vol < 0)
                flipped_tet.push_back(ti);
        }
        cerr << "# [info] flipped_tet number " << flipped_tet.size() << endl;
    }
    return 0;
}

int dump_out_normal_flipped_face(
        const zjucad::matrix::matrix<double> &node,
        const std::vector<zjucad::matrix::matrix<size_t> > &surface_patches,
        const std::vector<zjucad::matrix::matrix<double> > &patch_normal,
        const size_t i,
        const zjucad::matrix::matrix<double> &vtk_node)
{
    ostringstream os;
    os << i;
    string name = "flipped_face_iter_" + os.str() + ".vtk";

    vector<size_t> flipped_face;
    matrix<double> normal = zeros<double>(3,1), node_;
    for(size_t pi = 0; pi < surface_patches.size(); ++pi){
        const matrix<size_t> & one_patch = surface_patches[pi];
        for(size_t fi = 0; fi < one_patch.size(2); ++fi){
            node_ = node(colon(), one_patch(colon(),fi));
            calc_tri_normal_(&normal[0], &node_[0]);
            if(dot(normal, patch_normal[pi]) < 0.1)
                flipped_face.insert(flipped_face.end(),
                                    one_patch(colon(),fi).begin(),
                                    one_patch(colon(),fi).end());
        }
    }

    ofstream ofs(name.c_str());
    tri2vtk(ofs, &vtk_node[0], vtk_node.size(2), &flipped_face[0], flipped_face.size()/3);
    return 0;
}

double get_tet_arap_distortion(
        const zjucad::matrix::matrix<double> & orig_node,
        const zjucad::matrix::matrix<double> & new_node,
        const zjucad::matrix::matrix<size_t> & one_tet)
{
    matrix<double> orig_tet = orig_node(colon(), one_tet);
    matrix<double> new_tet = new_node(colon(), one_tet);
    matrix<double> dis = zeros<double>(3,3);
    calc_tet_arap_distortion(&orig_tet[0], &new_tet[0], &dis[0]);
    return norm(dis);
}

double tet_scaled_jacobian(const matrix<double> & tet_node)
{
    // tet_node is 3 * 4
    VERDICT_REAL tet_node_coor[4][3];
    for(size_t pi = 0; pi < 4; ++pi){
        for(size_t di = 0; di < 3; ++di){
            tet_node_coor[pi][di] = tet_node(di, pi);
        }
    }
    return v_tet_scaled_jacobian(tet_node.size(2), tet_node_coor);
}

int uniform_vertex_order(std::vector<size_t> &tri_a,
                         std::vector<size_t> &tri_b,
                         const zjucad::matrix::matrix<size_t> *cut_tet2tet)
{
    assert(tri_a.size() == 3);
    assert(tri_b.size() == 3);

    if(!cut_tet2tet){
        vector<size_t> same_points;
        for(size_t i = 0; i < tri_a.size(); ++i)
            for(size_t j = 0; j < tri_b.size(); ++j){
                if(tri_b[j] == tri_a[i]) {
                    same_points.push_back(tri_b[j]);
                    continue;
                }
            }
        if(same_points.size() != 2){
            cerr << "# [error] can not find common edge of two faces." << endl;
            return __LINE__;
        }
        size_t a_left_point = std::accumulate(tri_a.begin(), tri_a.end(), 0) -
                std::accumulate(same_points.begin(), same_points.end(),0);
        size_t b_left_point = std::accumulate(tri_b.begin(), tri_b.end(), 0) -
                std::accumulate(same_points.begin(), same_points.end(),0);
        tri_a[0] = same_points[0];
        tri_a[1] = same_points[1];
        tri_a[2] = a_left_point;

        tri_b[0] = same_points[1];
        tri_b[1] = same_points[0];
        tri_b[2] = b_left_point;
    }else{
        vector<size_t> same_points_a, same_points_b;
        for(size_t i = 0; i < tri_a.size(); ++i)
            for(size_t j = 0; j < tri_b.size(); ++j){
                if((*cut_tet2tet)[tri_b[j]] ==
                        (*cut_tet2tet)[tri_a[i]]) {
                    same_points_a.push_back(tri_a[i]);
                    same_points_b.push_back(tri_b[j]);
                    continue;
                }
            }
        if(same_points_a.size() != 2 || same_points_b.size() != 2){
            cerr << "# [error] can not find common edge of two faces." << endl;
            return __LINE__;
        }
        size_t a_left_point = std::accumulate(tri_a.begin(), tri_a.end(), 0) -
                std::accumulate(same_points_a.begin(), same_points_a.end(),0);
        size_t b_left_point = std::accumulate(tri_b.begin(), tri_b.end(), 0) -
                std::accumulate(same_points_b.begin(), same_points_b.end(),0);

        tri_a[0] = same_points_a[0];
        tri_a[1] = same_points_a[1];
        tri_a[2] = a_left_point;

        tri_b[0] = same_points_b[1];
        tri_b[1] = same_points_b[0];
        tri_b[2] = b_left_point;

        if(tri_a[0] == tri_b[1]){ // here means the two face share one point
            std::rotate(tri_b.begin(), tri_b.begin()+1, tri_b.end());
            assert(tri_a[0] == tri_b[0]);
        }else if(tri_a[1] == tri_b[0]){
            std::rotate(tri_a.begin(), tri_a.begin()+1, tri_a.end());
            assert(tri_a[0] == tri_b[0]);
        }
    }
    return 0;
}

int load_chain_file(const char * chain_file,
                    std::vector<std::vector<size_t> > & chains)
{
    string chain_name;
    size_t chain_index, chain_size, element;
    ifstream ifs(chain_file);
    if(ifs.fail()){
        cerr << "# [error] can not open chain file" << endl;
        return __LINE__;
    }

    chains.clear();
    while(!ifs.eof()){
        ifs >> chain_name >> chain_index >> chain_size;
        if(chain_name.empty()) break;
        vector<size_t> one_chain(chain_size);

        for(size_t i = 0; i < chain_size; ++i){
            ifs >> one_chain[i];
        }
        chains.push_back(one_chain);
        chain_name.clear();
    }
    return 0;
}



int dump_out_group_file(const char * group_file,
                        const std::vector<std::vector<size_t> > &groups,
                        const std::vector<size_t> &group_integer)
{
    ofstream ofs(group_file);
    if(ofs.fail()){
        cerr << "# [error] can not create group file." << endl;
        return __LINE__;
    }

    assert(groups.size() == group_integer.size());

    for(size_t gi = 0; gi < groups.size(); ++gi){
        const vector<size_t> & one_group = groups[gi];
        ofs << "g " << gi << " " << one_group.size() << " " << group_integer[gi] << endl;
        for(size_t ti = 0; ti < one_group.size(); ++ti)
            ofs << one_group[ti] << " ";
        ofs << endl;
    }

    return 0;
}

int dump_out_chain_file(const char * chain_file,
                        const vector<vector<size_t> > & chain)
{
    ofstream ofs(chain_file);
    if(ofs.fail()){
        cerr << "# [error] can not create chain_file." << endl;
        return __LINE__;
    }

    for(size_t ci = 0; ci  < chain.size(); ++ci){
        const vector<size_t> & one_chain = chain[ci];
        ofs << "chain " << ci << one_chain.size() << endl;
        for(size_t i = 0; i < one_chain.size(); ++i)
            ofs << one_chain[i] << " ";
        ofs << endl;
    }
    return 0;
}

static int split_surface_tets(matrix<size_t> & tet,
                              matrix<double> & node,
                              matrix<double> & def_node,
                              const jtf::mesh::face2tet_adjacent & fa,
                              const matrix<size_t> & outside_face,
                              const matrix<size_t> & outside_face_idx)
{
    set<size_t> surface_node(outside_face.begin(), outside_face.end());
    unique_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(outside_face));
    if(!ea.get()){
        cerr << "# [error] can not create edge2cell_adjacent." << endl;
        return __LINE__;
    }

    set<pair<size_t,size_t> > edges_need_to_split;
    pair<size_t,size_t> temp_edge;
    for(size_t fi = 0; fi < outside_face_idx.size(); ++fi){
        const pair<size_t,size_t> & face_pair = fa.face2tet_[outside_face_idx[fi]];
        assert(fa.is_outside_face(face_pair));
        const size_t tet_idx = face_pair.first==-1?face_pair.second:face_pair.first;
        const size_t other_point_idx =
                std::accumulate(tet(colon(),tet_idx).begin(),
                                tet(colon(),tet_idx).end(), 0)
                - std::accumulate(outside_face(colon(),fi).begin(),
                                  outside_face(colon(),fi).end(), 0);
        if(surface_node.find(other_point_idx) == surface_node.end()) continue;
        for(size_t pi = 0; pi < outside_face.size(1); ++pi){
            if(ea->get_edge_idx(outside_face(pi,fi), other_point_idx) == -1){
                temp_edge.first = outside_face(pi,fi);
                temp_edge.second = other_point_idx;
                break;
            }
        }
    }

    cerr << "# [info] edges need to split " << edges_need_to_split.size() << endl;

    if(edges_need_to_split.size() == 0) return 1;

    sxx::tet_mesh stm;
    stm.create_tetmesh(node, tet);

    map<size_t,pair<size_t,size_t> > new_node;
    for(set<pair<size_t,size_t> >::const_iterator cit = edges_need_to_split.begin();
        cit != edges_need_to_split.end(); ++cit){
        size_t rtn = stm.split_edge(*cit);
        if(rtn != -1){
            new_node[rtn] = *cit;
        }
    }

    stm.write_tetmesh_to_matrix(node, tet);
    assert(def_node.size(2) + new_node.size() == node.size(2));
    size_t i = def_node.size(2);

    def_node.resize(3, node.size(2));

    for(map<size_t,pair<size_t,size_t> >::const_iterator cit = new_node.begin();
        cit != new_node.end(); ++cit){
        const pair<size_t,size_t> & edge = cit->second;
        assert(cit->first == i);
        def_node(colon(),i) =
                (def_node(colon(), edge.first) + def_node(colon(), edge.second))/2.0;
    }

    return 0;
}

int remove_degenerated_face_of_tet(
        zjucad::matrix::matrix<size_t> &orig_tet,
        zjucad::matrix::matrix<double> &orig_node,
        zjucad::matrix::matrix<double> &def_node,
        zjucad::matrix::matrix<size_t> &tri_faces,
        const double min_angle_threshold,
        const double max_angle_threshold)
{
    // this function will:
    // 1. detect obstue angle > 170, split the largest edge
    // 2. detect angle < 5, collapse the other edge
    using namespace jtf::mesh;
    set<pair<size_t,size_t> > edges_need_to_collapse;
    set<pair<size_t,size_t> > edges_need_to_split;
    size_t min_angle_point = -1;
    size_t max_angle_point = -1;

    vector<double> angle(3);
    for(size_t fi = 0; fi < tri_faces.size(2); ++fi){

        cal_face_angle(tri_faces(colon(),fi), def_node, angle);
        if(*max_element(angle.begin(), angle.end()) > max_angle_threshold){
            max_angle_point = max_element(angle.begin(), angle.end()) - angle.begin();
            pair<size_t,size_t> edge(tri_faces((max_angle_point+1)%tri_faces.size(1),fi),
                                     tri_faces((max_angle_point-1+tri_faces.size(1))%tri_faces.size(1),fi));
            if(edge.first > edge.second)
                swap(edge.first, edge.second);
            edges_need_to_split.insert(edge);
        }
    }

    // detect angle > 170, split the large edge
    cerr << "# [info] edges_need_to_split " << edges_need_to_split.size() << endl;
    if(edges_need_to_split.size()){// do not split
        sxx::tet_mesh stm;
        vector<pair<size_t,size_t> > new_node;
        stm.create_tetmesh(orig_node, orig_tet);
        for(set<pair<size_t,size_t> >::const_iterator cit = edges_need_to_collapse.begin();
            cit != edges_need_to_collapse.end(); ++cit){
            size_t new_node_idx = stm.split_edge(*cit);
            if(new_node_idx != -1){
                new_node.push_back(*cit);
            }
        }

        stm.write_tetmesh_to_matrix(orig_node, orig_tet);

        size_t orig_node_number = def_node.size(2);
        def_node.resize(3, orig_node.size(2));
        for(size_t i = 0; i < new_node.size(); ++i){
            const pair<size_t,size_t> & node_pair = new_node[i];
            def_node(colon(),i+orig_node_number) =
                    (def_node(colon(),node_pair.first) + def_node(colon(), node_pair.second)) /2.0;
        }


        unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(orig_tet));
        if(!fa.get()){
            cerr << "# [error] can not create fa2tet_adjacent." << endl;
            return __LINE__;
        }
        get_outside_face(*fa, tri_faces, true);
    }

    for(size_t fi = 0; fi < tri_faces.size(2); ++fi){
        if(cal_min_angle_of_one_face(tri_faces(colon(),fi), def_node, min_angle_point) < min_angle_threshold){
            pair<size_t,size_t> edge(tri_faces((min_angle_point+1)%tri_faces.size(1),fi),
                                     tri_faces((min_angle_point-1+tri_faces.size(1))%tri_faces.size(1),fi));
            if(edge.first > edge.second)
                swap(edge.first, edge.second);
            edges_need_to_collapse.insert(edge);
        }
    }

    vector<pair<double,pair<size_t,size_t> > > edge_len_edges;
    edge_len_edges.reserve(edges_need_to_collapse.size());
    for(set<pair<size_t,size_t> >::const_iterator cit = edges_need_to_collapse.begin();
        cit != edges_need_to_collapse.end(); ++cit){
        const pair<size_t,size_t> & one_edge = *cit;
        edge_len_edges.push_back(make_pair(norm(def_node(colon(), one_edge.first)-
                                                def_node(colon(), one_edge.second)),
                                           one_edge));
    }

    sort(edge_len_edges.begin(), edge_len_edges.end());

    cerr << "# [info] edges_need_to_collapse " << edges_need_to_collapse.size() << endl;
    if(edges_need_to_collapse.empty()) return 1;

    sxx::tet_mesh stm, stm_def;
    matrix<size_t> def_tet = orig_tet;
    stm.create_tetmesh(orig_node, orig_tet);
    stm_def.create_tetmesh(def_node, def_tet);
    for(size_t ei = 0; ei < edge_len_edges.size(); ++ei){
        int rtn = stm.collapse_edge(edge_len_edges[ei].second, true);
        if(rtn != -1)
            stm_def.collapse_edge(edge_len_edges[ei].second, false);
    }

    stm.write_tetmesh_to_matrix(orig_node, orig_tet);
    stm_def.write_tetmesh_to_matrix(def_node, def_tet);

    unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(orig_tet));
    if(!fa.get()){
        cerr << "# [error] can not create fa2tet_adjacent." << endl;
        return __LINE__;
    }
    matrix<size_t> outside_face_idx;
    get_outside_face(*fa, tri_faces, true);
    get_outside_face_idx(*fa, outside_face_idx);

    int rtn = split_surface_tets(orig_tet, orig_node, def_node, *fa, tri_faces, outside_face_idx);
    if(rtn == 1) return 0;
    if(rtn == 0){
        unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(orig_tet));
        if(!fa.get()){
            cerr << "# [error] can not create fa2tet_adjacent." << endl;
            return __LINE__;
        }
        get_outside_face(*fa, tri_faces, true);
    }
    return 0;
}

