#include "singularity_adjustment.h"
#include "../common/vtk.h"
#include "../common/IO.h"
#include "common.h"
#include "topology_operation.h"
#include <iostream>
using namespace std;
using namespace zjucad::matrix;

//////////////////////////////////////////////////////////////
//  No used methods: just back up
//////////////////////////////////////////////////////////////

size_t get_new_begin_vertex(const vector<deque<size_t> > &singularity_edges,
                            const size_t begin_,
                            const size_t current_chain_idx)
{
    for(size_t t = 0 ;t < singularity_edges.size(); ++t) {
        if(t != current_chain_idx) {
            if(singularity_edges[t].front() == begin_) {
                return singularity_edges[t][1];
            }
            if(singularity_edges[t].back() == begin_) {
                return singularity_edges[t][singularity_edges[t].size() - 2];
            }
        }
    }
    return -1;
}



int find_around_tets(const pair<size_t,size_t> &edge,
                     vector<size_t> &around_tets,
                     const face2tet_adjacent &fa,
                     const matrix<size_t> &input_tet)
{
    //set<pair<size_t,size_t> > tet_segment;
    around_tets.clear();
    for(size_t t = 0; t < input_tet.size(2); ++t){
        if(find(input_tet(colon(),t).begin(), input_tet(colon(),t).end(),edge.first) != input_tet(colon(),t).end()
                && find(input_tet(colon(),t).begin(), input_tet(colon(),t).end(),edge.second) != input_tet(colon(),t).end())
        {
            around_tets.push_back(t);
            break;
        }
    }
    vector<size_t> other_vertex;//[2];
    for(size_t t = 0; t < 4; ++t)
    {
        if(input_tet(t,around_tets.front()) != edge.first
                && input_tet(t,around_tets.front()) != edge.second )
            other_vertex.push_back(input_tet(t,around_tets.front()));
    }
#if 0
    cerr << "# check other_vertex size = " << other_vertex.size() << endl;
#endif
    vector<size_t> face(3);
    face[0] = edge.first;
    face[1] = edge.second;
    face[2] = other_vertex[0];

    while(1)
    {
        const pair<size_t,size_t> & tet_pair = fa.face2tet_[fa.get_face_idx(&face[0])];
        if(tet_pair.first != around_tets.back()) around_tets.push_back(tet_pair.first);
        else if(tet_pair.second != around_tets.back()) around_tets.push_back(tet_pair.second);
        if(around_tets.back() == around_tets.front()) break;
        for(size_t t = 0; t < 4; ++t){
            if(input_tet(t,around_tets.back()) != face[0]
                    && input_tet(t,around_tets.back()) != face[1]
                    && input_tet(t,around_tets.back()) != face[2]) {
                face[2] = input_tet(t,around_tets.back());
                break;
            }
        }
    }
    return 0;
}

int find_around_tets_new(const pair<size_t,size_t> &edge,
                         const matrix<size_t> &tet,
                         const matrix<double> &node,
                         const face2tet_adjacent &fa,
                         const one_ring_tet_at_edge &ortae,
                         vector<size_t> &around_tets)
{
    typedef one_ring_tet_at_edge::e2tet_type::const_iterator oci;
    oci it = ortae.e2t_.find(edge);
    if(it == ortae.e2t_.end()) it = ortae.e2t_.find(make_pair(edge.second,edge.first));
    if(it == ortae.e2t_.end()) {
        cerr << "# strange, can not find edge." << endl;
        return __LINE__;
    }
    const vector<size_t> &tet_loop = it->second;
    if(tet_loop.front() != tet_loop.back()) {
        cerr << "# strange error may be sharp edge." << endl;
        return __LINE__;
    }
    else
    {
        around_tets.clear();
        size_t t = 0;
        if(tet_loop.front() == -1) t = 1;
        else t = 0;
        around_tets.push_back(tet_loop[t]);
        for(; t < tet_loop.size(); ++t){
            if(tet_loop[t] != around_tets.back())
                around_tets.push_back(tet_loop[t]);
        }
        if(around_tets.back() == -1) around_tets.pop_back();
    }
    return 0;
}

// no use
int reorder_singularity_chain(vector<deque<size_t> > &singularity_edges,
                              vector<deque<size_t> > &singularity_type,
                              const matrix<size_t> &outside_face)
{
    map<size_t,set<size_t> > chain_end_count;

    for(size_t t = 0; t < singularity_edges.size(); ++t){
        chain_end_count[singularity_edges[t].front()].insert(t);
        chain_end_count[singularity_edges[t].back()].insert(t);
    }

    //set<size_t> need_to_merge_chain;
    vector<bool> left_to_be_new_singularity_chain(singularity_edges.size(),true);
    const size_t black_type = 9;
    for(map<size_t,set<size_t> >::const_iterator mci = chain_end_count.begin();
        mci != chain_end_count.end(); ++mci)
    {
        if(mci->second.size() == 2 && !is_outside_vertex(mci->first,outside_face))
        {
            vector<size_t> merged_chain_idx(2);
            copy(mci->second.begin(),mci->second.end(),merged_chain_idx.begin());
            size_t from,to;
            if(find(singularity_type[merged_chain_idx[0]].begin(),
                    singularity_type[merged_chain_idx[0]].end(),
                    black_type) != singularity_type[merged_chain_idx[0]].end()) // merged_chain_idx[0] should be merged into merged_chain_idx[1]
            {
                from = merged_chain_idx[0];
                to = merged_chain_idx[1];
            }
            else if(find(singularity_type[merged_chain_idx[1]].begin(),
                         singularity_type[merged_chain_idx[1]].end(),
                         black_type) != singularity_type[merged_chain_idx[1]].end()) // merged_chain_idx[0] should be merged into merged_chain_idx[1]
            {
                from = merged_chain_idx[1];
                to = merged_chain_idx[0];
            }else
            {
                cerr << "# error singularity chain." << endl;
                continue;
            }

            left_to_be_new_singularity_chain[from] = false;
            if(singularity_edges[to].back() == singularity_edges[from].front())
            {
                for(size_t t = 0; t < singularity_edges[from].size(); ++t){
                    singularity_edges[to].push_back(singularity_edges[from][t]);
                    if(t % 2 == 0)
                        singularity_type[to].push_back(singularity_type[to].back());
                }
                \
            }
            if(singularity_edges[to].back() == singularity_edges[from].back())
            {
                for(size_t t = 0; t < singularity_edges[from].size(); ++t){
                    singularity_edges[to].push_back(singularity_edges[from][singularity_edges[from].size() - 1 - t]);
                    if(t % 2 == 0)
                        singularity_type[to].push_back(singularity_type[to].back());
                }
            }
            if(singularity_edges[to].front() == singularity_edges[from].front())
            {
                for(size_t t = 0; t < singularity_edges[from].size(); ++t){
                    singularity_edges[to].push_front(singularity_edges[from][singularity_edges[from].size() - 1 - t]);
                    if(t % 2 == 0)
                        singularity_type[to].push_front(singularity_type[to].front());
                }
            }
            if(singularity_edges[to].front() == singularity_edges[from].back())
            {
                for(size_t t = 0; t < singularity_edges[from].size(); ++t){
                    singularity_edges[to].push_front(singularity_edges[from][t]);
                    if(t % 2 == 0)
                        singularity_type[to].push_front(singularity_type[to].front());
                }
            }
        }
    }

    vector<deque<size_t> > tmp_singularity_chain_edge;
    vector<deque<size_t> > tmp_singularity_chain_type;
    for(size_t t = 0; t < left_to_be_new_singularity_chain.size(); ++t)
    {
        if(left_to_be_new_singularity_chain[t])
        {
            tmp_singularity_chain_edge.push_back(singularity_edges[t]);
            tmp_singularity_chain_type.push_back(singularity_type[t]);
        }
    }

    swap(tmp_singularity_chain_edge,singularity_edges);
    swap(tmp_singularity_chain_type,singularity_type);
    return 0;
}

bool is_vertex_in_set(const size_t t_idx,
                      const set<size_t> &s)
{
    if(find(s.begin(),s.end(),t_idx) == s.end())
        return false;
    return true;
}

//no use
int updata_singularity_chain(vector<deque<size_t> > &singularity_edges,
                             vector<deque<size_t> > &singularity_type,
                             const map<pair<size_t,size_t>,bool> &map_changed_edge,
                             const matrix<size_t> &outside_face)
{
    typedef map<pair<size_t,size_t>,bool >::const_iterator mci;
    typedef deque<size_t>::iterator dit;
    set<pair<size_t,size_t> > need_to_remove_edge;
    set<pair<size_t,size_t> > need_to_add_edge;
    for(mci it = map_changed_edge.begin(); it != map_changed_edge.end(); ++it){
        if(it->second)
            need_to_add_edge.insert(it->first);
        else
            need_to_remove_edge.insert(it->first);
    }

    // remove the non-singularity edges
    while(!need_to_remove_edge.empty())
    {
        for(set<pair<size_t,size_t> >::iterator sit = need_to_remove_edge.begin();
            sit != need_to_remove_edge.end(); ++sit){
            for(size_t t = 0; t < singularity_edges.size(); ++t){
                deque<size_t> &chain = singularity_edges[t];
                deque<size_t> &type_chain = singularity_type[t];
                dit it_first = find(chain.begin(),chain.end(),sit->first);
                dit it_second= find(chain.begin(),chain.end(),sit->second);
                if(it_first != chain.end() && it_second != chain.end()) {
                    if(it_first == chain.begin() || it_second == chain.begin()) {
                        //if(it_second != it_first + 2) cerr << "# strang singularity edge." << endl;
                        chain.pop_front(); // pop out the first edge(two points)
                        chain.pop_front();
                        type_chain.pop_front();
                        need_to_remove_edge.erase(sit);
                    }else if(it_second == chain.end() - 1 || it_first == chain.end() - 1) {
                        //if(it_first != chain.end() - 3) cerr << "# strang singularity edge." << endl;
                        chain.pop_back(); // pop out the last edge
                        chain.pop_back();
                        type_chain.pop_back();
                        need_to_remove_edge.erase(sit);
                    }
                }
            }
        }
    }
    const size_t black_type = 9;
    // collect all vertex on black lines, the any othen singularity chain if it touch one black_line, the touch end can be added new edges
    set<size_t> black_line_vertex;
    for(size_t t = 0; t < singularity_type.size(); ++t){
        for(size_t i = 0; i < singularity_type[t].size(); ++i){
            black_line_vertex.insert(singularity_type[t][i]);
        }
    }
    while(!need_to_add_edge.empty())
    {
        for(set<pair<size_t,size_t> >::iterator sit = need_to_add_edge.begin();
            sit != need_to_add_edge.end(); ++sit){
            for(size_t t = 0; t < singularity_edges.size(); ++t){
                deque<size_t> &chain = singularity_edges[t];
                deque<size_t> &type_chain = singularity_type[t];
                if(find(type_chain.begin(),type_chain.end(),black_type) != type_chain.end())
                    continue;
                dit it_first = find(chain.begin(),chain.end(),sit->first);
                dit it_second= find(chain.begin(),chain.end(),sit->second );
                if(it_first != chain.end() && it_second == chain.end())
                {
                    if(it_first == chain.end() - 1 && !is_vertex_in_set)
                    {
                        chain.push_back(sit->first);
                        chain.push_back(sit->second);
                        type_chain.push_back(type_chain.back());
                        need_to_add_edge.erase(sit);
                        break;
                    }else if(it_first == chain.begin())
                    {
                        chain.push_front(sit->first);
                        chain.push_front(sit->second);
                        type_chain.push_front(type_chain.front());
                        need_to_add_edge.erase(sit);
                        break;
                    }else
                        cerr << "# strang edge." << endl; // the edge order is wrong
                }
                if(it_second != chain.end() && it_first == chain.end())
                {
                    if(it_second == chain.begin() )
                    {
                        chain.push_front(sit->second);
                        chain.push_front(sit->first);
                        type_chain.push_front(type_chain.front());
                        need_to_add_edge.erase(sit);
                        break;
                    }else if(it_second == chain.end() - 1)
                    {
                        chain.push_back(sit->second);
                        chain.push_back(sit->first);
                        type_chain.push_back(type_chain.back());
                        need_to_add_edge.erase(sit);
                        break;
                    }else
                        cerr << "# strang edge." << endl; // the edge order is wrong
                }
                //            if(it_first != chain.end() && it_second != chain.end())
                //                 cerr << "# strange, edge is already on chain." << endl;
            }
        }
    }

#if 1
    cerr << "# check whether the singularity edge is right." << endl;
    size_t sum = 0;
    for(size_t t = 0; t < singularity_edges.size(); ++t){
        for(size_t i = 1; i < singularity_edges[t].size() -1 ; ++i)
        {
            if(i % 2 == 1)
                sum += singularity_edges[t][i];
            else
                sum -= singularity_edges[t][i];
        }
        if(sum != 0) cerr << "# error singularity chain." << sum << endl;
    }
#endif
    //reorder_singularity_chain(singularity_edges,singularity_type,outside_face);

    return 0;
}

//bool is_outside_edge(const size_t v0,
//                     const size_t v1,
//                     const matrix<size_t> &outside_face)
//{
//    for(size_t t = 0; t < outside_face.size(2); ++t){
//        const matrix<size_t> & face = outside_face(colon(),t);
//        if(find(face.begin(),face.end(),v0) != face.end()
//                && find(face.begin(),face.end(),v1) != face.end() )
//            return true;
//    }
//    return false;
//}


int  update_edges(map<pair<size_t,size_t>,bool> &map_changed_edge,
                  map<pair<size_t,size_t>,size_t> &edges,
                  const matrix<size_t> &outside_face)
{
    typedef map<pair<size_t,size_t>,bool>::iterator mci;
    typedef map<pair<size_t,size_t>,size_t>::iterator mi;
    for(mci ci = map_changed_edge.begin(); ci != map_changed_edge.end(); ++ci){

        //TODO: strange this function will result in finding a wrong edge which is actually inner
        if(is_outside_edge(ci->first.first,ci->first.second,outside_face)){
            map_changed_edge.erase(ci);
            continue;
        }

        if(ci->second) // is singularity edges
        {
            mi it = edges.find(ci->first);
            if(it == edges.end()) it = edges.find(make_pair(ci->first.second,ci->first.first));
            if(it == edges.end()) { // if it's a new singularity edge
                edges.insert(make_pair(ci->first,-1));
            }else{ // this singularity edge is relabelled as new type
                it->second = -1;
            }

            //            if(edges.find(ci->first) == edges.end()
            //                    && edges.find(make_pair(ci->first.second,ci->first.first)) == edges.end())
            //            {
            //                edges.insert(make_pair(ci->first,-1)); // -1 means unknow type
            //            }
        } else{
            mi it = edges.find(ci->first);
            if(it == edges.end()) it = edges.find(make_pair(ci->first.second,ci->first.first));
            if(it == edges.end())
                cerr << "# strange can not find edge." << ci->first.first << " " << ci->first.second << endl;
            else
                edges.erase(it);
        }
    }
    return 0;
}


int relabel_remove_black_lines(const matrix<size_t> &tet,
                               const matrix<double> &node,
                               const matrix<size_t> &outside_face,
                               const vector<deque<size_t> > &singularity_edges,
                               const vector<deque<size_t> > &singularity_type,
                               map<pair<size_t,size_t>,size_t> &edges,
                               map<pair<size_t,size_t>,bool> &map_changed_edge)
{
    //map_changed_edge.clear();
    auto_ptr<face2tet_adjacent> fa(face2tet_adjacent::create(tet));
    const size_t black_type = 9;
    size_t begin_,end_;
    //map<pair<size_t,size_t>,bool> map_changed_edge;
    for(size_t t = 0; t < singularity_type.size(); ++t)
    {
        if(find(singularity_type[t].begin(), singularity_type[t].end(), black_type)
                !=  singularity_type[t].end()) {
            if(is_outside_vertex(singularity_edges[t].front(),outside_face)
                    && is_outside_vertex(singularity_edges[t].back(),outside_face)){
                //cerr << "# strange black line." << endl;
                cerr << "# both of this black line reach the surface." << endl;
                {// remove the whole black line
                    for(size_t i = 0; i < singularity_edges[t].size(); i += 2){
                        map_changed_edge.insert(make_pair(make_pair(singularity_edges[t][i],singularity_edges[t][i+1]),false));
                    }
                }
                continue;//return __LINE__;
            }
            if(is_outside_vertex(singularity_edges[t].front(),outside_face)) {
                end_ = singularity_edges[t].front();
                begin_ = singularity_edges[t].back();
            } else {
                end_ = singularity_edges[t].back();
                begin_ = singularity_edges[t].front();
            }
            const deque<size_t> &black_chain = singularity_edges[t];
            size_t new_begin_ = get_new_begin_vertex(singularity_edges,begin_,t);
            if(new_begin_ == -1) {
                cerr << "# strange: can not find other chain which link this black chain." << endl;
                return __LINE__;
            }
            find_shortest_path(new_begin_,end_,tet,node,black_chain,map_changed_edge);
            //map_changed_edge.insert(make_pair(make_pair(new_begin_,begin_),false)); // remove the no need singularity edge
        }
    }
    size_t t = 0;
    for(map<pair<size_t,size_t>,bool>::const_iterator mci = map_changed_edge.begin();
        mci != map_changed_edge.end(); ++mci){
        if(mci->second == false){
            ++t;
        }
    }
    cerr << "# compound singularity edges num = " << t << endl;
    update_edges(map_changed_edge,edges,outside_face);
}

int reorder_the_ortmap_by_singularity_chain(const vector<deque<pair<size_t,size_t> > > &chain_list,
                                            one_ring_tet_at_edge &ort)
{
    typedef one_ring_tet_at_edge::e2tet_type::iterator oeit;
    for(size_t t = 0; t < chain_list.size(); ++t){
        for(size_t i = 0; i < chain_list[t].size(); ++i){
            oeit it = ort.e2t_.find(chain_list[t][i]);
            if(it != ort.e2t_.end()) continue;
            else{
                it = ort.e2t_.find(make_pair(chain_list[t][i].second,chain_list[t][i].first));
                if(it == ort.e2t_.end()){
                    cerr << "# strange can not find the singularity edge." << endl;
                    return __LINE__;
                }
                vector<size_t> loop = it->second;
                reverse(loop.begin(),loop.end());
                ort.e2t_.insert(make_pair(chain_list[t][i],loop));
                ort.e2t_.erase(it);
            }
        }
    }
    return 0;
}


// extract the chain may result in the edge up and down, so need to modify the ort map
int assemble_singularity_chain(const map<pair<size_t,size_t>,size_t > &edges,
                               const matrix<size_t> &outside_face,
                               vector<deque<size_t> > &singularity_edges,
                               vector<deque<size_t> > &singularity_type,
                               one_ring_tet_at_edge &ort)
{
    const size_t black_type = 9;
    set<pair<size_t,size_t> > edge_segment;
    typedef map<pair<size_t,size_t>,size_t>::const_iterator mci;
    {// copy edges
        for(mci it = edges.begin(); it != edges.end(); ++it)
            edge_segment.insert(it->first);
    }

    vector<deque<pair<size_t,size_t> > > chain_list;
    //    one_ring_tet_at_edge ort;
    ort.cut_singularity_edges_into_chain(edge_segment,outside_face,chain_list);

    reorder_the_ortmap_by_singularity_chain(chain_list,ort);

    // output
    singularity_edges.clear();
    singularity_type.clear();
    singularity_edges.resize(chain_list.size());
    for(size_t t = 0; t < chain_list.size(); ++t){
        //singularity_edges[t].resize(2 * chain_list[t].size());
        for(size_t i = 0; i < chain_list[t].size();++i){
            singularity_edges[t].push_back(chain_list[t][i].first);
            singularity_edges[t].push_back(chain_list[t][i].second);
        }
    }
    //cerr << singularity_edges[5].size() << endl;
    singularity_type.resize(chain_list.size());
    for(size_t t = 0; t < singularity_type.size(); ++t){
        //singularity_type[t].resize(chain_list[t].size());
        for(size_t i = 0; i < chain_list[t].size(); ++i)
        {
            mci it = edges.find(chain_list[t][i]);
            if(it == edges.end()) it = edges.find(make_pair(chain_list[t][i].second,chain_list[t][i].first));
            if(it == edges.end()) {
                cerr << "# strange edge segment" << endl;
                return __LINE__;
            }
            if(it->second == -1) {
                singularity_type[t].push_back(it->second);
                continue;
            }
            if(it->first == chain_list[t][i])
                singularity_type[t].push_back(it->second);
            else
            {
                size_t type = it->second;
                if(type != black_type)
                {
                    if(type - type/3 * 3 == 0) type += 2;
                    else if(type - type/3 * 3 == 2) type -= 2;
                }
                singularity_type[t].push_back(type);
            }
        }
#if 0
        //        size_t begin_idx = -1;
        //        for(size_t i = 0; i < singularity_type[t].size(); ++i){
        //            if(singularity_type[t][i] != black_type && singularity_type[t][i] != -1){
        //                begin_idx = i;
        //                break;
        //            }
        //        }
        //        if(begin_idx == -1) cerr << "# strange whole chain is black or unknown." << endl;

        //        size_t type = singularity_type[t][begin_idx];
        //        for(size_t i = begin_idx; i != -1; i--){
        //            if(singularity_type[t][i] == black_type || singularity_type[t][i] == -1){
        //                singularity_type[t][i] = type;
        //            }
        //            else type = singularity_type[t][i];
        //        }

        //        type = singularity_type[t][begin_idx];
        //        for(size_t i = begin_idx; i < singularity_type[t].size(); ++i){
        //            if(singularity_type[t][i] == black_type || singularity_type[t][i] == -1){
        //                singularity_type[t][i] = type;
        //            }else type = singularity_type[t][i];
        //        }
#endif
        for(size_t i = 0; i < singularity_type[t].size(); ++i){
            if(singularity_type[t][i] == black_type) singularity_type[t][i] = -1; // unknown type
        }
    }

    return 0;
}








bool points_is_in_a_tet(size_t index1, size_t index2, size_t index3,
                        const matrix<size_t> &tet)
{
    size_t fitNum;
    for(size_t i = 0; i < tet.size(2); ++i)
    {
        fitNum = 0;
        for(size_t j = 0; j < tet.size(1); ++j)
        {
            if(tet(j, i) == index1) ++fitNum;
            if(tet(j, i) == index2) ++fitNum;
            if(tet(j, i) == index3) ++fitNum;
        }
        if(fitNum == 3) return true;
    }
    return false;
}


int smooth_singularity_chain(const matrix<size_t> &tet,
                             const vector<deque<size_t> > &singualrity_chain,
                             map<pair<size_t,size_t>,bool> & edges)
{
    vector<bool> isExistEdge;
    for(size_t i = 0; i < singualrity_chain.size(); ++i)
    {
        const deque<size_t>& chain = singualrity_chain[i];
        size_t chainNum = chain.size() / 2 + 1;
        if(chainNum <= 2) continue;
        isExistEdge.resize(0);
        isExistEdge.resize(chainNum, true);
        size_t pre, next;
        for(size_t j = 1; j < chainNum - 1; ++j)
        {
            for(pre = j - 1; !isExistEdge[pre]; --pre);
            next = j + 1;
            if(points_is_in_a_tet(chain[j * 2], chain[pre * 2], chain[next * 2 - 1], tet)) isExistEdge[j] = false;
        }
        for(size_t j = 0; j < chainNum - 1; ++j)
        {
            if(isExistEdge[j] && isExistEdge[j + 1]) continue;
            pair<size_t, size_t> edge(chain[j * 2], chain[j * 2 + 1]);
            edges[edge] = false;
        }
        size_t prePoint = 0;
        bool addEdge = false;
        for(size_t j = 1; j < chainNum; ++j)
        {
            if(!isExistEdge[j]) addEdge = true;
            else
            {
                if(addEdge)
                {
                    pair<size_t, size_t> edge(chain[prePoint * 2], chain[j * 2 - 1]);
                    edges[edge] = true;
                    addEdge = false;
                }
                prePoint = j;
            }
        }
    }
    return 0;
}


int relabel_smooth_singulairty_edges(const matrix<size_t> &tet,
                                     const matrix<double> &node,
                                     const matrix<size_t> &outside_face,
                                     const vector<deque<size_t> > &singularity_edges,
                                     map<pair<size_t,size_t>,size_t> &edges,
                                     map<pair<size_t,size_t>,bool> & map_changed_edge)
{
    //map_changed_edge.clear();
    smooth_singularity_chain(tet, singularity_edges,map_changed_edge);
    size_t t = 0;
    for(map<pair<size_t,size_t>,bool>::const_iterator mci = map_changed_edge.begin();
        mci != map_changed_edge.end(); ++mci){
        if(mci->second == false){
            ++t;
        }
    }
    cerr << "# zigzag singularity edges num = " << t << endl;
    update_edges(map_changed_edge,edges,outside_face);
}

int init_singularity_edges(const vector<deque<size_t> > &singularity_edges,
                           const vector<deque<size_t> > &singularity_type,
                           map<pair<size_t,size_t>,size_t> &edges)
{
    edges.clear();
    for(size_t t = 0; t < singularity_type.size(); ++t){
        for(size_t i = 0; i < singularity_type[t].size(); ++i){
            edges.insert(make_pair(make_pair(
                                       singularity_edges[t][i * 2],
                                       singularity_edges[t][i * 2 + 1]),
                                   singularity_type[t][i]));
        }
    }
    return 0;
}



bool is_near_surface_chain(const deque<size_t> &chain,
                           const matrix<size_t> &outside_face,
                           const matrix<size_t> &tet)
{
    if(is_outside_vertex(chain.front(),outside_face)
            && is_outside_vertex(chain.back(),outside_face))
    {
        set<size_t> chain_vertex;
        for(size_t t = 0; t < chain.size(); ++t)
            chain_vertex.insert(chain[t]);

        if(chain_vertex.size() < 5) return true; // 0 1 2 3 ,0 and 3 is on surface, then 1 2 must be near surface
        size_t one_ring_near_surface_num = 0;
        set<size_t>::const_iterator sci_begin = chain_vertex.begin();

        set<size_t>::const_iterator sci_end = chain_vertex.end();

        for(set<size_t>::const_iterator sci = sci_begin; sci != sci_end; ++sci){
            if(is_one_ring_near_surface(*sci,tet,outside_face))
                ++one_ring_near_surface_num;
        }
        //assert(one_ring_near_surface_num >= 2); // TODO strange
        if(one_ring_near_surface_num * 1.0 / chain_vertex.size() > 0.5)
            return true;
        else
            return false;
    }else
        return false;

}

int relabel_remove_short_lines(const matrix<size_t> &tet,
                               const matrix<double> &node,
                               const matrix<size_t> &outside_face,
                               const vector<deque<size_t> > &singularity_edges,
                               map<pair<size_t,size_t>,size_t> &edges,
                               map<pair<size_t,size_t>,bool> & map_changed_edge)
{
    vector<bool> removed(singularity_edges.size(),false);
    for(size_t t = 0; t <singularity_edges.size(); ++t)
    {
        const deque<size_t> & chain = singularity_edges[t];
        if(is_near_surface_chain(chain,outside_face,tet))
            removed[t] = true;
    }
    for(size_t t = 0; t < removed.size(); ++t)
    {
        if(removed[t]){
            for(size_t i = 0; i < singularity_edges[t].size(); i+=2){
                map_changed_edge.insert(make_pair(make_pair(singularity_edges[t][i],singularity_edges[t][i+1]),false));
            }
        }
    }

    size_t t = 0;
    for(map<pair<size_t,size_t>,bool>::const_iterator mci = map_changed_edge.begin();
        mci != map_changed_edge.end(); ++mci){
        if(mci->second == false){
            ++t;
        }
    }
    cerr << "# snap to surface singularity edges num = " << t << endl;
    update_edges(map_changed_edge,edges,outside_face);
}

int relabel_triangle_singualrity_edge(const matrix<size_t> &tet,
                                      const matrix<double> &node,
                                      const matrix<size_t> &outside_face,
                                      const vector<deque<size_t> > &singularity_edges,
                                      map<pair<size_t,size_t>,size_t> &edges,
                                      map<pair<size_t,size_t>,bool> & map_changed_edge)
{
    vector<size_t> need_to_remove_triangle_singularity_chain_candidate;
    map<size_t,vector<size_t> > vertex_link_chain; // store these chains get together at vertex p. first,vertex_idx; second, chain idx
    map<size_t,bool> vertex_link_inner_chain;
    for(size_t t = 0; t < singularity_edges.size(); ++t){
        const size_t &begin = singularity_edges[t].front();
        const size_t &end = singularity_edges[t].back();
        bool is_outside_vertex_begin = is_outside_vertex(begin,outside_face);
        bool is_outside_vertex_end = is_outside_vertex(end,outside_face);
        if(is_outside_vertex_begin ^ is_outside_vertex_end){ //only one vertex is on the surface
            need_to_remove_triangle_singularity_chain_candidate.push_back(t);
            if(is_outside_vertex_begin){
                vertex_link_chain[end].push_back(t);
            }else{
                vertex_link_chain[begin].push_back(t);
            }
        }else if(!is_outside_vertex_begin && !is_outside_vertex_end){
            vertex_link_inner_chain[begin] = true;
            vertex_link_inner_chain[end] = true;
        }
    }

    typedef map<size_t,bool>::const_iterator mci;
    for(map<size_t,vector<size_t> >::const_iterator mcit = vertex_link_chain.begin();
        mcit != vertex_link_chain.end(); ++mcit){
        if(mcit->second.size() > 1){ // if there is only one chain which reach the surface, we do not need to remove it
            mci cit = vertex_link_inner_chain.find(mcit->first);
            if(cit == vertex_link_inner_chain.end()){ // this vertex has no inner sinuglarity chain, need to remove all associated chains
                for(size_t t = 0; t < mcit->second.size(); ++t){
                    const size_t &chain_idx = mcit->second[t];
                    for(size_t j = 0; j < singularity_edges[chain_idx].size(); j += 2){
                        map_changed_edge[make_pair(singularity_edges[chain_idx][j],
                                                   singularity_edges[chain_idx][j+1])] = false;
                    }
                }
            }else{// at least one inner singularity chain link at this vertex
                vector<pair<size_t,size_t> > chain_length_idx;
                for(size_t t = 0; t < mcit->second.size(); ++t){
                    const size_t &chain_idx  = mcit->second[t];
                    chain_length_idx.push_back(make_pair(singularity_edges[chain_idx].size(),chain_idx));
                    sort(chain_length_idx.begin(),chain_length_idx.end()); // sort the chain accroding the length, leave the shortest, cut off all long chain
                    for(size_t i = 1; i < chain_length_idx.size(); ++i){
                        const size_t & chain_idx_ = chain_length_idx[i].second;
                        for(size_t j = 0; j < singularity_edges[chain_idx_].size(); j +=2){
                            map_changed_edge[make_pair(singularity_edges[chain_idx_][j],
                                                       singularity_edges[chain_idx_][j+1])] = false; // remove the longer singulairty edges
                        }
                    }
                    // label the shortest singularity chain as unknown type
                    for(size_t j = 0; j < singularity_edges[chain_length_idx[0].second].size(); j += 2){
                        map_changed_edge[make_pair(singularity_edges[chain_length_idx[0].second][j],
                                                   singularity_edges[chain_length_idx[0].second][j+1])] = true;
                    }
                }
            }
        }
    }
    size_t t = 0;
    for(map<pair<size_t,size_t>,bool>::const_iterator mci = map_changed_edge.begin();
        mci != map_changed_edge.end(); ++mci){
        if(mci->second == false){
            ++t;
        }
    }
    cerr << "# U fan singularity edges num = " << t << endl;
    update_edges(map_changed_edge,edges,outside_face);
    return 0;
}


#if 0// no used
int relabel_singularity_chain_new(const matrix<size_t> &tet,
                                  const matrix<double> &node,
                                  const matrix<size_t> &outside_face,
                                  one_ring_tet_at_edge &ortae,
                                  vector<deque<size_t> > &singularity_edges,
                                  vector<deque<size_t> > &singularity_type
                                  //map<pair<size_t,size_t>,bool> &map_changed_tetpair_state
                                  )
{
    map<pair<size_t,size_t>,size_t> edges;
    map<pair<size_t,size_t>,bool> map_changed_edge;

    //remove black lines
#if 1
    init_singularity_edges(singularity_edges,singularity_type,edges);
    relabel_remove_black_lines(tet,node,outside_face,singularity_edges,singularity_type,edges,map_changed_edge);
    assemble_singularity_chain(edges,outside_face,singularity_edges,singularity_type,ortae);
#endif

#if 0 // not sure whether it's right in theory
    map_changed_edge.clear();
    init_singularity_edges(singularity_edges,singularity_type,edges);
    relabel_triangle_singualrity_edge(tet,node,outside_face,singularity_edges,edges,map_changed_edge);
    assemble_singularity_chain(edges,outside_face,singularity_edges,singularity_type,ortae);
#endif

    // clean up the short lines
#if 1
    map_changed_edge.clear();
    init_singularity_edges(singularity_edges,singularity_type,edges);
    relabel_remove_short_lines(tet,node,outside_face,singularity_edges,edges,map_changed_edge);
    assemble_singularity_chain(edges,outside_face,singularity_edges,singularity_type,ortae);
#endif
    //    //smooth lines
#if 1
    map_changed_edge.clear();
    init_singularity_edges(singularity_edges,singularity_type,edges);
    relabel_smooth_singulairty_edges(tet,node,outside_face,singularity_edges,edges,map_changed_edge);
    assemble_singularity_chain(edges,outside_face,singularity_edges,singularity_type,ortae);
#endif

    //dump_singularity_chain_to_vtk("test_haha3.vtk",node,singularity_edges,singularity_type);

#if 0
    {
        vector<size_t> lines;//(map_changed_edge.size() * 2);
        vector<double> color;//(map_changed_edge.size());
        typedef map<pair<size_t,size_t>,bool>::const_iterator mcit;
        for(mcit ci = map_changed_edge.begin(); ci != map_changed_edge.end(); ++ci)
        {
            lines.push_back(ci->first.first);
            lines.push_back(ci->first.second);
            if(ci->second)
                color.push_back(1);
            else
                color.push_back(0.5);
        }
        ofstream test_ofs("changed_line.vtk");
        line2vtk(test_ofs,&node[0],node.size(2),&lines[0],lines.size()/2);
        cell_data(test_ofs,&color[0],color.size(),"different");
    }
#endif

    //vector<size_t> tmp;
    //    {// output for cut face
    //        auto_ptr<face2tet_adjacent> fa(face2tet_adjacent::create(tet));
    //        typedef map<pair<size_t,size_t>,bool>::const_iterator mci;
    //        vector<size_t> around_tets;



    //        for(mci it = map_changed_edge.begin(); it != map_changed_edge.end(); ++it){
    //            //find_around_tets(it->first,around_tets,*fa,tet);
    //            find_around_tets_new(it->first,tet,node,*fa,ortae,around_tets);
    //            for(size_t t = 0; t <around_tets.size(); ++t) tmp.push_back(around_tets[t]);
    //            for(size_t t = 0; t < around_tets.size() - 1; ++t){
    //                map_changed_tetpair_state.insert(make_pair(make_pair(around_tets[t],around_tets[t + 1]),it->second));
    //            }
    //        }
#if 0
    ofstream test_tet_ofs("changed_tet.vtk");
    matrix<size_t> tets(4,tmp.size());
    for(size_t t = 0; t < tmp.size(); ++t){
        tets(colon(),t) = tet(colon(),tmp[t]);
    }
    tet2vtk(test_tet_ofs,&node[0],node.size(2),&tets[0],tmp.size());
#endif

    return 0;
}
#endif

int relabel_singularity_chain_new_test(const matrix<size_t> &tet,
                                       const matrix<double> &node,
                                       const matrix<size_t> &outside_face,
                                       one_ring_tet_at_edge &ortae,
                                       vector<deque<size_t> > &singularity_edges,
                                       vector<deque<size_t> > &singularity_type
                                       //map<pair<size_t,size_t>,bool> &map_changed_tetpair_state
                                       )
{
    map<pair<size_t,size_t>,size_t> edges;
    map<pair<size_t,size_t>,bool> map_changed_edge;

    //    //remove black lines
    init_singularity_edges(singularity_edges,singularity_type,edges);
    relabel_remove_black_lines(tet,node,outside_face,singularity_edges,singularity_type,edges,map_changed_edge);
    assemble_singularity_chain(edges,outside_face,singularity_edges,singularity_type,ortae);

    map_changed_edge.clear();
    init_singularity_edges(singularity_edges,singularity_type,edges);
    relabel_remove_short_lines(tet,node,outside_face,singularity_edges,edges,map_changed_edge);
    assemble_singularity_chain(edges,outside_face,singularity_edges,singularity_type,ortae);
    return 0;
}

#if 0// no use
int relabel_singularity_chain(const matrix<size_t> &input_tet,
                              const matrix<double> &node,
                              const matrix<size_t> &outside_face,
                              vector<deque<size_t> > &singularity_edges,
                              vector<deque<size_t> > &singularity_type,
                              map<pair<size_t,size_t>,bool> &map_changed_tetpair_state) // store the tet pair and whether is the around the singularity edge
{
    auto_ptr<face2tet_adjacent> fa(face2tet_adjacent::create(input_tet));
    const size_t black_type = 9;
    size_t begin_,end_;
    map<pair<size_t,size_t>,bool> map_changed_edge;
    for(size_t t = 0; t < singularity_type.size(); ++t)
    {
        if(find(singularity_type[t].begin(), singularity_type[t].end(), black_type)
                !=  singularity_type[t].end()) {
            if(is_outside_vertex(singularity_edges[t].front(),outside_face)
                    && is_outside_vertex(singularity_edges[t].back(),outside_face)){
                cerr << "# strange black line." << endl;
                return __LINE__;
            }
            if(is_outside_vertex(singularity_edges[t].front(),outside_face)) {
                end_ = singularity_edges[t].front();
                begin_ = singularity_edges[t].back();
            } else {
                end_ = singularity_edges[t].back();
                begin_ = singularity_edges[t].front();
            }
            const deque<size_t> &black_chain = singularity_edges[t];
            begin_ = get_new_begin_vertex(singularity_edges,begin_,t);
            if(begin_ == -1) {
                cerr << "# strange: can not find other chain which link this black chain." << endl;
                return __LINE__;
            }
            find_shortest_path(begin_,end_,input_tet,node,black_chain,map_changed_edge);
        }
    }

    // reassign the tetpair state
    {
        typedef map<pair<size_t,size_t>,bool>::const_iterator mci;
        vector<size_t> around_tets;
        for(mci it = map_changed_edge.begin(); it != map_changed_edge.end(); ++it){
            find_around_tets(it->first,around_tets,*fa,input_tet);
            for(size_t t = 0; t < around_tets.size() - 1; ++t){
                map_changed_tetpair_state.insert(make_pair(make_pair(around_tets[t],around_tets[t + 1]),it->second));
            }
        }
    }

    updata_singularity_chain(singularity_edges,singularity_type,map_changed_edge,outside_face);

    { // for vis 2
        map<size_t,vector<size_t> > chain_end_count;

        for(size_t t = 0; t < singularity_edges.size(); ++t){
            chain_end_count[singularity_edges[t].front()].push_back(t);
            chain_end_count[singularity_edges[t].back()].push_back(t);
        }
        set<size_t> nee_to_draw;
        for(map<size_t,vector<size_t> >::const_iterator mit = chain_end_count.begin();
            mit != chain_end_count.end(); ++mit)
        {
            if(mit->second.size() == 2){
                nee_to_draw.insert(mit->second[0]);
                nee_to_draw.insert(mit->second[1]);
            }
        }

    }
    {// for vis
        cerr << "# singularity chain num: " << singularity_edges.size() << endl;

        map<size_t,vector<size_t> > chain_end_count;

        for(size_t t = 0; t < singularity_edges.size(); ++t){
            chain_end_count[singularity_edges[t].front()].push_back(t);
            chain_end_count[singularity_edges[t].back()].push_back(t);
        }
        set<size_t> nee_to_draw;
        vector<size_t> points2;
        for(map<size_t,vector<size_t> >::const_iterator mit = chain_end_count.begin();
            mit != chain_end_count.end(); ++mit)
        {
            if(mit->second.size() == 2){
                nee_to_draw.insert(mit->second[0]);
                nee_to_draw.insert(mit->second[1]);
                points2.push_back(mit->first);
            }
        }

        vector<size_t> lines;
        vector<double> rgba_v;
        double color_table[30][4] = {{255,0,0,255},
                                     {255,80,80,255},
                                     {255,160,160,255},

                                     {0,255,0,255},
                                     {80,255,80,255},
                                     {160,255,160,255},

                                     {0,0,255,255},
                                     {80,80,255,255},
                                     {160,160,255,255},

                                     {0,0,0,255},

                                     {255,0,0,200},
                                     {255,80,80,200},
                                     {255,160,160,200},

                                     {0,255,0,200},
                                     {80,255,80,200},
                                     {160,255,160,200},

                                     {0,0,255,200},
                                     {80,80,255,200},
                                     {160,160,255,200},

                                     {0,0,0,200},

                                     {255,0,0,100},
                                     {255,80,80,100},
                                     {255,160,160,100},

                                     {0,255,0,100},
                                     {80,255,80,100},
                                     {160,255,160,100},

                                     {0,0,255,100},
                                     {80,80,255,100},
                                     {160,160,255,100},

                                     {0,0,0,100}};
        itr_matrix<double*> rgba_m(4,30,&color_table[0][0]);
        rgba_m /= 255.0;

        //        for(size_t t = 0 ;t < singularity_edges.size(); ++t){
        //            for(size_t i = 0; i < singularity_edges[t].size() ; ++i){
        //                lines.push_back(singularity_edges[t][i]);
        //            }
        //        }
        //        for(size_t t = 0; t < singularity_type.size(); ++t)
        //            for(size_t i = 0; i < singularity_type[t].size(); ++i)
        //                for(size_t j =0; j < 4; ++j)
        //                    //rgba_v.push_back(rgba_m(j,singularity_type[t][i]));
        //                    rgba_v.push_back(rgba_m(j,t));

        for(set<size_t>::const_iterator t = nee_to_draw.begin() ;
            t != nee_to_draw.end(); ++t){
            for(size_t i = 0; i < singularity_edges[*t].size() ; ++i){
                lines.push_back(singularity_edges[*t][i]);
            }
            for(size_t i = 0; i < singularity_type[*t].size(); ++i)
                for(size_t j =0; j < 4; ++j)
                    //rgba_v.push_back(rgba_m(j,singularity_type[t][i]));
                    rgba_v.push_back(rgba_m(j,*t ));
        }

        ofstream ofs("new_singularity2.vtk");

        line2vtk(ofs,&node[0],node.size(2),&lines[0],lines.size()/2);
        cell_data_rgba(ofs,&rgba_v[0],rgba_v.size() / 4, "new_singularity");
        ofstream ofs_p("new_points2.vtk");
        point2vtk(ofs_p,&node[0],node.size(2),&points2[0],points2.size());
    }


    return 0;
}
#endif
