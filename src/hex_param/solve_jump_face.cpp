#include "hex_param.h"
#include "find_singularities.h"
#include "../common/graph.h"
#include "../common/zyz.h"
#include "../common/transition.h"
#include "../common/vtk.h"
#include "../common/IO.h"
#include "common.h"

extern "C" {
#include "../spherical_harmonics/rot_cubic_f_SH.h"
}

#include <iostream>
#include <stack>
#include <fstream>
#include <ctime>
#include <limits>
#include <list>

#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/ptree/ptree.h>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/kruskal_min_spanning_tree.hpp>
//#include <boost/graph/dijkstra_shortest_paths.hpp>
//#include <boost/graph/graph_traits.hpp>

using namespace std;
using namespace zjucad::matrix;
using boost::property_tree::ptree;

int get_singularity_edge_type_loop(const size_t from,
                                   const size_t to,
                                   const map<pair<size_t,size_t>,size_t> &singularity_segment,
                                   vector<size_t> &loop)
{
    typedef map<pair<size_t,size_t>,size_t>::const_iterator mci;
    mci ci = singularity_segment.find(make_pair(from,to));
    if(ci != singularity_segment.end())
        return ci->second;
    else{
        ci = singularity_segment.find(make_pair(to,from));
        if(ci != singularity_segment.end()) {
            //reverse(loop.begin(),loop.end());
            size_t type = ci->second;
            if(type == -1) return type; // unknown type

            matrixd rot = type_transition2(type);
            inv(rot);
#if 0
            if(type - type/3 * 3 == 0 ) type += 2;
            else if(type - type/3 * 3 == 2) type -= 2;
            if(type != type_transition1(rot))
                cerr << "# strange the inv rot type is not correct." << endl;
#endif
            return type_transition1(rot);
            //            size_t type = ci->second;
            //            if(type - type/3 * 3 == 0) type += 2;
            //            else if(type - type/3 * 3 == 2) type -= 2;
            //            return type;

        }else
            return 9; // not singularity edge
    }
}

int get_singularity_edge_type(const size_t from,
                              const size_t to,
                              const map<pair<size_t,size_t>,size_t> &singularity_segment)
{
    typedef map<pair<size_t,size_t>,size_t>::const_iterator mci;
    mci ci = singularity_segment.find(make_pair(from,to));
    if(ci != singularity_segment.end())
        return ci->second;
    else{
        ci = singularity_segment.find(make_pair(to,from));
        if(ci != singularity_segment.end()) {
            cerr << "# can not touch here. < " << from << "," << to << ">." << endl;
            size_t type = ci->second;

            if(type == -1) return type; // unknown type

            matrixd rot = type_transition2(type);
            //inv(rot);
#if 0
            if(type - type/3 * 3 == 0 ) type += 2;
            else if(type - type/3 * 3 == 2) type -= 2;
            if(type != type_transition1(rot))
                cerr << "# strange the inv rot type is not correct." << endl;
#endif
            return type_transition1(trans(rot));
            //            size_t type = ci->second;
            //            if(type - type/3 * 3 == 0) type += 2;
            //            else if(type - type/3 * 3 == 2) type -= 2;
            //            return type;

        }else
            return 9; // not singularity edge
    }
}

int solve_jump_face(const matrixst &cut_tet,
                    const matrixd &cut_node,
                    const matrixst &face_pair,
                    matrixst &jump_face_type,
                    const matrixst &outside_face_cut,
                    const matrixst &outside_face_cut_idx,
                    const matrixst &cut_tet2tet,
                    const jtf::mesh::face2tet_adjacent &fa_cut,
                    const vector<deque<pair<size_t,size_t> > > &singularity_chain,
                    const vector<deque<size_t> > &singularity_type)
{
    std::map<pair<size_t,size_t>,size_t> tet_pair_type; // store the <beg tet,end tet idx> and the type
    set<pair<list<pair<size_t,size_t> >,size_t> > unknow_equations;
    enum JUMP_TYPE{
        JUMP_X_1 = 0, JUMP_X_2, JUMP_X_3, JUMP_Y_1, JUMP_Y_2, JUMP_Y_3, JUMP_Z_1, JUMP_Z_2, JUMP_Z_3, JUMP_Identity};

    map<pair<size_t,size_t>,vector<size_t> > outside_edge_face_pair; // store the face_pair <left,right> (right hand) beside the edge <e0,e1>

    for(size_t t = 0; t < outside_face_cut.size(2); ++t){
        for(size_t i = 0; i < 3; ++i){
            pair<size_t,size_t> edge;
            edge.first = outside_face_cut(i, t);
            edge.second = outside_face_cut((i+1)%3, t);
            if(edge.first > edge.second) swap(edge.first,edge.second);
            outside_edge_face_pair[edge].push_back(t);
        }
    }
    typedef map<pair<size_t,size_t>,vector<size_t> >::iterator mit;

#if 1 // check the valid
    for(mit ci = outside_edge_face_pair.begin(); ci != outside_edge_face_pair.end(); ++ci){
        if(ci->second.size() != 2) {
            cerr << "# strange outside edge "
                 << ci->first.first << " "
                 << ci->first.second << " "
                 << ci->second.size() << endl;
        }
    }
#endif
    for(mit it = outside_edge_face_pair.begin(); it != outside_edge_face_pair.end(); ++it){
        const matrixst &face0 = outside_face_cut(colon(),it->second[0]);
        size_t other_vertex;
        for(size_t i = 0; i < 3; ++i){
            if(face0[i] != it->first.first && face0[i] != it->first.second){
                other_vertex = face0[i];
            }
        }
        const pair<size_t,size_t>& tet_pair = fa_cut.face2tet_[outside_face_cut_idx[it->second[0]]];
        const size_t tet_idx = (tet_pair.first == -1)?tet_pair.second:tet_pair.first;
        matrixd normal = cross(cut_node(colon(),it->first.second) - cut_node(colon(),it->first.first),
                                      cut_node(colon(),other_vertex) - cut_node(colon(),it->first.second));
        matrixd face_bary = (cut_node(colon(),face0[0]) + cut_node(colon(),face0[1]) + cut_node(colon(),face0[2]))/3;
        matrixd tet_bary = (  cut_node(colon(),cut_tet(0,tet_idx))
                                   + cut_node(colon(),cut_tet(1,tet_idx))
                                   + cut_node(colon(),cut_tet(2,tet_idx))
                                   + cut_node(colon(),cut_tet(3,tet_idx)) ) /4.0;
        if(dot(normal,face_bary - tet_bary) < -1e-8) {
            swap(it->second[0],it->second[1]); // revert the face order
        }
    }

    map<pair<size_t,size_t>,size_t> singularity_segment; // stores the singulariy edge and type
    {// convert the singularity edges into segments
        for(size_t t = 0; t < singularity_chain.size(); ++t){
            for(size_t i = 0; i < singularity_chain[t].size(); ++i){
                singularity_segment.insert(make_pair(singularity_chain[t][i],singularity_type[t][i]));
            }
        }
    }

    for(mit it = outside_edge_face_pair.begin();
        it != outside_edge_face_pair.end(); ++it){
        if(face_pair[it->second[0]] != -1
                && face_pair[it->second[1]] != -1){ // both of the two face beside the edge are jump faces, we need to set equation
            int rtn = get_singularity_edge_type(cut_tet2tet[it->first.first],cut_tet2tet[it->first.second],singularity_segment);
            if(rtn == -1) // not singularity edge
            {
                list<pair<size_t,size_t> > list_;
                size_t left_tet_idx,left_tet_pair_idx;
                size_t right_tet_idx,right_tet_pair_tdx;

                const pair<size_t,size_t> &tet_pair_left = fa_cut.face2tet_[outside_face_cut_idx[it->second[0]]];
                left_tet_idx = (tet_pair_left.first == -1)?tet_pair_left.second:tet_pair_left.first;
                const pair<size_t,size_t> &tet_pair_right = fa_cut.face2tet_[outside_face_cut_idx[it->second[1]]];
                right_tet_idx = (tet_pair_right.first == -1)?tet_pair_right.second:tet_pair_right.first;

                const pair<size_t,size_t> &tet_pair_left_pair = fa_cut.face2tet_[outside_face_cut_idx[face_pair[it->second[0]]]];
                left_tet_pair_idx = (tet_pair_left_pair.first == -1)?tet_pair_left_pair.second:tet_pair_left_pair.first;
                const pair<size_t,size_t> &tet_pair_right_pair = fa_cut.face2tet_[outside_face_cut_idx[face_pair[it->second[1]]]];
                right_tet_pair_tdx = (tet_pair_right_pair.first == -1)?tet_pair_right_pair.second:tet_pair_right_pair.first;

                list_.push_back(make_pair(left_tet_idx,left_tet_pair_idx));
                list_.push_back(make_pair(right_tet_pair_tdx,right_tet_idx));
                unknow_equations.insert(make_pair(list_, static_cast<size_t>(JUMP_Identity)));
            }else{ // the singularity edge is that the left face and right is also the face pair
                list<pair<size_t,size_t> > list_;
                size_t left_tet_idx,left_tet_pair_idx;
                size_t right_tet_idx,right_tet_pair_idx;

                const pair<size_t,size_t> &tet_pair_left = fa_cut.face2tet_[outside_face_cut_idx[it->second[0]]];
                left_tet_idx = (tet_pair_left.first == -1)?tet_pair_left.second:tet_pair_left.first;
                const pair<size_t,size_t> &tet_pair_right = fa_cut.face2tet_[outside_face_cut_idx[it->second[1]]];
                right_tet_idx = (tet_pair_right.first == -1)?tet_pair_right.second:tet_pair_right.first;

                const pair<size_t,size_t> &tet_pair_left_pair = fa_cut.face2tet_[outside_face_cut_idx[face_pair[it->second[0]]]];
                left_tet_pair_idx = (tet_pair_left_pair.first == -1)?tet_pair_left_pair.second:tet_pair_left_pair.first;
                const pair<size_t,size_t> &tet_pair_right_pair = fa_cut.face2tet_[outside_face_cut_idx[face_pair[it->second[1]]]];
                right_tet_pair_idx = (tet_pair_right_pair.first == -1)?tet_pair_right_pair.second:tet_pair_right_pair.first;

                //assert(left_tet_pair_idx == right_tet_idx && right_tet_pair_tdx == left_tet_idx);
                //                cerr << "left_idx " << left_tet_idx << " right_idx " << right_tet_idx << endl;
                //                cerr << "left_pair_idx " << left_tet_pair_idx << " right_pair_idx " << right_tet_pair_idx << endl;

                tet_pair_type[make_pair(left_tet_idx,left_tet_pair_idx)] = rtn;
                size_t type = rtn;
                if(type -type/3*3  == 0) type += 2;
                else if(type - type/3 * 3 == 2) type -= 2;
                tet_pair_type[make_pair(left_tet_pair_idx,left_tet_idx)] = type;
            }
        }
    }// finish set up the equations

    // to solve the equations
    typedef set<pair<list<pair<size_t,size_t> >,size_t> >::iterator sit;
    typedef map<pair<size_t,size_t>,size_t >::const_iterator mpcit;
    typedef list<pair<size_t,size_t> >::iterator lit;
    while(!unknow_equations.empty())
    {
        for(sit it = unknow_equations.begin(); it != unknow_equations.end(); ++it){
            list<pair<size_t,size_t> > list_ = it->first;
            for(lit lit_ = list_.begin(); lit_ != list_.end(); ++lit_){
                mpcit ci = tet_pair_type.find(*lit_);
                if(ci == tet_pair_type.end()) continue;
                if(ci->second == static_cast<size_t>(JUMP_Identity)){
                    list_.erase(lit_);
                    lit_ = list_.begin(); // need to find the other tet_pair in this list
                    assert(list_.size() == 1);
                    tet_pair_type[*lit_] = static_cast<size_t>(JUMP_Identity); // add the left tet pair
                    tet_pair_type[make_pair(lit_->second,lit_->first)] = static_cast<size_t>(JUMP_Identity);
                    unknow_equations.erase(it); // delete this equation
                    break;
                }else{
                    size_t type = ci->second;
                    if(type - type/3 * 3 == 0) type += 2;
                    else if(type - type/3 * 3 == 2) type -= 2;
                    list_.erase(lit_);
                    lit_ = list_.begin();
                    tet_pair_type[*lit_] = type;
                    tet_pair_type[make_pair(lit_->second,lit_->first)] = ci->second;
                    unknow_equations.erase(it);
                    break;
                }
            }
        }
    }// end solve all the jump face

    jump_face_type = -1 * ones<size_t>(face_pair.size());
    for(size_t t = 0; t < face_pair.size(); ++t){
        if(face_pair[t] != -1){
            const pair<size_t,size_t> &tet_pair = fa_cut.face2tet_[outside_face_cut_idx[t]];
            const size_t tet_idx = (tet_pair.first == -1)?tet_pair.second:tet_pair.first;

            const pair<size_t,size_t> &tet_pair_other = fa_cut.face2tet_[outside_face_cut_idx[face_pair[t]]];
            const size_t tet_idx_other = (tet_pair_other.first == -1)?tet_pair_other.second:tet_pair_other.first;
            jump_face_type[t] = tet_pair_type[make_pair(tet_idx,tet_idx_other)];
        }
    }
    return 0;
}

int solve_jump_face_using_zyz(const matrixst &cut_tet,
                              const matrixd &cut_node,
                              const matrixst &face_pair,
                              matrixst &jump_face_type,
                              const matrixst &outside_face_cut,
                              const matrixst &outside_face_cut_idx,
                              const matrixst &cut_tet2tet,
                              const jtf::mesh::face2tet_adjacent &fa_cut,
                              const std::vector<std::deque<std::pair<size_t, size_t> > > &singularity_chain,
                              const std::vector<std::deque<size_t> > &singularity_type,
                              const zjucad::matrix::matrix<matrixd > &frame_inner)
{
    std::map<pair<size_t,size_t>,size_t> tet_pair_type; // store the <beg tet,end tet idx> and the type
    set<pair<list<pair<size_t,size_t> >,size_t> > unknow_equations;
    enum JUMP_TYPE{
        JUMP_X_1 = 0, JUMP_X_2, JUMP_X_3, JUMP_Y_1, JUMP_Y_2, JUMP_Y_3, JUMP_Z_1, JUMP_Z_2, JUMP_Z_3, JUMP_Identity};

    map<pair<size_t,size_t>,vector<size_t> > outside_edge_face_pair; // store the face_pair <left,right> (right hand) beside the edge <e0,e1>

    for(size_t t = 0; t < outside_face_cut.size(2); ++t){
        for(size_t i = 0; i < 3; ++i){
            pair<size_t,size_t> edge;
            edge.first = outside_face_cut(i, t);
            edge.second = outside_face_cut((i+1)%3, t);
            if(edge.first > edge.second) swap(edge.first,edge.second);
            outside_edge_face_pair[edge].push_back(t);
        }
    }
    typedef map<pair<size_t,size_t>,vector<size_t> >::iterator mit;

#if 1 // check the valid
    for(mit ci = outside_edge_face_pair.begin(); ci != outside_edge_face_pair.end(); ++ci){
        if(ci->second.size() != 2) {
            cerr << "# strange outside edge "
                 << ci->first.first << " "
                 << ci->first.second << " "
                 << ci->second.size() << endl;
        }
    }
#endif
    for(mit it = outside_edge_face_pair.begin(); it != outside_edge_face_pair.end(); ++it){
        const matrixst &face0 = outside_face_cut(colon(),it->second[0]);
        size_t other_vertex;
        for(size_t i = 0; i < 3; ++i){
            if(face0[i] != it->first.first && face0[i] != it->first.second){
                other_vertex = face0[i];
            }
        }
        const pair<size_t,size_t>& tet_pair = fa_cut.face2tet_[outside_face_cut_idx[it->second[0]]];
        const size_t tet_idx = (tet_pair.first == -1)?tet_pair.second:tet_pair.first;
        matrixd normal = cross(cut_node(colon(),it->first.second) - cut_node(colon(),it->first.first),
                                      cut_node(colon(),other_vertex) - cut_node(colon(),it->first.second));
        matrixd face_bary = (cut_node(colon(),face0[0]) + cut_node(colon(),face0[1]) + cut_node(colon(),face0[2]))/3;
        matrixd tet_bary = (  cut_node(colon(),cut_tet(0,tet_idx))
                                   + cut_node(colon(),cut_tet(1,tet_idx))
                                   + cut_node(colon(),cut_tet(2,tet_idx))
                                   + cut_node(colon(),cut_tet(3,tet_idx)) ) /4.0;
        if(dot(normal,face_bary - tet_bary) < -1e-8) {
            swap(it->second[0],it->second[1]); // revert the face order
        }
    }

    map<pair<size_t,size_t>,size_t> singularity_segment; // stores the singulariy edge and type
    {// convert the singularity edges into segments
        for(size_t t = 0; t < singularity_chain.size(); ++t){
            for(size_t i = 0; i < singularity_chain[t].size(); ++i){
                singularity_segment.insert(make_pair(singularity_chain[t][i],singularity_type[t][i]));
            }
        }
    }

    for(mit it = outside_edge_face_pair.begin();
        it != outside_edge_face_pair.end(); ++it){
        if(face_pair[it->second[0]] != -1
                && face_pair[it->second[1]] != -1){ // both of the two face beside the edge are jump faces, we need to set equation
            int rtn = get_singularity_edge_type(cut_tet2tet[it->first.first],cut_tet2tet[it->first.second],singularity_segment);
            if(rtn == -1) // not singularity edge
            {
                list<pair<size_t,size_t> > list_;
                size_t left_tet_idx,left_tet_pair_idx;
                size_t right_tet_idx,right_tet_pair_tdx;

                const pair<size_t,size_t> &tet_pair_left = fa_cut.face2tet_[outside_face_cut_idx[it->second[0]]];
                left_tet_idx = (tet_pair_left.first == -1)?tet_pair_left.second:tet_pair_left.first;
                const pair<size_t,size_t> &tet_pair_right = fa_cut.face2tet_[outside_face_cut_idx[it->second[1]]];
                right_tet_idx = (tet_pair_right.first == -1)?tet_pair_right.second:tet_pair_right.first;

                const pair<size_t,size_t> &tet_pair_left_pair = fa_cut.face2tet_[outside_face_cut_idx[face_pair[it->second[0]]]];
                left_tet_pair_idx = (tet_pair_left_pair.first == -1)?tet_pair_left_pair.second:tet_pair_left_pair.first;
                const pair<size_t,size_t> &tet_pair_right_pair = fa_cut.face2tet_[outside_face_cut_idx[face_pair[it->second[1]]]];
                right_tet_pair_tdx = (tet_pair_right_pair.first == -1)?tet_pair_right_pair.second:tet_pair_right_pair.first;

                list_.push_back(make_pair(left_tet_idx,left_tet_pair_idx));
                list_.push_back(make_pair(right_tet_pair_tdx,right_tet_idx));
                unknow_equations.insert(make_pair(list_, static_cast<size_t>(JUMP_Identity)));
            }else{ // the singularity edge is that the left face and right is also the face pair
                list<pair<size_t,size_t> > list_;
                size_t left_tet_idx,left_tet_pair_idx;
                size_t right_tet_idx,right_tet_pair_idx;

                const pair<size_t,size_t> &tet_pair_left = fa_cut.face2tet_[outside_face_cut_idx[it->second[0]]];
                left_tet_idx = (tet_pair_left.first == -1)?tet_pair_left.second:tet_pair_left.first;
                const pair<size_t,size_t> &tet_pair_right = fa_cut.face2tet_[outside_face_cut_idx[it->second[1]]];
                right_tet_idx = (tet_pair_right.first == -1)?tet_pair_right.second:tet_pair_right.first;

                const pair<size_t,size_t> &tet_pair_left_pair = fa_cut.face2tet_[outside_face_cut_idx[face_pair[it->second[0]]]];
                left_tet_pair_idx = (tet_pair_left_pair.first == -1)?tet_pair_left_pair.second:tet_pair_left_pair.first;
                const pair<size_t,size_t> &tet_pair_right_pair = fa_cut.face2tet_[outside_face_cut_idx[face_pair[it->second[1]]]];
                right_tet_pair_idx = (tet_pair_right_pair.first == -1)?tet_pair_right_pair.second:tet_pair_right_pair.first;


                matrixd rot_left(3,3),rot_right(3,3);
                get_best_alignment(&frame_inner[left_tet_idx][0],&frame_inner[left_tet_pair_idx][0],&rot_left[0]);
                get_best_alignment(&frame_inner[right_tet_pair_idx][0],&frame_inner[right_tet_idx][0],&rot_right[0]);

                if((norm(rot_left - eye<double>(3)) < 1e-8
                    && norm(rot_right - eye<double>(3)) < 1e-8)
                        || (norm(rot_left - eye<double>(3)) > 1e-8
                            && norm(rot_right - eye<double>(3)) > 1e-8)){ // this edge is labelled as singularity edges, the original

                    list<pair<size_t,size_t> > list_;
                    list_.push_back(make_pair(left_tet_idx,left_tet_pair_idx));
                    list_.push_back(make_pair(right_tet_pair_idx,right_tet_idx));
                    unknow_equations.insert(make_pair(list_,rtn));
                }else{
                    if(norm(rot_left - eye<double>(3)) < 1e-8 && norm(rot_right - eye<double>(3)) > 1e-8){
                        assert(rtn != 9 && rtn!= -1);
                        tet_pair_type[make_pair(left_tet_idx,left_tet_pair_idx)] = static_cast<size_t>(JUMP_Identity);
                        tet_pair_type[make_pair(left_tet_pair_idx,left_tet_idx)] = static_cast<size_t>(JUMP_Identity);
                        tet_pair_type[make_pair(right_tet_pair_idx,right_tet_idx)] = rtn;
                        size_t type = rtn;
                        if(type - type/3 * 3 == 0) type += 2;
                        else if(type - type/3 * 3 == 2) type -= 2;
                        tet_pair_type[make_pair(right_tet_idx,right_tet_pair_idx)] = type;
                        assert(type != 9 && type != -1);
                    }else if(norm(rot_right - eye<double>(3)) < 1e-8 && norm(rot_left - eye<double>(3)) > 1e-8){
                        assert(rtn != 9 && rtn!= -1);
                        tet_pair_type[make_pair(right_tet_idx,right_tet_pair_idx)] = static_cast<size_t>(JUMP_Identity);
                        tet_pair_type[make_pair(right_tet_pair_idx,right_tet_idx)] = static_cast<size_t>(JUMP_Identity);
                        tet_pair_type[make_pair(left_tet_idx,left_tet_pair_idx)] = rtn;
                        size_t type = rtn;
                        if(type - type/3 * 3 == 0) type += 2;
                        else if(type - type/3 * 3 == 2) type -= 2;
                        assert(type != 9 && type != -1);
                        tet_pair_type[make_pair(left_tet_pair_idx,left_tet_idx)] = type;
                    }
                }
                //assert(left_tet_pair_idx == right_tet_idx && right_tet_pair_tdx == left_tet_idx);
                //cerr << "left_idx " << left_tet_idx << " right_idx " << right_tet_idx << endl;
                //cerr << "left_pair_idx " << left_tet_pair_idx << " right_pair_idx " << right_tet_pair_idx << endl;

                //                assert(rtn != 9 && rtn != -1);
                //                tet_pair_type[make_pair(left_tet_idx,left_tet_pair_idx)] = rtn;
                //                size_t type = rtn;
                //                if(type -type/3*3  == 0) type += 2;
                //                else if(type - type/3 * 3 == 2) type -= 2;
                //                assert(type != 9 && type != -1);
                //                tet_pair_type[make_pair(left_tet_pair_idx,left_tet_idx)] = type;
            }
        }
    }// finish set up the equations

    // to solve the equations
    typedef set<pair<list<pair<size_t,size_t> >,size_t> >::iterator sit;
    typedef map<pair<size_t,size_t>,size_t >::const_iterator mpcit;
    typedef list<pair<size_t,size_t> >::iterator lit;
    while(!unknow_equations.empty())
    {
        for(sit it = unknow_equations.begin(); it != unknow_equations.end(); ++it){
            list<pair<size_t,size_t> > list_ = it->first;
            for(lit lit_ = list_.begin(); lit_ != list_.end(); ++lit_){
                mpcit ci = tet_pair_type.find(*lit_);
                if(ci == tet_pair_type.end()) continue;
                if(ci->second == static_cast<size_t>(JUMP_Identity)){
                    list_.erase(lit_);
                    lit_ = list_.begin(); // need to find the other tet_pair in this list
                    assert(list_.size() == 1);
                    tet_pair_type[*lit_] = static_cast<size_t>(JUMP_Identity); // add the left tet pair
                    tet_pair_type[make_pair(lit_->second,lit_->first)] = static_cast<size_t>(JUMP_Identity);
                    unknow_equations.erase(it); // delete this equation
                    break;
                }else{
                    size_t type = ci->second;
                    if(type - type/3 * 3 == 0) type += 2;
                    else if(type - type/3 * 3 == 2) type -= 2;
                    list_.erase(lit_);
                    lit_ = list_.begin();
                    tet_pair_type[*lit_] = type;
                    tet_pair_type[make_pair(lit_->second,lit_->first)] = ci->second;
                    unknow_equations.erase(it);
                    break;
                }
            }
        }
        cerr << "# unknown equations size " << unknow_equations.size() << endl;
    }// end solve all the jump face

#if 1 // check the jump face type
    for(mpcit ci = tet_pair_type.begin(); ci != tet_pair_type.end(); ++ci){
        if(ci->second == 9 || ci->second == -1){
            cerr << "# strange jump face, tet pair: <" << ci->first.first << "," << ci->first.second << ">" << " type " << ci->second << endl;
        }
    }
#endif

    jump_face_type = -1 * ones<size_t>(face_pair.size());
    for(size_t t = 0; t < face_pair.size(); ++t){
        if(face_pair[t] != -1){
            const pair<size_t,size_t> &tet_pair = fa_cut.face2tet_[outside_face_cut_idx[t]];
            const size_t tet_idx = (tet_pair.first == -1)?tet_pair.second:tet_pair.first;

            const pair<size_t,size_t> &tet_pair_other = fa_cut.face2tet_[outside_face_cut_idx[face_pair[t]]];
            const size_t tet_idx_other = (tet_pair_other.first == -1)?tet_pair_other.second:tet_pair_other.first;
            jump_face_type[t] = tet_pair_type[make_pair(tet_idx,tet_idx_other)];
            assert(jump_face_type[t] != 9 && jump_face_type[t] != -1);
        }
    }

    return 0;
}

int order_try_list(vector<size_t> &try_list,
                   const size_t MAX_JUMP_TYPE,
                   const size_t begin_try) // frome the begin to end if each item > 9 , then plus 1 to the behind and cut itself to 0
{
    for(size_t t = 0; t < try_list.size() - 1; ++t){
        if(try_list[t] > MAX_JUMP_TYPE) {
            try_list[t + 1] += try_list[t] / (MAX_JUMP_TYPE + 1);
            try_list[t] = try_list[t] % (MAX_JUMP_TYPE+1) + begin_try;
        }
    }
    return 0;
}

// returns: 0 resolve all equations , 1 will introduce conflict equations, 2 can not resolve all equations
int solve_equation(vector<vector<size_t> > &equation_variable_tmp,
                   vector<pair<list<pair<size_t,size_t> >,size_t> > &unknow_equations_tmp,
                   std::vector<bool> &equation_state_tmp,
                   std::map<pair<size_t,size_t>,size_t> &tet_pair_type_tmp,
                   vector<size_t> &try_list,
                   vector<pair<size_t,size_t> > &possible_variable)
{
    matrixd rot(3,3);
    for(size_t t = 0; t < possible_variable.size(); ++t){
        tet_pair_type_tmp[possible_variable[t]] = try_list[t];
        rot = type_transition2(try_list[t]);
        inv(rot);
        tet_pair_type_tmp[make_pair(possible_variable[t].second,possible_variable[t].first)] = type_transition1(rot);
    }

    typedef set<pair<list<pair<size_t,size_t> >,size_t> >::iterator sit;
    typedef list<pair<size_t,size_t> >::iterator lit;
    typedef map<pair<size_t,size_t>,size_t>::const_iterator mcit;
    typedef map<pair<size_t,size_t>,size_t>::iterator mit;

    size_t unresolved_equations = 0;
    while(find(equation_state_tmp.begin(),equation_state_tmp.end(),false) != equation_state_tmp.end()){ // still exist unresolved equation
        for(size_t t = 0; t < unknow_equations_tmp.size(); ++t){
            if(equation_state_tmp[t]) continue; // if the equation has been resolved
            list<pair<size_t,size_t> > &list_ = unknow_equations_tmp[t].first;

            for(size_t i = 0; i < equation_variable_tmp[t].size(); ++i) // for unresolved equation variable
            {
                if(equation_variable_tmp[t][i] != -1) continue; // !=-1 means this variable type is known
                lit lit_ = list_.begin();
                for(size_t j = 0; j < i; ++j) lit_++; // lit_ = list_.begin() + i
                mcit mcit_ = tet_pair_type_tmp.find(*lit_);
                if(mcit_ == tet_pair_type_tmp.end()) mcit_ = tet_pair_type_tmp.find(make_pair(lit_->second,lit_->first));
                if(mcit_ == tet_pair_type_tmp.end()) continue;

                size_t type = mcit_->second;
                if(mcit_->first == *lit_)
                    equation_variable_tmp[t][i] = type;
                else{
                    cerr << "# can not touch here." << endl;
                    matrixd rot = type_transition2(type);
                    inv(rot);
                    equation_variable_tmp[t][i] = type_transition1(rot);
                }
            }

            size_t count_unknown = count(equation_variable_tmp[t].begin(),equation_variable_tmp[t].end(),-1);

            if(count_unknown == 0) {
                // check whther the equation is correct
                matrixd rot = eye<double>(3);
                for(size_t i = 0; i < equation_variable_tmp[t].size(); ++i){
                    rot = temp(rot * type_transition2(equation_variable_tmp[t][i]));
                }
                if(norm(rot - type_transition2(unknow_equations_tmp[t].second)) < 1e-8){
                    equation_state_tmp[t] = true;
                    continue;
                }else{
                    cerr << "# error it's a conflict equation." << endl;
                    return 1;
                }
            }
            if(count_unknown  == 1){ // only one unknown variable, can be resolved
                matrixd left = eye<double>(3);
                matrixd right = eye<double>(3);
                size_t i = 0;
                for(; i < equation_variable_tmp[t].size(); ++i){
                    if(equation_variable_tmp[t][i] == -1) break; // find the unknown variable
                    left = temp(left * type_transition2(equation_variable_tmp[t][i]));
                }
                const size_t variable_idx = i;
                for(i = variable_idx + 1; i < equation_variable_tmp[t].size(); ++i){
                    right = temp( right * type_transition2(equation_variable_tmp[t][i]));
                }

                // to solve the left * X * right = Type ==> X = (left)^-1 * TYPE * (right)^-1
                matrixd Type = type_transition2(unknow_equations_tmp[t].second);
                inv(left);
                inv(right);
                matrixd itef_matrix = left * Type;
                itef_matrix *= right;
                size_t type_ = type_transition1(itef_matrix);
#if 0 // test
                matrixd test = eye<double>(3);
                inv(left);
                inv(right);
                test = temp(left * itef_matrix) * right ;
                if(type_transition1(test) != unknow_equations_tmp[t].second)
                    cerr << "# calculate error." << endl;
#endif
                lit lit_ = list_.begin();
                for(size_t i = 0; i < variable_idx; ++i) ++lit_;
                tet_pair_type_tmp[*lit_] = type_;
                inv(itef_matrix);
                tet_pair_type_tmp[make_pair(lit_->second,lit_->first)] = type_transition1(itef_matrix);
                equation_state_tmp[t] = true;
            }
        }
        if(unresolved_equations == count(equation_state_tmp.begin(),equation_state_tmp.end(),false)) return 2; // can not resolve all equations
        unresolved_equations = count(equation_state_tmp.begin(),equation_state_tmp.end(),false);
        cerr << "unresolved equations num " << unresolved_equations << endl;
    }
    return 0;
}

struct equation_state_str{
    vector<vector<size_t> > equation_variable_tmp;
    vector<pair<list<pair<size_t,size_t> >,size_t> > unknow_equations_tmp;
    std::vector<bool> equation_state_tmp;
    std::map<pair<size_t,size_t>,size_t> tet_pair_type_tmp;
    vector<pair<size_t,size_t> > possible_variable_tmp;
    vector<bool> possible_variable_resolved_tmp;
    vector<size_t> possible_variable_tried;
    size_t possible_variable_type;
    //size_t current_possible_variable_idx;
};

int solve_equation_state(equation_state_str & tmp)
{
    typedef set<pair<list<pair<size_t,size_t> >,size_t> >::iterator sit;
    typedef list<pair<size_t,size_t> >::iterator lit;
    typedef map<pair<size_t,size_t>,size_t>::const_iterator mcit;
    typedef map<pair<size_t,size_t>,size_t>::iterator mit;

    const size_t idx = tmp.possible_variable_tried.back();
    tmp.tet_pair_type_tmp[tmp.possible_variable_tmp[idx]] = tmp.possible_variable_type;
    matrixd rot;
    rot = type_transition2(tmp.possible_variable_type);
    inv(rot);
    tmp.tet_pair_type_tmp[make_pair(tmp.possible_variable_tmp[idx].second,tmp.possible_variable_tmp[idx].first)]
            = type_transition1(rot);


    //vector<size_t> has_been_resolved_possible_variable;
    size_t unresolved_equations = 0;
    while(find(tmp.equation_state_tmp.begin(),tmp.equation_state_tmp.end(),false) != tmp.equation_state_tmp.end()){ // still exist unresolved equation
        for(size_t t = 0; t < tmp.unknow_equations_tmp.size(); ++t){
            if(tmp.equation_state_tmp[t]) continue; // if the equation has been resolved
            list<pair<size_t,size_t> > &list_ = tmp.unknow_equations_tmp[t].first;

            for(size_t i = 0; i < tmp.equation_variable_tmp[t].size(); ++i) // for unresolved equation variable
            {
                if(tmp.equation_variable_tmp[t][i] != -1) continue; // !=-1 means this variable type is known
                lit lit_ = list_.begin();
                for(size_t j = 0; j < i; ++j) lit_++; // lit_ = list_.begin() + i
                mcit mcit_ = tmp.tet_pair_type_tmp.find(*lit_);
                if(mcit_ == tmp.tet_pair_type_tmp.end()) continue;

                size_t type = mcit_->second;
                tmp.equation_variable_tmp[t][i] = type;
            }

            size_t count_unknown = count(tmp.equation_variable_tmp[t].begin(),tmp.equation_variable_tmp[t].end(),-1);

            if(count_unknown == 0) {
                // check whther the equation is correct
                matrixd rot = eye<double>(3);
                for(size_t i = 0; i < tmp.equation_variable_tmp[t].size(); ++i){
                    rot = temp(rot * type_transition2(tmp.equation_variable_tmp[t][i]));
                }
                if(norm(rot - type_transition2(tmp.unknow_equations_tmp[t].second)) < 1e-8){
                    tmp.equation_state_tmp[t] = true;
                    continue;
                }else{
                    cerr << "# error it's a conflict equation." << endl;
                    return 1;
                }
            }
            if(count_unknown  == 1){ // only one unknown variable, can be resolved
                matrixd left = eye<double>(3);
                matrixd right = eye<double>(3);

                const vector<size_t>::const_iterator unknown_variable_it
                        = find(tmp.equation_variable_tmp[t].begin(),tmp.equation_variable_tmp[t].end(),-1); // find the unknown variable iterator
                assert(unknown_variable_it != tmp.equation_variable_tmp[t].end());

                lit lit_ = list_.begin();
                for(vector<size_t>::const_iterator it = tmp.equation_variable_tmp[t].begin();
                    it != unknown_variable_it; ++it,++lit_){
                    //if(equation_variable[t][i] == -1) break;
                    left = temp( left * type_transition2(*it));
                }
                //const size_t variable_idx = i;
                for(vector<size_t>::const_iterator it = unknown_variable_it + 1;
                    it != tmp.equation_variable_tmp[t].end(); ++it){
                    right = temp( right * type_transition2(*it));
                }

                // to solve the left * X * right = Type ==> X = (left)^-1 * TYPE * (right)^-1
                matrixd Type = type_transition2(tmp.unknow_equations_tmp[t].second);
                //inv(left);
                //inv(right);
                matrixd itef_matrix = trans(left) * Type;
                itef_matrix = temp(itef_matrix * trans(right));
                size_t type_ = type_transition1(itef_matrix);
#if 0 // test
                matrixd test = eye<double>(3);
                inv(left);
                inv(right);
                test = temp(left * itef_matrix) * right ;
                if(type_transition1(test) != unknow_equations_tmp[t].second)
                    cerr << "# calculate error." << endl;
#endif

                mcit cit_ = tmp.tet_pair_type_tmp.find(*lit_);
                if(cit_ == tmp.tet_pair_type_tmp.end()){
                    tmp.tet_pair_type_tmp[*lit_] = type_;
                    //inv(itef_matrix);
                    tmp.tet_pair_type_tmp[make_pair(lit_->second,lit_->first)] = type_transition1(trans(itef_matrix));
                    tmp.equation_state_tmp[t] = true;
                }else
                    if(cit_->second != type_) {
                        cerr << "# error it's a conflict equation. the exist variable does not equals this one." << endl;
                        return 1; // if the resolved variable is not equals that in map, means meat a conflict
                    }
            }
        }
        if(unresolved_equations == count(tmp.equation_state_tmp.begin(),tmp.equation_state_tmp.end(),false)) {
            {//updata the possible variable, those have been resolbed need to be remove
                typedef map<pair<size_t,size_t>,size_t>::const_iterator mci;
                for(size_t t = 0; t < tmp.possible_variable_tmp.size();++t){
                    mci it = tmp.tet_pair_type_tmp.find(tmp.possible_variable_tmp[t]);
                    if(it != tmp.tet_pair_type_tmp.end()){ // it's already unknown
                        tmp.possible_variable_resolved_tmp[t] = true;
                    }
                }
            }
            return 2; // can not resolve all equations
        }

        unresolved_equations = count(tmp.equation_state_tmp.begin(),tmp.equation_state_tmp.end(),false);
        cerr << "unresolved equations num " << unresolved_equations << endl;
    }


    {//updata the possible variable, those have been resolbed need to be remove
        // vector<pair<size_t,size_t> > left_possible_variable;
        typedef map<pair<size_t,size_t>,size_t>::const_iterator mci;
        for(size_t t = 0; t < tmp.possible_variable_tmp.size();++t){
            mci it = tmp.tet_pair_type_tmp.find(tmp.possible_variable_tmp[t]);
            if(it != tmp.tet_pair_type_tmp.end()){ // still unknown
                tmp.possible_variable_resolved_tmp[t] = true;
            }
        }
        //swap(left_possible_variable,tmp.possible_variable_tmp);
    }
    return 0;
}

int solve_jump_face_using_zyz_new(const matrixst &cut_tet,
                                  const matrixd &cut_node,
                                  const matrixst &tet,
                                  const matrixd &node,
                                  const matrixst &face_pair,
                                  matrixst &jump_face_type,
                                  const matrixst &outside_face_cut,
                                  const matrixst &outside_face_cut_idx,
                                  const matrixst &cut_tet2tet,
                                  const jtf::mesh::face2tet_adjacent &fa_cut,
                                  const jtf::mesh::face2tet_adjacent &fa,
                                  map<pair<size_t,size_t>,size_t> &singularity_segment,
                                  //                                  const std::vector<std::deque<std::pair<size_t, size_t> > > &singularity_chain,
                                  //                                  const std::vector<std::deque<size_t> > &singularity_type,
                                  const zjucad::matrix::matrix<matrixd > &frame_inner,
                                  const jtf::mesh::one_ring_tet_at_edge &ortae,
                                  const ptree& pt,
                                  std::map<pair<size_t,size_t>,size_t> &tet_pair_type)
{

#if 0 // check the singularity segment is right or wrong
    typedef map<pair<size_t,size_t>,size_t>::const_iterator mci;
    typedef one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
    for(mci it = singularity_segment.begin(); it != singularity_segment.end(); ++it){
        if(it->second != -1) // type is known, need to check whether is right
        {
            oecit cit = ortae.e2t_.find(it->first);
            if(cit == ortae.e2t_.end()) {
                cerr << "# strange. the order should be the same." << endl;
            }
            const vector<size_t> &loop = cit->second;
            matrixd rot_ = eye<double>(3);
            matrixd transit = eye<double>(3);
            for(size_t t = 0; t < loop.size() - 1; ++t){
                get_best_alignment(&frame_inner[loop[t]](0,0),&frame_inner[loop[t+1]](0,0),&rot_[0]);
                transit = temp(transit * rot_);
            }
            if(type_transition1(transit) != it->second){
                cerr << " # the input singularity type is wrong." << endl;
                cerr << " # the original singularity type is  " << it->second << " the calculated type is " << type_transition1(transit) << endl;
            }
        }
    }
#endif

    //std::map<pair<size_t,size_t>,size_t> tet_pair_type; // store the <beg tet,end tet idx> and the type
    vector<pair<list<pair<size_t,size_t> >,size_t> > unknow_equations;
    enum JUMP_TYPE{
        JUMP_X_1 = 0, JUMP_X_2, JUMP_X_3, JUMP_Y_1, JUMP_Y_2, JUMP_Y_3, JUMP_Z_1, JUMP_Z_2, JUMP_Z_3, JUMP_Identity};

    matrixst common_face(3);
    std::vector<vector<size_t> > equation_variable;
    vector<pair<list<pair<size_t,size_t> >,size_t> > uncertain_singularity_equations;
    std::vector<vector<size_t> > uncertain_equation_variable;
    typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
    for(oecit it = ortae.e2t_.begin(); it != ortae.e2t_.end(); ++it){
        const vector<size_t> &loop = it->second;
        const pair<size_t,size_t> &edge = it->first;
        if(loop.front() != loop.back()) continue;
        if(loop.front() == -1) continue;

        int rtn = get_singularity_edge_type(edge.first,edge.second,singularity_segment);
        if(rtn == 9){ // not singularity edge
            list<pair<size_t,size_t> > list_;
            for(size_t t = 0; t < loop.size() - 1; ++t){
                list_.push_back(make_pair(loop[t],loop[t+1]));
                if(jtf::mesh::find_common_face(tet(colon(), loop[t]),tet(colon(),loop[t+1]),&common_face[0])) {
                    cerr << "# strange can not find common face." << endl;
                    cerr << "# tet pair < " << loop[t] << "," << loop[t+1] << endl;
                }
            }
            unknow_equations.push_back(make_pair(list_,static_cast<size_t>(JUMP_Identity))); // jump identity
            vector<size_t> equation_variable_each(loop.size() - 1,-1); // set the equation variable as -1 unknown, record the jump type
            equation_variable.push_back(equation_variable_each);
        }else{
            if(rtn == -1) // uncertain singularity edge
            {
                list<pair<size_t,size_t> > list_;
                for(size_t t = 0; t < loop.size() - 1; ++t){
                    list_.push_back(make_pair(loop[t],loop[t+1]));
#if 0 // test
                    if(loop[t] == loop[t+1]){
                        cerr << "# strange this loop is error" << endl;
                        for(size_t i = 0; i < loop.size() - 1; ++i){
                            cerr << loop[i] << " " ;
                        }
                        cerr << endl;
                    }
#endif
                }
                uncertain_singularity_equations.push_back(make_pair(list_,rtn));
                //unknow_equations.push_back(make_pair(list_,static_cast<size_t>(JUMP_Identity))); // test
                vector<size_t> equation_variable_each(loop.size() - 1 ,-1); // set the equation variable as -1 unknown
                uncertain_equation_variable.push_back(equation_variable_each);
            }else{
                list<pair<size_t,size_t> > list_;
                for(size_t t = 0; t < loop.size() - 1; ++t){
                    list_.push_back(make_pair(loop[t],loop[t+1]));
#if 0 //test
                    if(loop[t] == loop[t+1]){
                        cerr << "# strange this loop is error" << endl;
                        for(size_t i = 0; i < loop.size() - 1; ++i){
                            cerr << loop[i] << " " ;
                        }
                        cerr << endl;
                    }
#endif
                }
                unknow_equations.push_back(make_pair(list_,rtn));
                //unknow_equations.push_back(make_pair(list_,static_cast<size_t>(JUMP_Identity))); // test
                vector<size_t> equation_variable_each(loop.size() - 1 ,-1); // set the equation variable as -1 unknown
                equation_variable.push_back(equation_variable_each);
            }
        }
    }// finish setting up the equations

    cerr << "# unknonw equations num " << unknow_equations.size() << endl;
    cerr << "# start to solve the equations." << endl;

    {// solve the equations
        std::vector<bool> equation_state(unknow_equations.size(),false); // stores whether the equation is resolved, true for yes, false for no
        assert(equation_variable.size() == unknow_equations.size());

        {// accroding to the cut_mesh, remove part of unknown variable // for each inner tet_pair in cut_tet mesh, assign the type as JUMP_Identity
            for(size_t t = 0; t < fa_cut.face2tet_.size(); ++t){
                if(fa_cut.is_outside_face(fa_cut.face2tet_[t])) continue;
                const pair<size_t,size_t> &tet_pair = fa_cut.face2tet_[t];
                assert(tet_pair.first != -1 && tet_pair.second != -1);
                tet_pair_type[make_pair(tet_pair.first,tet_pair.second)] = static_cast<size_t>(JUMP_Identity);
                tet_pair_type[make_pair(tet_pair.second,tet_pair.first)] = static_cast<size_t>(JUMP_Identity);
            }
        }
#if 0// test the original
        matrixd rot_01(3,3);
        matrixd rot_10(3,3);
        for(size_t t = 0; t < fa.face2tet_.size(); ++t){
            if(fa.is_outside_face(fa.face2tet_[t])) continue;
            const pair<size_t,size_t> &tet_pair = fa.face2tet_[t];
            const size_t tet_idx_0 = tet_pair.first;
            const size_t tet_idx_1 = tet_pair.second;
            get_best_alignment(&frame_inner[tet_idx_0][0],&frame_inner[tet_idx_1][0],&rot_01[0]);
            get_best_alignment(&frame_inner[tet_idx_1][0],&frame_inner[tet_idx_0][0],&rot_10[0]);
            tet_pair_type[make_pair(tet_idx_0,tet_idx_1)] = type_transition1(rot_01);
            tet_pair_type[make_pair(tet_idx_1,tet_idx_0)] = type_transition1(rot_10);
        }
#endif

        typedef set<pair<list<pair<size_t,size_t> >,size_t> >::iterator sit;
        typedef list<pair<size_t,size_t> >::iterator lit;
        typedef map<pair<size_t,size_t>,size_t>::const_iterator mcit;
        typedef map<pair<size_t,size_t>,size_t>::iterator mit;
        size_t unresolved_equations = 0;

        while(find(equation_state.begin(),equation_state.end(),false) != equation_state.end()){ // still exist unresolved equation

            cerr << "# unknown equations num " << count(equation_state.begin(),equation_state.end(),false) << endl;
            for(size_t t = 0; t < unknow_equations.size(); ++t){
                if(equation_state[t]) continue; // if the equation has been resolved
                list<pair<size_t,size_t> > &list_ = unknow_equations[t].first;

                lit lit_ = list_.begin();
                for(size_t i = 0; i < equation_variable[t].size(); ++i,++lit_) // for unresolved equation variable in equation t
                {
                    if(equation_variable[t][i] != -1) continue; // !=-1 means this variable type is known

                    mcit mcit_ = tet_pair_type.find(*lit_);
                    if(mcit_ == tet_pair_type.end()) continue;

                    size_t type = mcit_->second;
                    equation_variable[t][i] = type;
                }

                size_t count_unknown = count(equation_variable[t].begin(),equation_variable[t].end(),-1);

                if(count_unknown == 0) {
                    // check whther the equation is correct
                    matrixd rot = eye<double>(3);
                    for(size_t i = 0; i < equation_variable[t].size(); ++i){
                        rot = temp(rot * type_transition2(equation_variable[t][i]));
                    }
                    if(norm(rot - type_transition2(unknow_equations[t].second)) < 1e-8){
                        equation_state[t] = true; // means this equation is resolved
                        break;
                    }else{
                        cerr << "# error it's a conflict equation." << endl;
                        cerr << "# equation type shoule be " << unknow_equations[t].second << endl;
                        cerr << "# equation type actually is " << type_transition1(rot) << endl;
                        cerr << "# equation type list : " ;
                        for(size_t j = 0; j < equation_variable[t].size(); ++j){
                            cerr << equation_variable[t][j] << " ";
                        }
                        cerr << endl;
                    }
                }
                if(count_unknown == 1){ // has only one unknown variable, can be resolved
                    matrixd left = eye<double>(3);
                    matrixd right = eye<double>(3);

                    const vector<size_t>::const_iterator unknown_variable_it
                            = find(equation_variable[t].begin(),equation_variable[t].end(),-1); // find the unknown variable iterator
                    assert(unknown_variable_it != equation_variable[t].end());

                    lit lit_ = list_.begin();
                    for(vector<size_t>::const_iterator it = equation_variable[t].begin();
                        it != unknown_variable_it; ++it,++lit_){
                        left = temp( left * type_transition2(*it));
                    }

                    for(vector<size_t>::const_iterator it = unknown_variable_it + 1;
                        it != equation_variable[t].end(); ++it){
                        right = temp( right * type_transition2(*it));
                    }

                    const matrixd left_ = left;
                    const matrixd right_ = right;
                    matrixd Type = type_transition2(unknow_equations[t].second);
                    //inv(left);
                    //inv(right);
                    matrixd itef_matrix = trans(left) * Type;
                    itef_matrix = temp(itef_matrix * trans(right));
                    size_t type_ = type_transition1(itef_matrix);
#if 0 // the rotation type supposed to be 0~23, 25 means error.
                    if(type_ == 25){
                        cerr << "#  strange type " << "equation type is " << unknow_equations[t].second << endl;
                        cerr << "#  equation variable " ;
                        for(size_t j = 0; j < equation_variable[t].size(); ++j) cerr << equation_variable[t][j] << " ";
                        cerr << endl;

                        cerr << "#  left type " << type_transition1(left_) << left_ << endl;
                        cerr << "# right type " << type_transition1(right_) << right_ << endl;
                        cerr << "# unknown type " << type_transition1(itef_matrix) << endl;
                    }
#endif
#if 0 // test
                    matrixd test = eye<double>(3);
                    inv(left);
                    inv(right);
                    test = temp(left * itef_matrix) * right ;
                    if(type_transition1(test) != unknow_equations[t].second)
                        cerr << "# calculate error." << endl;
#endif

                    tet_pair_type[*lit_] = type_;
                    //inv(itef_matrix);
                    tet_pair_type[make_pair(lit_->second,lit_->first)] = type_transition1(trans(itef_matrix));
                    equation_state[t] = true;
                    // then the problem is left * X * right = TYPE, X = (left)^-1 * TYPE * (right)^-1
                }
            }
            if(unresolved_equations == count(equation_state.begin(),equation_state.end(),false)) break;
            unresolved_equations = count(equation_state.begin(),equation_state.end(),false);
            //cerr << "unresolved equations num " << unresolved_equations << endl;
        }

#if 0 // check the tet_pair validate
        matrixst common_face(3);
        for(map<pair<size_t,size_t>,size_t>::const_iterator it = tet_pair_type.begin();
            it != tet_pair_type.end(); ++it){
            if(find_common_face(tet,it->first.first,it->first.second,&common_face[0]))
            {
                cerr << "# error can not find common face for " << it->first.first << "," << it->first.second << endl;
            }
        }
#endif

#if 1 // check whether the resolution is correct first time
        set<size_t> error_tet;
        vector<size_t> error_lines;
        cerr << "## start the first check." << endl;
        {
            matrixd Type(3,3);
            matrixd rot;
            for(oecit it = ortae.e2t_.begin(); it != ortae.e2t_.end(); ++it){
                const vector<size_t> &loop = it->second;
                const pair<size_t,size_t> &edge = it->first;
                if(loop.front() != loop.back()) continue;
                if(loop.front() == -1) continue;

                int rtn = get_singularity_edge_type(edge.first,edge.second,singularity_segment);

                bool is_this_edge_need_to_be_checked = true;
                for(size_t t = 0; t < loop.size() - 1; ++t){
                    if(tet_pair_type.find(make_pair(loop[t],loop[t+1])) == tet_pair_type.end())
                    {
                        is_this_edge_need_to_be_checked = false;
                        break;
                    }
                }
                if(!is_this_edge_need_to_be_checked) continue;

                cerr << "# checking edge <" << edge.first << "," << edge.second << ">" << endl;
                if(rtn == 9) // not singularity edge
                {
                    rot = eye<double>(3);
                    for(size_t t = 0;t < loop.size() - 1; ++t){
                        //cerr << tet_pair_type[make_pair(loop[t],loop[t+1])] << " ";
                        mcit cit = tet_pair_type.find(make_pair(loop[t],loop[t+1]));
                        if(cit != tet_pair_type.end())
                            rot = temp(rot * type_transition2(cit->second));
                    }
                    //cerr << endl;
                    if(norm(rot - eye<double>(3)) > 1e-8){
                        cerr << " # error. this edge suppose to be non-singulairty edge, it's actual type is "
                             << type_transition1(rot) << endl;
                        error_lines.push_back(edge.first);
                        error_lines.push_back(edge.second);
                        for(size_t t = 0; t < loop.size() - 1; ++t){
                            error_tet.insert(loop[t]);
                        }
                    }
                }else
                {
                    if(rtn != -1)
                    {
                        Type = type_transition2(rtn);
                        rot = eye<double>(3);
                        for(size_t t = 0 ; t < loop.size() - 1; ++t){
                            //cerr << tet_pair_type[make_pair(loop[t],loop[t+1])] << " ";
                            mcit cit = tet_pair_type.find(make_pair(loop[t],loop[t+1]));
                            if(cit != tet_pair_type.end())
                                rot = temp(rot * type_transition2(cit->second));
                        }
                        //cerr << endl;
                        if(norm(Type - rot) > 1e-8){
                            cerr << " # error. this edge suppose to be type " << rtn << ", the actual type is " << type_transition1(rot) << endl;
                            error_lines.push_back(edge.first);
                            error_lines.push_back(edge.second);
                            for(size_t t = 0; t < loop.size() - 1; ++t){
                                error_tet.insert(loop[t]);
                            }
                        }
                    }else
                    {
                        rot = eye<double>(3);
                        for(size_t t = 0 ; t < loop.size() - 1; ++t){
                            //cerr << tet_pair_type[make_pair(loop[t],loop[t+1])] << " ";
                            mcit cit = tet_pair_type.find(make_pair(loop[t],loop[t+1]));
                            if(cit != tet_pair_type.end())
                                rot = temp(rot * type_transition2(cit->second));
                        }
                        if(type_transition1(rot) > 8){
                            cerr << " # error. this edge suppose to be singularity, the actual type is " << type_transition1(rot) << endl;
                            error_lines.push_back(edge.first);
                            error_lines.push_back(edge.second);
                            for(size_t t = 0; t < loop.size() - 1; ++t){
                                error_tet.insert(loop[t]);
                            }
                        }
                    }
                }
            }
        }

        //dump out the error tets and lines
        {
            if(!error_tet.empty()){
                vector<size_t> error_tets_vec(error_tet.size());
                copy(error_tet.begin(),error_tet.end(),error_tets_vec.begin());
                matrixst error_tet_ (4,error_tets_vec.size());
                for(size_t t = 0 ;t < error_tet_.size(2); ++t){
                    error_tet_(colon(),t) = tet(colon(),error_tets_vec[t]);
                }
                string error_tet_str = pt.get<string>("zyz.value");
                error_tet_str += ".error_tet.vtk";
                string error_line_str = pt.get<string>("zyz.value");
                error_line_str += ".error_line.vtk";
                ofstream err_tet_ofs(error_tet_str.c_str());
                ofstream err_line_ofs(error_line_str.c_str());

                tet2vtk(err_tet_ofs,&node[0],node.size(2),&error_tet_[0],error_tet_.size(2));
                line2vtk(err_line_ofs,&node[0],node.size(2),&error_lines[0],error_lines.size()/2);
            }
        }
#endif
        cerr << "unresolved equations num " << count(equation_state.begin(),equation_state.end(),false) << endl;
        //stack<pair<pair<size_t,size_t>,size_t> > try_stack;
        matrixst common_face(3);
        if(count(equation_state.begin(),equation_state.end(),false) != 0)
        {
            cerr << "# begin to recurisve searching..." << endl;
            vector<pair<size_t,size_t> > possible_variable;
            {// list all the posibile variable, each can choose the 24 types of rotation
                for(size_t t = 0; t < equation_state.size(); ++t){
                    if(equation_state[t]) continue;
                    list<pair<size_t,size_t> > &list_ = unknow_equations[t].first;
                    size_t idx = 0;
                    for(list<pair<size_t,size_t> >::const_iterator lci = list_.begin();
                        lci != list_.end(); ++lci,++idx){

                        if(equation_variable[t][idx] != -1) continue;

                        if(!possible_variable.empty()
                                && find(possible_variable.begin(),possible_variable.end(),*lci) == possible_variable.end()
                                && find(possible_variable.begin(),possible_variable.end(),make_pair(lci->second,lci->first)) == possible_variable.end())
                            possible_variable.push_back(*lci);
                        else if(possible_variable.empty()) possible_variable.push_back(*lci);
                    }
                }

#if 1 //test
                for(size_t t = 0; t <possible_variable.size(); ++t )
                  if(jtf::mesh::find_common_face(tet(colon(),possible_variable[t].first),tet(colon(),possible_variable[t].second),&common_face[0])) {
                        cerr << "# [in left unknown equation] strange can not find a common face." << endl;
                        cerr << "# tet " << possible_variable[t].first << tet(colon(),possible_variable[t].first) << endl;
                        cerr << "# tet " << possible_variable[t].second << tet(colon(),possible_variable[t].second) << endl;
                        cerr << __LINE__ << endl;
                        continue;
                    }
#endif
            }

            assert(possible_variable.size());
            if(possible_variable.size() == 0) {
                cerr << "# error no possible variable, but still have " << count(equation_state.begin(),equation_state.end(),false) << "unknown equations." << endl;
                return __LINE__;
            }
            const size_t begin_try_type = 0;
            const size_t MAX_TRY_TYPE = 24;
            equation_state_str begin;


            begin.equation_variable_tmp = equation_variable;
            begin.unknow_equations_tmp = unknow_equations;
            begin.equation_state_tmp = equation_state;
            begin.tet_pair_type_tmp = tet_pair_type;
            begin.possible_variable_tmp = possible_variable;

            vector<bool> pvrt(begin.possible_variable_tmp.size(),false);
            begin.possible_variable_resolved_tmp = pvrt;
            begin.possible_variable_tried.push_back(0);

            stack<equation_state_str> try_stack;
            try_stack.push(begin);
            bool is_solved = false;
            equation_state_str finial_solution;
            while(!try_stack.empty()){
                size_t j = begin_try_type ;
                for(; j < MAX_TRY_TYPE; ++j){
                    equation_state_str temp = try_stack.top();
                    temp.possible_variable_type = j;
                    int rtn = solve_equation_state(temp);
                    if(rtn == 0) {
                        is_solved = true;
                        cerr << "# finish resolving the problem." << endl;
                        finial_solution = temp;
                        break;
                    }
                    if(rtn == 1) {
                        cerr << "# meet a conflict equations." << endl;
                        continue;
                    }if(rtn == 2){// if can not resolve all equations, then should try some new possible_variable
                        try_stack.push(temp);
                        cerr << "# can not resolve all equations." << endl;
                        break;
                    }
                }
                if(is_solved) break;
                else{
                    if(j == MAX_TRY_TYPE){ // meet conflict equations on the last tryes MAX_TRY_TYPE - begin_try_type

                        while(!try_stack.empty()) { // maybe need to go back
                            size_t i = 0;
                            equation_state_str &temp = try_stack.top();
                            for(; i < temp.possible_variable_tmp.size(); ++i)
                                if(temp.possible_variable_resolved_tmp[i] == false
                                        && find(temp.possible_variable_tried.begin(),temp.possible_variable_tried.end(),i) == temp.possible_variable_tried.end()
                                        && i > temp.possible_variable_tried.back()) { // the last tried possible_variable is not correct, need to try a new
                                    temp.possible_variable_tried.back() = i;
                                    break;
                                }
                            if(i == temp.possible_variable_tmp.size()) // means can not find a new item to try, need to back again
                            {
                                try_stack.pop();
                            }else break;
                        }
                        if(try_stack.empty()) {
                            cerr << "# can not resolve equations, there are conflicts." << endl;
                            break;
                        }
                    }else{ // means can not resolve all equations, need to try some new possible_variable
                        equation_state_str &temp = try_stack.top();
                        size_t i = 0;
                        for(; i < temp.possible_variable_tmp.size(); ++i)
                            if(temp.possible_variable_resolved_tmp[i] == false
                                    &&find(temp.possible_variable_tried.begin(),temp.possible_variable_tried.end(),i) == temp.possible_variable_tried.end()) {
                                temp.possible_variable_tried.push_back(i); // add a new possible_variable_tried to try
                                break;
                            }
                        if(i == temp.possible_variable_tmp.size()){ // means can not resolve all equations, strange!
                            if(count(temp.possible_variable_resolved_tmp.begin(),temp.possible_variable_resolved_tmp.end(),false) == 0){
                                is_solved = true;
                                break;
                            }

                            cerr << "# strange, all possible_variable are tried, but still can not resolve all equations." << endl;
                            cerr << "# tried possible_variable size " << temp.possible_variable_tried.size() << endl;
                            cerr << "# possible variable size " << temp.possible_variable_tmp.size() << endl;

                            for(size_t t = 0; t < temp.unknow_equations_tmp.size(); ++t){
                                if(!temp.equation_state_tmp[t]){
                                    // output the unknown strange equation
                                    cerr << "# equation state: " << temp.unknow_equations_tmp[t].second << endl;
                                    const list<pair<size_t,size_t> > & list_ = temp.unknow_equations_tmp[t].first;
                                    for(list<pair<size_t,size_t> >::const_iterator lci = list_.begin();
                                        lci != list_.end(); ++lci){
                                        mcit mcit_ = temp.tet_pair_type_tmp.find(*lci);
                                        if(mcit_ != temp.tet_pair_type_tmp.end()){
                                            cerr << temp.tet_pair_type_tmp[*lci] << " " ;
                                        }else{
                                            cerr << -1 << " ";
                                        }
                                    }
                                    cerr << endl;
                                }
                            }
                            break;
                        }
                    }
                }
            }
            if(is_solved){
                tet_pair_type = finial_solution.tet_pair_type_tmp;
            }else {
                cerr << "# can not find jump type for each inner face." << endl;
                return __LINE__;
            }
        }

#if 0 // check whether the resolution can satisify the uncertain_singularity_edge
        vector<bool> uncertain_singularity_equations_state(uncertain_singularity_equations.size(),false);
        while(find(uncertain_singularity_equations_state.begin(),uncertain_singularity_equations_state.end(),false)
              != uncertain_singularity_equations_state.end()){ // still exist unresolved equation

            cerr << "# unknown uncertain singularity equations num " << count(uncertain_singularity_equations_state.begin(),
                                                                              uncertain_singularity_equations_state.end(),false) << endl;
            for(size_t t = 0; t < uncertain_singularity_equations.size(); ++t){
                if(uncertain_singularity_equations_state[t]) continue; // if the equation has been resolved
                list<pair<size_t,size_t> > &list_ = uncertain_singularity_equations[t].first;

                lit lit_ = list_.begin();
                for(size_t i = 0; i < uncertain_equation_variable[t].size(); ++i,++lit_) // for unresolved equation variable in equation t
                {
                    if(uncertain_equation_variable[t][i] != -1) continue; // !=-1 means this variable type is known

                    mcit mcit_ = tet_pair_type.find(*lit_);
                    if(mcit_ == tet_pair_type.end()) continue;

                    size_t type = mcit_->second;
                    uncertain_equation_variable[t][i] = type;
                }

                size_t count_unknown = count(uncertain_equation_variable[t].begin(),uncertain_equation_variable[t].end(),-1);

                if(count_unknown == 0) {
                    // check whther the equation is correct
                    matrixd rot = eye<double>(3);
                    for(size_t i = 0; i < uncertain_equation_variable[t].size(); ++i){
                        rot = temp(rot * type_transition2(uncertain_equation_variable[t][i]));
                    }
                    if(type_transition1(rot) < 9) {
                        uncertain_singularity_equations_state[t] = true;
                    }else{
                        cerr << "# error it's a conflict equation." << endl;
                        cerr << "# equation type actually is " << type_transition1(rot) << endl;
                        cerr << "# equation type list : " ;
                        for(size_t j = 0; j < uncertain_equation_variable[t].size(); ++j){
                            cerr << uncertain_equation_variable[t][j] << " ";
                        }
                        cerr << endl;
                    }
                }
            }
            if(unresolved_equations == count(equation_state.begin(),equation_state.end(),false)) break;
            unresolved_equations = count(equation_state.begin(),equation_state.end(),false);
            //cerr << "unresolved equations num " << unresolved_equations << endl;
        }
#endif
        typedef set<pair<list<pair<size_t,size_t> >,size_t> >::iterator sit;
        typedef map<pair<size_t,size_t>,size_t >::const_iterator mpcit;
        typedef list<pair<size_t,size_t> >::iterator lit;


#if 0 // check the jump face type
        for(mpcit ci = tet_pair_type.begin(); ci != tet_pair_type.end(); ++ci){
            if(ci->second == 9 || ci->second == -1){
                cerr << "# strange jump face, tet pair: <" << ci->first.first << "," << ci->first.second << ">" << " type " << ci->second << endl;
            }
        }
#endif

        jump_face_type = -1 * ones<size_t>(face_pair.size(),1);
        for(size_t t = 0; t < face_pair.size(); ++t){
            if(face_pair[t] != -1){
                const pair<size_t,size_t> &tet_pair = fa_cut.face2tet_[outside_face_cut_idx[t]];
                const size_t tet_idx = (tet_pair.first == -1)?tet_pair.second:tet_pair.first;

                const pair<size_t,size_t> &tet_pair_other = fa_cut.face2tet_[outside_face_cut_idx[face_pair[t]]];
                const size_t tet_idx_other = (tet_pair_other.first == -1)?tet_pair_other.second:tet_pair_other.first;
                //if()
                mpcit it = tet_pair_type.find(make_pair(tet_idx,tet_idx_other));
                if(it == tet_pair_type.end()) jump_face_type[t] = 9; // identity
                else
                    jump_face_type[t] = it->second;

                //assert(jump_face_type[t] != 9 && jump_face_type[t] != -1);
            }
        }
#if 1 // dump out the jump type
        if( zjucad::has("jump_type_dump.value",pt)){
            string dump_file = pt.get<string>("jump_type_dump.value");
            ofstream ofs(dump_file.c_str());
            if(ofs.fail()){
                cerr << "# can not open jump_type_dump file." << endl;
                return __LINE__;
            }
            for(map<pair<size_t,size_t>,size_t>::const_iterator mci = tet_pair_type.begin();
                mci != tet_pair_type.end(); ++mci){
                ofs << mci->first.first << " " << mci->first.second << "  " << mci->second << endl;
            }

            string dump_file_cut = dump_file + "_no_identity";
            ofstream ofs_cut(dump_file_cut.c_str());
            if(ofs_cut.fail()){
                cerr << "# can not open jump_type_dump file." << endl;
                return __LINE__;
            }
            for(map<pair<size_t,size_t>,size_t>::const_iterator mci = tet_pair_type.begin();
                mci != tet_pair_type.end(); ++mci){
                if(mci->second != 9) // not identity
                    ofs_cut << mci->first.first << " " << mci->first.second << "  " << mci->second << endl;
            }
            cerr << "# success dump jump_type file." << endl;
        }
#endif
#if 0 // check the jump type
        matrixd rot(3,3);
        for(mpcit it = tet_pair_type.begin(); it != tet_pair_type.end(); ++it){
            get_best_alignment(&frame_inner[it->first.first][0], &frame_inner[it->first.second][0], &rot[0]);
            const size_t type = type_transition(rot);
            if(type != it->second){
                cerr << "# jump type : original " << type << " new " << it->second << endl;
            }
        }
#endif

#if 1 // check whether the new jump face is right for singuolarity edge
        cerr << "# finial check the jump face type." << endl;
        vector<size_t> error_line_2;
        set<size_t> error_tet_2;
        {
            matrixd Type(3,3);
            matrixd rot;
            for(oecit it = ortae.e2t_.begin(); it != ortae.e2t_.end(); ++it){

                const vector<size_t> &loop = it->second;
                const pair<size_t,size_t> &edge = it->first;
                if(loop.front() != loop.back()) continue;
                if(loop.front() == -1) continue;
                int rtn = get_singularity_edge_type(edge.first,edge.second,singularity_segment);
                cerr << "# checking edge <" << edge.first << "," << edge.second << "> " << endl;
                if(rtn == 9) // not singularity edge
                {
                    rot = eye<double>(3);
                    for(size_t t = 0;t < loop.size() - 1; ++t){
                        //cerr << tet_pair_type[make_pair(loop[t],loop[t+1])] << " ";
                        mcit cit = tet_pair_type.find(make_pair(loop[t],loop[t+1]));
                        if(cit != tet_pair_type.end())
                            rot = temp(rot * type_transition2(cit->second));
                        //cerr << "t = " << t << rot << endl;
                    }
                    //cerr << endl;
                    if(norm(rot - eye<double>(3)) > 1e-8){
                        cerr << "\n #  error. this edge suppose to be non-singulairty edge, it's actual type is " << type_transition1(rot) << endl;
                        error_line_2.push_back(it->first.first);
                        error_line_2.push_back(it->first.second);
                        for(size_t t = 0; t < loop.size() - 1; ++t){
                            error_tet_2.insert(loop[t]);
                        }
                    }
                    //                    else
                    //                        cerr << " correct!" << endl;
                }else
                {
                    cerr << " singularity edge " << endl;
                    if (rtn != -1)
                    {
                        Type = type_transition2(rtn);
                        rot = eye<double>(3);
                        for(size_t t = 0 ; t < loop.size() - 1; ++t){
                            //cerr << tet_pair_type[make_pair(loop[t],loop[t+1])] << " ";
                            mcit cit = tet_pair_type.find(make_pair(loop[t],loop[t+1]));
                            if(cit != tet_pair_type.end())
                                rot = temp(rot * type_transition2(cit->second));
                        }
                        //cerr << endl;
                        if(norm(Type - rot) > 1e-8){
                            cerr << "\n #  error. this edge suppose to be type " << rtn << ", the actual type is " << type_transition1(rot) << endl;
                            error_line_2.push_back(it->first.first);
                            error_line_2.push_back(it->first.second);
                            for(size_t t = 0; t < loop.size() - 1; ++t){
                                error_tet_2.insert(loop[t]);
                            }
                        }/*else
                            cerr << " correct!" << endl;*/
                    }else{
                        typedef map<pair<size_t,size_t>,size_t>::iterator mit;
                        mit it = singularity_segment.find(edge);
                        if(it == singularity_segment.end()) it = singularity_segment.find(make_pair(edge.second,edge.first));

                        rot = eye<double>(3);
                        for(size_t t = 0 ; t < loop.size() - 1; ++t){
                            mcit cit = tet_pair_type.find(make_pair(loop[t],loop[t+1]));
                            if(cit != tet_pair_type.end())
                                rot = temp(rot * type_transition2(cit->second));
                        }
                        it->second = type_transition1(rot);
                        cerr << "\n # this edge suppose to be unknown , the actual type is " << type_transition1(rot) << endl;
                    }
                }
            }
        }

        {
            if(!error_tet_2.empty()){
                vector<size_t> error_tet_vec_2(error_tet_2.size());
                copy(error_tet_2.begin(),error_tet_2.end(),error_tet_vec_2.begin());
                matrixst error_tet_2_ (4,error_tet_vec_2.size());
                for(size_t t = 0 ;t < error_tet_2_.size(2); ++t){
                    error_tet_2_(colon(),t) = tet(colon(),error_tet_vec_2[t]);
                }
                string error_tet_str = pt.get<string>("zyz.value");
                error_tet_str += ".error_tet_2.vtk";
                string error_line_str = pt.get<string>("zyz.value");
                error_line_str += ".error_line_2.vtk";
                ofstream err_tet_ofs(error_tet_str.c_str());
                ofstream err_line_ofs(error_line_str.c_str());

                tet2vtk(err_tet_ofs,&node[0],node.size(2),&error_tet_2_[0],error_tet_2_.size(2));
                line2vtk(err_line_ofs,&node[0],node.size(2),&error_line_2[0],error_line_2.size()/2);
            }
        }
#endif

#if 0   // check the singularity edge info should be no unknown
        typedef map<pair<size_t,size_t>,size_t>::const_iterator mcit;
        for(mcit it = singularity_segment.begin(); it != singularity_segment.end(); ++it){
            cerr << "# this singularity edge type is " << it->second << endl;
            if(it->second == -1){
                cerr << "# error this singularity segment is sitll unknown." << endl;
                cerr << "# this edge <" << it->first.first << "," << it->first.second << "> is sitll unknown." << endl;
            }
        }
#endif

#if 1 // dump out the jump face to check

        map<size_t,size_t> jump_type_cluster; // first store the jump type, second store the render type
        {// set the jump face type
            set<pair<size_t,size_t> > jump_type_collect;
            matrixd rot;
            for(size_t t = 0; t < 24; ++t){
                rot = type_transition2(t);
                inv(rot);
                size_t type =type_transition1(rot);
                if(t > type) {
                    jump_type_collect.insert(make_pair(type,t));
                }else{
                    jump_type_collect.insert(make_pair(t,type));
                }
            }
            size_t render_type = 0;
            for(set<pair<size_t,size_t> >::const_iterator it = jump_type_collect.begin();
                it != jump_type_collect.end(); ++it){
                jump_type_cluster[it->first] = render_type;
                jump_type_cluster[it->second] = render_type;
                render_type++;
            }
        }

        {
            set<pair<pair<size_t,size_t>,size_t> > jump_tet_pair;
            matrixst common_face(3);
            vector<size_t> jump_face;
            vector<size_t> jump_face_render_type;
            for(map<pair<size_t,size_t>,size_t>::const_iterator ci = tet_pair_type.begin(); ci != tet_pair_type.end(); ++ci){
                //                if(find_common_face(tet,ci->first.first,ci->first.second,&common_face[0])) {
                //                    cerr << "# strange can not find a common face." << endl;
                //                    cerr << "# tet " << ci->first.first << tet(colon(),ci->first.first) << endl;
                //                    cerr << "# tet " << ci->first.second << tet(colon(),ci->first.second) << endl;
                //                    cerr << __LINE__ << endl;
                //                    continue;
                //                }
                if(ci->second != 9) // not identity
                {
                    if(ci->first.first > ci->first.second)
                        jump_tet_pair.insert(make_pair(make_pair(ci->first.second,ci->first.first),jump_type_cluster[ci->second]));
                    else
                        jump_tet_pair.insert(make_pair(ci->first,jump_type_cluster[ci->second]));
                }
            }
            cerr << "# step 0 . " << endl;
            for(set<pair<pair<size_t,size_t>,size_t> >::const_iterator sci = jump_tet_pair.begin();
                sci != jump_tet_pair.end(); ++sci){
                if(jtf::mesh::find_common_face(tet(colon(),sci->first.first), tet(colon(),sci->first.second),&common_face[0])) {
                    cerr << "# strange can not find a common face." << endl;
                    cerr << "# tet " << sci->first.first << tet(colon(),sci->first.first) << endl;
                    cerr << "# tet " << sci->first.second << tet(colon(),sci->first.second) << endl;
                    cerr << __LINE__ << endl;
                    continue;
                }
                for(size_t i = 0; i < 3; ++i)
                    jump_face.push_back(common_face[i]);
                jump_face_render_type.push_back(sci->second);
            }
            cerr << "# step 1 . " << endl;
            cerr << "# jump face size = " << jump_face.size()/3 << endl;
            {// for vis
                string jump_face_str = pt.get<string>("output_zyz.value");
                jump_face_str += ".jump_face.vtk";
                ofstream ofs(jump_face_str.c_str());
                tri2vtk(ofs,&node[0],node.size(2),&jump_face[0],jump_face.size()/3);
                cell_data(ofs,&jump_face_render_type[0],jump_face_render_type.size(),"jump_face_type");
            }
        }
#endif
        return 0;
    }
}

