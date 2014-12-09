#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <boost/smart_ptr.hpp>
#include <algorithm>
#include <boost/tuple/tuple_comparison.hpp>
#include <jtflib/util/util.h>
#include <fstream>
#include <boost/unordered_set.hpp>
#include <omp.h>
#include <zjucad/matrix/io.h>
#include "extractor.h"
#include "../common/vtk.h"
#include "feature.h"

using namespace std;
using namespace zjucad::matrix;
using namespace jtf::mesh;

#define __EPS  0.00001


namespace extractor {

void hex_extractor::calculate_tet_base_inverse(const zjucad::matrix::matrix<size_t>    &param_mesh,
                                               const zjucad::matrix::matrix<double>    &param_nodes,
                                               zjucad::matrix::matrix<double>          &inv_base,
                                               zjucad::matrix::matrix<bool>            &invertible)
{
    cout << "[INFO]the program is calculating tet base inverse\n";
    const size_t tet_num = param_mesh.size(2);
    zjucad::matrix::matrix<double> M(3, 3);
    inv_base.resize(3, 3 * tet_num);
    invertible.resize(1, tet_num);
    for (size_t id = 0; id < tet_num; ++id) {
        for (size_t offset = 0; offset < 3; ++offset)
            M(colon(), offset) = param_nodes(colon(), param_mesh(1 + offset , id)) - param_nodes(colon(), param_mesh(0, id));

        if ( inv(M) ) {
            cerr << "[INFO] inv fail\n";
            invertible(0, id) = false;
        }
        else {
            inv_base(colon(), colon(id * 3, id * 3 + 2)) = M;
            invertible(0, id) = true;
        }
    }
}

bool hex_extractor::point_in_triangle(const zjucad::matrix::matrix<double>  &p,
                                      const zjucad::matrix::matrix<double>  &v0,
                                      const zjucad::matrix::matrix<double>  &v1,
                                      const zjucad::matrix::matrix<double>  &v2)
{
    zjucad::matrix::matrix<double> c0(3, 1), c1(3, 1), c2(3, 1);
    c0 = zjucad::matrix::cross(p - v0, v1 - v0);
    c1 = zjucad::matrix::cross(p - v1, v2 - v1);
    c2 = zjucad::matrix::cross(p - v2, v0 - v2);
    return cast(zjucad::matrix::dot(c0, c1), __EPS) >= 0 && cast(zjucad::matrix::dot(c1, c2), __EPS) >= 0;
}

bool hex_extractor::contains(const zjucad::matrix::matrix<double>     &P,
                             const int                                id,
                             const zjucad::matrix::matrix<double>     &inv_base,
                             const zjucad::matrix::matrix<bool>       &invertible,
                             const zjucad::matrix::matrix<size_t>     &param_mesh,
                             const zjucad::matrix::matrix<double>     &param_nodes,
                             zjucad::matrix::matrix<double>           &_coeff)
{
    zjucad::matrix::matrix<double> a(3, 1);
    if ( !invertible(0, id) )
        return false;
    const double eps = 0.0001;
    a = inv_base(colon(), colon(id * 3, id * 3 + 2)) * ( P - param_nodes(colon(), param_mesh(0, id)));
    _coeff(0, 0) = cast(1 - a(0, 0) - a(1, 0) - a(2, 0), eps);
    _coeff(1, 0) = cast(a(0, 0), eps);
    _coeff(2, 0) = cast(a(1, 0), eps);
    _coeff(3, 0) = cast(a(2, 0), eps);
    for (size_t i = 0; i < 4; ++i) {
        if ( _coeff(i, 0) < 0 )
            return false;
    }
    return true;
}

int hex_extractor::extract_points(const zjucad::matrix::matrix<size_t>               &param_mesh,
                                  const zjucad::matrix::matrix<double>               &param_nodes,
                                  const double                                       offset,
                                  const zjucad::matrix::matrix<double>               &inv_base,
                                  const zjucad::matrix::matrix<bool>                 &invertible,
                                  std::vector<point_info>                            &point_mesh,
                                  std::vector<boost::tuple<double, double, double> > &point_nodes)
{
    if ( offset == 0 )
        cout << "[INFO]the program is extracting integer points\n";
    else
        cout << "[INFO]the program is extracting half integer points\n";

    zjucad::matrix::matrix<double> I(3, 1), coeff(4, 1);
    bounding_box BBox;

//#pragma omp parallel for
    for (size_t id = 0; id < param_mesh.size(2); ++id) {
        if ( !invertible(0, id) )
            continue;
        bound_ith_tet(id, param_mesh, param_nodes, BBox);

        for (int x = (int)floor(BBox.xmin); x <= (int)ceil(BBox.xmax); ++x) {
        for (int y = (int)floor(BBox.ymin); y <= (int)ceil(BBox.ymax); ++y) {
        for (int z = (int)floor(BBox.zmin); z <= (int)ceil(BBox.zmax); ++z) {
            I(0, 0) = x + offset; I(1, 0) = y + offset; I(2, 0) = z + offset;
            if ( contains(I, id, inv_base, invertible, param_mesh, param_nodes, coeff) ) {
                boost::tuple<double, double, double> II(I(0, 0), I(1, 0), I(2, 0));
                point_info mesh(point_nodes.size(), id, coeff);
//#pragma omp critical
                {
                point_mesh.push_back(mesh);
                point_nodes.push_back(II);
                }
            }
        }}}
    }
    cout << "[INFO]the num of extracted points is " << point_mesh.size() << endl;
    return 0;
}

void hex_extractor::merge_points_with_same_orig_coord(const zjucad::matrix::matrix<size_t>                     &param_mesh,
                                                      const std::vector<boost::tuple<double, double, double> > &point_nodes,
                                                      std::vector<point_info>                                  &point_mesh)
{
    cout << "[INFO]the program is merging integer points with identical coordinate in orginal space\n";

    father.resize(point_nodes.size());
    for (size_t i = 0; i < father.size(); ++i)
        father[i] = i;

    for (size_t p = 0; p < point_mesh.size() - 1; ++p) {
        for (size_t q =  p + 1; q < point_mesh.size(); ++q) {
            size_t pcid = point_mesh[p].coordinate_idx;
            size_t qcid = point_mesh[q].coordinate_idx;
            size_t ptid = point_mesh[p].tet_idx;
            size_t qtid = point_mesh[q].tet_idx;
            std::set<size_t> share;
            bool flag = false;
            if ( point_nodes[pcid] == point_nodes[qcid] ) {
                for (size_t k = 0; k < 4; ++k)
                    if ( point_mesh[p].b_coeff(k, 0) > 0 )
                        share.insert(param_mesh(k, ptid));
                for (size_t k = 0; k < 4; ++k)
                    if ( point_mesh[q].b_coeff(k, 0) > 0 )
                        if ( flag = share.insert(param_mesh(k, qtid)).second )
                            break;
                if ( !flag )
                    Union(pcid, qcid);
            }
        }
    }
    for (size_t id = 0; id < point_mesh.size(); ++id)
        point_mesh[id].coordinate_idx = father[point_mesh[id].coordinate_idx];
}

void hex_extractor::merging_vertices_via_cut_face(const zjucad::matrix::matrix<size_t>                      &param_mesh,
                                                  const zjucad::matrix::matrix<double>                      &param_nodes,
                                                  const std::vector<point_info>                             &point_mesh,
                                                  const std::vector<boost::tuple<double, double, double> >  &point_nodes,
                                                  const std::map<std::pair<size_t, size_t>, size_t>         &inner_type,
                                                  const boost::shared_ptr<jtf::mesh::face2tet_adjacent>     &orig_handle,
                                                  const zjucad::matrix::matrix<size_t>                      &cut_tet2tet)
{
    //cluster the duplicated points caused by numeric accuracy
    //aim : one vertex in orig space only has one copy in param space
    zjucad::matrix::matrix<double> R_st(3, 3), s_node(3, 1), t_node(3, 1), g_st(3, 1), s_vertex(3, 1), t_vertex(3, 1);

    std::multimap<size_t, size_t> tet2coord;
    std::multimap<size_t, size_t>::iterator ps, pt;
    std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator > seg_s, seg_t;
    for (size_t id = 0; id < point_mesh.size(); ++id)
        tet2coord.insert(std::make_pair(point_mesh[id].tet_idx, point_mesh[id].coordinate_idx));

    size_t order[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    std::map<std::pair<size_t, size_t>, size_t>::const_iterator it;
    std::set<std::pair<size_t, size_t> > vis;
    for (it = inner_type.begin(); it != inner_type.end(); ++it) {
        size_t s = it->first.first;
        size_t t = it->first.second;

        if ( vis.find(std::make_pair(s, t)) != vis.end() )
            continue;
        else {
            vis.insert(std::make_pair(s, t));
            vis.insert(std::make_pair(t, s));
        }

        for (size_t x = 0; x < 3; ++x)
            for (size_t y = 0; y < 3; ++y)
                R_st(x, y) = rot[it->second][x][y];

        for (size_t j = 0; j < 4; ++j) {
            std::pair<size_t, size_t> res = orig_handle->query(cut_tet2tet[param_mesh(order[j][0], s)],
                                                               cut_tet2tet[param_mesh(order[j][1], s)],
                                                               cut_tet2tet[param_mesh(order[j][2], s)]);
            if ( res.second == s )
                std::swap(res.first, res.second);
            if ( t == res.second ) {
                s_vertex = param_nodes(colon(), param_mesh(order[j][0], s));
                for (size_t k = 0; k < 4; ++k) {
                    if ( cut_tet2tet[param_mesh(k, t)] == cut_tet2tet[param_mesh(order[j][0], s)] ) {
                        t_vertex = param_nodes(colon(), param_mesh(k, t));
                        break;
                    }
                }
                g_st = t_vertex - R_st * s_vertex;
                break;
            }
        }
        ps = tet2coord.find(s);
        pt = tet2coord.find(t);
        if ( ps != tet2coord.end() && pt != tet2coord.end() ) {
            seg_s = tet2coord.equal_range(ps->first);
            seg_t = tet2coord.equal_range(pt->first);
            for (ps = seg_s.first; ps != seg_s.second; ++ps) {
                for (pt = seg_t.first; pt != seg_t.second; ++pt) {
                    s_node(0, 0) = point_nodes[ps->second].get<0>();
                    s_node(1, 0) = point_nodes[ps->second].get<1>();
                    s_node(2, 0) = point_nodes[ps->second].get<2>();
                    t_node(0, 0) = point_nodes[pt->second].get<0>();
                    t_node(1, 0) = point_nodes[pt->second].get<1>();
                    t_node(2, 0) = point_nodes[pt->second].get<2>();
                    s_node = temp(R_st * s_node) + g_st;
                    if ( norm(s_node - t_node) < 1.0 - EPSILON ) {
                        size_t fs = Find(ps->second);
                        size_t ft = Find(pt->second);
                        if ( fs != ft )
                            father[ft] = fs;
                    }
                }
            }
        }
    }
}

int hex_extractor::build_dual_mesh(const zjucad::matrix::matrix<size_t>               &param_mesh,
                                   const zjucad::matrix::matrix<double>               &param_nodes,
                                   const zjucad::matrix::matrix<size_t>               &orig_mesh,
                                   const std::map<std::pair<size_t, size_t>, size_t>  &inner_type,
                                   std::map<std::pair<size_t, size_t>, double>        &edge)
{
    cout << "[INFO]the program is building the dual mesh of cut mesh\n";
    int order[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3},{1, 2, 3}};

    zjucad::matrix::matrix<double> R_st(3, 3);
    zjucad::matrix::matrix<double> g_st(3, 1), f_s(3, 1), f_t(3, 1);
    zjucad::matrix::matrix<double> face(3, 1);

    zjucad::matrix::matrix<size_t> cut_tet2tet(zjucad::matrix::max(param_mesh) + 1);
    cut_tet2tet(param_mesh) = orig_mesh(colon());
    boost::shared_ptr<jtf::mesh::face2tet_adjacent> orig_handle(jtf::mesh::face2tet_adjacent::create(orig_mesh));

    for (size_t s = 0; s < param_mesh.size(2); ++s)
    {
        for (size_t i = 0; i < 4; ++i) {
            face(0, 0) = param_mesh(order[i][0], s);
            face(1, 0) = param_mesh(order[i][1], s);
            face(2, 0) = param_mesh(order[i][2], s);

            std::pair<size_t, size_t> res = orig_handle->query(cut_tet2tet[face(0, 0)],
                                                               cut_tet2tet[face(1, 0)],
                                                               cut_tet2tet[face(2, 0)]);
            if ( !orig_handle->is_outside_face(res) ) {
                if ( res.second == s )
                    std::swap(res.first, res.second);

                size_t t = res.second;

                auto it = inner_type.find(res);

                if ( it == inner_type.end() )
                    edge.insert(make_pair(res, 1));
                else
                {
                    for (size_t x = 0; x < 3; ++x)
                        for (size_t y = 0; y < 3; ++y)
                            R_st(x, y) = rot[it->second][x][y];

                    f_s = param_nodes(colon(), face(0, 0));
                    size_t correspond;
                    for (size_t p = 0; p < 4; p++) {
                        if ( cut_tet2tet[face(0, 0)] == cut_tet2tet[param_mesh(p, t)] ) {
                            correspond = param_mesh(p, t);
                            break;
                        }
                    }
                    f_t = param_nodes(colon(), correspond);
                    g_st = f_t - R_st * f_s;

                    if ( norm(g_st) < EPSILON )
                        edge.insert(make_pair(res, 1));
                }
            }
        }
    }
    return 0;
}

void hex_extractor::see_dual_graph(const zjucad::matrix::matrix<size_t>               &orig_mesh,
                                   const zjucad::matrix::matrix<double>               &orig_nodes,
                                   const std::map<std::pair<size_t, size_t>, double>  &edge,
                                   zjucad::matrix::matrix<size_t>                     &dual_mesh,
                                   zjucad::matrix::matrix<double>                     &dual_nodes)
{
    dual_mesh.resize(2, edge.size());
    dual_nodes.resize(3, orig_mesh.size(2));

    for (size_t i = 0; i < orig_mesh.size(2); ++i) {
        dual_nodes(colon(), i) = zjucad::matrix::zeros(3, 1);
        for (size_t j = 0; j < 4; ++j) {
            dual_nodes(colon(), i) += 0.25 * orig_nodes(colon(), orig_mesh(j, i));
        }
    }
    int ptr = 0;
    std::map<std::pair<size_t, size_t>, double>::const_iterator it;
    for (it = edge.begin(); it != edge.end(); ++it) {
        dual_mesh(0, ptr) = it->first.first;
        dual_mesh(1, ptr) = it->first.second;
        ++ptr;
    }
}

void hex_extractor::walk(const zjucad::matrix::matrix<double>                      &curr_src,
                         const size_t                                              curr_tet,
                         const zjucad::matrix::matrix<double>                      &curr_dir,
                         const double                                              curr_rest,
                         const boost::shared_ptr<jtf::mesh::face2tet_adjacent>     &orig_handle,
                         const zjucad::matrix::matrix<size_t>                      &cut_tet2tet,
                         const std::map<std::pair<size_t, size_t>, size_t>         &inner_type,
                         const zjucad::matrix::matrix<size_t>                      &param_mesh,
                         const zjucad::matrix::matrix<double>                      &param_nodes,
                         const zjucad::matrix::matrix<double>                      &inv_base,
                         const zjucad::matrix::matrix<bool>                        &invertible,
                         std::set<STATE>                                           &vis,
                         zjucad::matrix::matrix<double>                            &des,
                         size_t                                                    &des_tet,
                         size_t                                                    depth,
                         size_t                                                    dir_changes,
                         bool                                                      &arrive)
{
    zjucad::matrix::matrix<double> R_st(3, 3);
    zjucad::matrix::matrix<double> g_st(3, 1);
    zjucad::matrix::matrix<double> _coeff(4, 1);

    zjucad::matrix::matrix<double> P(3, 4);
    zjucad::matrix::matrix<size_t> face(3, 4);
    zjucad::matrix::matrix<int>    flag(4, 1);

    size_t next_tet;
    double next_rest;
    zjucad::matrix::matrix<double> next_src(3, 1), next_dir(3, 1);

    next_src = curr_src + curr_rest * curr_dir;

    if ( dir_changes % 2 == 0 &&
         contains(next_src, curr_tet, inv_base, invertible, param_mesh, param_nodes, _coeff) ) {


        for (int i = 0; i < 3; i++)
            des(i, 0) = cast(next_src(i, 0), 1.0);

        des_tet = curr_tet;
        arrive = true;
        return;
    }

    zjucad::matrix::matrix<double> end_port(3, 1);
    end_port = dir_changes % 2 == 0 ? next_src : curr_src + curr_dir;
    line_cross_tet(curr_tet, param_mesh, param_nodes, inv_base, invertible,
                   curr_src, end_port, dir_changes % 2, P, face, flag);
//    int cnt = 0;
//    for (int i = 0; i < 4; i++)
//        if ( flag(i, 0) == 1 )
//            cnt++;
//    cout << "cross face : " << cnt << endl;

    for (size_t f = 0; f < 4; ++f)
    {
         if ( flag(f, 0) == 1)
        {
            std::pair<size_t, size_t> res = orig_handle->query(cut_tet2tet[face(0, f)],
                                                               cut_tet2tet[face(1, f)],
                                                               cut_tet2tet[face(2, f)]);
            if ( res.second == curr_tet )
                std::swap(res.first, res.second);
            size_t s = res.first;
            size_t t = res.second;

            if ( !orig_handle->is_outside_face(res) )
            {
                auto it = inner_type.find(res);

                if ( it != inner_type.end() ) {
                    for (size_t i = 0; i < 3; ++i)
                        for (size_t j = 0; j < 3; ++j)
                            R_st(i, j) = rot[it->second][i][j];
                }
                else
                    R_st = zjucad::matrix::eye<double>(3);

                zjucad::matrix::matrix<double> f_s(3, 1);
                zjucad::matrix::matrix<double> f_t(3, 1);
                f_s = param_nodes(colon(), face(0, f));

                for (size_t i = 0; i < 4; ++i) {
                    size_t a_ = cut_tet2tet[param_mesh(i, t)];
                    size_t b_ = face(0, f);
                    size_t c_ = cut_tet2tet[b_];


                    if ( cut_tet2tet[param_mesh(i, t)] == cut_tet2tet[face(0, f)] ) {
                        f_t = param_nodes(colon(), param_mesh(i, t));
                        break;
                    }
                }
                g_st = f_t - R_st * f_s;

                size_t cnt;
                zjucad::matrix::matrix<double> tri_face(3, 3), normal(3, 1);
                for (size_t col = 0; col < 3; ++col)
                    tri_face(colon(), col) = R_st * param_nodes(colon(), face(col, f)) + g_st;
                std::map<size_t, size_t> orig2cut;
                for (size_t col = 0; col < 4; ++col) {
                    orig2cut.insert(std::make_pair(cut_tet2tet[param_mesh(col, s)], param_mesh(col, s)));
                    orig2cut.insert(std::make_pair(cut_tet2tet[param_mesh(col, t)], param_mesh(col, t)));
                }
                size_t face_idx_sum = cut_tet2tet[face(0, f)] + cut_tet2tet[face(1, f)] + cut_tet2tet[face(2, f)];
                size_t rest_s = orig2cut[cut_tet2tet[param_mesh(0, s)] + cut_tet2tet[param_mesh(1, s)] +
                                         cut_tet2tet[param_mesh(2, s)] + cut_tet2tet[param_mesh(3, s)] - face_idx_sum];
                size_t rest_t = orig2cut[cut_tet2tet[param_mesh(0, t)] + cut_tet2tet[param_mesh(1, t)] +
                                         cut_tet2tet[param_mesh(2, t)] + cut_tet2tet[param_mesh(3, t)] - face_idx_sum];

                normal = cross(tri_face(colon(), 1) - tri_face(colon(), 0), tri_face(colon(), 2) - tri_face(colon(), 0));

                double p1 = dot(R_st * param_nodes(colon(), rest_s) + g_st - tri_face(colon(), 0), normal);
                double p2 = dot(param_nodes(colon(), rest_t) - tri_face(colon(), 0), normal);

                if ( dir_changes % 2 == 0 )
                    next_rest = curr_rest - norm(P(colon(), f) - curr_src);
                else
                    next_rest = curr_rest + norm(P(colon(), f) - curr_src);

                next_src = R_st * P(colon(), f) + g_st;

                if ( p1 * p2 > 1e-15 ) {
                    cnt = dir_changes + 1;
                    next_dir = -R_st * curr_dir;
                }
                else {
                    cnt = dir_changes;
                    next_dir = R_st * curr_dir;
                }
                next_dir /= norm(next_dir);
                next_tet = res.second;

                boost::tuple<long long, long long, long long> next_coord((long long)(next_src(0, 0) * pow(10, 4)),
                                                                         (long long)(next_src(1, 0) * pow(10, 4)),
                                                                         (long long)(next_src(2, 0) * pow(10, 4)));
                if ( //vis.find(std::make_pair(next_coord, next_tet)) == vis.end())
                     find(vis.begin(), vis.end(), next_tet) == vis.end() )
                  {
                    //vis.insert(std::make_pair(next_coord, next_tet));
                    vis.insert(next_tet);
                    walk(next_src, next_tet, next_dir, next_rest,
                            orig_handle, cut_tet2tet, inner_type, param_mesh, param_nodes, inv_base, invertible,
                            vis, des, des_tet, depth + 1, cnt, arrive);
                    if ( arrive )
                        return;
                    //vis.erase(vis.find(std::make_pair(next_coord, next_tet)));
                    vis.erase(vis.find(next_tet));
                }
            }
        }
    }
    return;
}

int hex_extractor::extract_mesh(const zjucad::matrix::matrix<size_t>                      &orig_mesh,
                                const zjucad::matrix::matrix<double>                      &orig_nodes,
                                const zjucad::matrix::matrix<size_t>                      &param_mesh,
                                const zjucad::matrix::matrix<double>                      &param_nodes,
                                const zjucad::matrix::matrix<double>                      &inv_base,
                                const zjucad::matrix::matrix<bool>                        &invertible,
                                const std::vector<point_info>                             &int_mesh,
                                const std::vector<boost::tuple<double, double, double>>   &int_nodes,
                                const std::vector<point_info>                             &center_mesh,
                                const std::vector<boost::tuple<double, double, double>>   &center_nodes,
                                const std::map<std::pair<size_t, size_t>, size_t>         &inner_type,
                                zjucad::matrix::matrix<size_t>                            &hex_mesh,
                                zjucad::matrix::matrix<double>                            &hex_nodes )
{
    cout << "[INFO]the program is extracting the hex mesh\n";

    //********************PRE PROCESSING**************************
    std::multimap<boost::tuple<double, double, double>, point_info> int_set;
    RANGE seg;
    for (size_t i = 0; i < int_mesh.size(); ++i)
        int_set.insert(std::make_pair(int_nodes[int_mesh[i].coordinate_idx], int_mesh[i]));

    zjucad::matrix::matrix<size_t> cut_tet2tet(zjucad::matrix::max(param_mesh) + 1);
    cut_tet2tet(param_mesh) = orig_mesh(colon());

    ///build orig mesh handle
    boost::shared_ptr<jtf::mesh::face2tet_adjacent> orig_handle(jtf::mesh::face2tet_adjacent::create(orig_mesh));
    ///merging the points between cut faces
//    merging_vertices_via_cut_face(param_mesh, param_nodes, int_mesh, int_nodes, inner_type, orig_handle, cut_tet2tet);
    //*******************END PRE PROCESSING************************

    zjucad::matrix::matrix<double> dir(3, 8);
    dir(0, 0) = +0.5;  dir(0, 1) = +0.5;  dir(0, 2) = -0.5;  dir(0, 3) = -0.5;
    dir(1, 0) = +0.5;  dir(1, 1) = -0.5;  dir(1, 2) = +0.5;  dir(1, 3) = -0.5;
    dir(2, 0) = +0.5;  dir(2, 1) = +0.5;  dir(2, 2) = +0.5;  dir(2, 3) = +0.5;
    dir(0, 4) = +0.5;  dir(0, 5) = +0.5;  dir(0, 6) = -0.5;  dir(0, 7) = -0.5;
    dir(1, 4) = +0.5;  dir(1, 5) = -0.5;  dir(1, 6) = +0.5;  dir(1, 7) = -0.5;
    dir(2, 4) = -0.5;  dir(2, 5) = -0.5;  dir(2, 6) = -0.5;  dir(2, 7) = -0.5;

    std::vector<std::vector<size_t> > temp_hex_mesh;
    zjucad::matrix::matrix<double> C(3, 1), I(3, 1);

    std::vector<size_t> hex;
    std::set<STATE> vis;

    for (size_t id = 0; id < center_mesh.size(); ++id) {

//        cout << "This is " << id << endl;




        size_t cid = center_mesh[id].coordinate_idx;
        size_t tid = center_mesh[id].tet_idx;
        C(0, 0) = center_nodes[cid].get<0>();
        C(1, 0) = center_nodes[cid].get<1>();
        C(2, 0) = center_nodes[cid].get<2>();
        hex.clear();

//        if ( id == 2285 || id == 8610 || id == 9070 || id == 9071)
//            continue;

        size_t j;
        for (j = 0; j < 8; ++j) {
            vis.clear();
            size_t des_tet;
            boost::tuple<long long, long long, long long> begin((long long)(C(0, 0) * pow(10, 4)),
                                                                (long long)(C(1, 0) * pow(10, 4)),
                                                                (long long)(C(2, 0) * pow(10, 4)));
            //vis.insert(std::make_pair(begin, tid));
            vis.insert(tid);

            bool arrive = false;

            walk(C, tid, dir(colon(), j) / zjucad::matrix::norm(dir(colon(), j)), 0.5 * sqrt(3),
                      orig_handle, cut_tet2tet, inner_type, param_mesh, param_nodes,
                      inv_base, invertible, vis, I, des_tet, 0, 0, arrive);

            if ( arrive ) {
                auto it = int_set.find(boost::make_tuple(I(0, 0), I(1, 0), I(2, 0)));
                if ( it != int_set.end() ) {
                    seg = int_set.equal_range(it->first);
                    for (it = seg.first; it != seg.second; ++it) {
                        if ( it->second.tet_idx == des_tet) {
                            hex.push_back(it->second.coordinate_idx);
                            break;
                        }
                    }
                }
            }
//            else
//                break;
        }

//        if ( cid == 3914 ) {
//            see(C);
//            cout << hex.size() << endl;
//            for (size_t i = 0; i < hex.size(); ++i)
//                cout << hex[i] << " ";
//            cout << endl;
//        }


        if ( hex.size() == 8  ) {
            temp_hex_mesh.push_back(hex);
        }
    }





    //merge the vertices near the cut face
    for (size_t i = 0; i < temp_hex_mesh.size(); ++i)
        for (size_t j = 0; j < 8; ++j)
            temp_hex_mesh[i][j] = father[temp_hex_mesh[i][j]];

    //remove the dupliaceted hex
    std::vector<std::vector<size_t> > hex_set;
    boost::unordered_set<string> used;
    char hex_hash[20];
    for (size_t i = 0; i < temp_hex_mesh.size(); ++i) {
        std::vector<size_t> buffer(temp_hex_mesh[i]);
        sort(buffer.begin(), buffer.end());
        string tag;
        for (size_t j = 0; j < 8; ++j) {
            sprintf(hex_hash, "%lu", buffer[j]);
            tag += string(hex_hash);
            tag += ".";
        }
        if ( used.find(tag) == used.end() ) {
            used.insert(tag);
            hex_set.push_back(temp_hex_mesh[i]);
        }
    }

//    collapse

    generate_hex_feature_points(orig_mesh, orig_nodes, param_mesh, param_nodes, int_mesh, int_nodes, father, feature_points);

//    int cnt = 0;
    while ( true ) {
        if ( !collapse(hex_set, feature_points) )
            break;
    }

    //remove the degenerated hex
    temp_hex_mesh.clear();
    for (size_t i = 0; i < hex_set.size(); ++i) {
        boost::unordered_set<size_t> hash_set;
        for (size_t j = 0; j < 8; ++j)
            hash_set.insert(hex_set[i][j]);
        if ( hash_set.size() == 8 )
            temp_hex_mesh.push_back(hex_set[i]);
    }

    cout << "[INFO]the number of extracted hex is " << temp_hex_mesh.size() << endl;

    hex_mesh.resize(8, temp_hex_mesh.size());
    hex_nodes.resize(3, int_nodes.size());
    for (size_t i = 0; i < hex_mesh.size(2); ++i)
        for (size_t j = 0; j < 8; ++j)
            hex_mesh(j, i) = temp_hex_mesh[i][j];
    for (size_t j = 0; j < int_nodes.size(); ++j) {
        hex_nodes(0, j) = int_nodes[j].get<0>();
        hex_nodes(1, j) = int_nodes[j].get<1>();
        hex_nodes(2, j) = int_nodes[j].get<2>();
    }
    return 0;
}

bool hex_extractor::collapse(std::vector<std::vector<size_t> >     &hex_mesh,
                             const boost::unordered_set<size_t>    &feature_points)
{
    zjucad::matrix::matrix<size_t> enum_face(4, 6);
    enum_face(0, 0) = 0; enum_face(0, 1) = 4; enum_face(0, 2) = 0; enum_face(0, 3) = 2; enum_face(0, 4) = 1; enum_face(0, 5) = 0;
    enum_face(1, 0) = 1; enum_face(1, 1) = 5; enum_face(1, 2) = 1; enum_face(1, 3) = 3; enum_face(1, 4) = 3; enum_face(1, 5) = 2;
    enum_face(2, 0) = 2; enum_face(2, 1) = 6; enum_face(2, 2) = 4; enum_face(2, 3) = 6; enum_face(2, 4) = 5; enum_face(2, 5) = 4;
    enum_face(3, 0) = 3; enum_face(3, 1) = 7; enum_face(3, 2) = 5; enum_face(3, 3) = 7; enum_face(3, 4) = 7; enum_face(3, 5) = 6;

    boost::unordered_set<std::vector<size_t> > face_set;
    boost::unordered_set<std::vector<size_t> >::iterator fptr;
    bool update = false;

    for (size_t id = 0; id < hex_mesh.size(); ++id) {
        for (size_t j = 0; j < 6; ++j) {
            std::vector<size_t> tmp;
            for (size_t k = 0; k < 4; ++k)
                tmp.push_back(hex_mesh[id][enum_face(k, j)]);
            sort(tmp.begin(), tmp.end());
            face_set.insert(tmp);
        }
    }
    std::map<boost::tuple<size_t, size_t, size_t>, size_t> table;
    int arr[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    for ( fptr = face_set.begin(); fptr !=face_set.end(); ++fptr) {
        size_t sum = (*fptr)[0] + (*fptr)[1] + (*fptr)[2] + (*fptr)[3];
        for (size_t i = 0; i < 4; ++i) {
            boost::tuple<size_t, size_t, size_t> arg = boost::make_tuple((*fptr)[arr[i][0]], (*fptr)[arr[i][1]], (*fptr)[arr[i][2]]);
            size_t rest = sum - (*fptr)[arr[i][0]] - (*fptr)[arr[i][1]] - (*fptr)[arr[i][2]];
            auto it = table.find(arg);
            if ( it == table.end() )
                table.insert(std::make_pair(arg, rest));
            else {
                size_t f1 = Find(rest);
                size_t f2 = Find(it->second);
                if ( f1 != f2 ) {
                    update = true;
                    if ( feature_points.find(f1) != feature_points.end() )
                        father[f2] = f1;
                    else
                        father[f1] = f2;
                }
            }
        }
    }
    if ( !update )
        return false;
    for (size_t i = 0; i < hex_mesh.size(); ++i)
        for (size_t j = 0; j < 8; ++j)
            hex_mesh[i][j] = father[hex_mesh[i][j]];
    return true;
}

int hex_extractor::mesh_scaling(zjucad::matrix::matrix<double> &param_nodes, double coeff)
{
    param_nodes *= coeff;
    return 0;
}

int hex_extractor::bound_ith_tet(const size_t                           id,
					  			 const zjucad::matrix::matrix<size_t>   &param_mesh,
	                             const zjucad::matrix::matrix<double>   &param_nodes,
	                             bounding_box                           &box)
{
    zjucad::matrix::matrix<double> tet_nodes = param_nodes(colon(), param_mesh(colon(), id));
    box.xmin = zjucad::matrix::min(tet_nodes(0, colon()));
    box.ymin = zjucad::matrix::min(tet_nodes(1, colon()));
    box.zmin = zjucad::matrix::min(tet_nodes(2, colon()));
    box.xmax = zjucad::matrix::max(tet_nodes(0, colon()));
    box.ymax = zjucad::matrix::max(tet_nodes(1, colon()));
    box.zmax = zjucad::matrix::max(tet_nodes(2, colon()));
    return 0;
}

int hex_extractor::mapping_to_orig_space(const zjucad::matrix::matrix<size_t>    &orig_mesh,
                                         const zjucad::matrix::matrix<double>    &orig_node,
                                         const std::vector<point_info>           &int_mesh,
                                         const zjucad::matrix::matrix<size_t>    &media_mesh,
                                         const zjucad::matrix::matrix<double>    &media_nodes,
                                         zjucad::matrix::matrix<size_t>          &final_mesh,
                                         zjucad::matrix::matrix<double>          &final_nodes)
{
    cout << "[INFO]the program is mapping the mesh to original space\n";
    final_mesh.resize(8, media_mesh.size(2));
    final_nodes.resize(3, media_nodes.size(2));
    zjucad::matrix::matrix<double> orig_tet_coord(3, 4);
    zjucad::matrix::matrix<double> orig_int_coord(3, 1);
    for (size_t i = 0; i < int_mesh.size(); ++i) {
        size_t tet_id = int_mesh[i].tet_idx;
        size_t cor_id = int_mesh[i].coordinate_idx;
        orig_tet_coord = orig_node(colon(), orig_mesh(colon(), tet_id));
        orig_int_coord = orig_tet_coord * int_mesh[i].b_coeff;
        final_nodes(colon(), cor_id) = orig_int_coord;
    }
    final_mesh = media_mesh;
    mesh_orderring(final_nodes, final_mesh);

    const char* path = "/home/chenjiong/usr/WorkSpace/volume_frame/branches/for_cj/src/hex_extractor/dat/fea_point_hex.vtk";
    dump_feature_point(path, feature_points, final_nodes);
}

bool hex_extractor::line_cross_tet(const size_t                            id,
                                   const zjucad::matrix::matrix<size_t>    &param_mesh,
                                   const zjucad::matrix::matrix<double>    &param_nodes,
                                   const zjucad::matrix::matrix<double>    &inv_base,
                                   const zjucad::matrix::matrix<bool>      &invertible,
                                   const zjucad::matrix::matrix<double>    &src,
                                   const zjucad::matrix::matrix<double>    &des,
                                   const int                               type,
                                   zjucad::matrix::matrix<double>          &cross_point,
                                   zjucad::matrix::matrix<size_t>          &face_idx,
                                   zjucad::matrix::matrix<int>             &flag)
{
    bool mark = false;
    zjucad::matrix::matrix<double> vertex(3, 4);
    zjucad::matrix::matrix<size_t> v_id(4, 1);
    vertex = param_nodes(colon(), param_mesh(colon(), id));
    v_id = param_mesh(colon(), id);

    zjucad::matrix::matrix<double> arg(4, 1);
    int a[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};

    zjucad::matrix::matrix<double> v0(3, 1), v1(3, 1), v2(3, 1), x(3, 1);
    for (int i = 0; i < 4; i++) {
        v0 = vertex(colon(), a[i][0]);
        v1 = vertex(colon(), a[i][1]);
        v2 = vertex(colon(), a[i][2]);
        if ( line_cross_plane(src, des, v0, v1, v2, x) ) {
            bool expr0 = std::fabs(norm(x - src) + norm(des - x) - norm(des - src)) < EPSILON
                         && point_in_triangle(x, v0, v1, v2);
            bool expr1 = dot(des - src, x - src) > 0
                         && point_in_triangle(x, v0, v1, v2);
            bool expr = type == 0 ? expr0 : expr1;
            if ( expr ) {
                flag(i, 0) = 1;
                mark = true;
                cross_point(colon(), i) = x;
                face_idx(0, i) = v_id(a[i][0], 0);
                face_idx(1, i) = v_id(a[i][1], 0);
                face_idx(2, i) = v_id(a[i][2], 0);
            }
            else
                flag(i, 0) = 0;
        }
    }
    return mark;
}

int hex_extractor::line_cross_plane(const zjucad::matrix::matrix<double>     &p0,
                                    const zjucad::matrix::matrix<double>     &p1,
                                    const zjucad::matrix::matrix<double>     &v0,
                                    const zjucad::matrix::matrix<double>     &v1,
                                    const zjucad::matrix::matrix<double>     &v2,
                                    zjucad::matrix::matrix<double>           &x )
{
    zjucad::matrix::matrix<double> tan(3, 1);
    tan = p1 - p0;
    zjucad::matrix::matrix<double> nor(3, 1);
    nor = zjucad::matrix::cross(v1 - v0, v2 - v0);

    double A = zjucad::matrix::dot(tan, nor);
    if ( std::fabs(A) < EPSILON ) {
        if ( std::fabs(dot(cross(v1 - v0, v2 - v0), p0 - v0)) < EPSILON &&
             std::fabs(dot(cross(v1 - v0, v2 - v0), p1 - v0)) < EPSILON ) {
            x = p0;
            return 2;
        }
        return 0;
    }
    double b = -zjucad::matrix::dot(p0 - v0, nor);
    double t = b / A;
    x = tan * t + p0;
    return 1;
}

//int hex_extractor::mesh_orderring(const std::vector<boost::tuple<double, double, double> >  &int_nodes,
//                                  std::vector<std::vector<size_t> >                         &hex_mesh)
//{
//    for (size_t i = 0; i < hex_mesh.size(); ++i) {
//        std::swap(hex_mesh[i][2], hex_mesh[i][3]);
//        std::swap(hex_mesh[i][6], hex_mesh[i][7]);
//        zjucad::matrix::matrix<double> E01(3, 1), E02(3, 1), E04(3, 1);
//        E01(0, 0) = int_nodes[hex_mesh[i][1]].get<0>() - int_nodes[hex_mesh[i][0]].get<0>();
//        E01(1, 0) = int_nodes[hex_mesh[i][1]].get<1>() - int_nodes[hex_mesh[i][0]].get<1>();
//        E01(2, 0) = int_nodes[hex_mesh[i][1]].get<2>() - int_nodes[hex_mesh[i][0]].get<2>();

//        E02(0, 0) = int_nodes[hex_mesh[i][2]].get<0>() - int_nodes[hex_mesh[i][0]].get<0>();
//        E02(1, 0) = int_nodes[hex_mesh[i][2]].get<1>() - int_nodes[hex_mesh[i][0]].get<1>();
//        E02(2, 0) = int_nodes[hex_mesh[i][2]].get<2>() - int_nodes[hex_mesh[i][0]].get<2>();

//        E04(0, 0) = int_nodes[hex_mesh[i][4]].get<0>() - int_nodes[hex_mesh[i][0]].get<0>();
//        E04(1, 0) = int_nodes[hex_mesh[i][4]].get<1>() - int_nodes[hex_mesh[i][0]].get<1>();
//        E04(2, 0) = int_nodes[hex_mesh[i][4]].get<2>() - int_nodes[hex_mesh[i][0]].get<2>();

//        if ( zjucad::matrix::dot(cross(E01, E02), E04) < 0 ) {
//            std::swap(hex_mesh[i][1], hex_mesh[i][2]);
//            std::swap(hex_mesh[i][5], hex_mesh[i][6]);
//        }
//    }
//    return 0;
//}

int hex_extractor::mesh_orderring(const zjucad::matrix::matrix<double>  &final_nodes,
                                  zjucad::matrix::matrix<size_t>        &final_mesh)
{
    zjucad::matrix::matrix<double> E01(3, 1), E02(3, 1), E04(3, 1);
    for (size_t i = 0; i < final_mesh.size(2); ++i) {
        E01 = final_nodes(colon(), final_mesh(1, i)) - final_nodes(colon(), final_mesh(0, i));
        E02 = final_nodes(colon(), final_mesh(2, i)) - final_nodes(colon(), final_mesh(0, i));
        E04 = final_nodes(colon(), final_mesh(4, i)) - final_nodes(colon(), final_mesh(0, i));

        if ( zjucad::matrix::dot(cross(E01, E02), E04) < 0 ) {
            for (size_t j = 0; j < 4; j++)
                std::swap(final_mesh(j, i), final_mesh(j + 4, i));
        }

    }
    return 0;
}


double hex_extractor::oriented_area(const size_t                          id,
                                    const zjucad::matrix::matrix<size_t>  &param_mesh,
                                    const zjucad::matrix::matrix<double>  &param_nodes)
{
    zjucad::matrix::matrix<double> vertex(3, 4);
    vertex = param_nodes(colon(), param_mesh(colon(), id));
    return zjucad::matrix::dot(zjucad::matrix::cross(vertex(colon(), 1) - vertex(colon(), 0),
                                                     vertex(colon(), 2) - vertex(colon(), 0)),
                                                     vertex(colon(), 3) - vertex(colon(), 0));
}

size_t inline hex_extractor::Find(const size_t cid)
{
    return father[cid] == cid ? cid : father[cid] = Find(father[cid]);
}

void inline hex_extractor::Union(const size_t p, const size_t q)
{
    size_t fp = Find(p);
    size_t fq = Find(q);
    if ( fp != fq )
        father[fq] = fp;
}

double inline hex_extractor::cast(const double x, const double eps)
{
    return std::fabs(x - round(x)) < eps ? ( round(x) == -0 ? 0 : round(x) ) : x;
}

}




