#ifndef CJ_FEATURE_H
#define CJ_FEATURE_H

#include <fstream>
#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/mesh.h>
#include "../common/vtk.h"

namespace extractor {

using namespace zjucad::matrix;

struct point_info;

#define PI 3.141592653589793
#define FEATURE_ANGLE PI/4

typedef struct point3d {
    zjucad::matrix::matrix<double> coord;
    size_t tid;
    size_t cid;
//    bool operator<(const point3d &b) const {
//        return zjucad::matrix::norm(coord - b.coord) > 1 - 1e-8 && zjucad::matrix::norm(coord) < zjucad::matrix::norm(b.coord);
//    }
//    bool operator==(const point3d &b) const {
//        return zjucad::matrix::norm(coord - b.coord) < 1 - 1e-8 && tid == b.tid;
//    }
}point3d;

bool equal_point3d(const point3d &a, const point3d &b)
{
    return zjucad::matrix::norm(a.coord - b.coord) <= 1 - 1e-8 && a.tid == b.tid;
}

int extract_orig_feature_line(const zjucad::matrix::matrix<size_t>              &orig_mesh,
                              const zjucad::matrix::matrix<double>              &orig_nodes,
                              boost::unordered_set<std::pair<size_t, size_t>>   &f_line)
{
    boost::shared_ptr<jtf::mesh::face2tet_adjacent> f2t(jtf::mesh::face2tet_adjacent::create(orig_mesh));
    zjucad::matrix::matrix<size_t> tet_surface;
    jtf::mesh::get_outside_face(*f2t, tet_surface, true);

    boost::shared_ptr<jtf::mesh::edge2cell_adjacent> e2c(jtf::mesh::edge2cell_adjacent::create(tet_surface));
    boost::unordered_set<std::pair<size_t, size_t>> vis;

    for (size_t i = 0; i < tet_surface.size(2); ++i) {
        for (size_t j = 0; j < tet_surface.size(1); ++j) {
            size_t p = tet_surface(j, i);
            size_t q = tet_surface((j + 1) % tet_surface.size(1), i);
            if ( vis.find(std::make_pair(p, q)) != vis.end() )
                continue;
            std::pair<size_t, size_t> res = e2c->query(p, q);
            if ( !e2c->is_boundary_edge(res) ) {
                size_t s = res.first;
                size_t t = res.second;
                zjucad::matrix::matrix<double> nor_s(3, 1), nor_t(3, 1);
                nor_s = cross(orig_nodes(colon(), tet_surface(1, s)) - orig_nodes(colon(), tet_surface(0, s)),
                              orig_nodes(colon(), tet_surface(2, s)) - orig_nodes(colon(), tet_surface(1, s)));
                nor_t = cross(orig_nodes(colon(), tet_surface(1, t)) - orig_nodes(colon(), tet_surface(0, t)),
                              orig_nodes(colon(), tet_surface(2, t)) - orig_nodes(colon(), tet_surface(1, t)));
                double angle = dot(nor_s, nor_t) / norm(nor_s) / norm(nor_t);
                if ( angle < cos(FEATURE_ANGLE) ) {
                    f_line.insert(std::make_pair(p, q));
                    vis.insert(std::make_pair(p, q));
                    vis.insert(std::make_pair(q, p));
                }
            }
        }
    }

    std::cout << "[INFO]the number of feature lines : " << f_line.size() << std::endl;
    return 0;
}

//int get_a_chain(size_t                                curr,
//                zjucad::matrix::matrix<double>        curr_dir,
//                const zjucad::matrix::matrix<double>  &orig_nodes,
//                const std::vector<size_t>             &u,
//                const std::vector<size_t>             &v,
//                const std::vector<int>                &first,
//                const std::vector<int>                &next,
//                std::vector<size_t>                   &one_chain,
//                boost::unordered_set<size_t>          &vis,
//                std::vector<std::vector<size_t>>      &f_chain)
//{
//    vis.insert(curr);
//    one_chain.push_back(curr);
//    bool flag = false;
//    for (size_t e = first[curr]; e != -1; e = next[e]) {
//        size_t next = v[e];
//        zjucad::matrix::matrix<double> next_dir(3, 1);
//        next_dir = orig_nodes(colon(), next) - orig_nodes(colon(), curr);
//        if ( vis.find(next) == vis.end() && zjucad::matrix::dot(curr_dir, next_curr) < cos(PI / 6) ) {
//            get_a_chain(next, next_dir, orig_nodes, u, v, first, next, one_chain, vis, f_chain);
//            flag = true;
//        }
//    }
//    if ( !flag ) {
//        f_chain.push_back(one_chain);
//        one_chain.
//    }
//}

//int generate_feature_graph(const zjucad::matrix::matrix<double>                   &orig_nodes,
//                           const boost::unordered_set<std::pair<size_t, size_t>>  &f_line,
//                           std::vector<std::vector<size_t>>                       &f_chain)
//{
//    std::vector<size_t> u, v;
//    std::vector<int> first(orig_nodes.size(2), -1);
//    std::vector<int> next(f_line.size(), -1);

//    size_t i = 0;
//    for (auto it = f_line.first; it != f_line.end(); ++it, ++i) {
//        u.push_back(it->first);
//        v.push_back(it->second);
//        next[i] = first[u[i]];
//        first[u[i]] = i;
//    }
//}

void build_cut_points_set(const std::vector<point_info>                            &int_mesh,
                          const std::vector<boost::tuple<double, double, double>>  &int_nodes,
                          std::vector<point3d>                                     &table)
{
    point3d temp;
    temp.coord.resize(3, 1);
    for (size_t i = 0; i < int_mesh.size(); ++i) {
        temp.cid = int_mesh[i].coordinate_idx;
        temp.tid = int_mesh[i].tet_idx;
        temp.coord(0, 0) = int_nodes[temp.cid].get<0>();
        temp.coord(1, 0) = int_nodes[temp.cid].get<1>();
        temp.coord(2, 0) = int_nodes[temp.cid].get<2>();
        table.push_back(temp);
    }
}

//int binary_search(std::vector<size_t> &table, size_t l, size_t r, point3d *key)
//{
//    while ( l <= r ) {
//        size_t mid = (l + r) >> 1;
//        if ( table[mid] == key )
//            return mid;
//        else if ( key <  ) {
//            r  = mid - 1;
//        }
//        else
//            l = mid + 1;
//    }
//    return -1;
//}

bool point_on_line(const zjucad::matrix::matrix<double>  &p,
                   const zjucad::matrix::matrix<double>  &a,
                   const zjucad::matrix::matrix<double>  &b,
                   double                                &lamda)
{
    double PA = zjucad::matrix::norm(p - a);
    double PB = zjucad::matrix::norm(p - b);
    double AB = zjucad::matrix::norm(b - a);
    if ( AB < 1e-9 ) {
        lamda = 0;
        return PA < 1e-9;
    } else {
        lamda = PA / PB;
        return std::fabs(PA + PB - AB) < 1e-9;
    }
}

int extract_cut_feature_points(const zjucad::matrix::matrix<size_t>                     &orig_mesh,
                               const zjucad::matrix::matrix<size_t>                     &cut_mesh,
                               const zjucad::matrix::matrix<double>                     &cut_nodes,
                               const boost::unordered_set<std::pair<size_t, size_t>>    &f_line,
                               const std::vector<point_info>                            &i_mesh,
                               const std::vector<boost::tuple<double, double, double>>  &i_nodes,
                               const std::vector<size_t>                                &father,
                               boost::unordered_set<size_t>                             &feature_points)
{
    zjucad::matrix::matrix<size_t> cut_surface;
    boost::shared_ptr<jtf::mesh::face2tet_adjacent> f2c(jtf::mesh::face2tet_adjacent::create(cut_mesh));
    jtf::mesh::get_outside_face(*f2c, cut_surface, true);
    boost::shared_ptr<jtf::mesh::face2tet_adjacent> orig_f2c(jtf::mesh::face2tet_adjacent::create(orig_mesh));

    zjucad::matrix::matrix<size_t> cut_tet2tet(zjucad::matrix::max(cut_mesh) + 1);
    cut_tet2tet(cut_mesh) = orig_mesh(colon());

    std::vector<point3d> table;
    build_cut_points_set(i_mesh, i_nodes, table);

    for (size_t i = 0; i < cut_surface.size(2); ++i) {
        std::pair<size_t, size_t> tet_st = orig_f2c->query(cut_tet2tet(cut_surface(0, i)),
                                                           cut_tet2tet(cut_surface(1, i)),
                                                           cut_tet2tet(cut_surface(2, i)));
        size_t curr_tet = std::min(tet_st.first, tet_st.second);
        for (size_t j = 0; j < cut_surface.size(1); ++j) {
            size_t p = cut_surface(j, i);
            size_t q = cut_surface((j + 1) % cut_surface.size(1), i);
            size_t op = cut_tet2tet(p);
            size_t oq = cut_tet2tet(q);
            if ( f_line.find(std::make_pair(op, oq)) != f_line.end() ||
                 f_line.find(std::make_pair(oq, op)) != f_line.end() ) {
                zjucad::matrix::matrix<double> bbox(3, 2), I(3, 1);
                for (size_t k = 0; k < 3; ++k) {
                    bbox(k, 0) = std::min(cut_nodes(k, p), cut_nodes(k, q));
                    bbox(k, 1) = std::max(cut_nodes(k, p), cut_nodes(k, q));
                }
                for (int u = (int)floor(bbox(0, 0)); u <= (int)ceil(bbox(0, 1)); ++u) {
                for (int v = (int)floor(bbox(1, 0)); v <= (int)ceil(bbox(1, 1)); ++v) {
                for (int w = (int)floor(bbox(2, 0)); w <= (int)ceil(bbox(2, 1)); ++w) {
                    I(0, 0) = u; I(1, 0) = v; I(2, 0) = w;
                    double lamda;
                    if ( point_on_line(I, cut_nodes(colon(), p), cut_nodes(colon(), q), lamda) ) {
                        point3d curr{I, curr_tet, 0};
                        auto found =
                                std::find_if(table.begin(), table.end(), std::bind(&equal_point3d, std::placeholders::_1, curr));
                        if ( found != table.end() )
                            feature_points.insert(father[found->cid]);
//                        while ( found != table.end() ) {
//                            feature_points.insert(father[found->cid]);
//                            ++found;
//                            found = std::find_if(found, table.end(), std::bind(&equal_point3d, std::placeholders::_1, curr));
//                        }
//                        while ( found != table.end() ) {
//                            ++found;
//                            found = std::find_if(found, table.end(), std::bind(&equal_point3d, std::placeholders::_1, curr));
//                            if ( found != table.end() )
//                                feature_points.insert(father[found->cid]);
//                        }
                    }
                }}}
            }
        }
    }
    return 0;
}

void generate_hex_feature_points(const zjucad::matrix::matrix<size_t>                     &orig_mesh,
                                 const zjucad::matrix::matrix<double>                     &orig_nodes,
                                 const zjucad::matrix::matrix<size_t>                     &cut_mesh,
                                 const zjucad::matrix::matrix<double>                     &cut_nodes,
                                 const std::vector<point_info>                            &int_mesh,
                                 const std::vector<boost::tuple<double, double, double>>  &int_nodes,
                                 const std::vector<size_t>                                &father,
                                 boost::unordered_set<size_t>                             &feature_points)
{
    std::cout << "[INFO]the program is generating the feature points of hex mesh\n";
    boost::unordered_set<std::pair<size_t, size_t>> f_line;
    extract_orig_feature_line(orig_mesh, orig_nodes, f_line);
    extract_cut_feature_points(orig_mesh, cut_mesh, cut_nodes, f_line, int_mesh, int_nodes, father, feature_points);

    //treat points sharing the same param coord with extracted feature points above as new features
    boost::unordered_set<std::vector<int>> feature_coord;
    for (auto it = feature_points.begin(); it != feature_points.end(); ++it) {
        std::vector<int> coord(3);
        coord[0] = (int)int_nodes[*it].get<0>();
        coord[1] = (int)int_nodes[*it].get<1>();
        coord[2] = (int)int_nodes[*it].get<2>();
        feature_coord.insert(coord);
    }
    for (size_t id = 0; id < int_mesh.size(); ++id) {
        std::vector<int> coord(3);
        coord[0] = (int)int_nodes[father[int_mesh[id].coordinate_idx]].get<0>();
        coord[1] = (int)int_nodes[father[int_mesh[id].coordinate_idx]].get<1>();
        coord[2] = (int)int_nodes[father[int_mesh[id].coordinate_idx]].get<2>();
        if ( feature_coord.find(coord) != feature_coord.end() )
            feature_points.insert(father[int_mesh[id].coordinate_idx]);
    }
    std::cout << "[INFO]the number of feature points in hex mesh is " << feature_points.size() << std::endl;
}

void dump_feature_point(const char                            *path,
                        const boost::unordered_set<size_t>    &feature_points,
                        const zjucad::matrix::matrix<double>  &hex_nodes)
{
    std::ofstream out(path);
    zjucad::matrix::matrix<double> f_mesh;
    if ( feature_points.size() == 0 )
        return;
    f_mesh.resize(1, feature_points.size());
    size_t id = 0;
    for (auto it = feature_points.begin(); it != feature_points.end(); ++it)
        f_mesh(0, id++) = *it;
    point2vtk(out, &hex_nodes[0], hex_nodes.size(2), &f_mesh[0], f_mesh.size(2));

}



}

#endif
