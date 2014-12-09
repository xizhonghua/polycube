/* detail control param
/usr/bin/time ../../bin/polycube prog=polycube3 linear_solver/type=direct linear_solver/name=cholmod iter_w=20 output="a.tet" iter=100 tet=../../dat/kitty-4.8k.mesh-split-surface-tet.mesh normal_align_w=5e-2 L1_sqrt_eps=5e-1 epsg=1e-3 adj_normal_w=2e1
 */
#include <boost/property_tree/ptree.hpp>

#include <string>
#include <fstream>
#include <numeric>
#include <hjlib/function/func_aux.h>
#include <hjlib/math/polar.h>
#include <hjlib/sparse/sparse.h>


#include <zjucad/optimizer/optimizer.h>
#include <zjucad/ptree/ptree.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"
#include "../common/util.h"
#include "../common/IO.h"
#include <jtflib/mesh/mesh.h>
#include <jtflib/optimizer/optimizer.h>

#include "polycube_surface_func.h"
#include "smooth_L1.h"

#include "../numeric/util.h"
#include "../mesh_func/tri-area-normal-func.h"
#include "polycube_surface_func.h"

#include "../mesh_func/tri-area.h"
#include "tet_func.h"
#include "../mesh_func/tri-normal.h"
#include "../mesh_func/tri-area-normal.h"
#include "util.h"
#include "quality.h"

#include "adaptive.h"
#include "../tetmesh/subdivide_tet.h"

using namespace std;
using boost::property_tree::ptree;
using namespace hj::function;
using namespace zjucad::matrix;

double L1_sqrt_eps;
double arap_w;
double fix_zero_w;
double adj_normal_w;

//double trigger_div_L::k_ = 1.0;
//double trigger_div_L::L_ = 1.0;

//hj::function::function_t<double, int32_t> *
//build_polycube_edge(const  matrix<double> &node,
//                    const jtf::mesh::edge2cell_adjacent &ea,
//                    const matrixst &faces,
//                    const double normal_diff_w);

function_t<double, int32_t> *tetmesh_arap_for_diagnose = 0;

// move surface point to the center of one-ring neighborhood by a little
void pre_smooth(matrix<double> &node, const matrix<double> &faces)
{
    matrix<size_t> nb_num = zeros<size_t>(node.size(2), 1);
    matrix<double> ct = zeros<double>(3, node.size(2));
    for(size_t fi = 0; fi < faces.size(2); ++fi) {
        nb_num(faces(colon(), fi)) += 1;
        for(int i = 0; i < 3; ++i) {
            ct(colon(), faces(i, fi)) += node(colon(), faces((i+1)%3, fi)) + node(colon(), faces((i+2)%3, fi));
        }
    }
    const double alpha = 0.3;
    for(size_t ni = 0; ni < node.size(2); ++ni) {
        if(nb_num[ni] == 0) continue;
        node(colon(), ni) = (1-alpha)*node(colon(), ni) + alpha*ct(colon(), ni) / nb_num[ni] / 2;
    }
}

int polycube3(ptree &pt)
{
    jtf::mesh::meshes tm;
    if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("tet.value").c_str(), &tm.node_, &tm.mesh_))
        return __LINE__;

    pt.put("L1_sqrt_eps.desc", "L1 sqrt eps, default is 1e-1.");
    pt.put("normal_align_w.desc", "surface L1 weight, default is 1e2.");
    pt.put("adj_normal_w.desc", "surface smoothness weight, default is 0.");
    pt.put("anti_flip_w.desc", "surface anti-flip weight, default is 1e-2.");
    pt.put("arap_w.desc", "arap weight, default is 1.");
    pt.put("iter_w.desc", "iterative weight, default is 1");
    pt.put("fix_zero_w.desc", "fix first node to zero, default is 1");
    pt.put("div_L_w.desc", "div eps by L, default is sqrt2");
    pt.put("area_preserve_w.desc", "polycube area preserving, default is 1.0");
    pt.put("alpha_times.desc", "alpha times x, x default is 2");

    const double alpha_times = pt.get<double>("alpha_times.value",2);
    L1_sqrt_eps = pt.get<double>("L1_sqrt_eps.value", 0.5);
    arap_w = pt.get<double>("arap_w.value",1);
    double normal_align_w = pt.get<double>("normal_align_w.value", 5e-2);
    double anti_flip_w = pt.get<double>("anti_flip_w.value", 1e-2);
    adj_normal_w = pt.get<double>("adj_normal_w.value", 0.0);
    fix_zero_w = pt.get<double>("fix_zero_w.value", 1);

    const size_t iter_w = pt.get<size_t>("iter_w.value",1);

    const double div_L = pt.get<double>("div_L.value",sqrt(2.0));
    const double epsg = pt.get<double>("epsg.value", 1e-6);
    const matrix<double> ct = tm.node_*ones<double>(tm.node_.size(2), 1)/tm.node_.size(2);

    matrix<double> node = tm.node_;
    if(pt.get_child_optional("init.value")) {
        jtf::mesh::meshes init;
        if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("init.value").c_str(), &init.node_, &init.mesh_))
            return __LINE__;
        if(max(abs(init.mesh_ - tm.mesh_)) != 0) {
            cerr << "#incompatible init value" << endl;
            return __LINE__;
        }
        node = init.node_;
        cerr << "# use init value." << endl;
        // cerr << tm.node(colon(), colon(0, 5))
        //      << node(colon(), colon(0, 5)) << endl;
        // cerr << "# pre-smooth" << endl;
        // pre_smooth(node, faces);
    }
    if(pt.get_child_optional("subdivide.value")) {
        cerr << "# beg subdivide: " << tm.node_.size(2) << " " << tm.mesh_.size(2) << endl;
        matrix<double> *p_node[2] = {&tm.node_, &node};
        matrix<size_t> c_tet, node_parent;
        matrix<double> c_node;
        time_t beg = clock();
        boundary_uniform_subdivide(tm.node_, tm.mesh_, c_node, c_tet, node_parent);
        cout << "# subd time: " << (clock()-beg)/double(CLOCKS_PER_SEC) << endl;
        tm.mesh_ = c_tet;
        tm.node_ = c_node;

        subdivide_top2geo(node, node_parent, c_node);
        node = c_node;
        ofstream ofs("subdivide.vtk");
        tet2vtk(ofs, &tm.node_[0], tm.node_.size(2), &tm.mesh_[0], tm.mesh_.size(2));
        jtf::mesh::tet_mesh_write_to_zjumat("subdivide.tet", &tm.node_, &tm.mesh_);
        cerr << "# end subdivide: " << tm.node_.size(2) << " " << tm.mesh_.size(2) << endl;
    }

    unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
    if(!fa.get()){
        cerr  << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
        return __LINE__;
    }
    matrix<size_t> faces;
    get_outside_face(*fa, faces);

    cerr << "#tet size: " << tm.node_.size(2) << " " << tm.mesh_.size(2) << " " << faces.size(2) << endl;

    matrix<double> areas = zeros<double>(faces.size(2),1);
    for(size_t fi = 0; fi < faces.size(2); ++fi){
        areas[fi] = jtf::mesh::cal_face_area(faces(colon(),fi), tm.node_);
    }

    const double total_area = std::accumulate(areas.begin(), areas.end(), 0.0);
    cerr << "total_area: " << total_area << endl;
    matrix<double> R = eye<double>(3);
    ostringstream vtk_path_pref;
    vtk_path_pref << normal_align_w << "-" << L1_sqrt_eps << "-" << adj_normal_w << "-";

    double opt_time = 0;
    double remesh_time = 0;
    double assemble_time = 0;
    cerr << "# [info] initial polycube_L1_area_quality " << polycube_L1_area_quality(&node[0], node.size(2), faces) << endl;
    for(size_t i = 0; i < iter_w; ++i){
        clock_t beg = clock();
        cout << "begin iter: " << normal_align_w << " " << L1_sqrt_eps << endl;

        matrix<double> areas_k = zeros<double>(faces.size(2),1);
        for(size_t fi = 0; fi < faces.size(2); ++fi){
            areas_k[fi] = jtf::mesh::cal_face_area(faces(colon(),fi),node);
        }

        if(1) { // opt global R
            R = eye<double>(3);
            shared_ptr<jtf::function::functionN1_t<double,int32_t> > func(build_polycube_rot_func2(node, faces,areas_k));
            jtf::optimize(*func, R, pt, nullptr, nullptr, nullptr);
            cout << R << endl;
            hj::polar3d p;
            p(R);
            cout << R << endl;
            node = temp(R*node);
        }
        shared_ptr<function_t<double, int32_t> > func
                (build_polycube_function(tm.node_, node, node(colon(), 0), tm.mesh_, faces, tetmesh_arap_for_diagnose,
                                         nullptr, nullptr, adj_normal_w, 0, 0));
//            (build_polycube_function3(tm.node_, node, node(colon(), 0), tm.mesh_, faces));

        cerr << "# [info] add polycube surface normal L1 function, wegiht "
             << normal_align_w << endl;
        cerr << "# [info] L1_sqrt_eps = " << L1_sqrt_eps << endl;

        vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > sum;
        vector<pair<jtf::function::functionN1_t<double,int32_t> *, double> > wf;
        sum.push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(jtf::function::least_square_warpper(func)));

        shared_ptr<jtf::function::functionN1_t<double,int32_t> > arap_for_diagnose
                (jtf::function::least_square_warpper(*tetmesh_arap_for_diagnose));
        wf.push_back(make_pair(arap_for_diagnose.get(), 1));
        if(normal_align_w > 0) {
            sum.push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(build_smooth_L1_area_normal(node, faces, areas_k, normal_align_w)));
            cerr << "1: " << sum[1]->dim() << endl;;
            wf.push_back(make_pair(sum[1].get(), normal_align_w));
        }

        if(anti_flip_w > 0) {// use unflip
            matrix<double> weight;
            shared_ptr<function_t<double, int32_t> > adj_normal
                    (build_adj_normal_func(node, faces, weight));
            weight *= anti_flip_w;
            sum.push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                            jtf::function::neg_log_warpper(adj_normal, weight)));
            cerr << "# use anti-flip with weight: " << anti_flip_w << endl;
        }

        shared_ptr<jtf::function::functionN1_t<double,int32_t> > target(new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD>(sum));
        wf.push_back(make_pair(target.get(), (1+normal_align_w)));

        shared_ptr<jtf::function::functionN1_t<double,int32_t> > constraint(new area_sum(tm.node_.size(2), faces, total_area));
        unnormalized_normal_quality_checker cb(node, tm.mesh_, faces, wf, *constraint);

        double assemble_time_once = double(clock()-beg)/CLOCKS_PER_SEC ;
        assemble_time += assemble_time_once;

        clock_t opt_beg = clock();

        pt.put("epsg.value", epsg*(1+normal_align_w));

        int rtn = -1;
        if(normal_align_w > 0){
            vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > constraint_vec;
            constraint_vec.push_back(constraint);
          rtn = jtf::optimize(*target, node, pt, &constraint_vec,nullptr, &cb);
          }
        else
          rtn = jtf::optimize(*target, node, pt, nullptr, nullptr, nullptr);

        double opt_time_once = double(clock()-opt_beg)/CLOCKS_PER_SEC ;
        opt_time += opt_time_once;

        cout << "# a loop: " << double(clock()-beg)/CLOCKS_PER_SEC << endl;

        clock_t remesh_beg = clock();
       // remove_degenerated_face_of_tet(tm.mesh_, tm.node_, node, faces, 5);

        double remesh_time_once = double(clock()-remesh_beg)/CLOCKS_PER_SEC ;

        remesh_time += remesh_time_once;

        const matrix<double> cur_ct = node*ones<double>(tm.node_.size(2), 1)/(1.0*tm.node_.size(2));
        node += (ct-cur_ct)*ones<double>(1, node.size(2));

//         if(rtn == -1){ // interior fail, need to reoptimize the shape with current weighting
//            static vector<size_t> vis_log;

//            if(vis_log.empty() || vis_log.back() != i){ // meet interior point at new step
//                vis_log.clear();
//                vis_log.push_back(i);
//                --i;
//                continue;
//            }else{ // meet interior point many times,
//                vis_log.push_back(i);
//                if(vis_log.size() < 3){ // if stuck in this step for more than three times, skip on
//                    --i;
//                    continue;
//                }
//            }
//        }

        { // visualize
            ostringstream vtk_path;
            vtk_path << vtk_path_pref.str() << i << ".vtk";
            ofstream ofs(vtk_path.str().c_str());
            tet2vtk(ofs, &node[0], node.size(2), &tm.mesh_[0], tm.mesh_.size(2));
        }
        if(polycube_L1_area_quality(&node[0], node.size(2), faces) < 1e-3) { // the global condition
            cout << "global converge." << endl;
            break;
        }

        normal_align_w *= alpha_times;
        L1_sqrt_eps /= div_L;
        if(L1_sqrt_eps < 1e-2)
            L1_sqrt_eps = 1e-2;
    }

    cerr << "# [info] finial polycube_quality " << polycube_L1_area_quality(&node[0], node.size(2), faces) << endl;
    jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("output.value").c_str(), &node, &tm.mesh_);

    jtf::mesh::tet_mesh_write_to_zjumat("after_collapse_output.tet", &tm.node_, &tm.mesh_);

    cerr << "# [info] total opt time " << opt_time << endl;
    cerr << "# [info] total assemble time " << assemble_time << endl;
    cerr << "# [info] total remesh time " << remesh_time << endl;

    cerr << "success." << endl;
    return 0;
}
