#include <zjucad/ptree/ptree.h>
#include <zjucad/matrix/matrix.h>
#include "../common/zyz.h"
#include "../common/transition_type.h"
#include "../common/IO.h"
#include "../common/visualize_tool.h"
#include "../common/vtk.h"
#include "../tetmesh/tetmesh.h"
#include "../hex_param/global_alignment.h"
#include "../hex_param/find_singularities.h"
#include "../common/transition.h"
#include "../numeric/util.h"

using namespace std;
using boost::property_tree::ptree;
using namespace zjucad::matrix;

#define ASSERT_FUNC(f,str) if(f) throw std::logic_error(str);

int find_nearest_axis(const matrix<double> & frame,const matrix<double> & direction)
{
    vector<pair<double, size_t> > err;
    for(size_t di = 0; di < 3; ++di){
        err.push_back(make_pair(dot(frame(colon(),di),direction),2*di));
        err.push_back(make_pair(-1*dot(frame(colon(),di),direction),2*di+1));
    }
    auto it = max_element(err.begin(), err.end());
    return it->second;
}

matrix<double> get_axis_dir(const matrix<double> & frame, const size_t axis)
{
    assert(axis < 6);
    return (axis%2==1?-1:1)*frame(colon(),axis/2);
}

double calculate_dihedral_angle_degree_param(
        const pair<size_t,size_t> & tri_pair,
        const jtf::mesh::face2tet_adjacent &fa,
        const matrix<matrix<double> > & frame,
        const matrix<double> &n1, const matrix<double> & n2,
        const vector<size_t> & tet_loop)
{
    pair<size_t,size_t> tet_pair;
    const pair<size_t,size_t> &t0 = fa.face2tet_[tri_pair.first];
    const pair<size_t,size_t> &t1 = fa.face2tet_[tri_pair.second];
    assert(t0.first == -1 || t0.second == -1);
    assert(t1.first == -1 || t1.second == -1);
    tet_pair.first = (t0.first ==-1?t0.second:t0.first);
    tet_pair.second = (t1.first == -1?t1.second:t1.first);
    const size_t a0 = find_nearest_axis(frame[tet_pair.first], n1);
    const size_t a1 = find_nearest_axis(frame[tet_pair.second], n2);

    const matrix<double> I = eye<double>(3);
    matrix<double> n00 = get_axis_dir(I, a0);
    matrix<double> n11 = get_axis_dir(I, a1);

    vector<size_t> tet_loop_m = tet_loop;
    if(tet_loop_m.front() == -1) tet_loop_m.pop_back();
    tet_loop_m.erase(tet_loop_m.begin(), tet_loop_m.begin()+1);

    assert(tet_pair.first == tet_loop_m.front() ||
           tet_pair.second == tet_loop_m.front());
    if(tet_pair.second == tet_loop_m.front())
        reverse(tet_loop_m.begin(), tet_loop_m.end());
    assert(tet_pair.first == tet_loop_m.front());

    matrix<double> rot = I;
    for(size_t ti = 0; ti < tet_loop_m.size()-1; ++ti){
        get_best_alignment(&frame[tet_loop_m[ti]][0],
                &frame[tet_loop_m[ti+1]][0], &rot[0]);
        n00 = temp(trans(rot) * n00);
    }
    return calculate_dihedral_angle_degree(n00,n11);
}


int frame_field_check(ptree &pt)
{
    jtf::tet_mesh tm(pt.get<string>("input/tet.value").c_str());
    matrix<double> zyz;
    ASSERT_FUNC(read_zyz(pt.get<string>("input/zyz.value").c_str(), zyz), "invalid_zyz");
    matrix<double> new_zyz;

    hj_frame_alignemt(tm.tetmesh_.mesh_, *tm.fa_, zyz, new_zyz);

    {
      write_zyz("aligned.zyz", new_zyz);
    }
    vector<size_t> surface_edges;
    for(const auto & one_edge : tm.ea_outside_->edges_) {
        surface_edges.push_back(one_edge.first);
        surface_edges.push_back(one_edge.second);
    }
    {
        vector<double> surface_edge_m;
        surface_edge_m.resize(tm.ea_outside_->edge2cell_.size());
        for(size_t ei = 0; ei < tm.ea_outside_->edge2cell_.size(); ++ei){
            const pair<size_t,size_t> & tri_pair = tm.ea_outside_->edge2cell_[ei];
            assert(tm.ea_outside_->is_boundary_edge(tri_pair) == false);
            const double dihedral_angle =
                    calculate_dihedral_angle_degree(tm.outside_face_normal_(colon(), tri_pair.first),
                                                    tm.outside_face_normal_(colon(), tri_pair.second));
            surface_edge_m[ei] = dihedral_angle;
        }
        ofstream ofs("dihedral_edge_orig.vtk");
        line2vtk(ofs, &tm.tetmesh_.node_[0], tm.tetmesh_.node_.size(2), &surface_edges[0], surface_edges.size()/2);
        cell_data(ofs, &surface_edge_m[0], surface_edge_m.size(), "dihedral_angle");
    }

    {
        matrix<matrix<double> > new_frame;
        zyz2frame(new_zyz, new_frame);

        vector<double> surface_edge_frame_m(tm.ea_outside_->edge2cell_.size());
        for(size_t ei = 0; ei < tm.ea_outside_->edge2cell_.size(); ++ei){
            const pair<size_t,size_t> & edge = tm.ea_outside_->edges_[ei];
            if(edge.first == 45 && edge.second == 263)
              cerr << endl;
            if(edge.first == 263 && edge.second == 45)
              cerr << endl;
            const pair<size_t,size_t> & tri_pair = tm.ea_outside_->edge2cell_[ei];
            assert(tm.ea_outside_->is_boundary_edge(tri_pair) == false);
            auto edge_it = tm.ortae_.e2t_.find(tm.ea_outside_->edges_[ei]);
            assert(edge_it != tm.ortae_.e2t_.end());
            const double dihedral_angle_param =
                    calculate_dihedral_angle_degree_param(
                        make_pair(tm.outside_face_idx_[tri_pair.first],
                        tm.outside_face_idx_[tri_pair.second]), *tm.fa_,
                    new_frame, tm.outside_face_normal_(colon(), tri_pair.first),
                    tm.outside_face_normal_(colon(), tri_pair.second), edge_it->second);
            surface_edge_frame_m[ei] = dihedral_angle_param;
        }
        ofstream ofs("dihedral_edge_param.vtk");
        line2vtk(ofs, &tm.tetmesh_.node_[0], tm.tetmesh_.node_.size(2), &surface_edges[0], surface_edges.size()/2);
        cell_data(ofs, &surface_edge_frame_m[0], surface_edge_frame_m.size(), "dihedral_edge_param");
    }
    return 0;
}
