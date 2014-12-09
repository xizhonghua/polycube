#include <memory>
#include <jtflib/mesh/util.h>
#include "hex_frame_opt.h"
#include "../tetmesh/util.h"
#include "frame_function.h"

static int find_tet_adjacent_tet_idx(const matrixst &tet, const jtf::mesh::face2tet_adjacent &fa, vector<vector<size_t> > &tet_adjacent_idx)
{
    size_t tet_num = tet.size(2);
    tet_adjacent_idx.resize(tet_num);
    size_t valid_tet_num = 0;
    for(size_t i = 0; i < tet_num; ++i)
    {
        for(int j = 0; j < 4; ++j)
        {
            pair<size_t,size_t> nb_tet_idx = fa.query(tet(j,i), tet((j+1)%4,i), tet((j+2)%4,i));
            if(nb_tet_idx.first == -1 || nb_tet_idx.second == -1)
                continue;
            tet_adjacent_idx[i].push_back(nb_tet_idx.first == i?nb_tet_idx.second:nb_tet_idx.first);
        }
        if(tet_adjacent_idx[i].size() == 4)
            ++valid_tet_num;
    }
    return valid_tet_num;
}

static int calculate_element_tet_volume(const vector<vector<size_t> > &tet_adjacent_idx,
                                        const matrixst &tet,
                                        const matrixd &node,
                                        vector<double> &volume_of_tets)
{
    const int tet_num = tet.size(2);
    int valid_tet_num = 0;
    matrixd inner_position(3,tet_num);

    for(int i = 0; i < tet_num; ++i)
    {
        if(tet_adjacent_idx[i].size() == 4)
            ++valid_tet_num;
        inner_position(colon(),i)=(node(colon(),tet(0,i))
                                   +node(colon(),tet(1,i))
                                   +node(colon(),tet(2,i))
                                   +node(colon(),tet(3,i))
                                   )/4.0;
    }

    volume_of_tets.resize(valid_tet_num);

    for(int i = 0,idx = 0; i < tet_num; ++i)
    {
        const vector<size_t> & adjacent_idx = tet_adjacent_idx[i];
        if(adjacent_idx.size() == 4)
        {
            matrixd a = inner_position(colon(),adjacent_idx[0]) - inner_position(colon(),adjacent_idx[3]);
            matrixd b = inner_position(colon(),adjacent_idx[1]) - inner_position(colon(),adjacent_idx[3]);
            matrixd c = inner_position(colon(),adjacent_idx[2]) - inner_position(colon(),adjacent_idx[3]);
            double volume = (dot(a,cross(b,c)))/6.0;
            //cerr << "# " << i << " vol = " << volume << endl;
            volume_of_tets[idx] = fabs(volume);
            ++idx;
        }
    }
}

// if the face is on surface, return the tet_volume/4;
// if the face is inner surface, return the (tet_vi+tet_vj)/4

static double calculate_face_volume_weight(const size_t face_idx,const jtf::mesh::face2tet_adjacent& fa, const vector<double> & volume_tet)
{
    const size_t &vi = fa.face2tet_[face_idx].first;
    const size_t &vj = fa.face2tet_[face_idx].second;
    double volume = 0.0;
    if( vi != -1)
        volume += volume_tet[vi]/4.0;
    if( vj != -1)
        volume += volume_tet[vj] / 4.0;
    return volume;
}

int non_sym_frame_opt::setup_equations_new(
        const matrixst &tet,
        const matrixd &node,
        const matrixd &fixed_frame,
        const zjucad::matrix::matrix<zjucad::matrix::idx_type> &fixed_frame_idx,
        const boost::unordered_map<size_t,size_t> &surface_type,
        const matrixst &aligned_surface, // face
        const zjucad::matrix::matrix<zjucad::matrix::idx_type> &surf_aligned_idx,
        const double weight[2], //! [align, fix], default 1e6
        const matrixd &stiff, // tet weighting
        const boost::unordered_map<std::pair<size_t,size_t>,size_t> &jump_face_type,
        const double LP_surface,
        const double LP_smooth)
{
    unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
    const size_t nn  = tet.size(2), // elements
            inner_smooth_en = jump_face_type.size(),
            //fn = fixed_frame.size(2),
            surf_an = surf_aligned_idx.size();

    funcs_.reserve(inner_smooth_en + surf_an);

    std::vector<double> volume_of_tets(tet.size(2)),area_of_surf_tris(surf_aligned_idx.size());
    for(size_t ti = 0; ti < volume_of_tets.size(); ++ti)
        volume_of_tets[ti] = fabs(jtf::mesh::cal_tet_vol(node(colon(), tet(colon(), ti))));
    for(size_t ti = 0; ti < area_of_surf_tris.size(); ++ti)
        area_of_surf_tris[ti] = fabs(jtf::mesh::cal_face_area(&fa->faces_[surf_aligned_idx[ti]][0],fa->faces_[surf_aligned_idx[ti]].size(),node));
    const double area_of_surface = accumulate(area_of_surf_tris.begin(), area_of_surf_tris.end(), 0.0);
    const double total_volume = accumulate(volume_of_tets.begin(), volume_of_tets.end(), 0.0);

    // inner smoothing, including jump face
    typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mcit;

    cut_jump_smooth_function_.reserve(inner_smooth_en);
    for(mcit it = jump_face_type.begin(); it != jump_face_type.end(); ++it){
        matrixd bary_i = zeros<double>(3,1);
        matrixd bary_j = zeros<double>(3,1);

        const pair<size_t,size_t> &two_tets = it->first;
        for(size_t v = 0; v < 4; ++v) {
            bary_i += node(colon(),tet(v,two_tets.first) );
            bary_j += node(colon(),tet(v,two_tets.second));
        }
        bary_i /= 4;
        bary_j /= 4;
        const double dij = norm(bary_i - bary_j);
        if(dij > 1e-8)
        {
            const double weight_ = sqrt((volume_of_tets[two_tets.first] + volume_of_tets[two_tets.second] )/4) * pow(total_volume,-1.0/6.0) /dij;
            cut_jump_smooth_function_.push_back(
                        std::shared_ptr<cut_jump_smooth_function>(
                            new cut_jump_smooth_function(it->second,two_tets.first,two_tets.second,tet.size(2),LP_smooth)));
            funcs_.push_back(std::shared_ptr<function_t<double,int32_t> >(
                                 new scalar_op<std::multiplies,double,int32_t, nil_ptr>(*cut_jump_smooth_function_.back() ,
                                                                                        weight_)));
        }
    }


#if 1
    //original surface aligned
    cut_surface_align_function_.reserve(surf_an);

    for(size_t i = 0; i < surf_an && weight[0] > 0; ++i)
    {
        const size_t &vi = fa->face2tet_[surf_aligned_idx[i]].first;
        const size_t &vj = fa->face2tet_[surf_aligned_idx[i]].second;

        assert(vi == -1 || vj == -1);
        const size_t tet_idx = (vi == -1)?vj:vi;

        matrixd normal = zeros<double>(3,1);
        if(jtf::tetmesh::calculate_face_normal(tet,node,tet_idx,aligned_surface(colon(),i),normal)) continue; // degenerate face

        boost::unordered_map<size_t,size_t>::const_iterator it = surface_type.find(surf_aligned_idx[i]);
        if(it == surface_type.end()){
            cerr << "# strange, can not find this surface." << endl;
            continue;
        }else{
            const size_t surface_type_ = it->second;
            cut_surface_align_function_.push_back(
                        std::shared_ptr<cut_surface_align_function>(
                            new cut_surface_align_function(tet_idx, surface_type_,tet.size(2),normal,LP_surface)));
            funcs_.push_back(std::shared_ptr<function_t<double,int32_t> >(
                                 new scalar_op<std::multiplies,double,int32_t, nil_ptr>(*cut_surface_align_function_.back(),
                                                                                        sqrt(area_of_surf_tris[i]/area_of_surface*weight[0]) )));
        }
    }
#endif

    cerr << "# surface align functions num " << cut_surface_align_function_.size() << endl;
    cerr << "# inner smooth functions num " << cut_jump_smooth_function_.size() << endl;
    cerr << "# number of functions:" << funcs_.size() << endl;
    all_f_.reset(new_catenated_function_omp<double, int32_t>(funcs_));
    cerr << "# x: " << all_f_->dim_of_x()
         << ", #f: " << all_f_->dim_of_f()
         << ", #jac: " << all_f_->jac_nnz()
         << endl;
    return 0;
}



int non_sym_frame_opt::setup_equations_linear_new(
        const matrixst &tet,
        const matrixd &node,
        const matrixd &fixed_frame,
        const zjucad::matrix::matrix<zjucad::matrix::idx_type> &fixed_frame_idx,
        const boost::unordered_map<size_t,size_t> &surface_type,
        const matrixst &aligned_surface, // face
        const zjucad::matrix::matrix<zjucad::matrix::idx_type> &surf_aligned_idx,
        const double weight[2], //! [align, fix], default 1e6
        const matrixd &stiff, // tet weighting
        const boost::unordered_map<std::pair<size_t,size_t>,size_t> &jump_face_type,
        const double LP_surface,
        const double LP_smooth,
        const double LP_rtr,
        const double RTR_w
        )
{
#if 1
    unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
    const size_t nn  = tet.size(2), // elements
            inner_smooth_en = jump_face_type.size(),
            surf_an = surf_aligned_idx.size();

    funcs_.reserve(inner_smooth_en + surf_an);

    std::vector<double> volume_of_tets(tet.size(2)),area_of_surf_tris(surf_aligned_idx.size());
    for(size_t ti = 0; ti < volume_of_tets.size(); ++ti)
        volume_of_tets[ti] = fabs(jtf::mesh::cal_tet_vol(node(colon(), tet(colon(), ti))));
    for(size_t ti = 0; ti < area_of_surf_tris.size(); ++ti)
        area_of_surf_tris[ti] = fabs(jtf::mesh::cal_face_area(&fa->faces_[surf_aligned_idx[ti]][0],
                                     fa->faces_[surf_aligned_idx[ti]].size(), node ));
    const double area_of_surface = accumulate(area_of_surf_tris.begin(), area_of_surf_tris.end(), 0.0);
    const double total_volume = accumulate(volume_of_tets.begin(), volume_of_tets.end(), 0.0);

    // inner smoothing, including jump face
    typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mcit;
#endif
#if 1
    cut_jump_smooth_linear_function_new_.reserve(inner_smooth_en);
    matrixd bary_i = zeros<double>(3,1);
    matrixd bary_j = zeros<double>(3,1);

    for(mcit it = jump_face_type.begin(); it != jump_face_type.end(); ++it){
        bary_i = zeros<double>(3,1);
        bary_j = zeros<double>(3,1);

        const pair<size_t,size_t> &two_tets = it->first;
        for(size_t v = 0; v < 4; ++v) {
            bary_i += node(colon(),tet(v,two_tets.first) );
            bary_j += node(colon(),tet(v,two_tets.second) );
        }
        bary_i /= 4;
        bary_j /= 4;
        const double dij = norm(bary_i - bary_j);
        if(dij > 1e-8)
        {
            const double weight_ = sqrt((volume_of_tets[two_tets.first] + volume_of_tets[two_tets.second] )/4) * pow(total_volume,-1.0/6.0) /dij;
            cut_jump_smooth_linear_function_new_.push_back(
                        std::shared_ptr<cut_jump_smooth_linear_function_new>(
                            new cut_jump_smooth_linear_function_new(it->second,two_tets.first,two_tets.second,tet.size(2),LP_smooth)));
            funcs_.push_back(std::shared_ptr<function_t<double,int32_t> >(
                                 new scalar_op<std::multiplies,double,int32_t, nil_ptr>(*cut_jump_smooth_linear_function_new_.back() ,
                                                                                        weight_)));
        }
    }
#endif

#if 1
    //original surface aligned
    cut_surface_align_linear_function_new_.reserve(surf_an);

    for(size_t i = 0; i < surf_an && weight[0] > 0; ++i)
    {
        const size_t &vi = fa->face2tet_[surf_aligned_idx[i]].first;
        const size_t &vj = fa->face2tet_[surf_aligned_idx[i]].second;

        assert(vi == -1 || vj == -1);
        const size_t tet_idx = (vi == -1)?vj:vi;

        matrixd normal = zeros<double>(3,1);
        if(jtf::tetmesh::calculate_face_normal(tet,node,tet_idx,aligned_surface(colon(),i),normal)) continue; // degenerate face

        boost::unordered_map<size_t,size_t>::const_iterator it = surface_type.find(surf_aligned_idx[i]);
        if(it == surface_type.end()){
            cerr << "# strange, can not find this surface." << endl;
            continue;
        }else{
            const size_t surface_type_ = it->second;
            cut_surface_align_linear_function_new_.push_back(
                        std::shared_ptr<cut_surface_align_linear_function_new>(
                            new cut_surface_align_linear_function_new(tet_idx, surface_type_,tet.size(2),normal,LP_surface)));
            funcs_.push_back(std::shared_ptr<function_t<double,int32_t> >(
                                 new scalar_op<std::multiplies,double,int32_t, nil_ptr>(*cut_surface_align_linear_function_new_.back(),
                                                                                        sqrt(area_of_surf_tris[i]/area_of_surface*weight[0]) )));
        }
    }
#endif

#if 1
    cut_RTR_linear_function_.reserve(tet.size(2));
    for(size_t t = 0; t < tet.size(2) && fabs(RTR_w) > 1e-8;++t){
        cut_RTR_linear_function_.push_back(
                    std::shared_ptr<cut_RTR_linear_function>(
                        new cut_RTR_linear_function(tet.size(2),t,LP_rtr)));
        funcs_.push_back(std::shared_ptr<function_t<double,int32_t> >(
                             new scalar_op<std::multiplies,double,int32_t, nil_ptr>(*cut_RTR_linear_function_.back() ,
                                                                                    RTR_w)));
        break;
    }
#endif
#if 1
    cerr << "# surface align functions num " << cut_surface_align_linear_function_new_.size() << endl;
    cerr << "# inner smooth functions num " << cut_jump_smooth_linear_function_new_.size() << endl;
    cerr << "# rtr linear functions num " << cut_RTR_linear_function_.size() << endl;
    cerr << "# number of functions:" << funcs_.size() << endl;
    all_f_.reset(new_catenated_function_omp<double, int32_t>(funcs_));
    cerr << "# x: " << all_f_->dim_of_x()
         << ", #f: " << all_f_->dim_of_f()
         << ", #jac: " << all_f_->jac_nnz()
         << endl;
#endif
    return 0;
}
