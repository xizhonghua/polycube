#ifndef OPTIMAL_SH_GENERATOR_H
#define OPTIMAL_SH_GENERATOR_H

#include <zjucad/ptree/ptree.h>
#include <hjlib/sparse/sparse.h>
#include <fstream>
#include <jtflib/mesh/trimesh.h>
#include "../tetmesh/tetmesh.h"
#include <complex>
#include "angle_defect.h"
#include "../hex_param/global_alignment.h"
class optimal_sh_generator
{
public:
  optimal_sh_generator(jtf::tet_mesh &tm): tm_(tm){}
  virtual ~optimal_sh_generator(){}
  virtual void opt(zjucad::matrix::matrix<double> &field,
                   boost::property_tree::ptree &pt);

  virtual void project2zyz(const zjucad::matrix::matrix<double> &sh,
                           zjucad::matrix::matrix<double> &zyz)const;

  int write_zyz(const char *file, const zjucad::matrix::matrix<double> &field,
                bool global_align = false)const{
    std::ofstream ofs(file, std::ofstream::binary);
    if(ofs.fail()) return __LINE__;
    if(field.size(1) == 3)
      jtf::mesh::write_matrix(ofs, field);
    if(field.size(1) == 9){
        zjucad::matrix::matrix<double> zyz;
        project2zyz(field, zyz);
        if(global_align){
            zjucad::matrix::matrix<double> zyz_bkp = zyz;
            hj_frame_alignemt(tm_.tetmesh_.mesh_, *tm_.fa_, zyz_bkp, zyz);
          }
        jtf::mesh::write_matrix(ofs, zyz);
      }
    return 0;
  }

  int write_zyz(const zjucad::matrix::matrix<double> &field,
                zjucad::matrix::matrix<double> &zyz,
                bool global_align = false)const{
    assert(field.size(1) == 9);
    project2zyz(field, zyz);
    if(global_align){
        zjucad::matrix::matrix<double> zyz_bkp = zyz;
        hj_frame_alignemt(tm_.tetmesh_.mesh_, *tm_.fa_, zyz_bkp, zyz);
      }
    return 0;
  }

  int zyz2frame(const zjucad::matrix::matrix<double> & zyz,
                zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & frame);

  int write_sh(const char *file, zjucad::matrix::matrix<double> &field ) const
  {
    std::ofstream ofs(file, std::ofstream::binary);
    if(ofs.fail()) return __LINE__;
    if(field.size(1) != 9) return __LINE__;
    jtf::mesh::write_matrix(ofs, field);
    return __LINE__;
  }

private:
  void opt3(zjucad::matrix::matrix<double> &field,
            boost::property_tree::ptree &pt);

  void opt5(zjucad::matrix::matrix<double> &field,
            boost::property_tree::ptree &pt);

  void opt6(zjucad::matrix::matrix<double> &field,
            boost::property_tree::ptree &pt);

  void opt7(zjucad::matrix::matrix<double> &field,
            boost::property_tree::ptree &pt);

  void opt8(zjucad::matrix::matrix<double> &field,
            boost::property_tree::ptree &pt);

  void opt9(zjucad::matrix::matrix<double> &field,
            boost::property_tree::ptree &pt);

  void opt10(zjucad::matrix::matrix<double> &field,
             boost::property_tree::ptree &pt);

  void opt11(zjucad::matrix::matrix<double> &field,
             boost::property_tree::ptree &pt);

  void opt12(zjucad::matrix::matrix<double> &field,
             boost::property_tree::ptree &pt);

  void opt13(zjucad::matrix::matrix<double> &field,
             boost::property_tree::ptree &pt);

  void extract_singularity_from_connection(
      const std::map<std::pair<size_t,size_t>,
      zjucad::matrix::matrix<double> > & connection_matrix);

protected:
  optimal_sh_generator(){}
  optimal_sh_generator & operator() (const optimal_sh_generator &){}

  double get_rij_angle(const std::pair<size_t,size_t> & face_pair,
                       const zjucad::matrix::matrix<size_t> & faces,
                       const zjucad::matrix::matrix<double> & node,
                       const zjucad::matrix::matrix<double> & face_normal);

  double get_rij_angle2(const std::pair<size_t,size_t> & face_pair,
                        const zjucad::matrix::matrix<size_t> &faces,
                        const zjucad::matrix::matrix<double> & node,
                        const zjucad::matrix::matrix<double> & face_normal);

  ///
  /// @brief reorder_edge_according_face_vertex
  /// @param common_edge
  /// @param one_face
  /// @return 1, order is changed, or 0
  ///
  int reorder_edge_according_face_vertex(
      size_t * common_edge,
      const zjucad::matrix::matrix<size_t> &one_face);

  void init(boost::property_tree::ptree &pt);

  void init_connection(boost::property_tree::ptree &pt);

protected:
  void get_jacobian_matrix(hj::sparse::csc<double,int32_t> &A);

  void get_jacobian_matrix2(hj::sparse::csc<double,int32_t> &A,
                            boost::property_tree::ptree &pt);

  void smooth_surface_normal(zjucad::matrix::matrix<double> &normal,
                             const zjucad::matrix::matrix<size_t> &face,
                             const zjucad::matrix::matrix<double> & node,
                             std::shared_ptr<geometry_aware_angle_defect> ga,
                             boost::property_tree::ptree &pt,
                             zjucad::matrix::matrix<double> & reference_dir);

  void compute_dis2surface(zjucad::matrix::matrix<double> &sf,
                           boost::property_tree::ptree &pt);

  void compute_disrand(zjucad::matrix::matrix<double> &sh,
                       boost::property_tree::ptree &pt);

  void compute_dis2x(zjucad::matrix::matrix<double> &sf,
                     boost::property_tree::ptree &pt);

  void compute_dis2orig(zjucad::matrix::matrix<double> &sf,
                        boost::property_tree::ptree &pt);

  void compute_dis2param_vol(zjucad::matrix::matrix<double> &sf,
                             boost::property_tree::ptree &pt);

  void compute_dis2param_edge(zjucad::matrix::matrix<double> &sf,
                              boost::property_tree::ptree &pt);

  void compute_size_field(zjucad::matrix::matrix<double> &sf,
                          boost::property_tree::ptree &pt);

  void compute_inner_connection(zjucad::matrix::matrix<double> &sf,
                                std::map<std::pair<size_t,size_t>,
                                zjucad::matrix::matrix<double> > &inner_face_connection);

  void compute_inner_connection_jtf(zjucad::matrix::matrix<double> &sf,
                                    std::map<std::pair<size_t,size_t>, double> &inner_face_connection);

  void compute_inner_connection_jtf2(
      boost::property_tree::ptree &pt,
      std::map<std::pair<size_t,size_t>, zjucad::matrix::matrix<double> > &inner_face_connection);

  void compute_inner_connection_jtf3(
      boost::property_tree::ptree &pt,
      zjucad::matrix::matrix<double> & point_size,
      std::map<std::pair<size_t,size_t>, zjucad::matrix::matrix<double> > & connection_matrix);

  void load_edge_dense_field(std::map<std::pair<size_t,size_t>, double> & edge_dense,
                             boost::property_tree::ptree &pt)const;

  void load_guiding_field(zjucad::matrix::matrix<double> &gv,
                          boost::property_tree::ptree &pt);

  void load_surface_charts(const char * file,
                           std::map<std::pair<size_t,bool>, std::set<size_t> > & charts);

private:
  jtf::tet_mesh tm_;
  zjucad::matrix::matrix<double> rotation_angle_; // for connection
  zjucad::matrix::matrix<double> rij_;
  zjucad::matrix::matrix<double> edge_area_weight_;
  zjucad::matrix::matrix<double> local2global_rot_;
  zjucad::matrix::matrix<double> reference_dir_;
  zjucad::matrix::matrix<double> surface_target_dir_angle_;
  zjucad::matrix::matrix<size_t> surface_target_idx_, surface_target_idx_on_boundary_;
  std::map<std::pair<size_t,size_t>,double> edge_dense_field_;

  std::map<std::pair<size_t,size_t>,double> inner_face_connection_;
  std::map<std::pair<size_t,size_t>, zjucad::matrix::matrix<double> > connection_matrix_;
  std::vector<std::vector<size_t> > feature_lines_;
};

class optimal_sh_generator_seq : public optimal_sh_generator
{
public:
  optimal_sh_generator_seq(std::vector<jtf::tet_mesh> &tm_seq,
                           std::vector<double> &time_step)
    :tm_seq_(tm_seq), time_step_(time_step), optimal_sh_generator(tm_seq.front()){
    tm_ = tm_seq_.front();
  }
  virtual ~optimal_sh_generator_seq(){}
  virtual void opt(zjucad::matrix::matrix<double> &field,
                   boost::property_tree::ptree &pt){
    opt0(field, pt);
  }

  void opt0(zjucad::matrix::matrix<double> &field,
            boost::property_tree::ptree &pt);

  void deformation_project_sh(
      const jtf::tet_mesh &tm_current,
      const jtf::tet_mesh & tm_target,
      zjucad::matrix::matrix<double> &sh);
private:
  std::vector<jtf::tet_mesh> tm_seq_;
  std::vector<double> time_step_;
  jtf::tet_mesh tm_;

private:
  zjucad::matrix::matrix<double> rotation_angle_; // for connection
  zjucad::matrix::matrix<double> rij_;
  zjucad::matrix::matrix<double> edge_area_weight_;
  zjucad::matrix::matrix<double> local2global_rot_;
  zjucad::matrix::matrix<double> reference_dir_;
};

class optimal_surface_sh_generator : public optimal_sh_generator
{
public:
  optimal_surface_sh_generator(const jtf::mesh::tri_mesh &tm):tm_(tm){}
  virtual ~optimal_surface_sh_generator(){}

  virtual void opt(zjucad::matrix::matrix<double>& sh,
                   boost::property_tree::ptree &pt);

  void opt_cos_sin(zjucad::matrix::matrix<double> &tradition,
                   boost::property_tree::ptree &pt);

  void projection_z_axis_rotaiton(zjucad::matrix::matrix<double> &sh);

  void projection_z_axis_rotaiton_each_cell(
      zjucad::matrix::matrix<double> &sh);
private:
  optimal_surface_sh_generator(){}
  optimal_surface_sh_generator & operator() (const optimal_surface_sh_generator &){}

private:
  void init(boost::property_tree::ptree &pt);

  //! @brief init_connection, and flatten neigubour triangles.
  void init_connection(boost::property_tree::ptree &pt);

  void get_jacobian_matrix(hj::sparse::csc<double,int32_t>&A) const;

  void get_tradition_jacobian_matrix(hj::sparse::csc<double, int32_t> &A) const;

  //! @brief project the sh coefficients to the space where
  void projection(zjucad::matrix::matrix<double> &sh) const;

  double get_edge_adj_tri_area(
      const std::pair<size_t,size_t> & one_edge,
      const zjucad::matrix::matrix<double> & node,
      const zjucad::matrix::matrix<size_t> & face);

private:
  const jtf::mesh::tri_mesh tm_;
  zjucad::matrix::matrix<double> rotation_angle_; // for connection
  zjucad::matrix::matrix<double> rij_;
  zjucad::matrix::matrix<double> edge_area_weight_;
};

class optimal_sh_generator_polycube : public optimal_sh_generator
{
public:
  optimal_sh_generator_polycube(jtf::tet_mesh & orig_tet, jtf::tet_mesh & polycube_tet)
    :optimal_sh_generator(orig_tet), orig_tm_(orig_tet), polycube_tm_(polycube_tet){}
  virtual ~optimal_sh_generator_polycube(){}
  virtual void opt(zjucad::matrix::matrix<double> &field, boost::property_tree::ptree &pt);
private:
  void opt0(zjucad::matrix::matrix<double> & field, boost::property_tree::ptree &pt);
  void opt1(zjucad::matrix::matrix<double> & field, boost::property_tree::ptree &pt);
  void opt2(zjucad::matrix::matrix<double> & field, boost::property_tree::ptree &pt);
  void init(boost::property_tree::ptree &pt);
  void update_normal(const std::set<size_t> & patch_to_use_real_noraml,
                     const std::vector<std::vector<size_t> > & all_patches);
  bool valid_field(const zjucad::matrix::matrix<double> & field);
  bool rule0(const std::deque<std::pair<size_t,size_t> > & one_chain,
             const std::set<size_t> &surface_points,
             const zjucad::matrix::matrix<double> &point_normal);

  void search_nearest_k_chain(
      const std::deque<std::pair<size_t,size_t> > &one_chain,
      const size_t k_nearest, std::set<size_t> &nearest_chains);

  double dis_of_two_chains(
      const std::deque<std::pair<size_t,size_t> > & chain0,
      const std::deque<std::pair<size_t,size_t> > & chain1,
      const zjucad::matrix::matrix<double> & node);

  void cal_weighted_average_normal(const std::vector<size_t> & one_patch,
                                   const zjucad::matrix::matrix<double> & normal,
                                   const zjucad::matrix::matrix<double> & area,
                                   zjucad::matrix::matrix<double> & avg_normal);

  ///
  /// \brief adjust_surface_normal
  /// \param surface_normal
  /// \param frame
  /// \param ps
  /// \return 0: adjusted, 1: do not need adjusted, it's ok, -1: not able to find a solution
  ///
  int adjust_surface_normal( zjucad::matrix::matrix<double> & surface_normal,
                             const zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & frame,
                             const jtf::mesh::patch_separater &ps,
                             const zjucad::matrix::matrix<double> &polycube_face_normal);

  ///
  /// \brief load_picked_patches
  /// \param picked_patch_idx
  /// \param pt
  ///
  int load_picked_patches(std::vector<size_t> & picked_patch_idx, boost::property_tree::ptree &pt);
  ///
  /// \brief get_point_convex_info
  /// \param pi
  /// \param orfap
  /// \param ps
  /// \param normal
  /// \return 1: convex, 0: not convex, -1: invalid
  ///
  size_t get_point_convex_info(
      const size_t pi, const jtf::mesh::one_ring_face_at_point &orfap,
      const jtf::mesh::patch_separater & ps, const zjucad::matrix::matrix<double> & normal);
private:
  jtf::tet_mesh orig_tm_, polycube_tm_;
  std::vector<std::deque<std::pair<std::size_t,std::size_t> > > chain_;
  std::vector<bool> is_sharp_chain_;
};

#endif // OPTIMAL_SH_GENERATOR_H
