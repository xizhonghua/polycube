#ifndef JTF_angle_DEFECT_H
#define JTF_angle_DEFECT_H

#include <jtflib/mesh/mesh.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/ptree/ptree.h>

class angle_defect
{
public:
  angle_defect(){}
  virtual ~angle_defect(){}

  virtual void opt(const jtf::mesh::edge2cell_adjacent &ea,
                   const zjucad::matrix::matrix<size_t> &face,
                   const zjucad::matrix::matrix<double> &node,
                   zjucad::matrix::matrix<double> &angle,
                   boost::property_tree::ptree &pt){}
};

class naive_angle_defect : public angle_defect
{
public:
  naive_angle_defect(){}
  virtual ~naive_angle_defect(){}

  virtual void opt(const jtf::mesh::edge2cell_adjacent &ea,
                   const zjucad::matrix::matrix<size_t> &face,
                   const zjucad::matrix::matrix<double> &node,
                   zjucad::matrix::matrix<double> &angle,
                   boost::property_tree::ptree &pt);

  void cal_angle_defect(
      const jtf::mesh::one_ring_face_at_point & orfap,
      const zjucad::matrix::matrix<size_t> & face,
      const zjucad::matrix::matrix<double> & node,
      std::map<size_t,double> & angle_defect_map);
};

class manually_set_angle_defect : public angle_defect
{
public:
  manually_set_angle_defect(){}
  virtual ~manually_set_angle_defect(){}

  virtual void opt(const jtf::mesh::edge2cell_adjacent &ea,
                   const zjucad::matrix::matrix<size_t> &face,
                   const zjucad::matrix::matrix<double> &node,
                   zjucad::matrix::matrix<double> &angle,
                   boost::property_tree::ptree &pt);

  double calc_next_angle_along_circle(
      const size_t fi, const size_t fj,
      const zjucad::matrix::matrix<size_t> & face,
      const zjucad::matrix::matrix<double> & normal,
      const zjucad::matrix::matrix<double> & node,
      const double init_angle); // init_angle is based on the first edge of fi

  double calc_angle_from_d0_to_d1(
      const zjucad::matrix::matrix<double> & d0,
      const zjucad::matrix::matrix<double> & d1,
      const zjucad::matrix::matrix<double> & normal);

  void orient_edge( size_t * common_edge,
                    const zjucad::matrix::matrix<size_t> & one_face);

  //! @param one reference in tri which has the same topology order according to one_edge
  void find_safe_reference(
      const zjucad::matrix::matrix<size_t> & one_face,
      const zjucad::matrix::matrix<double> & node,
      const size_t *one_edge,
      zjucad::matrix::matrix<double> & ri);

protected:
  int load_angle_defect_set_file(const char * file,
                                 std::map<size_t,double> & angle_defect_map);

  int cal_angle_defect2(const jtf::mesh::one_ring_face_at_point & orfap,
                        const zjucad::matrix::matrix<size_t> & face,
                        const zjucad::matrix::matrix<double> & node,
                        const zjucad::matrix::matrix<double> & norma,
                        std::map<size_t,double> & angle_defect);

  double cal_cota_angle(
      const size_t face_idx,
      const zjucad::matrix::matrix<size_t> &face,
      const std::pair<size_t,size_t> & one_edge,
      const zjucad::matrix::matrix<double> &node);

  //! @brief this is calculated according to trivial connection paper, but fails
  //         when the mesh has flip elements
  double cal_angle_defect_along_circle(
      const std::vector<size_t> & face_loop,
      const zjucad::matrix::matrix<double> & node,
      const zjucad::matrix::matrix<size_t> & face,
      const zjucad::matrix::matrix<double> & normal);

  double cal_angle_defect_directly(
      const std::vector<size_t> & face_loop,
      const zjucad::matrix::matrix<double> & node,
      const zjucad::matrix::matrix<size_t> & face,
      const size_t this_point);

};


class even_angle_defect : public manually_set_angle_defect
{
public:
  even_angle_defect(){}
  virtual ~even_angle_defect(){}

  virtual void opt(const jtf::mesh::edge2cell_adjacent &ea,
                   const zjucad::matrix::matrix<size_t> &face,
                   const zjucad::matrix::matrix<double> &node,
                   zjucad::matrix::matrix<double> &angle,
                   boost::property_tree::ptree &pt);
};

class geometry_aware_angle_defect : public manually_set_angle_defect
{
public:
  geometry_aware_angle_defect(){}
  virtual ~geometry_aware_angle_defect(){}

  virtual void opt(const jtf::mesh::edge2cell_adjacent &ea,
                   const zjucad::matrix::matrix<size_t> &face,
                   const zjucad::matrix::matrix<double> &node,
                   zjucad::matrix::matrix<double> &angle,
                   boost::property_tree::ptree &pt);

  void cal_geodesic_distance(
      std::vector<std::vector<std::pair<size_t,double> > > & distance,
      const zjucad::matrix::matrix<double> &node,
      std::shared_ptr<const jtf::mesh::one_ring_point_at_point> orpap,
      const double max_dis,
      std::set<size_t> * feature_points = nullptr);

private:
  void filter_point_curvature(
      const double max_dis,
      const zjucad::matrix::matrix<size_t> &face,
      const zjucad::matrix::matrix<double> &node,
      const std::map<size_t,double> &k2,
      std::map<size_t,double> &kcorr,
      std::set<size_t> *feature_points = nullptr);
};
#endif // angle_DEFECT_H
