#ifndef CJ_EXTRACTOR_H
#define CJ_EXTRACTOR_H

#include <map>
#include <iostream>
#include <cstdio>
#include <vector>
#include <utility>
#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/util/vertex_connection.h>

namespace extractor {
  
#define EPSILON 1e-8
#define PI 3.141592653589793
//#define FEATURE_ANGLE PI/4


  typedef struct {
    double xmin, xmax, ymin, ymax, zmin, zmax;
  }bounding_box;

  typedef struct point_info {
    size_t coordinate_idx;
    size_t tet_idx;
    zjucad::matrix::matrix<double> b_coeff;
    
    point_info(size_t arg1, size_t arg2, const zjucad::matrix::matrix<double>& arg3) {
      this->coordinate_idx = arg1;
      this->tet_idx = arg2;
      this->b_coeff = arg3;
    }
  }point_info;

  //typedef std::pair<boost::tuple<long long, long long, long long>, size_t> STATE;
  typedef size_t STATE;

  typedef std::pair<std::multimap<boost::tuple<double, double, double>, point_info>::iterator,
                    std::multimap<boost::tuple<double, double, double>, point_info>::iterator> RANGE;

  const double rot [24][3][3]={
     {{1,0,0},{0,0,1},{0,-1,0}},  //  u1
     {{1,0,0},{0,-1,0},{0,0,-1}}, //  u2
     {{1,0,0},{0,0,-1},{0,1,0}},  //  u3
     {{0,0,-1},{0,1,0},{1,0,0}},  //  v1
     {{-1,0,0},{0,1,0},{0,0,-1}}, //  v2
     {{0,0,1},{0,1,0},{-1,0,0}},  //  v3
     {{0,1,0},{-1,0,0},{0,0,1}},  //  w1
     {{-1,0,0},{0,-1,0},{0,0,1}}, //  w2
     {{0,-1,0},{1,0,0},{0,0,1}},  //  w3
     {{1,0,0},{0,1,0},{0,0,1}},   //  I
     {{-1,0,0},{0,0,1},{0,1,0}},   //  u1 * v2, u3 * w2, w2 * u1
     {{-1,0,0},{0,0,-1},{0,-1,0}}, //  v2 * u1, u3 * v2
     {{0,0,1},{1,0,0},{0,1,0}},    //  v3 * u3, w3 * v3,  u3*w3
     {{0,0,-1},{1,0,0},{0,-1,0}},  //  v1 * u1, w3 * v1,
     {{0,0,-1},{-1,0,0},{0,1,0}},  //  v1 * u3, u3 * w1,
     {{0,0,1},{-1,0,0},{0,-1,0}},  //  v3 * u1, u1 * w1
     {{0,1,0},{0,0,1},{1,0,0}},    //  u1 * v1, v1 * w1, w1 * u1
     {{0,-1,0},{0,0,-1},{1,0,0}},  //  u3 * v1, w3 * u3
     {{0,0,1},{0,-1,0},{1,0,0}},   //  u2 * v1, v3 * u2
     {{0,1,0},{0,0,-1},{-1,0,0}},  //  u3 * v3, v3 * w1
     {{0,-1,0},{0,0,1},{-1,0,0}},  //  u1 * v3, w3 * u1, v3 * w3,
     {{0,0,-1},{0,-1,0},{-1,0,0}}, //  v1 * u2, w2 * v1, v3 * w2
     {{0,1,0},{1,0,0},{0,0,-1}},   //  u2 * w3, v2 * w1, w3 * v2
     {{0,-1,0},{-1,0,0},{0,0,-1}}  //  u2 * w1, w3 * u2
  };

  class extractor
  {
  public:
    virtual ~extractor(){}
    ///
    /// @brief extract      to extract quad/hex from a tri/tet mesh
    /// @param input_mesh   input mesh, either will be tri/tet
    /// @param param_node   mesh node in parameterizaton space
    /// @param orig_node    mesh node in original space
    /// @param inner_type   inner transition type defined on each two cells
    /// @param output_mesh  output mesh quad/hex
    /// @param output_node  node of output mesh
    /// @return             0 if works fine, or non-zeros
    ///
    virtual int extract(const zjucad::matrix::matrix<size_t> & input_mesh,
                        const zjucad::matrix::matrix<double> &param_node,
                        const zjucad::matrix::matrix<double> & orig_node,
                        const std::map<std::pair<size_t,size_t>,size_t> & inner_type,
                        zjucad::matrix::matrix<size_t> & output_mesh,
                        zjucad::matrix::matrix<double> & output_node) = 0;
  public:
    ///
    /// @brief extract_integer_points extract integer points, each should be
    ///                               associated with its coordinates and
    ///                               its tet/tri idx, and its barycentric coordinate
    /// @param param_mesh             input param mesh
    /// @param param_nodes            input param nodes
    /// @return                       0 if works fine, or non-zeros
    ///
    virtual int extract_integer_points(const zjucad::matrix::matrix<size_t> &param_mesh,
                                       const zjucad::matrix::matrix<double> & param_nodes) = 0;

    ///
    /// @brief extract_integer_lines extract integer lines according to integer
    ///                              points and transition on each faces (tet)
    /// @param param_mesh            input mesh
    /// @param param_nodes           input param node
    /// @param inner_type            transition type defined on each two cells (cells should be ordered)
    /// @return                      0 if works fine, or non-zeros
    ///
    virtual int extract_integer_lines(const zjucad::matrix::matrix<size_t> & param_mesh,
                                      const zjucad::matrix::matrix<double> &param_nodes,
                                      const std::map<std::pair<size_t,size_t>, size_t> & inner_type) = 0;

    ///
    /// @brief extract_mesh         extract mesh from integer lines
    /// @param param_mesh           input param mesh
    /// @param param_nodes          input nodes
    /// @return                     0 if works fine, or non-zeros
    ///
    virtual int extract_mesh(const zjucad::matrix::matrix<size_t> & param_mesh,
                             const zjucad::matrix::matrix<double> &param_nodes) = 0;

    ///
    /// @brief mapping_to_orig_mesh  map extracted mesh to original space
    /// @param orig_mesh             input origin mesh
    /// @param orig_node             input origin node
    /// @return
    ///
    virtual int mapping_to_orig_space(const zjucad::matrix::matrix<size_t> & orig_mesh,
                                      const zjucad::matrix::matrix<double> & orig_node) = 0;
  };


  class hex_extractor : public extractor
  {
  public :
      size_t SCALE;
      std::vector<size_t> father;
      boost::unordered_set<size_t> feature_points;

  public:
    virtual ~hex_extractor() {}

    ///
    /// @brief extract      to extract quad/hex from a tri/tet mesh
    /// @param input_mesh   input mesh, either will be tri/tet
    /// @param param_node   mesh node in parameterizaton space
    /// @param orig_node    mesh node in original space
    /// @param output_mesh  output mesh quad/hex
    /// @param output_node  node of output mesh
    /// @return             0 if works fine, or non-zeros

    virtual int extract(const zjucad::matrix::matrix<size_t> & input_mesh,
                        const zjucad::matrix::matrix<double> &param_node,
                        const zjucad::matrix::matrix<double> & orig_node,
                        const std::map<std::pair<size_t,size_t>,size_t> & inner_type,
                        zjucad::matrix::matrix<size_t> & output_mesh,
                        zjucad::matrix::matrix<double> & output_node)
    {/*
      if(extract_integer_points(input_mesh, param_node)) {
          std::cerr << "# [error] extract inerger_points fail." << std::endl;
          return __LINE__;
        }
      if(extract_integer_lines(input_mesh, param_node,inner_type)) {
          std::cerr << "# [error] extract inerger_lines fail." << std::endl;
          return __LINE__;
        }
      if(extract_mesh(input_mesh, param_node)) {
          std::cerr << "# [error] extract mesh fail." << std::endl;
          return __LINE__;
        }
      if(mapping_to_orig_space(input_mesh, param_node)) {
          std::cerr << "# [error] mappint to orig space fail." << std::endl;
          return __LINE__;
        }

      return 0;*/
    }



    ///
    /// @brief mesh_scaling           scale the coordinate of nodes by coefficient ceoff
    ///
    /// @return                       0 if works fine
    int mesh_scaling(zjucad::matrix::matrix<double> &param_nodes, double coeff);


  public:

    virtual int extract_integer_points(const zjucad::matrix::matrix<size_t> &param_mesh,
                                       const zjucad::matrix::matrix<double> & param_nodes) {}

    virtual int mapping_to_orig_space(const zjucad::matrix::matrix<size_t> & orig_mesh,
                                      const zjucad::matrix::matrix<double> & orig_node) {}

    ///
    /// \brief hex_extractor::extract_points
    /// \param param_mesh
    /// \param param_nodes
    /// \param offset
    /// \param inv_base
    /// \param invertible
    /// \param point_mesh
    /// \param point_nodes
    /// \return
    ///
    int extract_points(const zjucad::matrix::matrix<size_t>               &param_mesh,
                       const zjucad::matrix::matrix<double>               &param_nodes,
                       const double                                       offset,
                       const zjucad::matrix::matrix<double>               &inv_base,
                       const zjucad::matrix::matrix<bool>                 &invertible,
                       std::vector<point_info>                            &point_mesh,
                       std::vector<boost::tuple<double, double, double> > &point_nodes);


    ///
    /// \brief hex_extractor::merge_points_with_same_orig_coord
    /// \param param_mesh
    /// \param point_nodes
    /// \param point_mesh
    ///
    void merge_points_with_same_orig_coord(const zjucad::matrix::matrix<size_t>                     &param_mesh,
                                           const std::vector<boost::tuple<double, double, double> > &point_nodes,
                                           std::vector<point_info>                                  &point_mesh);

    ///
    /// @brief extract_integer_points extract integer points, each should be
    ///                               associated with its coordinates and
    ///                               its tet/tri idx, and its barycentric coordinate
    /// @param param_mesh             input param mesh
    /// @param param_nodes            input param nodes
    /// @return                       0 if works fine, or non-zeros
    ///
    virtual int extract_integer_points(const zjucad::matrix::matrix<size_t>                 &param_mesh,
                                       const zjucad::matrix::matrix<double>                 &param_nodes,
                                       const zjucad::matrix::matrix<double>                 &inv_base,
                                       const zjucad::matrix::matrix<bool>                   &invertible,
                                       std::vector<point_info>                              &int_mesh,
                                       std::vector<boost::tuple<double, double, double> >   &int_nodes) { }

    ///
    /// \brief hex_extractor::merge_int_points
    /// \param int_mesh
    /// \param orig_mesh
    /// \param orig_nodes
    ///
    void merge_int_points(std::vector<point_info>                &int_mesh,
                          const zjucad::matrix::matrix<size_t>   &orig_mesh,
                          const zjucad::matrix::matrix<double>   &orig_nodes);

    ///
    /// @brief extract_integer_lines extract integer lines according to integer
    ///                              points and transition on each faces (tet)
    /// @param param_mesh            input mesh
    /// @param param_nodes           input param node
    /// @param inner_type            transition type defined on each two cells (cells should be ordered)
    /// @return                      0 if works fine, or non-zeros
    ///
    virtual int extract_integer_lines(const zjucad::matrix::matrix<size_t> & param_mesh,
                                      const zjucad::matrix::matrix<double> &param_nodes,
                                      const std::map<std::pair<size_t,size_t>, size_t> & inner_type) {}

    ///
    /// @brief extract_mesh         extract mesh from integer lines
    /// @param param_mesh           input param mesh
    /// @param param_nodes          input nodes
    /// @return                     0 if works fine, or non-zeros
    ///
    virtual int extract_mesh(const zjucad::matrix::matrix<size_t> & param_mesh,
                             const zjucad::matrix::matrix<double> &param_nodes) { }


    ///
    ///
    ///
    ///
    ///
    virtual int mapping_to_orig_space(const zjucad::matrix::matrix<size_t>    &orig_mesh,
                                      const zjucad::matrix::matrix<double>    &orig_node,
                                      const std::vector<point_info>           &int_mesh,
                                      const zjucad::matrix::matrix<size_t>    &media_mesh,
                                      const zjucad::matrix::matrix<double>    &media_nodes,
                                      zjucad::matrix::matrix<size_t>          &final_mesh,
                                      zjucad::matrix::matrix<double>          &final_nodes);


   

    ///
    /// @brief bound_ith_tet          find the bounding cube of a tet
    /// @param 
    int bound_ith_tet(size_t id, 
                      const zjucad::matrix::matrix<size_t>  &param_mesh,
                      const zjucad::matrix::matrix<double>  &param_nodes,
                      bounding_box                          &box);

    ///
    /// @brief
    ///
    ///
    /// @return                       0 if works fine, otherwise non-zero
    bool line_cross_tet(const size_t                            id,
                        const zjucad::matrix::matrix<size_t>    &param_mesh,
                        const zjucad::matrix::matrix<double>    &param_nodes,
                        const zjucad::matrix::matrix<double>    &inv_base,
                        const zjucad::matrix::matrix<bool>      &invertible,
                        const zjucad::matrix::matrix<double>    &src,
                        const zjucad::matrix::matrix<double>    &des,
                        const int                               type,
                        zjucad::matrix::matrix<double>          &cross,
                        zjucad::matrix::matrix<size_t>          &face_idx,
                        zjucad::matrix::matrix<int>             &flag);

    ///
    /// @brief                    find the crossover point of a line through p0 p1 and a plane through v0 v1 v2
    ///
    ///
    ///
    /// @return                   0 if the line dosen't cross the plane,1 if there is only one cross point,
    ///                           2 if the line is on the plane
    int line_cross_plane(const zjucad::matrix::matrix<double>      &p0,
                          const zjucad::matrix::matrix<double>     &p1,
                          const zjucad::matrix::matrix<double>     &v0,
                          const zjucad::matrix::matrix<double>     &v1,
                          const zjucad::matrix::matrix<double>     &v2,
                          zjucad::matrix::matrix<double>           &ans);


    ///
    /// @brief                         function used for debugging, remove later
    ///
    void see(const zjucad::matrix::matrix<double> v)
    {
        //std::cout << v(0, 0) << " " << v(1, 0) << " " << v(2, 0) << std::endl;
        printf("%.12lf %.12lf %.12lf\n", v(0, 0), v(1, 0), v(2, 0));
    }

   ///
   ///
   ///
   ///
   ///
   ///
   int mesh_orderring(const std::vector<boost::tuple<double, double, double> >   &int_nodes,
                      std::vector<std::vector<size_t> >                          &hex_mesh );

   ///
   /// \brief oriented_area
   /// \param id
   /// \param param_mesh
   /// \param param_nodes
   /// \return
   ///
   double oriented_area(const size_t                          id,
                        const zjucad::matrix::matrix<size_t>  &param_mesh,
                        const zjucad::matrix::matrix<double>  &param_nodes);


   ///
   /// \brief cast
   /// \param x
   /// \return
   ///
   double cast(const double x, const double eps);


   ///
   /// \brief extract_half_integer_points
   /// \param param_mesh
   /// \param param_nodes
   /// \param half_mesh
   /// \param half_nodes
   /// \return
   ///
   int extract_half_integer_points(const zjucad::matrix::matrix<size_t>                &param_mesh,
                                   const zjucad::matrix::matrix<double>                &param_nodes,
                                   const zjucad::matrix::matrix<double>                &inv_base,
                                   const zjucad::matrix::matrix<bool>                  &invertible,
                                   std::vector<boost::tuple<size_t,size_t> >           &half_mesh,
                                   std::vector<boost::tuple<double, double, double> >  &half_nodes);

   ///
   /// \brief generate_hex
   /// \param param_mesh
   /// \param param_nodes
   /// \param int_mesh
   /// \param int_nodes
   /// \param half_mesh
   /// \param half_nodes
   /// \param hex_mesh
   /// \return
   ///
   virtual int extract_mesh(const zjucad::matrix::matrix<size_t>                      &orig_mesh,
                            const zjucad::matrix::matrix<double>                      &orig_nodes,
                            const zjucad::matrix::matrix<size_t>                      &param_mesh,
                            const zjucad::matrix::matrix<double>                      &param_nodes,
                            const zjucad::matrix::matrix<double>                      &inv_base,
                            const zjucad::matrix::matrix<bool>                        &invertible,
                            const std::vector<point_info>                             &int_mesh,
                            const std::vector<boost::tuple<double, double, double> >  &int_nodes,
                            const std::vector<point_info>                             &center_mesh,
                            const std::vector<boost::tuple<double, double, double> >  &center_nodes,
                            const std::map<std::pair<size_t, size_t>, size_t>         &inner_type,
                            zjucad::matrix::matrix<size_t>                            &hex_mesh,
                            zjucad::matrix::matrix<double>                            &hex_nodes);

   ///
   /// \brief walk
   /// \param src
   /// \param id
   /// \param dir
   /// \param dis
   /// \param orig_handle
   /// \param cut_tet2tet
   /// \param inner_type
   /// \param param_mesh
   /// \param param_nodes
   /// \param des
   /// \return
   ///
   void walk(const zjucad::matrix::matrix<double>                      &curr_src,
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
             size_t                                                    direction_change,
             bool                                                      &arrive);

   ///
   /// \brief extract_center_points
   /// \param param_mesh
   /// \param param_nodes
   /// \param center_mesh
   /// \param center_nodes
   /// \return
   ///
   int extract_center_points(const zjucad::matrix::matrix<size_t>                &param_mesh,
                             const zjucad::matrix::matrix<double>                &param_nodes,
                             const zjucad::matrix::matrix<double>                &inv_base,
                             const zjucad::matrix::matrix<bool>                  &invertible,
                             std::vector<point_info>                             &center_mesh,
                             std::vector<boost::tuple<double, double, double> >  &center_nodes);

   ///
   /// \brief calculate_tet_base_inverse
   /// \param param_mesh
   /// \param param_nodes
   /// \param inv_base
   /// \param invertible
   ///
   void calculate_tet_base_inverse(const zjucad::matrix::matrix<size_t>    &param_mesh,
                                   const zjucad::matrix::matrix<double>    &param_nodes,
                                   zjucad::matrix::matrix<double>          &inv_base,
                                   zjucad::matrix::matrix<bool>            &invertible);

   ///
   /// \brief contains
   /// \param P
   /// \param id
   /// \param inv_base
   /// \param invertible
   /// \param param_mesh
   /// \param param_nodes
   /// \param _coeff
   /// \return
   ///
   bool contains(const zjucad::matrix::matrix<double>     &P,
                 const int                                id,
                 const zjucad::matrix::matrix<double>     &inv_base,
                 const zjucad::matrix::matrix<bool>       &invertible,
                 const zjucad::matrix::matrix<size_t>     &param_mesh,
                 const zjucad::matrix::matrix<double>     &param_nodes,
                 zjucad::matrix::matrix<double>           &_coeff);

   ///
   /// \brief hex_extractor::build_dual_mesh
   /// \param param_mesh
   /// \param param_nodes
   /// \param orig_mesh
   /// \param inner_type
   /// \param edge
   /// \return
   ///
   int build_dual_mesh(const zjucad::matrix::matrix<size_t>               &param_mesh,
                       const zjucad::matrix::matrix<double>               &param_nodes,
                       const zjucad::matrix::matrix<size_t>               &orig_mesh,
                       const std::map<std::pair<size_t, size_t>, size_t>  &inner_type,
                       std::map<std::pair<size_t, size_t>, double>        &edge);

   ///
   /// \brief see_dual_graph
   /// \param orig_nodes
   /// \param edge
   /// \param dual_mesh
   /// \param dual_nodes
   ///
   void see_dual_graph(const zjucad::matrix::matrix<size_t>               &orig_mesh,
                       const zjucad::matrix::matrix<double>               &orig_nodes,
                       const std::map<std::pair<size_t, size_t>, double>  &edge,
                       zjucad::matrix::matrix<size_t>                     &dual_mesh,
                       zjucad::matrix::matrix<double>                     &dual_nodes);

   ///
   /// \brief choose_closest_point
   /// \param center_id
   /// \param center_mesh
   /// \param center_nodes
   /// \param orig_mesh
   /// \param orig_nodes
   /// \param seg
   /// \param graph
   /// \return
   ///
   std::pair<std::pair<size_t, size_t>, double> choose_closest_point(const size_t                                              &center_id,
                                                                     const std::vector<point_info>                             &center_mesh,
                                                                     const std::vector<boost::tuple<double, double, double> >  &center_nodes,
                                                                     const zjucad::matrix::matrix<size_t>                      &orig_mesh,
                                                                     const zjucad::matrix::matrix<double>                      &orig_nodes,
                                                                     const RANGE                                               &seg,
                                                                     const boost::shared_ptr<vertex_connection<DIRECT> >       &graph);

   ///
   /// \brief has_common_part
   /// \param a
   /// \param b
   /// \param param_mesh
   /// \return
   ///
   bool has_common_part(const size_t                          a,
                        const size_t                          b,
                        const zjucad::matrix::matrix<size_t>  &param_mesh);

   ///
   /// \brief mesh_orderring
   /// \param final_nodes
   /// \param final_mesh
   /// \return
   ///
   int mesh_orderring(const zjucad::matrix::matrix<double>  &final_nodes,
                      zjucad::matrix::matrix<size_t>        &final_mesh);

   ///
   /// \brief orig_coord
   /// \param id
   /// \param param_mesh
   /// \param param_nodes
   /// \param int_mesh
   /// \param o_coord
   ///
   void orig_coord(const size_t                                              id,
                   const zjucad::matrix::matrix<size_t>                      &param_mesh,
                   const zjucad::matrix::matrix<double>                      &param_nodes,
                   const std::vector<point_info>                             &int_mesh,
                   zjucad::matrix::matrix<double>                            &o_coord);

   ///
   /// \brief Find
   /// \param cid
   /// \param father
   /// \return
   ///
   size_t Find(const size_t cid);

   ///
   /// \brief Union
   /// \param p
   /// \param q
   /// \param father
   ///
   void Union(const size_t p, const size_t q);

   ///
   /// \brief merging_vertices_via_cut_face
   /// \param param_mesh
   /// \param param_nodes
   /// \param int_mesh
   /// \param int_nodes
   /// \param inner_type
   /// \param orig_handle
   /// \param cut_tet2tet
   /// \param father
   ///
   void merging_vertices_via_cut_face(const zjucad::matrix::matrix<size_t>                    &param_mesh,
                                    const zjucad::matrix::matrix<double>                      &param_nodes,
                                    const std::vector<point_info>                             &int_mesh,
                                    const std::vector<boost::tuple<double, double, double> >  &int_nodes,
                                    const std::map<std::pair<size_t, size_t>, size_t>         &inner_type,
                                    const boost::shared_ptr<jtf::mesh::face2tet_adjacent>     &orig_handle,
                                    const zjucad::matrix::matrix<size_t>                      &cut_tet2tet);

   ///
   /// \brief point_in_tri_face
   /// \param P
   /// \param v0
   /// \param v1
   /// \param v2
   /// \return
   ///
   bool point_in_triangle(const zjucad::matrix::matrix<double>  &P,
                          const zjucad::matrix::matrix<double>  &v0,
                          const zjucad::matrix::matrix<double>  &v1,
                          const zjucad::matrix::matrix<double>  &v2);

   ///
   /// \brief hex_extractor::collapse
   /// \param hex_nodes
   /// \param hex_mesh
   ///
   bool collapse(std::vector<std::vector<size_t> >     &hex_mesh,
                 const boost::unordered_set<size_t>    &feature_points);


  };
}

#endif
