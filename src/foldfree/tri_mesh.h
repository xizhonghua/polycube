#ifndef TRI_MESH_H
#define TRI_MESH_H

#include <vector>
#include <math.h>

#include<zjucad/matrix/matrix.h>

namespace jy{

struct tri_mesh
{
    size_t num_verts_;
    size_t num_faces_;

    double (*verts_)[3];
    size_t (*faces_)[3];


    tri_mesh()
        :verts_(NULL),faces_(NULL),
          num_verts_(0),num_faces_(0){
    }

    ~tri_mesh(){
        if(verts_){
            delete[] verts_;
        }
        if(faces_){
            delete[] faces_;
        }
    }
};

struct half_edge //prototype of the edge; used for mesh construction
{
    size_t face_id;
    size_t vertex_0; //adjacent vertices sorted by id value
    size_t vertex_1; //they are sorted, vertex_0 < vertex_1
};

////////////////////////////////////////////
bool read_obj_mesh(const std::string &file_name,tri_mesh &mesh);

bool write_obj_mesh(const std::string &file_name,tri_mesh &mesh,double *textures = NULL);


bool read_fixed_vertices(const std::string &file_name,
                         std::vector<size_t> &vert_ids,
                         std::vector<double> &vert_coords);

bool read_init_x(const std::string &file_name,std::vector<double> &x);

bool read_fix_and_init_vertices(const std::string &file_name,
                                std::vector<size_t> &vert_ids,
                                std::vector<double> &vert_coords,
                                std::vector<double> &x);

size_t get_vertex_id(const tri_mesh &mesh,const double point[3]);

size_t count_mesh_edges(const tri_mesh &mesh);

///////////////////////////////////////////////////////////////////////
inline bool operator < (const half_edge &x, const half_edge &y)
{
    return (x.vertex_0 == y.vertex_0) ? x.vertex_1 < y.vertex_1 : x.vertex_0 < y.vertex_0;
}

inline bool operator != (const half_edge &x, const half_edge &y)
{
    return x.vertex_0 != y.vertex_0 || x.vertex_1 != y.vertex_1;
}

inline bool operator == (const half_edge &x, const half_edge &y)
{
    return x.vertex_0 == y.vertex_0 && x.vertex_1 == y.vertex_1;
}

/////////////////////////////////////////////////////////////////////////
inline double dot(const double v1[3],const double v2[3])
{
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

inline double dist(const double v1[3],const double v2[3])
{
    double x = v1[0] - v2[0];
    double y = v1[1] - v2[1];
    double z = v1[2] - v2[2];

    return sqrt(x*x + y*y + z*z);
}

inline double dist2(const double v1[3],const double v2[3])
{
    double x = v1[0] - v2[0];
    double y = v1[1] - v2[1];
    double z = v1[2] - v2[2];

    return x*x + y*y + z*z;
}

}

#endif // TRI_MESH_H
