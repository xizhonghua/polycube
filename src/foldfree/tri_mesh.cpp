#include "tri_mesh.h"

#include <cstddef>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
//LOCAL
#include "inline.h"

namespace jy{

bool read_obj_mesh(const std::string &file_name, tri_mesh &mesh)
{
    //only OBJ format supported
    if ("obj" != szy::get_filename_suffix(file_name)) {
        std::cerr << __FILE__ << " " << __LINE__ << " "
                  << "only .obj supported" << std::endl;
        return false;
    }

    //build file object
    std::ifstream infile(file_name.c_str());
    if (!infile) {
        std::cerr << __FILE__ << " " << __LINE__ << " "
                  << "load file " << file_name << " for read failed" << std::endl;
        return false;
    }

    //get vertex and face number
    size_t num_vertices = 0, num_faces = 0;
    std::string line, word;
    while (getline(infile, line)) {
        if (szy::is_line_invalid(line))
            continue;

        std::istringstream instream(line);
        instream >> word;

        if ("v" == word || "V" == word) {
            ++num_vertices;
        } else if ('f' == word[0] || 'F' == word[0]) {
            ++num_faces;
        }
        word.clear();
    }
    infile.clear();
    infile.seekg(0, std::ios::beg);

    //allocate space for triangle mesh
    if (0 == num_vertices || 0 == num_faces) {
        std::cerr << __FILE__ << " " << __LINE__ << " "
                  << "no vertex or face data find" << std::endl;
        return false;
    }


    mesh.num_faces_ = num_faces;
    mesh.num_verts_ = num_vertices;
    if(mesh.verts_){
        delete[] mesh.verts_;
    }
    if(mesh.faces_){
        delete[] mesh.faces_;
    }
    mesh.verts_ = new double[num_vertices][3];
    mesh.faces_ = new size_t[num_faces][3];

    bool process_state = true;
    //push vertex and face data into triangle mesh
    size_t vertex_count = 0;
    size_t face_count = 0;
    while (getline(infile, line)) {
        if (szy::is_line_invalid(line))
            continue;

        std::istringstream instream(line);
        instream >> word;

        if ("v" == word || "V" == word) {
            double coordinate[3];
            double w = 1.0;
            instream >> coordinate[0] >> coordinate[1] >> coordinate[2];
            if (instream >> w) {
                if (w < 1e-6) {
                    std::cerr << __FILE__ << " " << __LINE__ << " "
                              << "error occured when read vertex coordinates" << std::endl;
                    process_state = false;
                    break;
                }
                coordinate[0] /= w;
                coordinate[1] /= w;
                coordinate[2] /= w;
            }

            double *v = mesh.verts_[vertex_count];

            v[0] = coordinate[0];
            v[1] = coordinate[1];
            v[2] = coordinate[2];

            ++vertex_count;
        } else if ('f' == word[0] || 'F' == word[0]) {
            //only triangle mesh supported
            std::string pair[3], test;
            instream >> pair[0] >> pair[1] >> pair[2] >> test;
            if (!test.empty() || !instream.eof()) {
                std::cerr << __FILE__ << " " << __LINE__ << " "
                          << "only triangle mesh supported" << std::endl;
                process_state = false;
                break;
            }

            //get vertex id in the this triangle face
            size_t vertex_ids[3];
            for (size_t i = 0; i < 3; ++i) {
                std::string vertex_str = pair[i].substr(0, pair[i].find('/'));

                sscanf(vertex_str.c_str(), "%lu", &vertex_ids[i]);
                --vertex_ids[i];
                if (vertex_ids[i] >= num_vertices) {
                    std::cerr << __FILE__ << " " << __LINE__ << " "
                              << "vertex index in face exceed limit" << std::endl;
                    process_state = false;
                    break;
                }
            }
            if (process_state != true)
                break;

            //push into faces
            size_t *face = mesh.faces_[face_count];
            face[0] = vertex_ids[0];
            face[1] = vertex_ids[1];
            face[2] = vertex_ids[2];

            ++face_count;
        }
        word.clear();
    }

    infile.close();

    if (process_state != true) {
        delete[] mesh.verts_;
        delete[] mesh.faces_;
        mesh.verts_ = NULL;
        mesh.faces_ = NULL;
        return process_state;
    }

    return true;
}


bool write_obj_mesh(const std::string &file_name, tri_mesh &mesh,double *textures)
{
    //only OBJ format supported
    if ("obj" != szy::get_filename_suffix(file_name)) {
        std::cerr << __FILE__ << " " << __LINE__ << " "
                  << "only .obj supported" << std::endl;
        return false;
    }

    //build file object
    std::ofstream outfile(file_name.c_str());
    if (!outfile) {
        std::cerr << __FILE__ << " " << __LINE__ << " "
                  << "load file " << file_name << " to write failed" << std::endl;
        return false;
    }

    outfile << "#Generated by QuadRangular" << std::endl;

    //write vertex coordinates infos
    const size_t num_vertices = mesh.num_verts_;
    for (size_t vid = 0; vid < num_vertices; ++vid) {
        //vertex &vert = vertices_[vid];
        //outfile << "v";

        double *v = mesh.verts_[vid];
        outfile << "v " << v[0]<< " " << v[1]<< " " << v[2]<<std::endl;
    }

    //write texture coordinates infos
    if(textures){
        for (size_t tid = 0; tid < num_vertices; ++tid) {
            outfile << "vt " << textures[2*tid] << " " << textures[2*tid+1]<< std::endl;
        }
    }

    //write face infos
    const size_t num_faces = mesh.num_faces_;
    for (size_t fid = 0; fid < num_faces; ++fid) {
        //faces_[fid].get_vertices_id(fv);
        outfile << "f";
        size_t *f = mesh.faces_[fid];
        for (size_t i = 0; i < 3; ++i)
            outfile << " " << f[i]+1 << "/" << f[i]+1;
        outfile << std::endl;
    }

    outfile.close();
    return true;
}

bool read_fixed_vertices(const std::string &file_name,
                         std::vector<size_t> &vert_ids,
                         std::vector<double> &vert_coords)
{
    size_t num_vert = 0;
    double coord1,coord2;
    std::string line, word;
    std::ifstream infile(file_name.c_str());
    if (!infile) {
        std::cerr <<"read fixed vertices failed!" << std::endl;
        return false;
    }
    while (getline(infile, line)) {
        if (szy::is_line_invalid(line))
            continue;

        std::istringstream instream(line);
        instream >> word;

        if ("#" != word) {
            ++num_vert;
        }

        word.clear();
    }
    infile.clear();
    infile.seekg(0, std::ios::beg);

    vert_ids.reserve(num_vert);
    vert_coords.reserve(num_vert*2);
    while (getline(infile, line)) {
        if (szy::is_line_invalid(line))
            continue;

        std::istringstream instream(line);
        instream >> word;

        if ("#" != word) {
            instream >> coord1 >> coord2;
            vert_ids.push_back(atoi(word.c_str())-1);
            vert_coords.push_back(coord1);
            vert_coords.push_back(coord2);
        }
        word.clear();
    }
    infile.close();

    return true;
}

bool read_init_x(const std::string &file_name, std::vector<double> &x)
{
    size_t index;
    double u,v;
    std::string line;
    std::ifstream infile(file_name.c_str());
    if (!infile) {
        std::cerr <<"read init x failed!" << std::endl;
        return false;
    }

    index = 0;
    while (getline(infile, line)) {
        if (szy::is_line_invalid(line))
            continue;

        std::istringstream instream(line);

        instream >> u >> v;

        x[2*index]   = u;
        x[2*index+1] = v;

        ++index;
    }
    infile.close();
}


bool read_fix_and_init_vertices(const std::string &file_name,
                                  std::vector<size_t> &vert_ids,
                                  std::vector<double> &vert_coords,
                                  std::vector<double> &x)
{
    size_t num_vert = 0,num_init = 0;
    double coord1,coord2;
    std::string line, word;
    std::ifstream infile(file_name.c_str());
    if (!infile) {
        std::cerr <<"read fixed vertices failed!" << std::endl;
        return false;
    }
    while (getline(infile, line)) {
        if (szy::is_line_invalid(line))
            continue;

        if("# init" == line){
            break;
        }

        std::istringstream instream(line);
        instream >> word;

        if ("#" != word) {
            ++num_vert;
        }
        word.clear();
    }

    while (getline(infile, line)) {
        if (szy::is_line_invalid(line))
            continue;

        std::istringstream instream(line);
        instream >> word;

        if ("#" != word) {
            ++num_init;
        }
        word.clear();
    }

    infile.clear();
    infile.seekg(0, std::ios::beg);

    vert_ids.reserve(num_vert);
    vert_coords.reserve(num_vert*2);
    while (getline(infile, line)) {
        if (szy::is_line_invalid(line))
            continue;

        if("# init" == line){
            break;
        }

        std::istringstream instream(line);
        instream >> word;

        if ("#" != word) {
            instream >> coord1 >> coord2;
            vert_ids.push_back(atoi(word.c_str())-1);
            vert_coords.push_back(coord1);
            vert_coords.push_back(coord2);
        }
        word.clear();
    }

    x.reserve(2*num_init);
    while (getline(infile, line)) {
        if (szy::is_line_invalid(line))
            continue;

        std::istringstream instream(line);
        instream >> word;

        if ("#" != word) {
            instream >> coord1 >> coord2;
            x.push_back(coord1);
            x.push_back(coord2);
        }
        word.clear();
    }
    infile.close();

    return true;
}

size_t get_vertex_id(const tri_mesh &mesh,const double point[3])
{
    size_t vertex_id = 0;
    double min_diff = 1.0;
    for (size_t i = 0; i < mesh.num_verts_; ++i) {
        double tmp = dist2(point,mesh.verts_[i]);
        if (tmp < min_diff) {
            min_diff = tmp;
            vertex_id = i;
        }
    }
    if (min_diff > 0.1)
        return 0;

    return vertex_id+1;
}

size_t count_mesh_edges(const tri_mesh &mesh)
{
    std::vector<half_edge> he(mesh.num_faces_*3);
    for(size_t i = 0,k = 0; i < mesh.num_faces_; ++i){
        size_t *f = mesh.faces_[i];
        for(size_t j = 0; j < 3; ++j,++k){
            he[k].face_id = i;
            size_t vertex_id_1 = f[j];
            size_t vertex_id_2 = f[(j+1) % 3];
            he[k].vertex_0 = std::min(vertex_id_1, vertex_id_2);
            he[k].vertex_1 = std::max(vertex_id_1, vertex_id_2);
        }
    }

    std::sort(he.begin(), he.end());

    size_t num_edges = 1;
    size_t num_hes = he.size();
    for(size_t i = 1; i< num_hes; ++i){
        if(he[i] != he[i-1]){
            ++num_edges;
        }
        else{
            if(i<num_hes-1){ //sanity check: there should be at most two equal half-edges
                assert(he[i] != he[i+1]);
            }
        }
    }

    return num_edges;
}


}
