#include <iostream>
#include <fstream>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>

#include "util.h"
using namespace std;
using namespace zjucad::matrix;

namespace jtf
{
namespace quadmesh
{

using namespace jtf::mesh;

int extract_quadmesh_feature_line(
        const matrixst & quad,
        const matrixd & node,
        matrixst & feature_line,
        const double cosin_theta)
{
    // TODO: left for shenxinxin
    static std::shared_ptr<edge2cell_adjacent> edge_quad_mesh(edge_quad_mesh -> create(quad));
    if(edge_quad_mesh == 0)
        return __LINE__;

    size_t face1_id, face2_id;
    matrixd face1_normal, face2_normal;
    size_t feature_line_num ;
    vector<size_t> feature_vec;
    for(size_t i = 0; i < edge_quad_mesh -> edges_.size(); ++i){
        if(! edge_quad_mesh -> is_boundary_edge(edge_quad_mesh -> edge2cell_[i])){
            face1_id = edge_quad_mesh -> edge2cell_[i].first;
            face2_id = edge_quad_mesh -> edge2cell_[i].second;
            face1_normal = cross( (node(colon(), quad(0,face1_id)) -
                                   node(colon(), quad(1,face1_id))),
                                  (node(colon(), quad(1,face1_id)) -
                                   node(colon(), quad(2,face1_id))) );
            if(norm(face1_normal) < 1e-8){
                face1_normal = cross( (node(colon(), quad(1,face1_id)) -
                                       node(colon(), quad(2,face1_id))),
                                      (node(colon(), quad(2,face1_id)) -
                                       node(colon(), quad(3,face1_id))) );
                if(norm(face1_normal) < 1e-8){
                    cerr << "# [info] meet degenerated quad: " << endl;
                    cerr << "# ------ ";
                    cerr << quad(0, face1_id) << " "
                         << quad(1, face1_id) << " "
                         << quad(2, face1_id) << " "
                         << quad(3, face1_id) << endl;
                    continue;
                }
            }

            face2_normal = cross( (node(colon(), quad(0,face2_id)) -
                                   node(colon(), quad(1,face2_id))),
                                  (node(colon(), quad(1,face2_id)) -
                                   node(colon(), quad(2,face2_id))) );
            if(norm(face2_normal) < 1e-8)
            {
                face2_normal = cross( (node(colon(), quad(1,face2_id)) -
                                       node(colon(), quad(2,face2_id))),
                                      (node(colon(), quad(2,face2_id)) -
                                       node(colon(), quad(3,face2_id))) );
                if(norm(face2_normal) < 1e-8){
                    cerr << "# [info] meet degenerated quad." << endl;
                    cerr << "# ------ ";
                    cerr << quad(0, face2_id) << " "
                         << quad(1, face2_id) << " "
                         << quad(2, face2_id) << " "
                         << quad(3, face2_id) << endl;
                    continue;
                }
            }
            face1_normal = face1_normal / norm(face1_normal);
            face2_normal = face2_normal / norm(face2_normal);
            if( fabs( dot(face1_normal, face2_normal)) <= cosin_theta )
            {
                feature_vec.push_back(edge_quad_mesh -> edges_[i].first);
                feature_vec.push_back(edge_quad_mesh -> edges_[i].second);
            }
        }
    }
    feature_line_num = feature_vec.size() / 2;
    feature_line.resize(2,feature_line_num);
    for(size_t j = 0; j < feature_vec.size(); ++j)
        feature_line[j] = feature_vec[j];

    return 0;
}

int extract_quadmesh_singularity_line(const char *path1,
                                      const char *path2)
{
  matrixd nodes;
  matrixst quad;
  jtf::mesh::load_obj(path1, quad, nodes);
   static std::shared_ptr<edge2cell_adjacent> edge_quad_mesh(edge_quad_mesh -> create(quad));
  if(!edge_quad_mesh)
    return __LINE__;

  //find adjacent nodes
  std::vector<vector<size_t> > adj_vertex(nodes.size(2));
  std::vector<vector<bool> > is_visited(nodes.size(2));
  const std::vector<std::pair<size_t, size_t> > edges = edge_quad_mesh->edges_;
  for(size_t i = 0; i < edges.size(); ++i)
    {
      adj_vertex[edges[i].first].push_back(edges[i].second);
      is_visited[edges[i].first].push_back(false);
      adj_vertex[edges[i].second].push_back(edges[i].first);
      is_visited[edges[i].second].push_back(false);
    }

  //std::sort(adj_vertex.begin(), adj_vertex.end(), cmp_vec_size);

  //find singularity line
  std::vector<vector<size_t> > singularity_lines;
//  for(size_t i = 0;  i < adj_vertex.size(); ++i)
//    if(adj_vertex[i].size() != 4)
//      std::cout << i <<": " << adj_vertex[i].size() << std::endl;
  for(size_t i = 0; i < adj_vertex.size(); ++i)
    {
      if(adj_vertex[i].size() == 4)
        continue;
      find_vertex_lines(i, singularity_lines, *edge_quad_mesh,
                        adj_vertex, quad, is_visited);
    }
  write_singularity_line(path2, singularity_lines);
  return 0;
}


// find the singularity lines around the singularity point
void find_vertex_lines(size_t index,
                       std::vector<std::vector<size_t> > &lines,
                       const edge2cell_adjacent &edge2quad,
                       const std::vector<std::vector<size_t> > &adj_ver,
                       const matrixst &quad,
                       std::vector<std::vector<bool> > &is_visited)
{
  size_t next_node, cur_node;
  std::vector<size_t> line;
  for(size_t i = 0; i < adj_ver[index].size(); ++i)
    {
      line.clear();
      if((is_visited[index])[i])
        continue;
      line.push_back(index);
      (is_visited[index])[i] = true;
      cur_node = index;
      next_node = (adj_ver[index])[i];
      while(adj_ver[next_node].size() == 4)
        {
          line.push_back(next_node);
          choose_next_node(cur_node, next_node,
                           edge2quad, adj_ver, quad);
        }
      line.push_back(next_node);
      lines.push_back(line);
      for(size_t j = 0; j < adj_ver[next_node].size(); ++j)
        {
          if((adj_ver[next_node])[j] == cur_node)
            (is_visited[next_node])[j] = true;
        }
    }
}

// find the next node accord the direction
void choose_next_node(size_t &cur_node, size_t &next_node,
                      const jtf::mesh::edge2cell_adjacent &edge2quad,
                      const vector<vector<size_t> > &adj_ver, const matrixst &quad)
{
  assert(adj_ver[next_node].size() == 4);
  size_t next;
  std::pair<size_t, size_t> adj_quad = edge2quad.query(cur_node, next_node);
  assert(adj_quad.first != -1 && adj_quad.second != -1);
  std::vector<size_t> adj_quad_index;
  for(size_t i = 0; i < quad.size(1); ++i)
    adj_quad_index.push_back(quad(i, adj_quad.first));
  for(size_t i = 0; i < quad.size(1); ++i)
    adj_quad_index.push_back(quad(i, adj_quad.second));

  for(size_t i = 0;  i < adj_ver[next_node].size(); ++i)
    {
      if( adj_ver[next_node][i] != cur_node &&
          !is_in_vec(adj_ver[next_node][i], adj_quad_index))
        {
          next = adj_ver[next_node][i];
          break;
        }
    }
  cur_node = next_node;
  next_node = next;
}

// write singularity lines to file
int write_singularity_line(const char* path2,
                           const std::vector<std::vector<size_t> > &singularity_lines)
{
  std::ofstream ofs(path2);
  if(!ofs)
    {
      std::cerr << "save singularity lines file error!" << std::endl;
      return 1;
    }

  size_t num = singularity_lines.size();
  ofs << num << std::endl;
  for(size_t i = 0; i < num; ++i)
    {
      ofs << singularity_lines[i].size() << std::endl;
      for(size_t j = 0; j < singularity_lines[i].size(); ++j)
        {
          ofs << (singularity_lines[i])[j] << " ";
        }
      ofs << std::endl;
    }
  std::cout << "save singularity lines file success!" << std::endl;
  return 0;
}

}
}
