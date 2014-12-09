#include <iostream>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/util/vertex_connection.h>
#include <jtflib/mesh/trimesh.h>

using namespace std;
using namespace zjucad::matrix;

size_t find_closest_point(const matrix<double> & p,
                          const matrix<double> & nodes)
{
  assert(nodes.size(1) == p.size(1));

  static vector<double> diff(nodes.size(2),1);

  size_t i = 0;
#pragma omp parallel for private(i)
  for(i = 0; i < diff.size(); ++i)
    diff[i] = norm(nodes(colon(),i)-p);
  return size_t(min_element(diff.begin(), diff.end())-diff.begin());
}

int transplant_feature_line(int argc, char *argv[])
{
  if(argc != 4){
      cerr << "# [usage] transplant_feature_line orig_obj orig_fl new_obj" << endl;
      return __LINE__;
    }

  jtf::mesh::tri_mesh trm(argv[1]);
  jtf::mesh::tri_mesh trm2(argv[3]);

  vector<vector<size_t> > feature_line;
  if(jtf::mesh::load_feature_line(argv[2], feature_line))
    throw std::invalid_argument("invalid feature line.");

  map<pair<size_t,size_t>,double> edge_map;
  for(size_t ei = 0; ei < trm2.ea_->edges_.size(); ++ei){
      edge_map[trm2.ea_->edges_[ei]] = norm(trm2.trimesh_.node_(colon(), trm2.ea_->edges_[ei].first)
                                            - trm2.trimesh_.node_(colon(), trm2.ea_->edges_[ei].second));
    }
  shared_ptr<vertex_connection<UNDIRECT> > vc(vertex_connection<UNDIRECT>::create(edge_map));
  vector<size_t> path;
  vector<vector<size_t> > new_feature_lines(feature_line.size());
  for(size_t i = 0; i  < feature_line.size(); ++i){
      for(size_t pi = 0; pi < feature_line[i].size(); ++pi){
          size_t closet_point = find_closest_point(trm.trimesh_.node_(colon(), feature_line[i][pi]), trm2.trimesh_.node_);
          if(new_feature_lines[i].empty()) new_feature_lines[i].push_back(closet_point);
          else{
              const vector<size_t> & one_feature_line = new_feature_lines[i];
              if(trm2.ea_->get_edge_idx(one_feature_line.back(), closet_point) != -1)
                new_feature_lines[i].push_back(closet_point);
               else{
                  if(vc->get_shortest_path(one_feature_line.back(), closet_point, path))
                    throw std::invalid_argument("strange: can not find path");
                  else{
                      for(size_t j = 1; j < path.size(); ++j)
                        new_feature_lines[i].push_back(path[j]);
                    }
                }
            }
        }
    }

  jtf::mesh::save_feature_line("new.fl", new_feature_lines);
  return 0;
}
