#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/util/util.h>
#include <iostream>
#include <boost/unordered_set.hpp>
#include "../common/util.h"

using namespace std;
using namespace zjucad::matrix;

template <typename T>
class group : public boost::unordered_set<T>
{
public:
  typedef typename boost::unordered_set<T>::const_iterator const_iterator;
public:

  void operator << (const T & a){
    this->insert(a);
  }
  void operator << (const group<T> & a){
    this->insert(a.begin(), a.end());
  }
  void operator >> (group<T> &a)const{
    a.insert(this->begin(), this->end());
  }
  void merge(group<T> &a){
    if(&a == this)
      return ;
    this->insert(a.begin(), a.end());
    a.clear();
  }
  friend std::ostream& operator << (std::ostream &output,
                                    const group<T> &g)
  {
    output << "# [info] group has: ";
    for(group<T>::const_iterator it = g.begin(); it != g.end(); ++it)
      output << *it << " ";
    output << std::endl;
    return output;
  }
};


int merge_obj(int argc, char * argv[])
{
  if(argc != 3 && argc != 4){
      cerr << "# [error] merge_obj input_obj output_obj [epsilon]" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes obj;
  if(jtf::mesh::load_obj(argv[1], obj.mesh_, obj.node_))
    return __LINE__;


  double avg_len = 0;
  size_t edge_num = 0;
  for(size_t fi = 0; fi < obj.mesh_.size(2); ++fi){
      for(size_t pi = 0 ; pi < obj.mesh_.size(1); ++pi){
          avg_len += norm(obj.node_(colon(), obj.mesh_(pi,fi)) -
                          obj.node_(colon(), obj.mesh_((pi+1)%obj.mesh_.size(1),fi)));
          ++edge_num;
        }
    }
  avg_len /= edge_num;

  double epsilon = 1e-6;
  if(argc == 4)
    epsilon = atof(argv[3]);
  set<pair<size_t,size_t> > merged_v;
  for(size_t pi = 0; pi < obj.node_.size(2); ++pi){

      double min_len = std::numeric_limits<double>::infinity();
      size_t picked_point_index = -1;
      for(size_t pj = 0; pj < obj.node_.size(2); ++pj){
          if(pi == pj) continue;
          const double len = norm(obj.node_(colon(), pi) - obj.node_(colon(), pj));
          if(len < min_len){
              min_len = len;
              picked_point_index = pj;
            }
          if(len < epsilon) break;
        }
      if(min_len > epsilon) continue;
      pair<size_t,size_t> one_pair(pi,picked_point_index);
      if(one_pair.first > one_pair.second) swap(one_pair.first, one_pair.second);
      merged_v.insert(one_pair);
    }

  cerr << "# [info] merged_vertex_pairs " << merged_v.size() << endl;

  vector<group<size_t> > groups(obj.node_.size(2));
  for(size_t i = 0; i < groups.size(); ++i) groups[i] << i;
  for(const auto & one_pair : merged_v){
      groups[one_pair.first].merge(groups[one_pair.second]);
    }

  for(const auto & one_group : groups){
      if(one_group.empty()) continue;
      const size_t target = *(one_group.begin());
      for(const auto & one_idx : one_group){
          if(one_idx == target) continue;
          for(auto & p : obj.mesh_){
              if(p == one_idx)
                p = target;
            }
        }
    }

  remove_extra_node(obj.mesh_, obj.node_);

  cerr << "# [info] vertex " << obj.node_.size(2) << " face " << obj.mesh_.size(2) << endl;
  if(jtf::mesh::save_obj(argv[2], obj.mesh_, obj.node_))
    return __LINE__;

  cerr << "# [info] succeed." << endl;
  return 0;
}


int test(int argc, char * argv[])
{
  if(argc != 4){
      return __LINE__;
    }
  jtf::mesh::meshes obj;
  if(jtf::mesh::load_obj(argv[1], obj.mesh_, obj.node_))
    return __LINE__;

  set<size_t> face_idx;
  ifstream ifs(argv[2]);
  size_t idx;
  while(!ifs.eof()){
      ifs >> idx;
      face_idx.insert(idx);
    }
  matrix<size_t> new_faces(obj.mesh_.size(1), face_idx.size());
  size_t face_i = 0;
  for(const auto & one_face_idx: face_idx){
      new_faces(colon(),face_i++) = obj.mesh_(colon(), one_face_idx);
    }

  remove_extra_node(new_faces, obj.node_);
  cerr << "# [info] vertex " << obj.node_.size(2) << " face " << new_faces.size(2) << endl;
  if(jtf::mesh::save_obj(argv[3], new_faces, obj.node_))
    return __LINE__;
  return 0;
}
