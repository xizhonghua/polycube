#include <string>
#include <zjucad/ptree/ptree.h>
#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/io.h>
#include "../hex_param/topology_operation.h"
#include "../hex_param/io.h"
#include "../common/vtk.h"
#include "../equation_graph/equation_graph.h"
#include <jtflib/util/util.h>

using namespace std;
using namespace jtf::mesh;
using namespace zjucad::matrix;
using boost::property_tree::ptree;

void pt_descriptor(ptree &pt)
{
  pt.put("input/cut_tet.desc", "input cut tetmesh");
}

template <typename T1>
void get_edges_from_tet(
    const matrix_expression<T1> & tet,
    vector<pair<size_t,size_t> > & edges)
{
  set<pair<size_t,size_t> > edges_set;
  for(size_t ti = 0; ti < tet().size(2); ++ti){
      for(size_t pi = 0; pi < tet().size(1); ++pi){
          for(size_t pj = pi + 1; pj < tet().size(1); ++pj){
              pair<size_t,size_t> one_edge(tet()(pi,ti), tet()(pj,ti));
              if(one_edge.first > one_edge.second)
                swap(one_edge.first, one_edge.second);
              edges_set.insert(one_edge);
            }
        }
    }
  edges.resize(edges_set.size());
  copy(edges_set.begin(), edges_set.end(), edges.begin());
}

class surface_type_cleaner_cut{
public:
  surface_type_cleaner_cut(
      const matrix<size_t> & cut_tet,
      const matrix<double> & cut_node,
      const matrix<size_t> & uncut_tet,
      const matrix<double> & uncut_node,
      const boost::unordered_map<size_t,size_t> & surface_type_uncut,
      const vector<vector<size_t> > & node_group)
    :uncut_tet_(uncut_tet),
      uncut_node_(uncut_node), surface_type_uncut_(surface_type_uncut),
      node_group_(node_group) {init(cut_tet, cut_node);}

public:
  ///
  /// @brief remove_degeneration, detect surface type degeneration, and update surface_type_uncut
  /// @param surface_type_uncut
  ///
  void remove_degeneration(boost::unordered_map<size_t,size_t> &surface_type_uncut)const;

  ///
  /// @brief remove_degeneration_with_equation_graph, HANGOUT
  /// @param surface_type_uncut
  ///
  void remove_degeneration_with_equation_graph(
      boost::unordered_map<size_t, size_t> &surface_type_uncut,
      vector<vector<size_t> > & node_groups) const;

  /// @brief remove degeneration, based on cut tet mesh geometey node, remove zigzag and turning point
  /// @param surface type uncut
  void remove_degeneration_with_cut_node(
      boost::unordered_map<size_t,size_t> & surface_type_uncut);

  ///
  /// \brief remove_type_patch_degree_degeneration, this function detect type patches
  /// of one cut patch, and modify those whose degree is less than 3
  /// \param face_type
  /// \param type_patch_vec
  /// \param patch_faces
  /// \param ea
  ///
  void remove_type_patch_degree_degeneration(
      matrix<size_t> & face_type,
      std::vector<vector<size_t> > &type_patch_vec,
      const matrix<size_t> & patch_faces,
      const jtf::mesh::edge2cell_adjacent & ea);

protected:
  void init(const matrix<size_t> & cut_tet,
            const matrix<double> & cut_node){
    merged_cut_tet_ = cut_tet;
    merged_cut_node_ = cut_node;
    {
      if(node_group_.size()){
          matrix<size_t> cut_node_merged(cut_node.size(1),cut_node.size(2));
          for(size_t i = 0 ; i < cut_node_merged.size(); ++i) cut_node_merged[i] = i;

          for(const auto & one_group: node_group_){
              if(one_group.size() < 2) continue;
              for(size_t vi = 1; vi < one_group.size(); ++vi){
                  cut_node_merged[one_group[vi]] = cut_node_merged[one_group.front()];
                }
            }
          map<vector<size_t>,vector<size_t> > coord2points;
          vector<size_t> one_coord(3);
          for(size_t pi = 0; pi < cut_node_merged.size(2); ++pi){
              std::copy(cut_node_merged(colon(),pi).begin(),
                        cut_node_merged(colon(),pi).end(), one_coord.begin());
              coord2points[one_coord].push_back(pi);
            }

          matrix<size_t> point_mapping(max(merged_cut_tet_)+1,1);
          size_t one_group_idx = 0;
          merged_cut_node_.resize(3, coord2points.size()); // real nodes
          for(const auto & one_group_point: coord2points ){
              const vector<size_t> & one_vec = one_group_point.second;
              for(const auto & idx : one_vec){
                  point_mapping[idx] = one_group_idx;
                }
              merged_cut_node_(colon(), one_group_idx) = cut_node(colon(),one_vec.front());
              ++one_group_idx;
            }
          for(size_t i = 0; i < merged_cut_tet_.size(); ++i){
              merged_cut_tet_[i] = point_mapping[merged_cut_tet_[i]];
            }
        }
    }
    // this fa_cut will contains non-manifold surface, but it will not affect following steps
    fa_cut_.reset(jtf::mesh::face2tet_adjacent::create(merged_cut_tet_));
    fa_uncut_.reset(jtf::mesh::face2tet_adjacent::create(uncut_tet_));
    if(!fa_cut_.get() || !fa_uncut_.get()){
        throw std::invalid_argument("can not build face2tet_adjacent.");
      }
    jtf::mesh::get_outside_face(*fa_cut_, outside_face_cut_);

    cut_tet2tet_.resize(max(merged_cut_tet_)+1);
    cut_tet2tet_(merged_cut_tet_) = uncut_tet_(colon());

    vector<size_t> orig_face_in_cut_vec;
    for(size_t fi = 0 ;fi < outside_face_cut_.size(2); ++fi){
        const size_t face_idx =
            fa_uncut_->get_face_idx(cut_tet2tet_[outside_face_cut_(0,fi)],
            cut_tet2tet_[outside_face_cut_(1,fi)],
            cut_tet2tet_[outside_face_cut_(2,fi)]);
        assert(face_idx != -1);
        const auto it = surface_type_uncut_.find(face_idx);
        if(it == surface_type_uncut_.end()) continue;

        orig_face_in_cut_vec.insert(orig_face_in_cut_vec.end(),
                                    outside_face_cut_(colon(),fi).begin(),
                                    outside_face_cut_(colon(),fi).end());
        const size_t face_cut_idx = fa_cut_->get_face_idx(&outside_face_cut_(0,fi));
        assert(face_cut_idx != -1);
        surface_type_cut_[face_cut_idx] = it->second;
        orig_face_idx_in_cut_.push_back(face_cut_idx);
      }
    itr_matrix<const size_t*>  orig_face_in_cut_mat(3, orig_face_in_cut_vec.size()/3, &orig_face_in_cut_vec[0]);
    orig_face_in_cut_ = orig_face_in_cut_mat;

    {
      ofstream debug_ofs("cut_out_side_face.vtk");
      tri2vtk(debug_ofs, &merged_cut_node_[0], merged_cut_node_.size(2),
          &orig_face_in_cut_[0], orig_face_in_cut_.size(2));
    }

    ea_.reset(jtf::mesh::edge2cell_adjacent::create(orig_face_in_cut_));
    if(!ea_.get())
      throw std::logic_error("can not build edge2cell adjacent.");
  }

  void convert_cut_surface_type2uncut_surface_type(
      const boost::unordered_map<size_t,size_t> & cut_surface_type,
      boost::unordered_map<size_t,size_t> & uncut_surface_type)const
  {
    uncut_surface_type.clear();
    for(const auto & one_face: cut_surface_type){
        const vector<size_t> & one_face_vec = fa_cut_->faces_[one_face.first];
        const size_t face_idx_original = fa_uncut_->get_face_idx(
              cut_tet2tet_[one_face_vec[0]], cut_tet2tet_[one_face_vec[1]],
            cut_tet2tet_[one_face_vec[2]]);
        assert(face_idx_original != -1);
        uncut_surface_type[face_idx_original] = one_face.second;
      }
  }

  ///
  /// @brief get_restricted_edges_from_one_piece
  /// @param one_patch it stores idx to orig_face_in_cut_
  /// @param restricted_edges output restricted edges
  ///
  void get_restricted_edges_from_one_piece(
      const vector<size_t> &one_patch, vector<pair<size_t,size_t> > & restricted_edges)const;

  ///
  /// @brief handle_each_turning_point
  /// @param idx
  ///
  void handle_each_turning_point(
      const size_t idx,
      const pair<size_t,size_t> & two_point_ends,
      const map<size_t,size_t> &type2number);

  ///
  /// \brief handle_each_turning_point_with_type_patch, it modify face type, depends its group number
  /// \param idx
  /// \param two_point_ends
  /// \param type_patch
  /// \param one_type_face2patches
  ///
  void handle_each_turning_point_with_type_patch(
      const size_t idx,
      const pair<size_t,size_t> & two_point_ends,
      vector<set<size_t> > &type_patch,
      map<size_t,size_t> & one_type_face2patches);

  void update_type2number_one_piece(
      const vector<size_t> & one_patch,
      map<size_t,size_t>  &type2number);
public:
  unique_ptr<jtf::mesh::face2tet_adjacent> fa_cut_, fa_uncut_;
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea_;
  matrix<size_t> outside_face_cut_;
  matrix<size_t> cut_tet2tet_;
  matrix<size_t> orig_face_in_cut_;
  vector<size_t> orig_face_idx_in_cut_;
  matrix<size_t> merged_cut_tet_;
  matrix<double> merged_cut_node_;
  const matrix<size_t> & uncut_tet_;
  const matrix<double> & uncut_node_;
  const boost::unordered_map<size_t,size_t> & surface_type_uncut_;
  boost::unordered_map<size_t,size_t> surface_type_cut_;
  const vector<vector<size_t> > &node_group_;
private:
  surface_type_cleaner_cut();
};

void surface_type_cleaner_cut::update_type2number_one_piece(
    const vector<size_t> & one_patch,
    map<size_t,size_t>  &type2number)
{
  for(size_t fi = 0; fi < one_patch.size(); ++fi){
      const auto it = surface_type_cut_.find(orig_face_idx_in_cut_[one_patch[fi]]);
      assert(it != surface_type_cut_.end());
      auto type_it = type2number.find(it->second);
      if(type_it == type2number.end()) type2number[it->second] = 1;
      else
        ++type_it->second;
    }
}

void surface_type_cleaner_cut::handle_each_turning_point(
    const size_t idx,
    const pair<size_t,size_t> & two_point_ends,
    const map<size_t,size_t> &type2number)
{
  size_t two_ends[2] = {two_point_ends.first, two_point_ends.second};
  for(size_t i = 0; i < 2; ++i){
      const size_t edge_idx = ea_->get_edge_idx(idx, two_ends[i]);
      assert(edge_idx != -1);
      const pair<size_t,size_t> & two_cell = ea_->edge2cell_[edge_idx];
      auto face_type_first = surface_type_cut_.find(orig_face_idx_in_cut_[two_cell.first]);
      auto face_type_second = surface_type_cut_.find(orig_face_idx_in_cut_[two_cell.second]);

      assert(face_type_first != surface_type_cut_.end());
      assert(face_type_second != surface_type_cut_.end());

      if(face_type_first->second == face_type_second->second) {// has already been handled

          continue;
        }
      const auto it_first = type2number.find(face_type_first->second);
      const auto it_second = type2number.find(face_type_second->second);
      const size_t face_type_number_first = (it_first != type2number.end()? it_first->second:0);
      const size_t face_type_number_second = (it_second != type2number.end()? it_second->second:0);

      if(face_type_number_first < face_type_number_second){
          face_type_first->second = face_type_second->second;
        }else{
          face_type_second->second = face_type_first->second;
        }
    }
}

void surface_type_cleaner_cut::handle_each_turning_point_with_type_patch(
    const size_t idx,
    const pair<size_t,size_t> & two_point_ends,
    vector<set<size_t> > &type_patch,
    map<size_t,size_t> & one_type_face2patches)
{
  size_t two_ends[2] = {two_point_ends.first, two_point_ends.second};
  for(size_t i = 0; i < 2; ++i){
      const size_t edge_idx = ea_->get_edge_idx(idx, two_ends[i]);
      assert(edge_idx != -1);
      const pair<size_t,size_t> & two_cell = ea_->edge2cell_[edge_idx];

      auto face2patch_first = one_type_face2patches.find(two_cell.first);
      auto face2patch_second = one_type_face2patches.find(two_cell.second);

      assert(face2patch_first != one_type_face2patches.end());
      assert(face2patch_second != one_type_face2patches.end());

      if(face2patch_first->second == face2patch_second->second) {// has already been handled
          continue;
        }

      auto face_type_first = surface_type_cut_.find(orig_face_idx_in_cut_[two_cell.first]);
      auto face_type_second = surface_type_cut_.find(orig_face_idx_in_cut_[two_cell.second]);

      if(face_type_first->second == face_type_second->second) continue;

      auto & first_patch = type_patch[face2patch_first->second];
      auto & second_patch = type_patch[face2patch_second->second];

      if(first_patch.size() > second_patch.size()){
          face2patch_second->second = face2patch_first->second;
          first_patch.insert(two_cell.second);
          auto it_second = second_patch.find(two_cell.second);
          assert(it_second != second_patch.end());
          second_patch.erase(it_second);
          face_type_second->second = face_type_first->second;
        }else{
          face2patch_first->second = face2patch_second->second;
          second_patch.insert(two_cell.first);
          auto it_first = first_patch.find(two_cell.first);
          assert(it_first != first_patch.end());
          first_patch.erase(it_first);
          face_type_first->second = face_type_second->second;
        }
    }
}

void surface_type_cleaner_cut::remove_degeneration(
    boost::unordered_map<size_t, size_t> &updated_surface_type_uncut) const
{
  boost::unordered_map<size_t, size_t> updated_surface_type_cut = surface_type_cut_;
  boost::unordered_set<pair<size_t,size_t> > boundary_edges;
  vector<vector<size_t> > patches;
  matrix<size_t> boundary_edge_mat;
  jtf::mesh::get_boundary_edge(*ea_, boundary_edge_mat);
  for(size_t ei = 0; ei < boundary_edge_mat.size(2); ++ei){
      if(boundary_edge_mat(0,ei) > boundary_edge_mat(1, ei))
        boundary_edges.insert(make_pair(boundary_edge_mat(1,ei), boundary_edge_mat(0,ei)));
      else
        boundary_edges.insert(make_pair(boundary_edge_mat(0,ei), boundary_edge_mat(1,ei)));
    }

  get_face_patches_according_to_boundary(orig_face_in_cut_, *ea_, boundary_edges, patches);

  vector<size_t> arounding_type;
  vector<size_t> arounding_idx;

  matrix<size_t> count_type_number = zeros<size_t>(3,1); // 3 faces

  for(size_t gi = 0; gi < patches.size(); ++gi){
      const vector<size_t> & one_patch = patches[gi];
      {//visualization
        matrix<size_t> one_patch_faces(3, one_patch.size());
        for(size_t fi = 0; fi < one_patch.size(); ++fi){
            for(size_t di = 0; di < 3; ++di){
                one_patch_faces(di,fi) = cut_tet2tet_[orig_face_in_cut_(di,one_patch[fi])];
              }
          }
        stringstream gi_ss;
        gi_ss << "uncut_patach_" << gi << ".vtk";
        ofstream ofs(gi_ss.str().c_str());
        tri2vtk(ofs, &uncut_node_[0], uncut_node_.size(2), &one_patch_faces[0], one_patch_faces.size(2));
      }

      {//visualization
        matrix<size_t> one_patch_faces(3, one_patch.size());
        for(size_t fi = 0; fi < one_patch.size(); ++fi){
            for(size_t di = 0; di < 3; ++di){
                one_patch_faces(di,fi) = orig_face_in_cut_(di,one_patch[fi]);
              }
          }
        stringstream gi_ss;
        gi_ss << "cut_patach_" << gi << ".vtk";
        ofstream ofs(gi_ss.str().c_str());
        tri2vtk(ofs, &merged_cut_node_[0], merged_cut_node_.size(2), &one_patch_faces[0], one_patch_faces.size(2));
      }
      while(1){
          bool find_zigzag_edge = false;
          for(size_t fi = 0; fi < one_patch.size(); ++fi){
              const size_t &face_idx_in_orig_face_in_cut_ = one_patch[fi];
              auto my_type_it =
                  updated_surface_type_cut.find(orig_face_idx_in_cut_[face_idx_in_orig_face_in_cut_]);

              arounding_idx.clear();
              arounding_type.clear();

              assert(my_type_it != updated_surface_type_cut.end());


              for(size_t pi = 0; pi < orig_face_in_cut_.size(1); ++pi){
                  const size_t edge_idx = ea_->get_edge_idx(
                        orig_face_in_cut_(pi,face_idx_in_orig_face_in_cut_),
                        orig_face_in_cut_((pi+1)%orig_face_in_cut_.size(1),face_idx_in_orig_face_in_cut_) );
                  const pair<size_t,size_t> & one_edge = ea_->edges_[edge_idx];

                  assert(edge_idx != -1);
                  const pair<size_t,size_t> & edge2cell = ea_->edge2cell_[edge_idx];
                  if(ea_->is_boundary_edge(edge2cell)) continue;
                  assert(edge2cell.first == face_idx_in_orig_face_in_cut_ ||
                         edge2cell.second == face_idx_in_orig_face_in_cut_);
                  const size_t other_face_idx =
                      edge2cell.first + edge2cell.second -
                      face_idx_in_orig_face_in_cut_;
                  arounding_idx.push_back(other_face_idx);
                  const auto it =
                      updated_surface_type_cut.find(orig_face_idx_in_cut_[other_face_idx]);
                  assert(it != updated_surface_type_cut.end());

                  arounding_type.push_back(it->second);
                }

              if(arounding_type.size() <2) continue;
              count_type_number *= 0;
              for(size_t i = 0 ; i < arounding_type.size(); ++i){
                  ++count_type_number[arounding_type[i]];
                }
              const size_t max_arounding_type =
                  max_element(count_type_number.begin(), count_type_number.end())
                  - count_type_number.begin();
              if(count_type_number[max_arounding_type] == 1) continue; // three edges have three type
              if(max_arounding_type == my_type_it->second) continue;

              my_type_it->second = max_arounding_type;
              find_zigzag_edge = true;

              { // visualization
                boost::unordered_map<size_t,size_t> update_surface_type_uncut_i;
                convert_cut_surface_type2uncut_surface_type(updated_surface_type_cut, update_surface_type_uncut_i);

                static size_t count_i = 0;
                stringstream ss;
                ss << "surface_type_uncut_" << count_i++ << ".vtk";
                ofstream ofs_vtk(ss.str().c_str());
                static matrix<size_t> outside_face_uncut;
                static matrix<size_t> outside_face_idx_uncut;
                static matrix<size_t> outside_face_type_uncut;
                if(outside_face_uncut.size() == 0){
                    get_outside_face(*fa_uncut_, outside_face_uncut);
                  }
                if(outside_face_idx_uncut.size() == 0){
                    get_outside_face_idx(*fa_uncut_, outside_face_idx_uncut);
                  }
                if(outside_face_type_uncut.size() == 0)
                  outside_face_type_uncut.resize(outside_face_idx_uncut.size(),1);

                for(size_t fi = 0; fi < outside_face_idx_uncut.size(); ++fi){
                    outside_face_type_uncut[fi] = update_surface_type_uncut_i[outside_face_idx_uncut[fi]];
                  }
                tri2vtk(ofs_vtk, &uncut_node_[0], uncut_node_.size(2), &outside_face_uncut[0], outside_face_uncut.size(2));
                cell_data(ofs_vtk, &outside_face_type_uncut[0], outside_face_type_uncut.size(), "restricted_type");
              }
            }
          if(find_zigzag_edge == false) break;
        }
    }
  convert_cut_surface_type2uncut_surface_type(updated_surface_type_cut, updated_surface_type_uncut);

  {
    ofstream ofs("modified_surface_type_cut.vtk");
    tri2vtk(ofs, &merged_cut_node_[0], merged_cut_node_.size(2), &orig_face_in_cut_[0], orig_face_in_cut_.size(2));
    vector<size_t> surface_type_cut;
    for(size_t fi = 0 ; fi < orig_face_idx_in_cut_.size(); ++fi){
        const auto it = updated_surface_type_cut.find(orig_face_idx_in_cut_[fi]);
        surface_type_cut.push_back(it->second);
      }
    cell_data(ofs, &surface_type_cut[0], surface_type_cut.size(), "type");
  }
}

void surface_type_cleaner_cut::remove_degeneration_with_equation_graph(
    boost::unordered_map<size_t, size_t> &surface_type_uncut,
    vector<vector<size_t> > & node_groups) const
{
  throw std::invalid_argument("invalid function, not complete.");
  equation_graph eg(merged_cut_node_.size());
  vector<pair<size_t,size_t> > edges;
  get_edges_from_tet(merged_cut_tet_, edges);

  while(1){
      eg.assemble_equation_graph_generally(node_groups, edges);
      bool rtn = eg.check_valid(true);
    }
}

void surface_type_cleaner_cut::get_restricted_edges_from_one_piece(
    const vector<size_t> &one_patch, vector<pair<size_t,size_t> > & restricted_edges)const
{
  set<pair<size_t,size_t> > edges;
  for(size_t fi = 0; fi < one_patch.size(); ++fi){
      for(size_t pi = 0; pi < orig_face_in_cut_.size(1); ++pi){
          pair<size_t,size_t> one_edge(orig_face_in_cut_(pi,one_patch[fi]),
                                       orig_face_in_cut_((pi+1)%3, one_patch[fi]));

          const size_t edge_idx = ea_->get_edge_idx(one_edge.first, one_edge.second);
          assert(edge_idx != -1);
          const pair<size_t,size_t> & two_cells = ea_->edge2cell_[edge_idx];
          if(ea_->is_boundary_edge(two_cells)) continue;
          const auto type_it_0 = surface_type_cut_.find(orig_face_idx_in_cut_[two_cells.first]);
          const auto type_it_1 = surface_type_cut_.find(orig_face_idx_in_cut_[two_cells.second]);
          if(type_it_0 == surface_type_cut_.end() ||
             type_it_1 == surface_type_cut_.end()) {
              throw std::logic_error("can not find surface type ");
            }
          if(type_it_0->second != type_it_1->second){
              if(one_edge.first > one_edge.second)
                swap(one_edge.first, one_edge.second);
              edges.insert(one_edge);
            }
        }
    }

  restricted_edges.resize(edges.size());
  std::copy(edges.begin(), edges.end(), restricted_edges.begin());
}


void surface_type_cleaner_cut::remove_type_patch_degree_degeneration(
    matrix<size_t> & face_type,
    std::vector<vector<size_t> > &type_patch_vec,
    const matrix<size_t> & patch_faces,
    const jtf::mesh::edge2cell_adjacent & ea)
{
  // first build type patch connection degree
  map<size_t,size_t> face2groups;
  map<size_t, set<size_t> > p2p; // type_patch to type_patch
  {
    for(size_t pi = 0; pi < type_patch_vec.size(); ++pi){
        const vector<size_t> & one_type_patch = type_patch_vec[pi];
        for(const auto & face_idx : one_type_patch){
            face2groups[face_idx] = pi;
          }
      }
    face2groups[-1]=-1;
    for(size_t pi = 0; pi < type_patch_vec.size(); ++pi){
        const vector<size_t>& one_type_patch = type_patch_vec[pi];
        for(const auto & face_idx : one_type_patch){
            for(size_t pj = 0; pj < patch_faces.size(1); ++pj){
                const size_t edge_idx =
                    ea.get_edge_idx(patch_faces(pj,face_idx),
                                    patch_faces((pj+1)%patch_faces.size(1), face_idx));
                assert(edge_idx != -1);
                const pair<size_t,size_t> & face_pair = ea.edge2cell_[edge_idx];
                const size_t other_face_idx = face_pair.first + face_pair.second - face_idx;
                if(face2groups[other_face_idx] == pi) continue;
                p2p[pi].insert(face2groups[other_face_idx]);
              }
          }
      }
  }
  // for each type patch, if it's aroundded by one patches, merge it to that one
  {
    while(1){
        bool merged = false;
        for(auto it = p2p.begin(); it != p2p.end();){
            const auto & one_p2p = *it;
            if(one_p2p.second.size() == 1 && *(one_p2p.second.begin()) != -1){
                // this patch is arrounded by one patch, merge it.
                const size_t other_patch_idx = *(one_p2p.second.begin());
                const size_t merged_to_type = face_type[type_patch_vec[other_patch_idx].front()];
                type_patch_vec[other_patch_idx].insert(type_patch_vec[other_patch_idx].end(),
                                                       type_patch_vec[one_p2p.first].begin(),
                    type_patch_vec[one_p2p.first].end());

                for(const auto & face_idx : type_patch_vec[one_p2p.first]){
                    face_type[face_idx] = merged_to_type;
                  }
                { // remove one_p2p.first from other_patch
                  auto & other_patch = p2p[other_patch_idx];
                  auto it_this_patch  = other_patch.find(one_p2p.first);
                  assert(it_this_patch  != other_patch.end());
                  other_patch.erase(it_this_patch);
                }
                p2p.erase(it++);
                merged = true;
                type_patch_vec[one_p2p.first].clear();
                continue;
              }else if(one_p2p.second.size() == 2){
                const auto & to_other = one_p2p.second;
                if(to_other.find(-1) == to_other.end()){
                    size_t other_patch_idx = -1;

                    const size_t other_patch_idx_0 = *(one_p2p.second.begin());
                    auto other_patch_idx_it = one_p2p.second.begin();
                    ++other_patch_idx_it;
                    const size_t other_patch_idx_1 = *other_patch_idx_it;
                    if(type_patch_vec[other_patch_idx_0].size() > type_patch_vec[other_patch_idx_1].size())
                      other_patch_idx = other_patch_idx_0;
                    else
                      other_patch_idx = other_patch_idx_1;

                    const size_t merged_to_type = face_type[type_patch_vec[other_patch_idx].front()];
                    type_patch_vec[other_patch_idx].insert(type_patch_vec[other_patch_idx].end(),
                                                           type_patch_vec[one_p2p.first].begin(),
                        type_patch_vec[one_p2p.first].end());

                    for(const auto & face_idx : type_patch_vec[one_p2p.first]){
                        face_type[face_idx] = merged_to_type;
                      }
                    { // remove one_p2p.first from other_patch
                      auto & other_patch = p2p[other_patch_idx];
                      auto it_this_patch  = other_patch.find(one_p2p.first);
                      assert(it_this_patch  != other_patch.end());
                      other_patch.erase(it_this_patch);
                    }
                    p2p.erase(it++);
                    type_patch_vec[one_p2p.first].clear();
                    merged = true;
                    continue;
                  }
              }
            ++it;
          }
        if(merged == false) break;
      }

    std::vector<vector<size_t> > type_patch_vec_temp;
    for(size_t i = 0; i < type_patch_vec.size(); ++i){
        if(type_patch_vec[i].empty()) type_patch_vec_temp.push_back(type_patch_vec[i]);
      }
    type_patch_vec = type_patch_vec_temp;
  }
}


void surface_type_cleaner_cut::remove_degeneration_with_cut_node(
    boost::unordered_map<size_t,size_t> & updated_surface_type_uncut)
{
  //boost::unordered_map<size_t, size_t> updated_surface_type_cut = surface_type_cut_;
  boost::unordered_set<pair<size_t,size_t> > boundary_edges;
  vector<vector<size_t> > patches;
  matrix<size_t> boundary_edge_mat;
  jtf::mesh::get_boundary_edge(*ea_, boundary_edge_mat);
  for(size_t ei = 0; ei < boundary_edge_mat.size(2); ++ei){
      if(boundary_edge_mat(0,ei) > boundary_edge_mat(1, ei))
        boundary_edges.insert(make_pair(boundary_edge_mat(1,ei), boundary_edge_mat(0,ei)));
      else
        boundary_edges.insert(make_pair(boundary_edge_mat(0,ei), boundary_edge_mat(1,ei)));
    }

  get_face_patches_according_to_boundary(orig_face_in_cut_, *ea_, boundary_edges, patches);

  vector<size_t> arounding_idx, arounding_type;
  matrix<size_t> count_type_number(3,1);
  static int ppppi = 0;
  for(size_t pi = 0; pi < patches.size() ; ++pi){
      const vector<size_t> & one_patch = patches[pi]; // notice that  here idx is based on patches

      { //visualization
        matrix<size_t> patch_face(3, one_patch.size());
        matrix<size_t> patch_face_type(one_patch.size(),1);
        for(size_t fi = 0; fi < one_patch.size(); ++fi){
            patch_face(colon(),fi) = orig_face_in_cut_(colon(), one_patch[fi]);
            const auto it = surface_type_cut_.find(orig_face_idx_in_cut_[one_patch[fi]]);
            if(it == surface_type_cut_.end()){
                throw std::logic_error("strange can not find face type.");
              }
            patch_face_type[fi] = it->second;
          }
        stringstream ss;
        ss << "before_modified_patch_type_" << pi << ".vtk";
        ofstream ofs(ss.str().c_str());
        tri2vtk(ofs, &merged_cut_node_[0], merged_cut_node_.size(2),
            &patch_face[0], patch_face.size(2));
        cell_data(ofs, &patch_face_type[0], patch_face_type.size(), "type" );
      }

      matrix<size_t> patch_faces(3, one_patch.size());
      for(size_t fi = 0 ; fi < patch_faces.size(2); ++fi){
          patch_faces(colon(), fi) = orig_face_in_cut_(colon(), one_patch[fi]);
        }
      unique_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(patch_faces));
      if(!ea.get()){
          cerr << "# [error] strange can not build edge2cell_adjacent." << endl;
          return ;
        }
      jtf::mesh::one_ring_face_at_point orfap;
      orfap.add_all_faces(patch_faces, *ea);

      matrix<size_t> face_type(one_patch.size(),1);

      // handle turning points, and handle each non-straight restrcted edges
      while(1){
          bool remove_zigzag = false;
          bool remove_turing_point = false;
          /* while(1)*/{ // remove turing point
            vector<set<size_t> > type_patches;
            map<size_t,size_t> face2patch_idx;

            { // assemble type patch and face to it connection.
              vector<vector<size_t> > type_patch_vec;
              for(size_t fi = 0; fi < one_patch.size(); ++fi){
                  face_type[fi] = surface_type_cut_[orig_face_idx_in_cut_[one_patch[fi]]];
                }

              // here I can remove patch degeneration
              type_patches.clear();
              get_face_patches_according_to_type(patch_faces, face_type, orfap, type_patch_vec);
              type_patches.resize(type_patch_vec.size());

              for(size_t i = 0; i < type_patch_vec.size(); ++i){
                  const vector<size_t> & one_type_patch = type_patch_vec[i];
                  for(const auto & idx : one_type_patch){
                      face2patch_idx[one_patch[idx]] = i;
                      type_patches[i].insert(one_patch[idx]);
                    }
                }

              {
                vector<size_t> patch_group_idx(patch_faces.size(2));
                for(size_t tpi = 0; tpi < type_patch_vec.size(); ++tpi){
                    const vector<size_t> &one_type_patch = type_patch_vec[tpi];
                    for(const auto & face_idx : one_type_patch){
                        patch_group_idx[face_idx] = tpi;
                      }
                  }
                stringstream ss;
                ss << "before_surface_type_" << ppppi << ".vtk";
                ofstream ofs(ss.str().c_str());
                tri2vtk(ofs, &merged_cut_node_[0], merged_cut_node_.size(2), &patch_faces[0], patch_faces.size(2));
                cell_data(ofs, &face_type[0], face_type.size(), "type");
                vtk_data(ofs, &patch_group_idx[0], patch_group_idx.size(), "type_patch");

                stringstream ss_uncut;
                ss_uncut << "before_surface_type_uncut" << ppppi << ".vtk";
                ofstream ofs_uncut(ss_uncut.str().c_str());
                matrix<size_t> patch_faces_uncut = patch_faces;
                for(size_t i = 0; i < patch_faces_uncut.size(); ++i)
                  patch_faces_uncut[i] = cut_tet2tet_[patch_faces_uncut[i]];
                tri2vtk(ofs_uncut, &uncut_node_[0], uncut_node_.size(2), &patch_faces_uncut[0], patch_faces_uncut.size(2));
                cell_data(ofs_uncut, &face_type[0], face_type.size(), "type");
                vtk_data(ofs_uncut, &patch_group_idx[0], patch_group_idx.size(), "type_patch");
              }

              remove_type_patch_degree_degeneration(face_type, type_patch_vec, patch_faces, *ea);
              for(size_t fi = 0; fi < one_patch.size(); ++fi){
                  surface_type_cut_[orig_face_idx_in_cut_[one_patch[fi]]] = face_type[fi];
                }
              {
                stringstream ss;
                ss << "surface_type_" << ppppi << ".vtk";
                ofstream ofs(ss.str().c_str());
                tri2vtk(ofs, &merged_cut_node_[0], merged_cut_node_.size(2), &patch_faces[0], patch_faces.size(2));
                cell_data(ofs, &face_type[0], face_type.size(), "type");

                stringstream ss_uncut;
                ss_uncut << "surface_type_uncut" << ppppi << ".vtk";
                ofstream ofs_uncut(ss_uncut.str().c_str());
                matrix<size_t> patch_faces_uncut = patch_faces;
                for(size_t i = 0; i < patch_faces_uncut.size(); ++i)
                  patch_faces_uncut[i] = cut_tet2tet_[patch_faces_uncut[i]];
                tri2vtk(ofs_uncut, &uncut_node_[0], uncut_node_.size(2), &patch_faces_uncut[0], patch_faces_uncut.size(2));
                cell_data(ofs_uncut, &face_type[0], face_type.size(), "type");

                ++ppppi;
              }
            }

            //              vector<pair<size_t,size_t> > restricted_edges;
            //              get_restricted_edges_from_one_piece( one_patch, restricted_edges);

            //              if(restricted_edges.empty()) break;

            //              // these edges should be divided into different types
            //              vector<deque<pair<size_t,size_t> > > chains;
            //              {
            //                vector<deque<pair<size_t,size_t> > > one_chain;
            //                vector<vector<pair<size_t,size_t> > > uvw_edges(3);
            //                for(const auto & one_edge : restricted_edges){
            //                    const size_t edge_idx = ea->get_edge_idx(one_edge.first, one_edge.second);
            //                    assert(edge_idx != -1);
            //                    const pair<size_t,size_t> & two_cells = ea->edge2cell_[edge_idx];
            //                    assert(two_cells.first != -1 && two_cells.second != -1);
            //                    assert(face_type[two_cells.first] != face_type[two_cells.second]);
            //                    uvw_edges[3 - face_type[two_cells.first] - face_type[two_cells.second]].push_back(one_edge);
            //                  }
            //                for(size_t i = 0; i < uvw_edges.size(); ++i){
            //                    jtf::util::extract_chain_from_edges(uvw_edges[i], one_chain);
            //                  }
            //                chains.insert(chains.end(), one_chain.begin(), one_chain.end());
            //              }

            //              set<size_t> turning_points;
            //              map<size_t, pair<size_t,size_t> > turning_point2edge_ends;

            //              {
            //                stringstream ss;
            //                ss << "restricted_edges_" << ppppi << ".vtk";
            //                ofstream ofs(ss.str().c_str());
            //                vector<size_t> restricted_edges_vec;
            //                vector<size_t> edge_type;
            //                size_t edge_idx = 0;
            //                for(const auto & one_chain : chains){
            //                    for(const auto & one_edge : one_chain){
            //                        restricted_edges_vec.push_back(one_edge.first);
            //                        restricted_edges_vec.push_back(one_edge.second);
            //                        edge_type.push_back(edge_idx);
            //                      }
            //                    ++edge_idx;
            //                  }
            //                line2vtk(ofs, &merged_cut_node_[0],merged_cut_node_.size(2),
            //                    &restricted_edges_vec[0], restricted_edges_vec.size()/2);
            //                cell_data(ofs, &edge_type[0], edge_type.size(), "edge_idx");
            //              }
            //              for(size_t ci = 0; ci < chains.size(); ++ci){
            //                  const deque<pair<size_t,size_t> > & one_chain = chains[ci];

            //                  //                  const matrix<double> dir0 =
            //                  //                      merged_cut_node_(colon(), one_chain.front().second) - merged_cut_node_(colon(), one_chain.front().first);
            //                  //                  const matrix<double> dir1 =
            //                  //                      merged_cut_node_(colon(), one_chain.back().second) - merged_cut_node_(colon(), one_chain.back().first);
            //                  /*                  if(dot(dir0,dir1) < 0){ // if the chain is not straightful
            //                      //                      for(size_t ei = 0; ei < one_chain.size()-1; ++ei){
            //                      //                          turning_points.insert(one_chain[ei].second);
            //                      //                          turning_point2edge_ends[one_chain[ei].second] =
            //                      //                              make_pair(one_chain[ei].first, one_chain[ei+1].second);
            //                      //                        }
            //                      continue;
            //                    }else*/
            //                  { // if the chain is straightful, find turning points
            //                    for(size_t ei = 0; ei < one_chain.size()-1; ++ei){
            //                        const matrix<double> dir0 =
            //                            merged_cut_node_(colon(), one_chain[ei].second) - merged_cut_node_(colon(), one_chain[ei].first);
            //                        const matrix<double> dir1 =
            //                            merged_cut_node_(colon(), one_chain[ei+1].first) - merged_cut_node_(colon(), one_chain[ei+1].second);
            //                        if(dot(dir0,dir1) < 0) continue;
            //                        turning_points.insert(one_chain[ei].second);
            //                        turning_point2edge_ends[one_chain[ei].second] =
            //                            make_pair(one_chain[ei].first, one_chain[ei+1].second);
            //                      }
            //                  }
            //                }

            //              if(turning_points.empty()) break;

            //              //              remove_turing_point = true;
            //              //                        map<size_t, size_t> type2number;
            //              //                        update_type2number_one_piece(one_patch, type2number);

            //              for(const auto & idx : turning_points){
            //                  handle_each_turning_point_with_type_patch(
            //                        idx, turning_point2edge_ends[idx], type_patches, face2patch_idx);
            //                  //   handle_each_turning_point(idx, turning_point2edge_ends[idx], type2number);
            //                }


          }


          while(1){ // remove zigzag
              bool find_zigzag_edge = false;
              for(size_t fi = 0; fi < one_patch.size(); ++fi){
                  const size_t &face_idx_in_orig_face_in_cut_ = one_patch[fi];
                  auto my_type_it =
                      surface_type_cut_.find(orig_face_idx_in_cut_[face_idx_in_orig_face_in_cut_]);

                  arounding_idx.clear();
                  arounding_type.clear();

                  assert(my_type_it != surface_type_cut_.end());

                  for(size_t pi = 0; pi < orig_face_in_cut_.size(1); ++pi){
                      const size_t edge_idx = ea_->get_edge_idx(
                            orig_face_in_cut_(pi,face_idx_in_orig_face_in_cut_),
                            orig_face_in_cut_((pi+1)%orig_face_in_cut_.size(1),face_idx_in_orig_face_in_cut_) );
                      const pair<size_t,size_t> & one_edge = ea_->edges_[edge_idx];

                      assert(edge_idx != -1);
                      const pair<size_t,size_t> & edge2cell = ea_->edge2cell_[edge_idx];
                      if(ea_->is_boundary_edge(edge2cell)) continue;
                      assert(edge2cell.first == face_idx_in_orig_face_in_cut_ ||
                             edge2cell.second == face_idx_in_orig_face_in_cut_);
                      const size_t other_face_idx =
                          edge2cell.first + edge2cell.second -
                          face_idx_in_orig_face_in_cut_;
                      arounding_idx.push_back(other_face_idx);
                      const auto it =
                          surface_type_cut_.find(orig_face_idx_in_cut_[other_face_idx]);
                      assert(it != surface_type_cut_.end());

                      arounding_type.push_back(it->second);
                    }

                  if(arounding_type.size() <2) continue;
                  count_type_number *= 0;
                  for(size_t i = 0 ; i < arounding_type.size(); ++i){
                      ++count_type_number[arounding_type[i]];
                    }
                  const size_t max_arounding_type =
                      max_element(count_type_number.begin(), count_type_number.end())
                      - count_type_number.begin();
                  if(count_type_number[max_arounding_type] == 1) continue; // three edges have three type
                  if(max_arounding_type == my_type_it->second) continue;

                  my_type_it->second = max_arounding_type;
                  find_zigzag_edge = true;
                  remove_zigzag = true;
                }
              if(find_zigzag_edge == false) break;
            }

          { //visualization
            matrix<size_t> patch_face(3, one_patch.size());
            matrix<size_t> patch_face_type(one_patch.size(),1);
            for(size_t fi = 0; fi < one_patch.size(); ++fi){
                patch_face(colon(),fi) = orig_face_in_cut_(colon(), one_patch[fi]);
                const auto it = surface_type_cut_.find(orig_face_idx_in_cut_[one_patch[fi]]);
                if(it == surface_type_cut_.end()){
                    throw std::logic_error("strange can not find face type.");
                  }
                patch_face_type[fi] = it->second;
              }
            stringstream ss;
            ss << "after_modified_patch_type_" << pi << ".vtk";
            ofstream ofs(ss.str().c_str());
            tri2vtk(ofs, &merged_cut_node_[0], merged_cut_node_.size(2),
                &patch_face[0], patch_face.size(2));
            cell_data(ofs, &patch_face_type[0], patch_face_type.size(), "type" );
          }
          if(!remove_turing_point && !remove_zigzag) break;
        }
    }

  convert_cut_surface_type2uncut_surface_type(surface_type_cut_, updated_surface_type_uncut);

  {
    ofstream ofs("modified_surface_type_cut.vtk");
    tri2vtk(ofs, &merged_cut_node_[0], merged_cut_node_.size(2), &orig_face_in_cut_[0], orig_face_in_cut_.size(2));
    vector<size_t> surface_type_cut;
    for(size_t fi = 0 ; fi < orig_face_idx_in_cut_.size(); ++fi){
        const auto it = surface_type_cut_.find(orig_face_idx_in_cut_[fi]);
        surface_type_cut.push_back(it->second);
      }
    cell_data(ofs, &surface_type_cut[0], surface_type_cut.size(), "type");
  }

}

int remove_surface_degenerated_patch_in_cut(ptree &pt)
{
  pt_descriptor(pt);

  matrix<size_t> cut_tet;
  matrix<double> cut_node;

  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("input/cut_tet.value").c_str(),
                                          &cut_node, &cut_tet))
    return __LINE__;

  matrix<size_t> uncut_tet;
  matrix<double> uncut_node;
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("input/uncut_tet.value").c_str(),
                                          &uncut_node, &uncut_tet))
    return __LINE__;

  boost::unordered_map<size_t,size_t> surface_type;
  if(load_surface_restricted_type(pt.get<string>("input/surface_type_uncut.value").c_str(),
                                  surface_type))
    return __LINE__;

  vector<vector<size_t> > node_groups;
  if(zjucad::has("input/node_group.value",pt)){
      if(load_group_file(pt.get<string>("input/node_group.value").c_str(),
                         node_groups))
        return __LINE__;
    }
  surface_type_cleaner_cut stcc(cut_tet, cut_node, uncut_tet, uncut_node, surface_type,
                                node_groups);

  stcc.remove_degeneration_with_cut_node(surface_type);

  {
    ofstream ofs("modified_surface_type_uncut");
    for(const auto & one_face : surface_type){
        ofs << one_face.first << " " << one_face.second << endl;
      }

    ofstream ofs_vtk("modified_surface_type_uncut.vtk");
    matrix<size_t> outside_face_uncut, outside_face_idx_uncut;
    get_outside_face(*stcc.fa_uncut_, outside_face_uncut);
    get_outside_face_idx(*stcc.fa_uncut_, outside_face_idx_uncut);
    matrix<size_t> outside_face_type_uncut(outside_face_idx_uncut.size(),1);
    for(size_t fi = 0; fi < outside_face_idx_uncut.size(); ++fi){
        outside_face_type_uncut[fi] = surface_type[outside_face_idx_uncut[fi]];
      }
    tri2vtk(ofs_vtk, &uncut_node[0], uncut_node.size(2), &outside_face_uncut[0], outside_face_uncut.size(2));
    cell_data(ofs_vtk, &outside_face_type_uncut[0], outside_face_type_uncut.size(), "restricted_type");
  }

  return 0;
}

