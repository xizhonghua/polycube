#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/math/math.h>
#include <fstream>
#include "../common/vtk.h"
#include <jtflib/mesh/trimesh.h>
#include "../hex_frame_opt/angle_defect.h"
#include "../numeric/util.h"

using namespace std;
using namespace zjucad::matrix;

static int load_fv(const char * filename,
                    zjucad::matrix::matrix<double> &field)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      return __LINE__;
    }
  size_t num;
  ifs >> num;
  field.resize(6, num);
  for(size_t i = 0 ; i < num ; ++i){
      for(size_t j = 0; j < 6; ++j){
          ifs >> field(j,i);
        }
    }
  return 0;
}

double cal_dihedral_angle(const size_t fi,const size_t fj,
                          const zjucad::matrix::matrix<size_t> & mesh,
                          const zjucad::matrix::matrix<double> & node,
                          zjucad::matrix::matrix<double> &dir)
{
  matrix<double> ni,nj;
  jtf::mesh::cal_face_normal(mesh(colon(), fi), node, ni);
  jtf::mesh::cal_face_normal(mesh(colon(), fj), node, nj);

  double dihedral_angle = jtf::math::safe_acos(dot(ni,nj));
  // adjust angle sign
  size_t common_edge[2];
  jtf::mesh::find_common_edge(mesh(colon(), fi), mesh(colon(), fj), common_edge);
  auto it = std::find(mesh.begin(), mesh.end(), common_edge[0]);
  assert(it != mesh.end());
  const size_t i = it-mesh.begin();
  if(mesh[(i+1)%mesh.size(1)] != common_edge[1]){
      swap(common_edge[0],common_edge[1]);
    }

  dir = node(colon(), common_edge[1]) - node(colon(), common_edge[0]);

  const zjucad::matrix::matrix<double> cross_ni_nj = cross(ni,nj);
  if(norm(cross_ni_nj) > 1e-8) {// not parallel
      if(dot(dir, cross_ni_nj) < 0)
        dihedral_angle *= -1;
    }
  return dihedral_angle;
}

double get_angle_rot(const size_t fi, const size_t fj,  const jtf::mesh::tri_mesh &trim,
                     const zjucad::matrix::matrix<double> &field)
{
  itr_matrix<const double*> frame_i(3,2, &field(0,fi));
  itr_matrix<const double*> frame_j(3,2, &field(0,fj));

  matrix<double> dir;
  double angle = cal_dihedral_angle(fi,fj,trim.trimesh_.mesh_, trim.trimesh_.node_, dir);
  matrix<double> rot(3,3);
  from_angle_to_rotation_matrix(angle, dir, rot);

  matrix<double> new_frame_i = rot * frame_i;
  vector<double> dot_angle(4);
  for(size_t i = 0; i < 2; ++i){
      dot_angle[2*i+0] = dot(new_frame_i(colon(),0), frame_j(colon(),i));
      dot_angle[2*i+1] = dot(new_frame_i(colon(),0), -1*frame_j(colon(),i));
    }

  auto max_it = max_element(dot_angle.begin(), dot_angle.end());
  size_t max_idx = max_it-dot_angle.begin();

  angle = jtf::math::safe_acos(dot_angle[max_idx]);
  if(dot(trim.face_normal_(colon(),fj),
         cross(new_frame_i(colon(),0), (max_idx%2==0?1.0:-1.0)*frame_j(colon(), max_idx/2)))<0)
    angle *= -1;

  return angle;
}

double get_point_index(const size_t pi,
                       const jtf::mesh::tri_mesh &trim,
                       const zjucad::matrix::matrix<double> &field,
                       const jtf::mesh::one_ring_face_at_point &orfap,
                       const map<size_t, double> &angle_defect_map)
{
  const auto it = orfap.p2f_.find(pi);
  if(it == orfap.p2f_.end()){
      throw std::logic_error("input pi is not inside angle_defect_map");
    }

  const vector<size_t> & around_face = it->second;
  if(around_face.front()==-1) return 0; // vetex on boundary
  assert(around_face.front() == around_face.back());
  assert(around_face.front() != -1);
  double total_angle = 0;
  for(size_t i = 0; i < around_face.size()-1; ++i){
      total_angle += get_angle_rot(around_face[i], around_face[i+1], trim, field);
    }
  total_angle += angle_defect_map.at(pi);
  return  float_mod(total_angle, 2*My_PI())/(2*My_PI());
}

int extract_4sys_field_singularity(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] extract_4sys_field_singularity input_obj fv singularity_vtk" << endl;
      return __LINE__;
    }

  jtf::mesh::tri_mesh trim(argv[1]);

  matrix<double> field;
  load_fv(argv[2], field);

  jtf::mesh::one_ring_face_at_point orfap;
  orfap.add_all_faces(trim.trimesh_.mesh_,*trim.ea_);
  orfap.sort_int_loop_with_normal_info(trim.trimesh_.mesh_, trim.trimesh_.node_,
                                       *trim.ea_, trim.face_normal_);
  naive_angle_defect ad;
  map<size_t,double> angle_defect_map;
  ad.cal_angle_defect(orfap, trim.trimesh_.mesh_, trim.trimesh_.node_, angle_defect_map);

  set<size_t> face_used_nodes(trim.trimesh_.mesh_.begin(), trim.trimesh_.mesh_.end());
  matrix<double> index(face_used_nodes.size(),1);
  size_t ii = 0;
  for(const auto & one_idx : face_used_nodes){
      index[ii++] = get_point_index(one_idx, trim, field, orfap, angle_defect_map);
    }

  cerr << "# [info] total index " << std::accumulate(index.begin(),index.end(),0.0) << endl;
  cerr << "# [info] gauss bonent " << static_cast<int>(trim.trimesh_.node_.size(2)) + static_cast<int>(trim.trimesh_.mesh_.size(2))
          -static_cast<int>(trim.ea_->edges_.size()) << endl;
  ofstream ofs(argv[3]);
  tri2vtk(ofs, &trim.trimesh_.node_[0], trim.trimesh_.node_.size(2), &trim.trimesh_.mesh_[0], trim.trimesh_.mesh_.size(2));
  point_data(ofs, &index[0], index.size(), "index");

  {
    index *= 4;
    map<double, int> idx2point_num;
    int idx = 0;
    for(size_t i = 0; i < index.size(); ++i){
        idx = round_integer(index[i]);
        auto it = idx2point_num.find(idx*0.25);
        if(it == idx2point_num.end()){
            idx2point_num[idx*0.25] = 1;
          }else
          ++it->second;
      }
    for(const auto & one_idx : idx2point_num){
        if(std::fabs(one_idx.first) < 1e-8) continue;
        cerr << "# [info] singularity idx, number " << one_idx.first << " " << one_idx.second << endl;
      }
  }
  return 0;
}
