#include <iostream>
#include <fstream>
#include <numeric>
#include <sstream>
#include <omp.h>

#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <jtflib/mesh/mesh.h>

#include "visualize_tool.h"
#include "transition.h"
#include "vtk.h"
#include "util.h"

#include "LineToCylinder.h"
#include <jtflib/util/container_operation.h>
#include "../common/transition_type.h"

using namespace std;
using namespace zjucad::matrix;

static double color_table[30][4] =
{{255,0,0,255}, {255,80,80,255}, {255,160,160,255},// red for u axis
 {0,255,0,255}, {80,255,80,255}, {160,255,160,255},// green for v axis
 {0,0,255,255}, {80,80,255,255}, {160,160,255,255},// blue for w axis
 {0,0,0,255}, {255,0,0,200},{255,80,80,200},
 {255,160,160,200}, {0,255,0,200}, {80,255,80,200},
 {160,255,160,200}, {0,0,255,200}, {80,80,255,200},
 {160,160,255,200}, {0,0,0,200}, {255,0,0,100},
 {255,80,80,100}, {255,160,160,100}, {0,255,0,100},
 {80,255,80,100}, {160,255,160,100}, {0,0,255,100},
 {80,80,255,100}, {160,160,255,100},{255,255,255,255}};

int  dump_singularity_to_vtk(
    const char * file_name ,
    const matrixd &node,
    const vector<deque<pair<size_t,size_t> > > &singularity_edges)
{
  ofstream ofs(file_name);
  if(ofs.fail()){
      std::cerr << "# can not write file." << endl;
      return __LINE__;
    }

  vector<size_t> lines;
  vector<size_t> tag;
  for(size_t t = 0; t < singularity_edges.size(); ++t){
      for(size_t i = 0; i < singularity_edges[t].size(); ++i){
          lines.push_back(singularity_edges[t][i].first);
          lines.push_back(singularity_edges[t][i].second);
          tag.push_back(t);
        }
    }
  line2vtk(ofs,&node[0],node.size(2),&lines[0],lines.size()/2);
  cell_data(ofs, &tag[0], tag.size(), "chain_idx");
  return 0;
}


int dump_singularity_to_vtk(const char * vtk_file,
                            const matrixd &node,
                            const vector<deque<size_t> > &singularity_edges,
                            const vector<deque<size_t> > &singularity_type)
{

  cerr << "# singularity chain num: " << singularity_edges.size() << endl;

  vector<size_t> lines;
  vector<double> rgba_v;

  itr_matrix<const double*> rgba_m(4,30,&color_table[0][0]);
  //rgba_m /= 255.0;

  for(size_t t = 0 ;t < singularity_edges.size(); ++t){
      for(size_t i = 0; i < singularity_edges[t].size() ; ++i){
          lines.push_back(singularity_edges[t][i]);
        }
    }
  for(size_t t = 0; t < singularity_type.size(); ++t)
    for(size_t i = 0; i < singularity_type[t].size(); ++i)
      for(size_t j =0; j < 4; ++j)
        {
          if(singularity_type[t][i] != -1)
            rgba_v.push_back(rgba_m(j,singularity_type[t][i])/255.0);
          else
            rgba_v.push_back(rgba_m(j,29)/255.0);
        }
  //rgba_v.push_back(rgba_m(j,t));
  ofstream ofs(vtk_file);

  line2vtk(ofs,&node[0],node.size(2),&lines[0],lines.size()/2);
  cell_data_rgba(ofs,&rgba_v[0],rgba_v.size() / 4, "new_singularity");
}


int dump_singularity_to_vtk(const char * vtk_file,
                            const matrixd &node,
                            const map<pair<size_t,size_t>,size_t> &edge_info)
{

  vector<size_t> lines;
  vector<size_t> line_type;
  vector<double> rgba_v;

  itr_matrix<const double*> rgba_m(4,30,&color_table[0][0]);
  //rgba_m /= 255.0;

  for(map<pair<size_t,size_t>,size_t>::const_iterator mci = edge_info.begin();
      mci != edge_info.end(); ++mci){
      lines.push_back(mci->first.first);
      lines.push_back(mci->first.second);
      line_type.push_back((mci->second == -1)?29:mci->second);
      for(size_t j = 0; j < 4; ++j)
        {
          if(mci->second != -1)
            if(mci->second > 9)
              rgba_v.push_back(rgba_m(j,9)/255.0);
            else
              rgba_v.push_back(rgba_m(j,mci->second)/255.0);
          else
            rgba_v.push_back(rgba_m(j,29)/255.0);
        }
    }

  ofstream ofs(vtk_file);

  line2vtk(ofs,&node[0],node.size(2),&lines[0],lines.size()/2);
  ofs << "CELL_DATA " << rgba_v.size() / 4 << "\n";
  vtk_data_rgba(ofs, &rgba_v[0], rgba_v.size() / 4, "new_singularity_after_jump", "my_table0");
  vtk_data(ofs,&line_type[0],line_type.size(), "singularity_type","my_table1");
}


int dump_singularity_to_vtk(
    const char * vtk_file,
    const matrixd &node,
    const std::vector<std::deque<pair<size_t,size_t > > > &singularity_edges,
    const std::vector<size_t> &singularity_type)
{
  cerr << "# singularity chain num: " << singularity_edges.size() << endl;

  vector<size_t> lines;
  vector<double> rgba_v;

  itr_matrix<const double*> rgba_m(4,30,&color_table[0][0]);
  //rgba_m /= 255.0;

  for(size_t t = 0 ;t < singularity_edges.size(); ++t){
      for(size_t i = 0; i < singularity_edges[t].size() ; ++i){
          lines.push_back(singularity_edges[t][i].first);
          lines.push_back(singularity_edges[t][i].second);
        }
    }
  for(size_t t = 0; t < singularity_type.size(); ++t)
    //for(size_t i = 0; i < singularity_type[t].size(); ++i)
    for(size_t j =0; j < 4; ++j)
      {
        if(singularity_type[t] != -1)
          if(singularity_type[t] > 9)
            rgba_v.push_back(rgba_m(j,9)/255.0);
          else
            rgba_v.push_back(rgba_m(j,singularity_type[t])/255.0);
        else
          rgba_v.push_back(rgba_m(j,29)/255.0);
      }
  //rgba_v.push_back(rgba_m(j,t));
  ofstream ofs(vtk_file);

  line2vtk(ofs,&node[0],node.size(2),&lines[0],lines.size()/2);
  cell_data_rgba(ofs,&rgba_v[0],rgba_v.size() / 4, "new_singularity");
  return 0;
}


int dump_singularity_to_cylinder(
    const char * obj_file,
    const matrixd &node,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    const double radious_per_bounding_sphere_r)
{
  cerr << "# singularity chain num: " << singularity_edges.size() << endl;

  const double d = calc_bounding_sphere_size(node);

  const double radious = d * radious_per_bounding_sphere_r;

  vector<vector<line_type> >  line_segments;

  vector<line_type> group;
  for(size_t chi = 0; chi < singularity_edges.size(); ++chi){
      const deque<pair<size_t,size_t> > & one_chain = singularity_edges[chi];
      line_type one_chain_points;

      for(size_t pi = 0; pi < one_chain.size(); ++pi){
          one_chain_points.push_back(node(colon(), one_chain[pi].first));
        }
      one_chain_points.push_back(node(colon(), one_chain.back().second));
      group.push_back(one_chain_points);
    }

  line_segments.push_back(group);

  vector<vector<matrixd > > point_coord;
  vector<vector<matrix<matrix<int> > > > face_index;

  convertToCylinder(line_segments, point_coord, face_index, radious);
  saveCylinder(obj_file, point_coord, face_index);
  return 0;
}

int dump_singularity_chain_to_vtk_2(
    const char * vtk_file,
    const matrixd &node,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    const std::vector<std::deque<size_t> > &singularity_type)
{
  cerr << "# singularity chain num: " << singularity_edges.size() << endl;

  vector<size_t> lines;
  vector<double> rgba_v;

  itr_matrix<const double*> rgba_m(4,30,&color_table[0][0]);
  //rgba_m /= 255.0;

  for(size_t t = 0 ;t < singularity_edges.size(); ++t){
      for(size_t i = 0; i < singularity_edges[t].size() ; ++i){
          lines.push_back(singularity_edges[t][i].first);
          lines.push_back(singularity_edges[t][i].second);

          for(size_t j =0; j < 4; ++j)
            {
              if(singularity_type[t][i] != -1)
                if(singularity_type[t][i] > 9)
                  rgba_v.push_back(rgba_m(j,9)/255.0);
                else
                  rgba_v.push_back(rgba_m(j,singularity_type[t][i])/255.0);
              else
                rgba_v.push_back(rgba_m(j,29)/255.0);
            }
        }
    }

  //rgba_v.push_back(rgba_m(j,t));
  ofstream ofs(vtk_file);

  line2vtk(ofs,&node[0],node.size(2),&lines[0],lines.size()/2);
  cell_data_rgba(ofs,&rgba_v[0],rgba_v.size() / 4, "new_singularity");
  return 0;
}

int dump_singularity_chain_to_vtk_2_type(
    const char * vtk_file,
    const matrixd &node,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    const std::vector<std::deque<size_t> > &singularity_type)
{
  vector<size_t> lines;
  vector<size_t> type;

  for(size_t t = 0 ;t < singularity_edges.size(); ++t){
      for(size_t i = 0; i < singularity_edges[t].size() ; ++i){
          lines.push_back(singularity_edges[t][i].first);
          lines.push_back(singularity_edges[t][i].second);
          type.push_back(singularity_type[t][i]);
        }
    }

  //rgba_v.push_back(rgba_m(j,t));
  ofstream ofs(vtk_file);

  line2vtk(ofs,&node[0],node.size(2),&lines[0],lines.size()/2);
  cell_data(ofs,&type[0],type.size(), "free_axis");
  return 0;
}

int dump_singularity_chain_to_vtk_3(
    const char * vtk_file,
    const matrixd &node,
    const std::vector<deque<std::tuple<size_t,size_t,size_t> > > &singularity_edges)
{

  vector<size_t> lines;
  vector<size_t> type;

  for(size_t t = 0 ;t < singularity_edges.size(); ++t){
      for(size_t i = 0; i < singularity_edges[t].size() ; ++i){
          lines.push_back(get<0>(singularity_edges[t][i]));
          lines.push_back(get<1>(singularity_edges[t][i]));
          type.push_back(get<2>(singularity_edges[t][i]));
        }
    }

  //rgba_v.push_back(rgba_m(j,t));
  ofstream ofs(vtk_file);

  line2vtk(ofs,&node[0],node.size(2),&lines[0],lines.size()/2);
  cell_data(ofs,&type[0],type.size(), "free_axis");
  return 0;
}

int dump_singularity_to_vtk(
    const char *vtk_file,
    const matrixd & node,
    const set<std::tuple<size_t,size_t,size_t> > & edges)
{
  vector<size_t> lines;
  vector<size_t> type;

  for(const auto & one_edge:edges){
      lines.push_back(get<0>(one_edge));
      lines.push_back(get<1>(one_edge));
      type.push_back(get<2>(one_edge));
    }

  ofstream ofs(vtk_file);

  line2vtk(ofs,&node[0],node.size(2),&lines[0],lines.size()/2);
  cell_data(ofs,&type[0],type.size(), "free_axis");
  return 0;
}

int dump_uvw(const char* vtk_file,
             const matrixst &tet,
             const matrixd &node,
             const matrixd &uvw)
{
  const matrixd & u = uvw(0,colon());
  const matrixd & v = uvw(1,colon());
  const matrixd & w = uvw(2,colon());

  matrixd separate_node(3, tet.size(2) * 4);
  for(size_t t = 0; t < tet.size(2); ++t){
      for(size_t i = 0; i < 4; ++i)
        separate_node(colon(),t*4 + i) = node(colon(),tet(i,t));
    }

  vector<size_t> idx(separate_node.size(2));
  for(size_t t = 0; t < idx.size(); ++t) idx[t] = t;

  string str = vtk_file;
  str += "_uvw.vtk";
  ofstream ofs(str.c_str());
  point2vtk(ofs,&separate_node[0],separate_node.size(2),&idx[0],idx.size());
  ofs << "CELL_DATA " << separate_node.size(2) << "\n";
  vtk_data(ofs, &u[0], u.size(), "u", "my_table");
  vtk_data(ofs, &v[0], v.size(), "v", "my_table");
  vtk_data(ofs, &w[0], w.size(), "w", "my_table");

  return 0;
}

int dump_uvw_wave(const char * vtk_file,
                  const matrixst &tet,
                  const matrixd &node,
                  const matrixd &uvw,
                  const double freq)
{
  matrixd  wave_value = zeros<double>(tet.size(2),1);
  matrixd  u_wave_value = zeros<double>(tet.size(2),1);
  matrixd  v_wave_value = zeros<double>(tet.size(2),1);
  matrixd  w_wave_value = zeros<double>(tet.size(2),1);
  matrixd inner_uvw = zeros<double>(3,1);
  size_t t = 0;
  const double pow_n = 1;
  const double wave_f = freq * (1.0/pow_n);
#pragma omp parallel for private(t,inner_uvw)
  for(t= 0; t < tet.size(2); ++t){
      inner_uvw = zeros<double>(3,1);
      for(size_t i = 0; i < 4; ++i) {
          inner_uvw += uvw(colon(),t*4 + i);
        }
      inner_uvw /= 4.0;
      u_wave_value[t] = pow(cos(wave_f* inner_uvw[0] * My_PI()),pow_n);
      v_wave_value[t] = pow(cos(wave_f* inner_uvw[1] * My_PI()),pow_n);
      w_wave_value[t] = pow(cos(wave_f* inner_uvw[2] * My_PI()),pow_n);
      wave_value[t] = pow(cos(wave_f* inner_uvw[0] * My_PI()) * cos(wave_f*inner_uvw[1]* My_PI()) * cos(wave_f*inner_uvw[2]* My_PI()),pow_n);
    }


  string str = vtk_file;
  str += "_wave.vtk";
  ofstream ofs(str.c_str());
  tet2vtk(ofs,&node[0],node.size(2),&tet[0],tet.size(2));
  ofs << "CELL_DATA " << tet.size(2) << "\n";
  vtk_data(ofs, &wave_value[0], wave_value.size(), "wave", "my_table0");
  vtk_data(ofs, &u_wave_value[0], u_wave_value.size(), "u", "my_table");
  vtk_data(ofs, &v_wave_value[0], v_wave_value.size(), "v", "my_table");
  vtk_data(ofs, &w_wave_value[0], w_wave_value.size(), "w", "my_table");

  return 0;
}

int dump_uvw_vol(const char * vtk_file,
                 const matrixst &tet_,
                 const matrixd &node_,
                 const matrixd &uvw)
{
  matrixd u(tet_.size(2),1);
  matrixd v(tet_.size(2),1);
  matrixd w(tet_.size(2),1);

  size_t t = 0;
#pragma omp parallel for private(t)
  for(t = 0; t < tet_.size(2); ++t){
      u[t] = (uvw(0,t * 4) + uvw(0,t * 4 + 1) + uvw(0,t * 4 + 2) + uvw(0,t * 4 + 3))/4;
      v[t] = (uvw(1,t * 4) + uvw(1,t * 4 + 1) + uvw(1,t * 4 + 2) + uvw(1,t * 4 + 3))/4;
      w[t] = (uvw(2,t * 4) + uvw(2,t * 4 + 1) + uvw(2,t * 4 + 2) + uvw(2,t * 4 + 3))/4;
    }

  string str = vtk_file;
  str += "_uvw_vol.vtk";
  ofstream ofs(str.c_str());
  tet2vtk(ofs,&node_[0],node_.size(2),&tet_[0],tet_.size(2));

  ofs << "CELL_DATA " << tet_.size(2) << "\n";
  vtk_data(ofs, &u[0], u.size(), "u", "my_table");
  vtk_data(ofs, &v[0], v.size(), "v", "my_table");
  vtk_data(ofs, &w[0], w.size(), "w", "my_table");
  return 0;
}

int dump_frame_align_error(const char * vtk_file,
                           const matrix<matrixd > &frame_in_tet,
                           const matrixst &tet,
                           const matrixd &node,
                           const matrixd &uvw)
{
  matrixd P = ones<double>(4, 4);
  matrixd err = zeros<double>(tet.size(2),1);
  matrixd M_(4,3);
  //matrixd f = zeros<double>(3,4);
  for(size_t t = 0; t < tet.size(2); ++t){
      P = ones<double>(4, 4);
      P(colon(0, 2), colon(0, 3)) = node(colon(),tet(colon(),t));
      if(inv(P))
        std::cerr << "inverse fail." << std::endl;
      M_ = P(colon(), colon(0, 2));
      matrixd f0 =  uvw(colon(),colon(t * 4,t * 4 + 3));
      assert(f0.size(2) == M_.size(1));
      err[t] = norm(f0 * M_ - trans(frame_in_tet[t]));
    }

  string str = vtk_file;
  str += "_.vtk";
  ofstream ofs(str.c_str());
  //tet2vtk(ofs,&separate_node[0],separate_node.size(2),&idx[0],idx.size());
  tet2vtk(ofs,&node[0],node.size(2),&tet[0],tet.size(2));
  ofs << "CELL_DATA " << tet.size(2) << "\n";
  vtk_data(ofs, &err[0], err.size(), "frame_align_error", "my_table");

  return 0;
}

int dump_sh_difference_to_vtk(
    const char * vtk_file,
    const matrixst & tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const matrixd& sh)
{
  ofstream ofs(vtk_file);
  if(ofs.fail()){
      cerr << "# [error] can not open vtk file." << endl;
      return __LINE__;
    }
  vector<double> residual(tet.size(2));

  vector<vector<size_t> > tet_adjacent_tet;
  tet_adjacent_tet.resize(tet.size(2));
  const size_t face_num = fa.faces_.size();
  for(size_t fi = 0; fi < face_num; ++fi) {
      const size_t &vi = fa.face2tet_[fi].first;
      const size_t &vj = fa.face2tet_[fi].second;
      if(vi != -1 && vj != -1){
          tet_adjacent_tet[vi] << vj;
          tet_adjacent_tet[vj] << vi;
        }
    }

  for(size_t vi = 0; vi < tet_adjacent_tet.size(); ++vi) {
      double residual_each_tet = 0.0;
      for(size_t advi = 0; advi < tet_adjacent_tet[vi].size(); ++advi) {
          residual_each_tet += norm(sh(colon(), vi) -
                                    sh(colon(),tet_adjacent_tet[vi][advi]));
        }
      if(tet_adjacent_tet[vi].size())
        residual[vi] = residual_each_tet / tet_adjacent_tet[vi].size();
      else{
          std::cerr << "# [error] tet " << vi
                    << " is not beside with anyother tets." << std::endl;
          return __LINE__;
        }
    }
  tet2vtk(ofs, &node[0],node.size(2), &tet[0], tet.size(2));
  cell_data(ofs, &residual[0], residual.size(), "residual");
  return 0;
}

int dump_frame_difference_to_vtk(
    const char * vtk_file,
    const matrixst & tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const zjucad::matrix::matrix<matrixd >& frame)
{
  ofstream ofs(vtk_file);
  if(ofs.fail()){
      cerr << "# [error] can not open vtk file." << endl;
      return __LINE__;
    }
  vector<double> residual(tet.size(2));

  vector<vector<size_t> > tet_adjacent_tet;
  tet_adjacent_tet.resize(tet.size(2));
  const size_t face_num = fa.faces_.size();
  for(size_t fi = 0; fi < face_num; ++fi) {
      const size_t &vi = fa.face2tet_[fi].first;
      const size_t &vj = fa.face2tet_[fi].second;
      if(vi != -1 && vj != -1){
          tet_adjacent_tet[vi] << vj;
          tet_adjacent_tet[vj] << vi;
        }
    }

  for(size_t vi = 0; vi < tet_adjacent_tet.size(); ++vi) {
      double residual_each_tet = 0.0;
      for(size_t advi = 0; advi < tet_adjacent_tet[vi].size(); ++advi) {
          residual_each_tet += norm(frame[vi] - frame[tet_adjacent_tet[vi][advi]]);
        }
      if(tet_adjacent_tet[vi].size())
        residual[vi] = residual_each_tet / tet_adjacent_tet[vi].size();
      else{
          std::cerr << "# [error] tet " << vi
                    << " is not beside with anyother tets." << std::endl;
          return __LINE__;
        }
    }
  tet2vtk(ofs, &node[0],node.size(2), &tet[0], tet.size(2));
  cell_data(ofs, &residual[0], residual.size(), "residual");
  return 0;
}


int dump_out_singularity_chain(
    const char * file_name,
    const vector<deque<pair<size_t,size_t> > > &singularities_chain,
    const vector<size_t> & singularities_type
    )
{
  ofstream ofs(file_name);
  if(ofs.fail())
    {
      cerr << "# error: can not open " << file_name << endl;
      return __LINE__;
    }

  ofs << singularities_chain.size() << endl;
  for(size_t t = 0, j = 0; t < singularities_chain.size(); ++t){
      ofs << singularities_chain[t].size() << endl;
      for(size_t i = 0; i < singularities_chain[t].size(); ++i){
          ofs << singularities_chain[t][i].first << " "
              << singularities_chain[t][i].second << " ";
        }
      ofs << endl;
      for(size_t i = 0;i < singularities_chain[t].size(); ++i, ++j){
          ofs << singularities_type[j] << " ";
        }
      ofs << endl;
    }
  return 0;
}

int dump_out_singularity_chain_2(
    const char * file_name,
    const vector<deque<size_t> > &singularities_chain,
    const vector<deque<size_t> > & singularities_type
    )
{
  cerr << "# begin to dump out the singularity edge. " << file_name << endl;
  cerr << "# singularity chain num " << singularities_chain.size() << endl;
  ofstream ofs(file_name);
  if(ofs.fail())
    {
      cerr << "# error: can not open " << file_name << endl;
      return __LINE__;
    }
  assert(singularities_chain.size() == singularities_type.size());
  ofs << singularities_chain.size() << endl;
  for(size_t t = 0; t < singularities_chain.size(); ++t){
      ofs << singularities_chain[t].size() << endl;
      for(size_t i = 0; i < singularities_chain[t].size(); ++i){
          ofs << singularities_chain[t][i]<< " ";
        }
      ofs << endl;
      for(size_t i = 0;i < singularities_type[t].size(); ++i){
          ofs << singularities_type[t][i] << " ";
        }
      ofs << endl;
    }
  return 0;
}

int dump_out_singularity_chain_3(const char * file_name,
                                 const vector<deque<pair<size_t,size_t > > > &singularities_chain,
                                 const vector<deque<size_t> > & singularities_type)
{
  cerr << "# begin to dump out the singularity edge. " << file_name << endl;
  cerr << "# singularity chain num " << singularities_chain.size() << endl;
  ofstream ofs(file_name);
  if(ofs.fail())
    {
      cerr << "# error: can not open " << file_name << endl;
      return __LINE__;
    }
  assert(singularities_chain.size() == singularities_type.size());
  ofs << singularities_chain.size() << endl;
  for(size_t t = 0; t < singularities_chain.size(); ++t){
      ofs << singularities_chain[t].size() << endl;
      for(size_t i = 0; i < singularities_chain[t].size(); ++i){
          ofs << singularities_chain[t][i].first << " ";
          ofs << singularities_chain[t][i].second << " ";
        }
      ofs << endl;
      for(size_t i = 0;i < singularities_type[t].size(); ++i){
          ofs << singularities_type[t][i] << " ";
        }
      ofs << endl;
    }
  return 0;
}

int dump_out_singularity_chain_with_tet_ends(
    const char * file_name,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const std::vector<std::deque<std::pair<size_t,size_t > > > &singularities_chain,
    const std::vector<std::deque<size_t> > & singularities_type)
{
  cerr << "# begin to dump out the singularity edge. " << file_name << endl;
  cerr << "# singularity chain num " << singularities_chain.size() << endl;
  ofstream ofs(file_name);
  if(ofs.fail()) {
      cerr << "# error: can not open " << file_name << endl;
      return __LINE__;
    }

  assert(singularities_chain.size() == singularities_type.size());
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  ofs << singularities_chain.size() << endl;
  for(size_t t = 0; t < singularities_chain.size(); ++t){
      const deque<pair<size_t,size_t> > & one_chain = singularities_chain[t];
      ofs << one_chain.size() << endl;
      for(size_t i = 0; i < one_chain.size(); ++i){
          ofs << one_chain[i].first << " ";
          ofs << one_chain[i].second << " ";
        }
      ofs << endl;
      for(size_t i = 0;i < singularities_type[t].size(); ++i){
          ofs << singularities_type[t][i] << " ";
        }
      ofs << endl;
      for(size_t i = 0; i < one_chain.size(); ++i){
          oecit it = ortae.e2t_.find(one_chain[i]);
          if(it == ortae.e2t_.end())
            it = ortae.e2t_.find(make_pair(one_chain[i].second,one_chain[i].first));
          if(it == ortae.e2t_.end()) {
              cerr << "# [error] can not find edge <" << one_chain[i].first
                   << "," << one_chain[i].second << ">." << endl;
              return __LINE__;
            }
          const vector<size_t> & around_tets = it->second;
          if(find(around_tets.begin(), around_tets.end(),-1) != around_tets.end()){
              cerr << "# [error] strange inner edge meets outside tet." << endl;
              return __LINE__;
            }
          assert(around_tets.size() > 2);
          ofs << around_tets[0] << " " << around_tets[1] << " ";
        }
      ofs << endl;
    }
  return 0;
}

int dump_out_singularity_axis(const char *file_name,
                              const matrixd &node,
                              const vector<deque<pair<size_t,size_t> > > &singularities_chain,
                              const vector<vector<size_t> > &singularities_loop,
                              const vector<size_t> & singularities_type,
                              const matrix<matrixd > & frame_inner)
{
  ofstream ofs(file_name);
  if(ofs.fail())
    {
      cerr << "# error: can not open " << file_name << endl;
      return __LINE__;
    }

  ofs << singularities_chain.size() << endl;
  for(size_t t = 0,j = 0; t < singularities_chain.size(); ++t){
      ofs << singularities_chain[t].size() << endl;
      for(size_t i = 0; i < singularities_chain[t].size(); ++i,++j){

          ofs << singularities_chain[t][i].first <<  " "
              << singularities_chain[t][i].second << " ";

          if(is_black_line_new(singularities_type[j])){
              ofs << " " << 0 << " " << 0 << " " << 0 << endl;
              continue;
            }

          const size_t col = singularities_type[j]/3;

          matrixd axis = zeros<double>(3,1);
          for(size_t k = 0; k < singularities_loop[j].size(); ++k){
              const size_t &tet_idx = singularities_loop[j][k];
              axis += frame_inner[tet_idx](colon(),col);
            }
          axis /= singularities_loop[j].size();

          matrixd edge = node(colon(),singularities_chain[t][i].second) - node(colon(),singularities_chain[t][i].first);
          if(dot(edge,axis) < -1e-8)
            axis *= -1;
          ofs << axis[0] << " " << axis[1] << " " << axis[2] << endl; // output the axis(x,y,z)
        }
    }
  return 0;
}

int dump_out_faces_in_cut_tet(const char * file,
                              const std::vector<size_t> faces,
                              const jtf::mesh::face2tet_adjacent &fa,
                              const matrixst &cut_tet2tet,
                              const matrixd &node)
{
  ofstream ofs(file);
  if(ofs.fail()){
      cerr << "# [error] can not open file." << endl;
      return __LINE__;
    }
  matrixst out_faces(3,faces.size());
  for(size_t t = 0; t < faces.size(); ++t){
      const vector<size_t> & face_v = fa.faces_[faces[t]];
      for(size_t i = 0; i < 3; ++i){
          out_faces(i,t) = cut_tet2tet[face_v[i]];
        }
    }
  tri2vtk(ofs,&node[0],node.size(2),&out_faces[0],out_faces.size(2));
  return 0;
}

int load_singularity_chain(const char * file_name,
                           std::vector<std::deque<std::pair<size_t,size_t> > > &singularities_chain,
                           std::vector<size_t> & singularities_type)
{
  ifstream ifs(file_name);
  if(ifs.fail())
    {
      cerr << "# error: can not open " << file_name << endl;
      return __LINE__;
    }

  size_t chain_num = 0;
  size_t edge_num = 0;
  size_t type = -1;
  ifs >> chain_num;
  pair<size_t,size_t> edge;
  singularities_chain.reserve(chain_num);
  for(size_t t = 0; t < chain_num; ++t){
      ifs >> edge_num;
      deque<pair<size_t,size_t> > chain;

      for(size_t i = 0; i < edge_num; ++i){
          ifs >> edge.first >> edge.second;
          chain.push_back(edge);
        }
      singularities_chain.push_back(chain);

      for(size_t i = 0;i < edge_num; ++i){
          ifs >> type;
          singularities_type.push_back(type);
        }
    }
  return 0;
}

int load_singularity_chain_new(const char * file_name,
                               std::vector<std::deque<std::pair<size_t,size_t> > > &singularities_chain,
                               std::vector<deque<size_t> > & singularities_type)
{
  ifstream ifs(file_name);
  if(ifs.fail())
    {
      cerr << "# error: can not open " << file_name << endl;
      return __LINE__;
    }

  size_t chain_num = 0;
  size_t edge_num = 0;
  size_t type = -1;
  ifs >> chain_num;
  pair<size_t,size_t> edge;
  singularities_chain.reserve(chain_num);
  singularities_type.reserve(chain_num);
  for(size_t t = 0; t < chain_num; ++t){
      ifs >> edge_num;
      deque<pair<size_t,size_t> > chain;
      deque<size_t> chain_type;
      for(size_t i = 0; i < edge_num; ++i){
          ifs >> edge.first >> edge.second;
          chain.push_back(edge);
        }
      singularities_chain.push_back(chain);

      for(size_t i = 0;i < edge_num; ++i){
          ifs >> type;
          chain_type.push_back(type);
        }
      singularities_type.push_back(chain_type);
    }
  return 0;
}

std::string bitset_to_hex_string(const boost::dynamic_bitset<> & bit)
{
  boost::dynamic_bitset<> bit_ = bit;
  bit_.resize(bit_.size()/4 * 4);
  string hex;
  for(size_t n = 0; n < bit_.size()/4; ++n){
      int a = bit_[n*4] * 8 + bit_[n*4+1] * 4 + bit_[n*4+2] * 2 + bit_[n*4+3];
      if(a < 10)
        hex.push_back('0' + a);
      else
        hex.push_back('A' + a - 10);
    }
  return hex;
}
