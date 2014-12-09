#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <jtflib/mesh/io.h>

#include "io.h"
#include "topology_analysis.h"
#include "../numeric/util.h"
#include <jtflib/algorithm/gauss_elimination.h>
#include "../common/transition_type.h"
#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"
#include "../equation_graph/equation_graph.h"

using namespace std;
using namespace jtf::algorithm;

int dump_out_group_file(
    const char * group_file,
    const matrixst & fnode_idx,
    const singularity_graph & sg,
    const boost::unordered_set<size_t> & restricted_nodes,
    boost::unordered_map<size_t,size_t> * fnode2group_ptr = 0)
{
  if(!group_file) {
      cerr << "# [error] can not open group_file." << endl;
      return __LINE__;
    }

  ofstream ofs(group_file);

  for(size_t t = 0,gi = 0; t < sg.groups_.size(); ++t){
      if(!sg.groups_[t].empty()){
          if(sg.groups_[t].size() == 1) continue;
          boost::unordered_set<size_t>::const_iterator bucit =
              restricted_nodes.find(*(sg.groups_[t].begin()));
          if(bucit == restricted_nodes.end())
            ofs << "g " << gi << " " << sg.groups_[t].size() << " " <<  0 << endl;
          else
            ofs << "g " << gi << " " << sg.groups_[t].size() << " " << 1 << endl;

          for(boost::unordered_set<size_t>::const_iterator buscit =
              sg.groups_[t].begin(); buscit != sg.groups_[t].end(); ++buscit){
              ofs << fnode_idx[*buscit] << " ";
              if(fnode2group_ptr)
                (*fnode2group_ptr)[fnode_idx[*buscit]] = gi;
            }
          ofs << endl;
          ++gi;
        }
    }

  return 0;
}

//static int add_g_equation_to_f_equation(
//    const singularity_graph & sg,
//    const vector<size_t> & rot_type,
//    const jtf::mesh::face2tet_adjacent & fa_cut,
//    const matrixst & cut_tet2tet,
//    const boost::unordered_map<pair<size_t,size_t>,size_t> & jump_face_pair2_rot_idx,
//    const jtf::algorithm::equation<double> &g_eqn,
//    vector<jtf::algorithm::equation<double> > &f_eqns)
//{
//  f_eqns.clear();
//  f_eqns.resize(3);
//  vector<size_t> face_from(3) , face_to(3);
//  matrixd face_from_p(3,1), face_to_p(3,1);

//  for(jtf::algorithm::equation<double>::eq_const_iterator cit = g_eqn.begin();
//      cit != g_eqn.end(); ++cit){
//    const jtf::algorithm::expression<double> & exp = *cit;
//    const pair<size_t,size_t> face_pair = sg.get_jump_face_from_gnode_idx(exp.index);
//    boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit
//        = jump_face_pair2_rot_idx.find(face_pair);

//    const vector<size_t> & face_from_vec = fa_cut.faces_[face_pair.first];
//    const vector<size_t> & face_to_vec = fa_cut.faces_[face_pair.second];

//    for(size_t i = 0; i < 3; ++i){
//      for(size_t j = 0; j < 3; ++j){
//        if(cut_tet2tet[face_from_vec[i]] == cut_tet2tet[face_to_vec[j]]){
//          face_from[i] = face_from_vec[i];
//          face_to[i] = face_to_vec[j];
//          break;
//        }
//      }
//    }

//    for(size_t j = 0; j < 3; ++j){
//      face_from_p(j,0) = 3 * face_from[0] + j;
//      face_to_p(j,0) = 3 * face_to[0] + j;
//    }

//    const matrixd rot = trans(type_transition2(rot_type[cit->second]));
//    face_from_p = temp(rot * face_from_p);

//    face_from_p
//  }

//  return 0;
//}

int dump_out_r_equation_file(
    const char * equation_file,
    const matrixst & fnode_idx,
    const singularity_graph & sg,
    const matrixst & cut_tet2tet,
    const vector<size_t> & rot_type,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & jump_face_pair2_rot_idx,
    const vector<pair<size_t,size_t> > & g_unknown_face_pair,
    const boost::unordered_map<size_t,size_t> * fnode2group_ptr = 0)
{
  if(!equation_file){
      cerr << "# [error] can not open equation_file. " << endl;
      return __LINE__;
    }

  ofstream ofs(equation_file);
  if(!sg.ge_ptr.get()){
      cerr << "# [error] can not access guass_eliminator" << endl;
      return __LINE__;
    }

  vector<double> fnode(3 * cut_tet2tet.size());
  boost::dynamic_bitset<> fnode_flag(fnode.size());
  std::shared_ptr<jtf::algorithm::gauss_eliminator<double> > fe(
        new jtf::algorithm::gauss_eliminator<double>(fnode, fnode_flag));

  vector<size_t> face_from(3) , face_to(3);
  matrixd face_from_p(3,3), face_to_p(3,3);
  matrixd temp_to_recored_sign(3,1);

  for(size_t fi = 0; fi < g_unknown_face_pair.size(); ++fi){
      const pair<size_t,size_t> & face_pair = g_unknown_face_pair[fi];
      boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator cit
          = jump_face_pair2_rot_idx.find(face_pair);
      if(cit == jump_face_pair2_rot_idx.end()){
          cerr << "# [error] strange can not find face pair. " << __LINE__ << endl;
          return __LINE__;
        }

      const vector<size_t> & face_from_vec = fa_cut.faces_[face_pair.first];
      const vector<size_t> & face_to_vec = fa_cut.faces_[face_pair.second];

      for(size_t i = 0; i < 3; ++i){
          for(size_t j = 0; j < 3; ++j){
              if(cut_tet2tet[face_from_vec[i]] == cut_tet2tet[face_to_vec[j]]){
                  face_from[i] = face_from_vec[i];
                  face_to[i] = face_to_vec[j];
                  break;
                }
            }
        }

      for(size_t i = 0; i < 3; ++i){
          for(size_t j = 0; j < 3; ++j){
              face_from_p(j,i) = 3 * face_from[i] + j;
              face_to_p(j,i) = 3 * face_to[i] + j;
            }
        }

      const matrixd rot = trans(type_transition2(rot_type[cit->second]));

      temp_to_recored_sign = zjucad::matrix::ones<double>(3,1);
      face_from_p = temp(rot * face_from_p);
      temp_to_recored_sign = temp(rot * temp_to_recored_sign);

      // face_to_p(colon(), 0) - face_to_p(colon(), 1) =
      // R^T( face_from_p(colon(),0) - face_from_p(colon(),1))
      for(size_t i = 0; i < 3; ++i){
          jtf::algorithm::equation<double> eq0;
          eq0.add_expression(jtf::algorithm::make_expression(
                               number_rounding(fabs(face_to_p(i,0))),
                               1.0));
          eq0.add_expression(jtf::algorithm::make_expression(
                               number_rounding(fabs(face_to_p(i,1))),
                               -1.0));
          eq0.add_expression(jtf::algorithm::make_expression(
                               number_rounding(fabs(face_from_p(i,0))),
                               -1.0 * (temp_to_recored_sign[i]>0?1.0:-1.0)));

          eq0.add_expression(jtf::algorithm::make_expression(
                               number_rounding(fabs(face_from_p(i,1))),
                               (temp_to_recored_sign[i]>0?1.0:-1.0)));
          fe->add_equation(eq0);

          jtf::algorithm::equation<double> eq1;
          eq1.add_expression(jtf::algorithm::make_expression(
                               number_rounding(fabs(face_to_p(i,0))),
                               1.0));
          eq1.add_expression(jtf::algorithm::make_expression(
                               number_rounding(fabs(face_to_p(i,2))),
                               -1.0));
          eq1.add_expression(jtf::algorithm::make_expression(
                               number_rounding(fabs(face_from_p(i,0))),
                               -1.0 * (temp_to_recored_sign[i]>0?1:-1)));

          eq1.add_expression(jtf::algorithm::make_expression(
                               number_rounding(fabs(face_from_p(i,2))),
                               (temp_to_recored_sign[i]>0?1.0:-1.0)));
          fe->add_equation(eq1);
        }
    }

  //  { // add left g equation to fnode equation
  //    for(jtf::algorithm::gauss_eliminator<double>::const_equation_ptr epr
  //        = sg.ge_ptr->begin(); epr != sg.ge_ptr->end(); ++epr){
  //      if(epr->state() == 2){
  //        const jtf::algorithm::equation<double> & g_eqn = *epr;
  //        vector<jtf::algorithm::equation<double> > f_eqns;
  //        add_g_equation_to_f_equation(g_eqn, f_eqns);
  //        for(size_t i = 0; i < f_eqns.size(); ++i)
  //          fe->add_equation(f_eqns[i]);
  //      }
  //    }
  //  }

  size_t eqidx = 0;
  for(jtf::algorithm::gauss_eliminator<double>::const_equation_ptr epr
      = fe->begin(); epr != fe->end(); ++epr){
      if(epr->state() == 2){
          vector<pair<size_t,double> > equation_temp;
          const jtf::algorithm::equation<double> & eq = *epr;
          for(jtf::algorithm::equation<double>::eq_const_iterator eqcit = eq.begin();
              eqcit != eq.end(); ++eqcit){
              equation_temp.push_back(make_pair(eqcit->index, eqcit->coefficient));
            }
          ofs << "eq " << eqidx++ << " " << equation_temp.size() << endl;
          for(size_t i = 0; i < equation_temp.size(); ++i){
              ofs << equation_temp[i].first << " ";
            }
          ofs << endl;
          for(size_t i = 0; i < equation_temp.size(); ++i){
              ofs << equation_temp[i].second << " ";
            }
          ofs << endl;
        }
    }
  return 0;
}

int dump_out_equation_file(
    const char * equation_file,
    const matrixst & fnode_idx,
    const singularity_graph & sg,
    const matrixst & cut_tet2tet,
    const vector<size_t> & rot_type,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & jump_face_pair2_rot_idx,
    const boost::unordered_map<size_t,size_t> * fnode2group_ptr = 0)
{
  if(!equation_file){
      cerr << "# [error] can not open equation_file. " << endl;
      return __LINE__;
    }

  ofstream ofs(equation_file);
  if(!sg.ge_ptr.get()){
      cerr << "# [error] can not access guass_eliminator" << endl;
      return __LINE__;
    }

  std::vector<double> fnode_vec;
  boost::dynamic_bitset<> fnode_flag;
  std::shared_ptr<jtf::algorithm::gauss_eliminator<double> > fe;
  sg.convert_g_eq_to_f_eq(rot_type, cut_tet2tet, fa_cut, jump_face_pair2_rot_idx,
                          fe, fnode_vec, fnode_flag);

  if(!fe.get()){
      cerr << "conver to fe fail." << endl;
      return __LINE__;
    }
  const std::list<equation<double> > & equations = fe->get_equations();
  //size_t eqi = 0;
  cerr << "# LINE " << __LINE__ << endl;

  set<vector<pair<size_t,double> > > equations_set;
  for(std::list<equation<double> >::const_iterator lcit = equations.begin();
      lcit != equations.end(); ++lcit){
      if(lcit->state() == 0) continue;
      if(lcit->state() == -1) {
          cerr << "# [conflict]" << endl;
          continue;
        }

      const equation<double> & eq = *lcit;
      //cerr << eq << endl;

      if(eq.e_vec_.size() == 2){
          if(fabs(eq.e_vec_.front().coefficient
                  + eq.e_vec_.back().coefficient) < 1e-6)
            {
              const size_t fnode_idx_0 = fnode_idx[eq.e_vec_.front().index];
              const size_t fnode_idx_1 = fnode_idx[eq.e_vec_.back().index];
              if(fnode2group_ptr){
                  boost::unordered_map<size_t,size_t>::const_iterator cit_0 =
                      fnode2group_ptr->find(fnode_idx_0);
                  boost::unordered_map<size_t,size_t>::const_iterator cit_1 =
                      fnode2group_ptr->find(fnode_idx_1);
                  if(cit_0 != fnode2group_ptr->end() &&
                     cit_1 != fnode2group_ptr->end()){
                      if(cit_0->second == cit_1->second)
                        continue;
                      //            cerr << "# [info] " << fnode_idx_0 << " and " << fnode_idx_1
                      //                 << " should be glued " << endl;
                    }
                }
            }
        }
      vector<pair<size_t,double> > one_equation;
      for(equation<double>::eq_const_iterator eqcit = eq.begin();
          eqcit != eq.end();++eqcit){
          const expression<double> & exp = *eqcit;
          one_equation.push_back(make_pair(fnode_idx[exp.index],exp.coefficient));
        }
      equations_set.insert(one_equation);
      //    ofs << "eq " <<  eqi << " " <<  lcit->e_vec_.size() << endl;
      //    for(equation<double>::eq_const_iterator eqcit = eq.begin();
      //        eqcit != eq.end();++eqcit){
      //      const expression<double> & exp = *eqcit;
      //      ofs << fnode_idx[exp.index] << " ";
      //    }
      //    ofs << endl;
      //    for(equation<double>::eq_const_iterator eqcit = eq.begin();
      //        eqcit != eq.end();++eqcit){
      //      const expression<double> & exp = *eqcit;
      //      ofs << exp.coefficient << " ";
      //    }
      //    ofs << endl;
      //++eqi;
    }
  size_t eqi = 0;
  for(set<vector<pair<size_t,double> > >::const_iterator scit =
      equations_set.begin(); scit != equations_set.end(); ++scit, ++eqi){
      const vector<pair<size_t,double> > &one_eq = *scit;
      ofs << "eq " <<  eqi << " " <<  one_eq.size() << endl;
      for(size_t expi = 0; expi < one_eq.size(); ++expi){
          ofs << one_eq[expi].first << " ";
        }
      ofs << endl;
      for(size_t expi = 0; expi < one_eq.size(); ++expi){
          ofs << one_eq[expi].second << " ";
        }
      ofs << endl;
    }
  return 0;
}

int dump_out_chain_file(const char * chain_file,
                        const matrixst & fnode_idx,
                        const singularity_graph & sg)
{
  if(!chain_file){
      cerr << "# [error] can not open chain file" << endl;
      return __LINE__;
    }

  ofstream ofs(chain_file);
  cerr << "# [info] chain num " << sg.chains_.size() << endl;
  for(size_t ci = 0; ci < sg.chains_.size(); ++ci){
      ofs << "chain " << ci << " " << sg.chains_[ci].size()+1 << endl;
      for(size_t ei = 0; ei < sg.chains_[ci].size(); ++ei){
          ofs << fnode_idx[sg.chains_[ci][ei].first] << " ";
        }
      ofs << fnode_idx[sg.chains_[ci].back().second] << endl;
    }
  return 0;
}

int dump_out_param_config(
    const char * group_file,
    const char * equation_file,
    const char * chain_file,
    const char * cut_tet_file,
    const matrixst & cut_tet,
    const matrixd & cut_node,
    const matrixst & cut_tet2tet,
    const vector<size_t> & rot_type,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const singularity_graph & sg,
    const boost::unordered_set<size_t> & restricted_nodes,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & jump_face2rot_idx,
    const std::vector<std::pair<size_t,size_t> > &g_unknown_face_pair)
{
  matrixst fnode_to_cut_point_with_uvw(sg.fnode_.size(),1);
  for(size_t fi = 0; fi < fnode_to_cut_point_with_uvw.size(); ++fi)
    fnode_to_cut_point_with_uvw[fi] =
        3 * sg.get_point_idx_from_node_idx(fi).first +
        sg.get_point_idx_from_node_idx(fi).second;

  boost::unordered_map<size_t,size_t> fnode2group;
  // dump out group_file
  dump_out_group_file(group_file, fnode_to_cut_point_with_uvw, sg,
                      restricted_nodes, &fnode2group);
  //  dump_out_equation_file(equation_file, fnode_to_cut_point_with_uvw, sg,
  //                         cut_tet2tet, rot_type, fa_cut, jump_face2rot_idx,
  //                         &fnode2group);
  dump_out_r_equation_file(equation_file, fnode_to_cut_point_with_uvw, sg,
                           cut_tet2tet, rot_type, fa_cut, jump_face2rot_idx,
                           g_unknown_face_pair, &fnode2group);

  dump_out_chain_file(chain_file, fnode_to_cut_point_with_uvw, sg);
  //dump_out_node_mapping_file(node_to_fnode, sg);

  jtf::mesh::tet_mesh_write_to_zjumat(cut_tet_file, &cut_node, &cut_tet);
  return 0;
}

int load_inner_face_jump_type(
    const char * filename,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open inner face jump type" << endl;
      return __LINE__;
    }
  size_t t0 = -1,t1 = -1,type = -1;
  while(!ifs.eof()){
      ifs >> t0 >> t1 >> type;
      if(t0 == -1 && t1 == -1 && type == -1){
          cerr << "# [info] empty inner_jump_type file." << endl;
          return 0;
        }
      if(type > 23){
          cerr << "# [error] jump type should only be [0,23]." << endl;
          return __LINE__;
        }
      if(is_trivial_type(type)) continue;
      inner_face_jump_type[make_pair(t0,t1)] = type;
    }
  return 0;
}

int load_fnode_group(
    const char * filename,
    const matrixst & cut_tet,
    std::vector<size_t> & fnode)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open fnode_group file." << endl;
      return __LINE__;
    }

  map<size_t,size_t> fnode_to_group;
  size_t group_idx, group_size, type_trash;
  string trash;
  size_t f_node = -1;
  while(!ifs.eof()){
      ifs >> trash >> group_idx >> group_size >> type_trash;
      if(trash.empty()) break;
      for(size_t i = 0; i < group_size; ++i){
          ifs >> f_node;
          map<size_t,size_t>::const_iterator mcit = fnode_to_group.find(f_node);
          if(mcit != fnode_to_group.end())
            std::cerr << "# [error] strange " << f_node << " is in group "
                      << mcit->second << " and " << group_idx << std::endl;
          fnode_to_group[f_node] = group_idx;
        }
      trash.clear();
    }

  fnode.resize((max(cut_tet) + 1) * 3);
  for(size_t i = 0; i < fnode.size(); ++i)
    fnode[i] = i;
  for(map<size_t,size_t>::const_iterator mcit = fnode_to_group.begin();
      mcit != fnode_to_group.end(); ++mcit){
      fnode[mcit->first] = mcit->second;
    }
  return 0;
}

int dump_surface_restricted_type_to_vtk(
    const char * filename,
    const std::string & type_name,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const boost::unordered_map<size_t,size_t> & surface_type)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      cerr << "# [error] can not open file." << endl;
      return __LINE__;
    }

  vector<size_t> face_vec;
  face_vec.reserve(3 * surface_type.size());
  vector<size_t> face_type;
  face_type.reserve(surface_type.size());
  for(boost::unordered_map<size_t,size_t>::const_iterator cit = surface_type.begin();
      cit != surface_type.end(); ++cit){
      const vector<size_t> & one_face = fa.faces_.at(cit->first);
      face_vec.insert(face_vec.end(), one_face.begin(), one_face.end());
      face_type.push_back(cit->second);
    }

  tri2vtk(ofs, &node[0], node.size(2), &face_vec[0], face_vec.size()/3);
  cell_data(ofs, &face_type[0], face_type.size(), const_cast<char*>(type_name.c_str()));
  return 0;
}


int load_restricted_edges(
    const char * filename,
    boost::unordered_set<std::pair<size_t,size_t> > &degenerated_edges)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open restricted edge." << endl;
      return __LINE__;
    }

  degenerated_edges.clear();
  size_t edge_num = 0;
  ifs >> edge_num;
  pair<size_t,size_t> edge;
  for(size_t ei = 0; ei < edge_num; ++ei){
      ifs >> edge.first >> edge.second;
      degenerated_edges.insert(edge);
    }

  return 0;
}

int load_g_unknown_face(
    const char * filename,
    std::vector<std::pair<size_t,size_t> > & face_pairs)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open g_unknown face file." << endl;
      return __LINE__;
    }

  size_t face_num = 0;
  ifs >> face_num;
  face_pairs.resize(face_num);
  for(size_t i = 0; i < face_num; ++i){
      ifs >> face_pairs[i].first >> face_pairs[i].second;
    }

  return 0;
}

int dump_surface_restricted_type(
    const char * filename,
    const boost::unordered_map<size_t,size_t> & surface_type)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      cerr << "# [error] can not open surface restricted type." << endl;
      return __LINE__;
    }
  for(boost::unordered_map<size_t,size_t>::const_iterator cit =
      surface_type.begin(); cit != surface_type.end(); ++cit){
      ofs << cit->first << " " << cit->second << endl;
    }

  return 0;
}

int load_surface_restricted_type(
    const char * filename,
    boost::unordered_map<size_t,size_t> & surface_type)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open surface restricted type." << endl;
      return __LINE__;
    }
  bool is_restricted_type = true;
  size_t face_idx, face_type;
  while(!ifs.eof()){
      ifs >> face_idx >> face_type;
      if(face_type > 2)
        is_restricted_type = false;
      surface_type[face_idx] = face_type;
    }

  if(!is_restricted_type) return 1;
  return 0;
}



int load_surface_normal_align_type(
    const char * filename,
    const jtf::mesh::face2tet_adjacent &fa,
    boost::unordered_map<size_t,size_t> &outside_face_type)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open surface normal align type file" << endl;
      return __LINE__;
    }

  size_t f0,f1,f2,type;
  while(!ifs.eof()){
      ifs >> f0 >> f1 >> f2 >> type;
      const size_t face_real_idx = fa.get_face_idx(f0,f1,f2);
      if(face_real_idx == -1){
          cerr << "# [error] face " << f0 << " " << f1 << " " << f2
               << "do not exist." << endl;
          return __LINE__;
        }
      outside_face_type[face_real_idx] = type;
    }
  return 0;
}

int convert_surface_normal_type2_restricted_type(
    boost::unordered_map<size_t,size_t> & surface_type)
{
  using namespace zjucad::matrix;
  for(auto & one_face : surface_type){
      const size_t trans_type = get_trans_type(one_face.second);
      size_t rt = -1;
      matrix<double> ttm = trans(type_transition2(trans_type));
      for(size_t i = 0; i < 3; ++i)
        if(fabs(fabs(ttm(0,i)) -1) < 1e-6){
            rt = i; break;
          }
      if(rt == -1) return __LINE__;
      one_face.second = rt;
    }

  return 0;
}

int load_surface_type(const char * filename,
                      const jtf::mesh::face2tet_adjacent &fa,
                      const matrixst & outside_face_idx,
                      matrixst & outside_face_type)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open surface normal align type file" << endl;
      return __LINE__;
    }

  map<size_t,size_t> outside_face_real_idx2_outidx;
  for(size_t t = 0; t < outside_face_idx.size(); ++t){
      outside_face_real_idx2_outidx[outside_face_idx[t]] = t;
    }

  outside_face_type.resize(outside_face_idx.size());

  size_t f0,f1,f2,type;
  while(!ifs.eof()){
      ifs >> f0 >> f1 >> f2 >> type;
      const size_t face_real_idx = fa.get_face_idx(f0,f1,f2);
      if(face_real_idx == -1){
          cerr << "# [error] face " << f0 << " " << f1 << " " << f2
               << "do not exist." << endl;
          return __LINE__;
        }
      outside_face_type[outside_face_real_idx2_outidx[face_real_idx]] = type;
    }
  return 0;
}

int load_surface_type(
    const char * filename,
    boost::unordered_map<size_t,size_t> & surface_type,
    const jtf::mesh::face2tet_adjacent *fa)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open surface type file" << endl;
      return __LINE__;
    }

  string face_str;
  getline(ifs, face_str);
  vector<string> split_string;
  boost::split(split_string, face_str, boost::is_any_of(" "));
  if(split_string.size() == 4) {// f0 f1 f2 type
      if(fa == 0) {
          cerr << "# [error] I need face2tet_adjacent." << endl;
          return __LINE__;
        }
      int rtn = load_surface_normal_align_type(filename, *fa, surface_type);
      if(rtn == 0) // if read ok, return 1 to demonstrate that read surface normal alignment type
        return 1;
      else
        return rtn;
    }else if(split_string.size() == 2){ // f_idx type
      return load_surface_restricted_type(filename, surface_type);
    }else{
      cerr << "# [error] should not be here. " << __FILE__ << __LINE__ << endl;
      return __LINE__;
    }
}

int dump_inner_face_jump_type(
    const char * filename,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      cerr << "# [error] can not open filename" << endl;
      return __LINE__;
    }
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
  for(mpscit it = inner_face_jump_type.begin();
      it != inner_face_jump_type.end(); ++it){
      if(is_trivial_type(it->second)) continue;
      ofs << it->first.first << " "
          << it->first.second << " "
          << it->second << endl;
    }
  return 0;
}


int dump_surface_normal_align_type(
    const char * filename,
    const matrixst & outside_face,
    const matrixst & outside_face_type)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      cerr << "# [error] can not open surface normal align type file" << endl;
      return __LINE__;
    }
  for(size_t t = 0; t < outside_face.size(2); ++t){
      ofs << outside_face(0,t) << " "
          << outside_face(1,t) << " "
          << outside_face(2,t) << " "
          << outside_face_type[t] << endl;
    }
  return 0;
}


int dump_surface_normal_align_type_map(
    const char * filename,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const boost::unordered_map<size_t,size_t> & face_jump_type)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      std::cerr << "# [error] can not open surface normal align type file" << std::endl;
      return __LINE__;
    }
  for(size_t t = 0; t < outside_face.size(2); ++t){
      auto it = face_jump_type.find(outside_face_idx[t]);
      if(it == face_jump_type.end()) {
          cerr << "# [error] strange can not find outside_face_idx " << outside_face_idx[t] << endl;
          continue;
        }
      ofs << outside_face(0,t) << " "
          << outside_face(1,t) << " "
          << outside_face(2,t) << " "
          << it->second << std::endl;
    }
  return 0;
}

int load_loop_points(
    const char * filename,
    boost::unordered_set<size_t> & loop_points)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open loop_point file." << endl;
      return __LINE__;
    }

  loop_points.clear();
  string loop;
  size_t loop_idx, loop_point_num, point_idx;
  while(!ifs.eof()){
      ifs >> loop >> loop_idx >> loop_point_num;
      if(loop.empty()) break;
      for(size_t i = 0; i < loop_point_num; ++i){
          ifs >> point_idx;
          loop_points.insert(point_idx);
        }
      loop.clear();
    }

  return 0;
}

int load_restricted_edge(
    const char * filename,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & restrict_edge_type)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open " << filename << endl;
      return __LINE__;
    }

  size_t p0, p1, type, num = 0;

  restrict_edge_type.clear();
  ifs >> num;
  for(size_t i = 0; i < num; ++i){
      ifs >> p0 >> p1 >> type;
      if(p0 > p1) swap(p0,p1);
      restrict_edge_type[make_pair(p0, p1)] = type;
    }
  return 0;
}

int dump_restricted_edge(
    const char *filename,
    const std::vector<std::pair<size_t,size_t> > & restriced_edges,
    const std::vector<size_t> & restriced_edges_type)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      cerr << "# [error] can not open " << filename << endl;
      return __LINE__;
    }

  if(restriced_edges.size() != restriced_edges_type.size()){
      cerr << "# [error] not compatiable edges." << endl;
      return __LINE__;
    }

  ofs << restriced_edges.size() << endl;
  for(size_t ei = 0; ei < restriced_edges.size(); ++ei){
      ofs << restriced_edges[ei].first << " " << restriced_edges[ei].second << " ";
      ofs << restriced_edges_type[ei] << endl;
    }

  return 0;
}

int load_group_file(const char * group_file,
		    std::vector<std::vector<size_t> > & groups,
		    std::vector<bool> * integer_group_flag)
{
  string group_name;
  size_t group_index, group_size, element, integer_or_not;
  ifstream ifs(group_file);
  if(ifs.fail()){
      cerr << "# [error] can not open group file" << endl;
      return __LINE__;
    }

  groups.clear();
  while(!ifs.eof()){
      ifs >> group_name >> group_index >> group_size >> integer_or_not;
      if(group_name.empty()) break;
      vector<size_t> one_group(group_size);

      for(size_t i = 0; i < group_size; ++i){
          ifs >> one_group[i];
        }
      groups.push_back(one_group);
      group_name.clear();
      if(integer_group_flag){
          integer_group_flag->push_back(integer_or_not==0?false:true);
        }
    }
  return 0;
}

int load_chain_file(const char * chain_file,
                    std::vector<std::vector<size_t> > & chains)
{
  string chain_name;
  size_t chain_index, chain_size, element;
  ifstream ifs(chain_file);
  if(ifs.fail()){
      cerr << "# [error] can not open chain file" << endl;
      return __LINE__;
    }

  chains.clear();
  while(!ifs.eof()){
      ifs >> chain_name >> chain_index >> chain_size;
      if(chain_name.empty()) break;
      vector<size_t> one_chain(chain_size);

      for(size_t i = 0; i < chain_size; ++i){
          ifs >> one_chain[i];
        }
      chains.push_back(one_chain);
      chain_name.clear();
    }
  return 0;
}

int dump_surface_normal_align_type2vtk(
    const char * filename,
    const zjucad::matrix::matrix<double> &node,
    const jtf::mesh::face2tet_adjacent &fa,
    const boost::unordered_map<size_t,size_t> & outside_face_type)
{
  ofstream ofs(filename);
  if(ofs.fail()){
      cerr << "# [error] can not open surface normal align type file" << endl;
      return __LINE__;
    }

  vector<size_t> faces;
  vector<size_t> type;

  for(const auto & one_face: outside_face_type){
      const vector<size_t> &f = fa.faces_[one_face.first];
      faces.insert(faces.end(), f.begin(), f.end());
      const matrixd rot = trans(type_transition2(one_face.second));
      for(size_t i = 0; i < 3; ++i)
        if(fabs(fabs(rot(0,i))-1) < 1e-7) {
            type.push_back(i);
            break;
          }
    }

  tri2vtk(ofs, &node[0], node.size(2), &faces[0], faces.size()/3);
  cell_data(ofs, &type[0], type.size(), "surface_type");
  return 0;
}


int dump_inner_face_jump_type2vtk(
    const char *filename,
    const zjucad::matrix::matrix<double> & node,
    const zjucad::matrix::matrix<size_t> & tet,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_type)
{
  using namespace zjucad::matrix;
  vector<size_t> all_faces;
  vector<size_t> face(3);
  for(const auto & one_face : inner_type){
      if(one_face.second == TRIVIAL_TYPE) continue;
      jtf::mesh::find_common_face(tet(colon(), one_face.first.first),
                                  tet(colon(), one_face.first.second),
                                  &face[0]);
      all_faces.insert(all_faces.end(), face.begin(), face.end());
    }
  ofstream ofs(filename);
  tri2vtk(ofs, &node[0], node.size(2), &all_faces[0], all_faces.size()/3);
}
