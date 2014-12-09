#include "io.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <boost/dynamic_bitset.hpp>
#include <jtflib/algorithm/gauss_elimination.h>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace zjucad::matrix;

int load_equation_file(
    const char * file,
    vector<vector<size_t> > &variable_idx,
    vector<vector<double> > &coefficient,
    const zjucad::matrix::matrix<size_t> * node_mapping)
{
  ifstream ifs(file);
  if(ifs.fail()){
      cerr << "# [error] can not open equation_file." << endl;
      return __LINE__;
    }
  string temp;
  size_t variable_number = 0;

  variable_idx.clear();
  coefficient.clear();

  vector<size_t> one_variable_idx;
  vector<double> one_coefficients;

  vector<vector<size_t> > variable_idx_temp;
  vector<vector<double> > coefficient_temp;

  while(!ifs.eof()){
      ifs >> temp >> temp >> variable_number; // eq eq_idx eq_variable_number
      if(ifs.eof()) break;
      one_coefficients.resize(variable_number);
      one_variable_idx.resize(variable_number);
      for(size_t i = 0; i < variable_number; ++i){
          ifs >> one_variable_idx[i];
        }
      for(size_t i = 0; i < variable_number; ++i){
          ifs >> one_coefficients[i];
        }
      variable_idx_temp.push_back(one_variable_idx);
      coefficient_temp.push_back(one_coefficients);
    }

  if(!node_mapping){
      variable_idx = variable_idx_temp;
      coefficient = coefficient_temp;
    }else{
      vector<double> node_v(node_mapping->size());
      boost::dynamic_bitset<> node_v_flag(node_mapping->size());
      jtf::algorithm::gauss_eliminator<double> ge(node_v, node_v_flag);
      for(size_t eqi = 0; eqi < variable_idx_temp.size(); ++eqi){
          jtf::algorithm::equation<double> one_eqn;
          const vector<size_t> & one_eqn_idx = variable_idx_temp[eqi];
          const vector<double> & one_eqn_coeff = coefficient_temp[eqi];
          for(size_t vi = 0; vi < one_eqn_idx.size(); ++vi)
            one_eqn.add_expression(
                  jtf::algorithm::make_expression(
                    (*node_mapping)[one_eqn_idx[vi]], one_eqn_coeff[vi]));
          one_eqn.standardization();
          const size_t state = one_eqn.state();
          if(state == 0) continue; // cleared;
          if(state == 1 || state == -1)
            throw std::invalid_argument("group node mapping is not suitable with equation");
          ge.add_equation(one_eqn);
        }
      for(jtf::algorithm::gauss_eliminator<double>::const_equation_ptr it
          = ge.begin(); it != ge.end(); ++it){
          if(it->state() == 2){
              one_variable_idx.clear();
              one_coefficients.clear();
              const auto & one_eqn = *it;
              for(const auto & one_exp : one_eqn){
                  one_variable_idx.push_back(one_exp.index);
                  one_coefficients.push_back(one_exp.coefficient);
                }
              variable_idx.push_back(one_variable_idx);
              coefficient.push_back(one_coefficients);
            }
        }
    }
  return 0;
}


int dump_equation_file(
    const char * file,
    const std::vector<std::vector<size_t> > & variable_idx,
    const std::vector<std::vector<double> > & variable_coeff)
{
  assert(variable_coeff.size() == variable_coeff.size());
  ofstream ofs(file);
  if(ofs.fail()){
      cerr << "# [error] can not write equation file." << endl;
      return __LINE__;
    }

  for(size_t ei = 0; ei < variable_idx.size(); ++ei){
      ofs << "eqi " << ei << " " << variable_idx[ei].size() << endl;
      for(size_t vi = 0; vi < variable_idx[ei].size(); ++vi)
          ofs << variable_idx[ei][vi] << " " ;
      ofs << endl;
      for(size_t vi = 0; vi < variable_coeff[ei].size(); ++vi)
        ofs << variable_coeff[ei][vi] << " ";
      ofs << endl;
    }

  return 0;
}

int load_variable_group_file(const char * file,
                             vector<vector<size_t> > &groups)
{
  ifstream ifs(file);
  if(ifs.fail()){
      cerr << "# [error] can not open variable group file." << endl;
      return __LINE__;
    }
  groups.clear();

  string temp;
  vector<string> idx_str;
  vector<size_t> empty_g;
  size_t st_v, max_st = 0;
  while(!ifs.eof()){
      getline(ifs, temp); // group name or other things
      if(ifs.eof()) break;
      getline(ifs, temp); // variables in this group
      idx_str.clear();
      empty_g.clear();
      boost::split(idx_str, temp, boost::is_any_of(" "));
      for(const auto & idx: idx_str){
          stringstream ss;
          ss << idx; ss >> st_v;
          empty_g.push_back(st_v);
          max_st = max(max_st, st_v);
        }
      groups.push_back(empty_g);
    }

  return 0;
}

int load_integer_constraint(const char * filename,
                            std::vector<size_t> & integer_variables)

{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open integer constraint file." << endl;
      return __LINE__;
    }
  set<size_t> integer_set;
  size_t temp, v_size;
  int integer_or_not; // 1: integer, 0: not integer
  string str_temp;
  while(!ifs.eof()){
      ifs >> str_temp >> temp >> v_size >> integer_or_not;
      if(ifs.eof()) break;
      for(size_t i = 0; i < v_size; ++i){
          ifs >> temp;
          if(integer_or_not)
            integer_set.insert(temp);
        }
    }
  integer_variables.resize(integer_set.size());
  std::copy(integer_set.begin(),integer_set.end(), integer_variables.begin());
  return 0;
}

int load_restricted_path(
    const char * filename,
    std::vector<pair<size_t,size_t> > & restricted_path)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open restricted patch file." << endl;
      return __LINE__;
    }
  size_t path_number;
  ifs >> path_number;
  restricted_path.resize(path_number);
  size_t edge_number;
  std::tuple<size_t,size_t,size_t> one_edge;
  for(size_t i = 0; i < path_number; ++i){
      ifs >> edge_number;
      for(size_t ei = 0; ei < edge_number; ++ei){
          // only record first and last point of each pathl
          if(ei == 0)
            ifs >> std::get<0>(one_edge) >> std::get<1>(one_edge) >> std::get<2>(one_edge) ;
          else
            ifs >> std::get<1>(one_edge) >> std::get<1>(one_edge) >> std::get<2>(one_edge) ;
        }

      restricted_path[i].first = 3 * std::get<0>(one_edge) + std::get<2>(one_edge);
      restricted_path[i].second = 3 * std::get<1>(one_edge) + std::get<2>(one_edge);
    }
  return 0;
}
