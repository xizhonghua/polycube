#include <boost/algorithm/string.hpp>
#include <sstream>

#include <zjucad/matrix/matrix.h>
#include <zjucad/ptree/ptree.h>
#include "util.h"

using namespace std;
using namespace zjucad::matrix;
using boost::property_tree::ptree;

// assume cons_weight is allocated
void load_cons_weight(
    ptree &pt,
    const vector<tuple<string,size_t,size_t> > & cons_type,
    matrix<double> &cons_weight)
{

  if(!zjucad::has("cons/weight.value",pt))  return;
  const string weight = pt.get<string>("cons/weight.value");
  vector<string> temp_str;

  boost::split(temp_str, weight, boost::is_any_of(","));

#if 0 // this code can not work on g++ 4.6.3 and boost 1.46
  {
    // detect float number:integer:integer
    boost::regex range("[-+]?[0-9]*\.?[0-9]+:[^\\d+$]:[^\\d+$]");
    for(const auto & str : temp_str){
        if(boost::regex_match(str,range)){
            cerr << "match " << str << endl;
          }else
          cerr << "not match " << str << endl;
      }
  }
#endif

  if(temp_str.size() != cons_type.size()){
      throw std::logic_error("# [error] cons_type is not compatiable with weight.");
    }
  double w = 0;
  for(size_t ti = 0; ti < temp_str.size(); ++ti){
      stringstream ss;
      ss << temp_str[ti]; ss >> w;
      const tuple<string,size_t,size_t> & one_type = cons_type[ti];
      for(size_t i =  get<1>(one_type); i < get<2>(one_type);  ++i){
          cons_weight[i] = w;
        }
    }
}
