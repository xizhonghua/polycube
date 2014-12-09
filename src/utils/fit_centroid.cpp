#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <fstream>
#include <jtflib/mesh/util.h>
#include <zjucad/matrix/matrix.h>
#include <hjlib/function/function.h>
#include <zjucad/matrix/io.h>
#include <zjucad/optimizer/optimizer.h>
#include "../common/vtk.h"
using namespace std;
using namespace zjucad::matrix;

int read_tet_abq(const char * filename, matrix<double> & node, matrix<size_t> & tet)
{
  std::ifstream ifs(filename);
  if(ifs.fail()) return __LINE__;

  string temp;
  node.resize(3,2398);
  tet.resize(4,6658);
  while(!ifs.eof()){
      ifs >> temp;
      if(temp == "*NODE"){
          for(size_t i = 0; i < 2398; ++i){
              ifs >> temp;
              for(size_t di = 0; di < 3; ++di){
                  ifs >> node(di,i);
                }
            }
        }else if(temp == "*ELEMENT"){
          ifs >> temp;
          for(size_t i = 0; i < 6658; ++i){
              ifs >> temp;
              for(size_t di = 0 ; di < 4; ++di){
                  ifs >> tet(di,i);
                }
            }
        }

    }

  tet -= 1;
  return 0;
}

void get_centroid(const matrix<double> &node,
                  const matrix<size_t> &tet,
                  const matrix<double> &seq,
                  matrix<double> &centroid)
{
  centroid = zeros<double>(3, seq.size(2));

  matrix<double> vol(tet.size(2),1);
  matrix<double> node_temp;
  matrix<double> centeroid_one(3,1);
  matrix<double> one_tet_node;
  for(size_t ti = 0; ti < seq.size(2); ++ti){
      node_temp = node(colon()) + seq(colon(),ti);
      itr_matrix<const double*> node_temp_m(3, node.size(2), &node_temp[0]);

      for(size_t teti = 0; teti < tet.size(2); ++teti){
          one_tet_node =  node_temp_m(colon(), tet(colon(),teti));
          vol[teti] = fabs(jtf::mesh::cal_tet_vol(one_tet_node));
          centeroid_one *= 0;
          for(size_t pi = 0; pi < tet.size(1); ++pi){
              centeroid_one += one_tet_node(colon(),pi);
            }
          centeroid_one /= 4;

          centroid(colon(),ti) += vol[teti] * centeroid_one;
        }
      centroid(colon(),ti) /= std::accumulate(vol.begin(), vol.end(), 0.0);

      // centroid(colon(),ti) = node_temp_m(colon(),1974);
    }
}

typedef hj::function::function_t<double,int32_t> hj_func;
typedef std::shared_ptr<const hj_func> hj_func_cons_ptr;

class fit_func : public hj_func
{
public:
  fit_func(const size_t t,
           const double y,
           double w = 1):y_(y), t_(t), w_(w){}
  virtual size_t dim_of_f() const {return 1;}
  virtual size_t dim_of_x() const {return 3;}
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx) const
  {
    *f = (x[0] * t_*t_ + x[1]*t_+x[2] - y_)*w_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr, int32_t *idx, hj::function::func_ctx *ctx) const
  {
    ptr[1] = ptr[0]+3;
    idx[ptr[0]] = 0;
    val[ptr[0]] = w_ * t_*t_;

    idx[ptr[0]+1] = 1;
    val[ptr[0]+1] = w_*t_;

    idx[ptr[0]+2] = 2;
    val[ptr[0]+2] = w_;

    return 0;
  }
private:
  const double y_;
  const size_t t_;
  const double w_;
};



int fit_centroid(int argc , char *argv[])
{
  if(argc != 5){
      std::cerr << "# [error] fit_centroid tet.abq u.b  start end" << std::endl;
      return __LINE__;
    }

  matrix<double> node;
  matrix<size_t> tet;
  matrix<double> seq;

  if(read_tet_abq(argv[1], node, tet))
    return __LINE__;

  cerr << "tet number " << tet.size(2) << std::endl;
  cerr << "node number " << node.size(2) << std::endl;

  if(jtf::mesh::read_matrix(argv[2], seq))
    return __LINE__;

  cerr << "seq " << seq.size(1) << " " << seq.size(2) << endl;

  size_t start = atoi(argv[3]);
  size_t end = atoi(argv[4]);

  cerr << "range " << start << " " << end << endl;
  matrix<double> centroid(3, seq.size(2));
  get_centroid(node, tet, seq, centroid);

  {
    ofstream ofs("centroid.vtk");
    vector<size_t> lines;
    matrix<double> new_centroid = centroid;
    for(size_t qi = 0; qi < seq.size(2)-1; ++qi){
        new_centroid(0,qi) = qi*1.0/1000;
        new_centroid(0,qi+1) = (qi+1)*1.0/1000;
        lines.push_back(qi);
        lines.push_back(qi+1);
      }
    line2vtk(ofs, &new_centroid[0],new_centroid.size(2), &lines[0], lines.size()/2);
  }

  cerr << "finish centroid " << endl;
  // assume only y changes
  vector<hj_func_cons_ptr> funcs;
  double w = 1;
  for(size_t i = start; i < end+1; ++i){
      if(i == start || i == end) w = 1;
      else w = 1;
      funcs.push_back(hj_func_cons_ptr(new fit_func(i-start,centroid(1,i), w)));
    }

  hj_func_cons_ptr all_func(hj::function::new_catenated_function<double,int32_t>(funcs));
  matrix<double> abc(3,1);
  matrix<double> residual(all_func->dim_of_f(),1);

  boost::property_tree::ptree pt;
  pt.put("package.value","hj");
  pt.put("iter.value",20);

  zjucad::optimize(*all_func, abc, residual, pt);

  cerr << "abc " << abc << endl;

  matrix<double> difference = zeros<double>(3,seq.size(2));
  matrix<double> abc_center(1,end+1-start);
  for(size_t i = start; i < end+1; ++i){
      difference(1,i) = abc[0]*(i-start)*(i-start)+abc[1]*(i-start)+abc[2]-centroid(1,i);
      abc_center[i-start] = abc[0]*(i-start)*(i-start)+abc[1]*(i-start)+abc[2];
    }

  {
    ofstream ofs("abc_centroid.vtk");
    vector<size_t> lines;
    matrix<double> new_abc_center = zeros<double>(3,abc_center.size());
    new_abc_center(1,colon()) = abc_center;
    for(size_t qi = 0; qi < new_abc_center.size(2)-1; ++qi){
        new_abc_center(0,qi) = qi*1.0/1000;
        new_abc_center(0,qi+1) = (qi+1)*1.0/1000;
        lines.push_back(qi);
        lines.push_back(qi+1);
      }
    line2vtk(ofs, &new_abc_center[0],new_abc_center.size(2), &lines[0], lines.size()/2);
  }


  cerr << zjucad::matrix::max(difference) << std::endl;

  for(size_t qi = 0; qi < seq.size(2); ++qi){
      itr_matrix<double*> seq_m(3, node.size(2), &seq(0,qi));
      for(size_t pi = 0; pi < node.size(2); ++pi){
          seq_m(1,pi) += difference(1,qi);
        }
    }

  matrix<double> node_temp;

  for(size_t qi = 0; qi < seq.size(2); ++qi){
      node_temp = node(colon()) + seq(colon(),qi);
      itr_matrix<double*> node_temp_m(3, node.size(2), &node_temp[0]);
      stringstream ss;
      ss << "output_" << qi << ".vtk";
      ofstream ofs(ss.str().c_str());
      tet2vtk(ofs, &node_temp_m[0], node_temp_m.size(2), &tet[0], tet.size(2));
    }

  {
    get_centroid(node, tet, seq, centroid);
    ofstream ofs("new_centroid.vtk");
    vector<size_t> lines;
    matrix<double> new_centroid = centroid;
    for(size_t qi = 0; qi < seq.size(2)-1; ++qi){
        new_centroid(0,qi) = qi*1.0/1000;
        new_centroid(0,qi+1) = (qi+1)*1.0/1000;
        lines.push_back(qi);
        lines.push_back(qi+1);
      }
    line2vtk(ofs, &new_centroid[0],new_centroid.size(2), &lines[0], lines.size()/2);
  }

  if(jtf::mesh::write_matrix("new_Ub", seq))
    return __LINE__;

  return 0;
}
