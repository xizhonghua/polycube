#include "gauss_elimination_test.h"
#include <vector>

using namespace std;
using namespace jtf::algorithm;

//  2x1 + 2x0 + 2x0 = 2
//  x0 = 2;
//  5X0 + 2*X2 = 2
//  0.5x3+x4=2


void gauss_eliminator_test::setUp()
{

}

void gauss_eliminator_test::tearDown()
{

}

void gauss_eliminator_test::gauss_elimination2()
{
  double A[] = {
    0.769580 , 0.309191,
    0.249734 , 0.099414
  };

  double B [] = {
    0.119502,
    0.030924
  };
  double x [] = {
    -3.2734,
    8.5340
  };
  typedef double val_type;
  const size_t row = 2;
  vector<val_type> node(row,0);
  boost::dynamic_bitset<> node_flag(row);

  jtf::algorithm::gauss_eliminator<val_type> ge(node, node_flag);
  {
    for(size_t ei = 0; ei < row; ++ei){
      equation<val_type> eq;
      eq.value() = B[ei];
      for(size_t i = 0; i < row; ++i){
        eq.add_expression(jtf::algorithm::make_expression(
                            static_cast<size_t>(i), A[ei*row + i]));
      }
      ge.add_equation(eq);
    }
  }

  for(size_t xi = 0; xi < row; ++xi){

    if(fabs(node[xi] - x[xi])>1e-6)
      cerr << " error : " << fabs(node[xi] - x[xi]) << " "
           << xi  << " = " << node[xi] << endl;
  }
}

void gauss_eliminator_test::gauss_elimination4()
{
  double A[] = {
    0.525642 ,  0.290499 ,  0.543517 ,  0.544226,
    0.737130 ,  0.179899 ,  0.797661 ,  0.353023,
    0.193812 ,  0.369145 ,  0.068310 ,  0.177968,
    0.314097 ,  0.280042 ,  0.436713 ,  0.506834
  };

  double B [] = {
    0.71630,
    0.25728,
    0.50137,
    0.90972
  };
  double x [] = {
    -1.7389,
    1.5572,
    1.1117,
    1.0543
  };
  typedef double val_type;
  const size_t row = 4;
  vector<val_type> node(row,0);
  boost::dynamic_bitset<> node_flag(row);

  jtf::algorithm::gauss_eliminator<val_type> ge(node, node_flag);
  {
    for(size_t ei = 0; ei < row; ++ei){
      equation<val_type> eq;
      eq.value() = B[ei];
      for(size_t i = 0; i < row; ++i){
        eq.add_expression(jtf::algorithm::make_expression(
                            static_cast<size_t>(i), A[ei*row + i]));
      }
      ge.add_equation(eq);
    }
  }

  for(size_t xi = 0; xi < row; ++xi){

    if(fabs(node[xi] - x[xi])>1e-6)
      cerr << " error : " << fabs(node[xi] - x[xi]) << " "
           << xi  << " = " << node[xi] << endl;
  }
}

void gauss_eliminator_test::gauss_elimination8()
{
  double A[] = {
    0.0555267 ,  0.0598231 ,  0.8163864 ,  0.0450317 ,  0.7636139 ,  0.0166776 ,  0.4378259 ,  0.7238870,
    0.2101941 ,  0.7574916 ,  0.0471822 ,  0.7717254 ,  0.2535038 ,  0.0566050 ,  0.7066258 ,  0.3106741,
    0.2050745 ,  0.0875299 ,  0.4345902 ,  0.8384609 ,  0.3174239 ,  0.0143568 ,  0.2404941 ,  0.3915063,
    0.2992646 ,  0.6297010 ,  0.1205361 ,  0.9044948 ,  0.8548211 ,  0.0728032 ,  0.6400402 ,  0.2896103,
    0.7327456 ,  0.8003895 ,  0.2122077 ,  0.2271898 ,  0.3370035 ,  0.2145170 ,  0.0042949 ,  0.8020700,
    0.2070253 ,  0.1203956 ,  0.4986220 ,  0.1892236 ,  0.5626788 ,  0.4195724 ,  0.8800979 ,  0.7536519,
    0.8336110 ,  0.0109277 ,  0.4024928 ,  0.6122628 ,  0.3172605 ,  0.2051002 ,  0.1758442 ,  0.1625469,
    0.0935987 ,  0.0323862 ,  0.7350551 ,  0.6689585 ,  0.7509560 ,  0.5510155 ,  0.8769548 ,  0.6400437
  };

  double B [] = {
    0.890103,
    0.896679,
    0.328070,
    0.087875,
    0.263768,
    0.735404,
    0.903774,
    0.856485
  };
  double x [] = {
    0.67288,
    1.28147,
    2.90459,
    -0.55278,
    -1.30952,
    -0.28631,
    1.28196,
    -1.55663

  };
  typedef double val_type;
  vector<val_type> node(8,0);
  boost::dynamic_bitset<> node_flag(8);

  jtf::algorithm::gauss_eliminator<val_type> ge(node, node_flag);
  {
    for(size_t ei = 0; ei < 8; ++ei){
      equation<val_type> eq;
      eq.value() = B[ei];
      for(size_t i = 0; i < 8; ++i){
        eq.add_expression(jtf::algorithm::make_expression(
                            static_cast<size_t>(i), A[ei*8 + i]));
      }
      ge.add_equation(eq);
    }
  }

  for(size_t xi = 0; xi < 8; ++xi){

    if(fabs(node[xi] - x[xi])>1e-6)
      cerr << " error : " << fabs(node[xi] - x[xi]) << " "
           << xi  << " = " << node[xi] << endl;
  }
}

void gauss_eliminator_test::gauss_elimination3()
{
  double A[] = {
    0.689408 ,  0.367467 ,  0.657361,
    0.400340 ,  0.200086 ,  0.479519,
    0.178435 ,  0.994344 ,  0.091681
  };

  double B [] = {
    0.80140,
    0.88989,
    0.76432
  };
  double x [] = {
    -3.6462,
    1.0101,
    4.4785
  };
  typedef double val_type;
  const size_t row = 3;
  vector<val_type> node(row,0);
  boost::dynamic_bitset<> node_flag(row);

  jtf::algorithm::gauss_eliminator<val_type> ge(node, node_flag);
  {
    for(size_t ei = 0; ei < row; ++ei){
      equation<val_type> eq;
      eq.value() = B[ei];
      for(size_t i = 0; i < row; ++i){
        eq.add_expression(jtf::algorithm::make_expression(
                            static_cast<size_t>(i), A[ei*row + i]));
      }
      ge.add_equation(eq);
    }
  }

  for(size_t xi = 0; xi < row; ++xi){
    if(fabs(node[xi] - x[xi])>1e-6)
      cerr << " error : " << fabs(node[xi] - x[xi]) << " "
           << xi  << " = " << node[xi] << endl;
  }
}

void gauss_eliminator_test::gauss_elimination1()
{
  typedef double val_type;
  vector<val_type> node(5,0);
  boost::dynamic_bitset<> node_flag(5);

  jtf::algorithm::gauss_eliminator<val_type> ge(node, node_flag);

  {
    // x0  + x1 + x2 = 2
    equation<val_type> eq;
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(0), static_cast<val_type>(1)));
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(1), static_cast<val_type>(1)));
//    eq.add_expression(jtf::algorithm::make_expression(
//                        static_cast<size_t>(2), static_cast<val_type>(1)));
    eq.value() = 2;
    ge.add_equation(eq);
  }

  {
    // x2 + x3 = 4
    equation<val_type> eq;
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(2), static_cast<val_type>(1)));
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(3), static_cast<val_type>(1)));
    eq.value() = 4;
    ge.add_equation(eq);
  }


  {
    // x1 + x2 = 1
    equation<val_type> eq;
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(1), static_cast<val_type>(1)));
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(2), static_cast<val_type>(1)));
    eq.value() = 1;
    ge.add_equation(eq);
  }

  for(size_t t = 0; t < node_flag.size(); ++t){
    if(node_flag[t] == true){
      cerr << "# node " << t << " = " << node[t] << endl;
    }else
      cerr << "# node " << t << " unknown." << endl;
  }
  for(jtf::algorithm::gauss_eliminator<val_type>::const_equation_ptr
      it = ge.begin(); it != ge.end(); ++it){
    cerr << *it << endl;
  }
}


void gauss_eliminator_test::gauss_elimination()
{
  typedef double val_type;
  vector<val_type> node(5,0);
  boost::dynamic_bitset<> node_flag(5);

  jtf::algorithm::gauss_eliminator<val_type> ge(node, node_flag);

  {
    equation<val_type> eq;
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(1), static_cast<val_type>(2)));
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(0), static_cast<val_type>(2)));
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(0), static_cast<val_type>(2)));
    eq.value() = 2;
    ge.add_equation(eq);
  }
  {
    equation<val_type> eq;
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(0), static_cast<val_type>(1)));
    eq.value() = 2;
    ge.add_equation(eq);
  }
  {
    equation<val_type> eq;
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(0), static_cast<val_type>(5)));
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(2), static_cast<val_type>(2)));
    eq.value() = 2;
    ge.add_equation(eq);
  }
  {
    equation<val_type> eq;
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(3), static_cast<val_type>(0.5)));
    eq.add_expression(jtf::algorithm::make_expression(
                        static_cast<size_t>(4), static_cast<val_type>(1)));
    eq.value() = 2;
    ge.add_equation(eq);
  }


  for(size_t t = 0; t < node_flag.size(); ++t){
    if(node_flag[t] == true){
      cerr << "# node " << t << " = " << node[t] << endl;
    }else
      cerr << "# node " << t << " unknown." << endl;
  }

  for(jtf::algorithm::gauss_eliminator<val_type>::const_equation_ptr
      it = ge.begin(); it != ge.end(); ++it){
    if(it->state() == 0) continue;
    if(it->state() == -1){
      std::cerr << "# [conflict] " << endl;
    }
    cerr << *it << endl;
  }

  cerr << "# test file IO." << endl;
  {
    ofstream ofs("test_equations_io", ios::binary);
    ge.save_equations(ofs);
  }
  cerr << "# [io] finish write out" << endl;
  {
    ifstream ifs("test_equations_io", ios::binary);
    ge.load_equations(ifs);
  }
  //  {
  //    equation<val_type> eq;
  //    eq.add_expression(jtf::algorithm::make_expression(
  //                        static_cast<size_t>(0), static_cast<val_type>(1)));
  //    eq.value() = 1;
  //    ge.add_equation(eq);
  //  }

  for(size_t t = 0; t < node_flag.size(); ++t){
    if(node_flag[t] == true){
      cerr << "# node " << t << " = " << node[t] << endl;
    }else
      cerr << "# node " << t << " unknown." << endl;
  }
  for(jtf::algorithm::gauss_eliminator<val_type>::const_equation_ptr
      it = ge.begin(); it != ge.end(); ++it){
    cerr << *it << endl;
  }
}
