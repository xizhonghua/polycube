#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

#include <zjucad/matrix/matrix.h>
#include "../common/def.h"
#include "../numeric/util.h"
using namespace zjucad::matrix;

int show_sh_func(int argc, char *argv[])
{
  matrixd v(3), v2;
  while(!cin.eof()) {
      cin >> v[0] >> v[1] >> v[2];
      if(cin.fail())
        break;
      v2 = pow(v, 2);
      const double f = v2[0]*v2[1] + v2[1]*v2[2] + v2[2]*v2[0];
      cout << "f " << f << "\n";
      const double sh = -(15/(2*sqrt(My_PI()))*f+8)*sqrt(7.0);
      cout << "sh " << sh << "\n";
    }
  return 0;
}
