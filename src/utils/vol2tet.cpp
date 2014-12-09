#include <iostream>
#include <fstream>

#include <jtflib/mesh/io.h>

#include "../tetmesh/hex_io.h"
#include "../common/def.h"

using namespace std;

int vol2tet(const char * filename,
	    matrixst & tet,
	    matrixd & node)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open vol file." << endl;
      return __LINE__;
    }
  
  string line;
  while(!ifs.eof()){
      getline(ifs, line);
      if(line  == "dimension") {
          size_t d = 0 ;
          ifs >> d;
          if(d != 3){
              cerr << "# [error] dimesnsion is " << d << ", not 3." << endl;
              return __LINE__;
            }
        }

      if(line == "geomtype"){
          size_t geo_t = 0;
          ifs >> geo_t;
          if(geo_t != 11){
              cerr << "# [error] geomtye is " << geo_t << ", not 11." << endl;
              return __LINE__;
            }
        }


      if(line[0] == '#'){
          ifs >> line;
          if(line == "volumeelements"){
              size_t tet_num = 0;
              ifs >> tet_num;
              tet.resize(4, tet_num);
              size_t trash;
              for(size_t ti = 0; ti < tet_num; ++ti){
                  ifs >> trash >> trash;
                  for(size_t di = 0; di < 4; ++di)
                    ifs >> tet(di,ti);
                }
              tet -= 1;
            }

          if(line == "points"){
              size_t point_num;
              ifs >> point_num;
              node.resize(3, point_num);
              for(size_t pi = 0; pi < point_num; ++pi)
                for(size_t di = 0; di < 3; ++di)
                  ifs >> node(di,pi);
            }
        }
    }
  return 0;
}

int vol2tet(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [usage] vol2tet vol tet." << endl;
      return __LINE__;
    }

  matrixst tet;
  matrixd node;
  vol2tet(argv[1], tet, node);
  orient_tet(node, tet);
  if(jtf::mesh::tet_mesh_write_to_zjumat(argv[2], &node, &tet))
    return __LINE__;
  cerr << "# [info] success." << endl;
  return 0;
}
