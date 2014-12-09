#ifndef READ_TET_H
#define READ_TET_H

#include "../../common/def.h"

namespace lq {

  class read_tet {

 public:

    void load_data(const char *path, matrixd *node,
                   matrixst *tet, matrixst *tri);
    
  };
}

#endif
