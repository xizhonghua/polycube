#ifndef HJ_OBJ_LIKE_FORMAT_MAT_H_
#define HJ_OBJ_LIKE_FORMAT_MAT_H_

#include <zjucad/matrix/matrix.h>
#include <iostream>

#include "../common/obj_like_format.h"

template <typename T>
int read_section(std::istream &is, size_t ele_per_line,
				 zjucad::matrix::matrix<T> &con, const char *name = 0) {
	size_t num;
	if(read_head(is, num, name))
		return __LINE__;
	con.resize(ele_per_line, num);
	if(read_body(is, num, ele_per_line, &con[0]))
		return __LINE__;
	return 0;
}

template <typename T>
int write_section(std::ostream &os, const zjucad::matrix::matrix<T> &con,
				  const char *name = "unknown", const char *id = "unknown") {
	return write_section(os, con.size(2), con.size(1), &con[0], name, id);
}

#endif
