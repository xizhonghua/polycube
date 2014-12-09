#ifndef HJ_OBJ_LIKE_FORMAT_H_
#define HJ_OBJ_LIKE_FORMAT_H_

#include <iostream>
#include <string>

inline int read_head(std::istream &is, size_t &num, const char *name = 0) {
	std::string word;
	is >> word >> word >> num;
	return is.fail() || (name && word != name);
}

template <typename T>
int read_body(std::istream &is, size_t num, size_t ele_per_line, T *body) {
	std::string word;
	for(size_t i = 0; i < num; ++i) {
		is >> word;
		for(size_t ei = 0; ei < ele_per_line; ++ei, ++body)
			is >> *body;
	}
	return is.fail();
}

template <typename Con>
int read_section(std::istream &is, size_t ele_per_line, Con &con, const char *name = 0) {
	size_t num;
	if(read_head(is, num, name))
		return __LINE__;
	con.resize(ele_per_line*num);
	if(read_body(is, num, ele_per_line, &con[0]))
		return __LINE__;
	return 0;
}

template <typename T>
int write_section(std::ostream &os, size_t num, size_t ele_per_line, const T *body,
				  const char *name = "unknown", const char *id = "unknown") {
	os << "# " << name << ' ' << num << '\n';
	for(size_t i = 0; i < num; ++i) {
		os << id;
		for(size_t ei = 0; ei < ele_per_line; ++ei)
			os << ' ' << body[i*ele_per_line+ei];
		os << '\n';
	}
	return os.fail();
}

#endif
