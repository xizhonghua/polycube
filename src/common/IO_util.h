#ifndef HJ_IO_UTIL_H_
#define HJ_IO_UTIL_H_

#include <iostream>
#include <cassert>

template <typename T>
inline void write_val(std::ostream &os, const T &val) {
	os.write(reinterpret_cast<const char *>(&val), sizeof(T));
}

template <typename T>
inline void write_arr(std::ostream &os, const T *val, size_t size) {
	assert(val);
	os.write(reinterpret_cast<const char *>(val), sizeof(T)*size);
}

template <typename T>
inline void read_val(std::istream &is, T &val) {
	is.read(reinterpret_cast<char *>(&val), sizeof(T));
}

template <typename T>
inline void read_arr(std::istream &is, T *val, size_t size) {
	assert(val);
	is.read(reinterpret_cast<char *>(val), sizeof(T)*size);
}

#endif
