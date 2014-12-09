#ifndef HJ_HEXGEN_COMMON_CONFIG_H_
#define HJ_HEXGEN_COMMON_CONFIG_H_

#ifdef HEXGEN_COMMON_EXPORT
#define HEXGEN_COMMON_API __declspec(dllexport)
#else
#define HEXGEN_COMMON_API
#endif

#if defined(__GNUC__)
# define GCC_VERSION (__GNUC__ * 10000 + __GNUC__MINOR__ * 100)
# if (GCC_VERSION >= 40600)
#include <cstddef>
# endif
#endif

#endif
