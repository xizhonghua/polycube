#ifndef BYTE_CONFIG_H
#define BYTE_CONFIG_H
#include <vector>
#include <algorithm>
#include <assert.h>

namespace jtf{
class Endianness{
public:
    static bool is_little_endian(void) {
        unsigned char endian_test[2] = {1,0};
        short x = *(short*) endian_test;
        return (x == 1);
    }

    template<typename IDX_TYPE>
    static void swap_each_word(unsigned char * p, IDX_TYPE char_in_word) {
        assert(p);
        std::reverse(p,p + char_in_word);
    }

    template<typename VAL_TYPE, typename IDX_TYPE>
    static void swap_endian(VAL_TYPE *array, IDX_TYPE size) {
        assert(array);
        IDX_TYPE num = sizeof(VAL_TYPE);
        unsigned char *p = reinterpret_cast<unsigned char *>(array);
        for(IDX_TYPE t = 0; t < size; ++t){
            swap_each_word(p,num);
            p += num;
        }
    }
};
}
#endif // BYTE_CONFIG_H
