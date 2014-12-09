#include "endianness_test.h"
#include <iostream>
#include <vector>
void Endianness_test::swap_each_word_test(){
    double temp = 1;
    jtf::Endianness::swap_each_word((unsigned char*)&temp,sizeof(double));
    std::cout << temp << std::endl;
    CPPUNIT_ASSERT(temp != 1);
}

void Endianness_test::swap_endian_test(){
    std::vector<size_t> temp(10);
    for(size_t t = 0; t < temp.size(); ++t) temp[t] = t;
    jtf::Endianness::swap_endian(&temp[0],temp.size());
    for(size_t t = 0; t < temp.size(); ++t)
    std::cout << temp[t] << std::endl;
}
