#ifndef NUMERIC_CHECK_H
#define NUMERIC_CHECK_H

#include <boost/iterator/iterator_traits.hpp>
#include <boost/type_traits.hpp>
//#include <memory>
#include <stack>
template <typename T>
struct Iterator_Traits
{
  typedef typename T::value_type value_type;
};

template <typename type>
struct Iterator_Traits<type*>
{
  typedef type value_type;
};


//! @brief check whether the usinged add will result in overflow
template<typename input_iterator>
int is_unsinged_add_overflow(input_iterator begin,
                             input_iterator last)
{
  typedef typename Iterator_Traits<input_iterator>::value_type value_type;
  if(!(boost::is_same<value_type,size_t>::value)
      && !(boost::is_same<value_type,unsigned int>::value)
      && !(boost::is_same<value_type,unsigned long>::value)
      && !(boost::is_same<value_type,unsigned short>::value)){
    return -1;
  }
  std::stack<value_type> temp_sum;
  input_iterator it = begin;
  while(it != last){
    value_type a = *it++;
    if(it == last) {
      temp_sum.push(a);
      break;
    }
    value_type b = *it++;
    value_type sum = a + b;
    if(sum < a || sum < b) return 0; // overflow
    temp_sum.push(sum);
  }
  while(!temp_sum.empty()){
    if(temp_sum.size() == 1) return 1;// no overflow
    std::stack<value_type> forwad_sum;
    while(!temp_sum.empty()){
      value_type a = temp_sum.top();
      temp_sum.pop();
      if(temp_sum.empty()) forwad_sum.push(a);
      else{
        value_type b = temp_sum.top();
        temp_sum.pop();
        value_type sum = a + b;
        if(sum < a || sum < b) return 0; // overflow
        forwad_sum.push(sum);
      }
    }
   swap(temp_sum,forwad_sum);
  }
}

#endif
