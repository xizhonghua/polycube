#ifndef JTF_TIMER_H
#define JTF_TIMER_H

#include <time.h>

class timer{
public:
  timer():begin_(0), end_(0){}
  void start(){ begin_ = clock();}
  void finish(){end_ = clock();}
  long result()const{return (end_ - begin_)/1000;}
  long result_c()const{return (end_ - begin_);}
private:
  long begin_, end_;
};

#endif //  JTF_TIMER_H
