#ifndef RANDOMH
#define RANDOMH

#include <stdio.h> 
#include <iostream> 
#include <time.h>

inline double random_double() {
    
    //利用c++内置随机数生成器，除以随机数的最大值，生成[0,1]区间的浮点数
    return ((double)rand()/(double)RAND_MAX);
}
#endif