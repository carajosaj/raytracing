#ifndef RANDOMH
#define RANDOMH

#include <stdio.h> 
#include <iostream> 
#include <time.h>

inline double random_double() {
    
    //����c++�������������������������������ֵ������[0,1]����ĸ�����
    return ((double)rand()/(double)RAND_MAX);
}
#endif