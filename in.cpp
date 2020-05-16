#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <stdlib.h>

using namespace std;
int sign(int x){
    /*
    return 1 if x is positive, -1 negative and 0 if it's zero
    */
   if (x > 0)
   {
       return 1;
   }else if(x < 0){
       return -1;
   }else{
       return 0;
   } 
}
int vertex_sign(array P, array O){
    int result = sign(P[0] - O[0]) or sign(P[1] - O[1]) or sign(P[2] - O[2])
    if not result:
        
}