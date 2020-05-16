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
struct polyhedron{
    int vertex_num;
    int facet_num;
    float **vertices; 
    int **facets;
};

int main(int argc, char* argv[]) {

  struct polyhedron instPoly;
  std::ifstream infile("0.off");
  // Building array for facets
  instPoly.facets = (int**)malloc(276 * sizeof(int*));
  for(int i = 0; i < 276; i++)
      instPoly.facets[i] = (int*)malloc(3 * sizeof(int));
  //Building array vertex
  instPoly.vertices = (float**)malloc(149 * sizeof(int*));
  for(int i = 0; i < 141; i++)
      instPoly.vertices[i] = (float*)malloc(3 * sizeof(int));;
  int count = 0;
  string line;
  while(getline(infile, line))
  {
    count++;
    istringstream iss(line);
    float a, b,c;
    int ia,ib,ic,id;
    // assign vertices in array
    if(count > 1 && count < 142){
        if (!(infile >> a >> b >> c )) { 
            break; 
        } 
        instPoly.vertices[count-2][0] = a;
        instPoly.vertices[count-2][1] = b;
        instPoly.vertices[count-2][2] = c;
    // assign facets in array
    }else if(count > 142){
        if (!(infile >> ia >> ib >> ic >>id)) { 
            break; 
        } 
        instPoly.facets[count-143][0] = ib;
        instPoly.facets[count-143][1] = ic;
        instPoly.facets[count-143][2] = id;
    }
  }

  float x0 = 0;
  float x1 = 0;
  float x2 = 0;
  float y0 = 0;
  float y1 = 0;
  float y2 = 0;
  float z0 = 0;
  float z1 = 0;
  float z2 = 0;
  float check = 0;
  //get total Face
  for(int i = 0; i < 276; i++){
    x0 = instPoly.vertices[instPoly.facets[i][0]][0];
    x1 = instPoly.vertices[instPoly.facets[i][0]][1];
    x2 = instPoly.vertices[instPoly.facets[i][0]][2];
    y0 = instPoly.vertices[instPoly.facets[i][1]][0];
    y1 = instPoly.vertices[instPoly.facets[i][1]][1];
    y2 = instPoly.vertices[instPoly.facets[i][1]][2];
    z0 = instPoly.vertices[instPoly.facets[i][2]][0];
    z1 = instPoly.vertices[instPoly.facets[i][2]][1];
    z2 = instPoly.vertices[instPoly.facets[i][2]][2];
    // cout << i <<"   "<< instPoly.facets[i][0] <<"   "<< instPoly.facets[i][1] <<"      " <<instPoly.facets[i][2] <<"    "<< "wise lock" << endl;
    check = x0*(y1*z2 - z1*y2) + y0*(z1*x2 - x1*z2)+z0*(x1*y2 - y1*x2);
    if( check == 0){
        cout << i << "Same line" << endl;
    }else if(check > 0){
        cout << "counter-wise lock" << endl;
    }else{
        cout << i <<"   "<< instPoly.facets[i][0] <<"   "<< instPoly.facets[i][1] <<"      " << "wise lock" << endl;
    }
    // new_a.push_back(init);
    //cout << instPoly.facets[i][j] << "  "<<  "" <<" "<< init.x <<" "<< init.y<< " " << init.z << endl;
      
      
  }

  return 0;
}

