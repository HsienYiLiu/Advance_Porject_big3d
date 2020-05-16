/*
Filename: readFile.cpp
Maintainer: Hsien-Yi Liu
Description: Read 0.off file and save the data into array
Data: Feb 7, 2020
*/

#include <fstream>
#include <sstream>
#include <string>
#include<iostream>
#include <cstdlib>
#include <stdlib.h>

using namespace std;
struct polyhedron{
    int vertex_num;
    int facet_num;
    float **vertices; 
    int **facets;
};

int main()
{
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
    if(count == 1){
        cout<< infile << endl;
    }
    istringstream iss(line);
    float a, b,c,d;
    int ia,ib,ic,id;
    // assign vertices in array
    if(count > 1 && count < 143){
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
  //tesing output
  
  for(int i=0;i<140;i++)
  {
	    cout<< instPoly.vertices[i][0] << " "<< instPoly.vertices[i][1] <<" " << instPoly.vertices[i][2] << " " << endl;
  }
  
  
  cout << instPoly.vertices[instPoly.facets[7][2]][1]<< endl;
  free(instPoly.vertices);
  free(instPoly.facets);
}

