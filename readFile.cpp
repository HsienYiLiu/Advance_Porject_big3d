#include <fstream>
#include <sstream>
#include <string>
#include<iostream>
#include <cstdlib>

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
  // Build for facets
  instPoly.facets = (int**)malloc(276 * sizeof(int*));
  for(int i = 0; i < 276; i++)
      instPoly.facets[i] = (int*)malloc(3 * sizeof(int));
  //Build for vertex
  instPoly.vertices = (float**)malloc(149 * sizeof(int*));
  for(int i = 0; i < 141; i++)
      instPoly.vertices[i] = (float*)malloc(3 * sizeof(int));;
  int count = 0;
  string line;
  while(getline(infile, line))
  {
    count++;
    istringstream iss(line);
    float a, b,c,d;
    int ia,ib,ic,id;
    if(count > 1 && count < 143){
        if (!(infile >> a >> b >> c )) { 
            break; 
        } 
        //cout << a << endl;
        instPoly.vertices[count-2][0] = a;
        instPoly.vertices[count-2][1] = b;
        instPoly.vertices[count-2][2] = c;
    }else if(count > 142){
        if (!(infile >> ia >> ib >> ic >>id)) { 
            break; 
        } 
        //cout << b << endl;
        instPoly.facets[count-143][0] = ib;
        instPoly.facets[count-143][1] = ic;
        instPoly.facets[count-143][2] = id;
    }
  }

  for(int i=0;i<140;i++)
  {
	    cout<< instPoly.vertices[i][0] << " "<< instPoly.vertices[i][1] <<" " << instPoly.vertices[i][2] << " " << endl;
  }
}

