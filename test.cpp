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

struct Vector {
  double x, y, z;


  Vector operator-(Vector p) const {
    return Vector{x - p.x, y - p.y, z - p.z};
  }

  Vector cross(Vector p) const {
    return Vector{
      y * p.z - p.y * z,
      z * p.x - p.z * x,
      x * p.y - p.x * y
    };
  }

  double dot(Vector p) const {
    return x * p.x + y * p.y + z * p.z;
  }

  double norm() const {
    return sqrt(x*x + y*y + z*z);
  }
};
using Point = Vector;
struct Face {
  vector<Point> v;

  Vector normal() const {
    assert(v.size() > 2);
    Vector dir1 = v[1] - v[0];
    Vector dir2 = v[2] - v[0];
    Vector n  = dir1.cross(dir2);
    double d = n.norm();
    return Vector{n.x / d, n.y / d, n.z / d};
  }
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
    istringstream iss(line);
    float a, b,c;
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


  Point p;
  Point c;
  vector<Face> f;
  Face f2;
  vector<Face> cube{ // faces with 4 points, last point is ignored
  /*
    Face{{Point{0,0,0}, Point{1,0,0}, Point{1,0,1}, Point{0,0,1}}}, // front
    Face{{Point{0,1,0}, Point{0,1,1}, Point{1,1,1}, Point{1,1,0}}}, // back
    Face{{Point{0,0,0}, Point{0,0,1}, Point{0,1,1}, Point{0,1,0}}}, // left
    Face{{Point{1,0,0}, Point{1,1,0}, Point{1,1,1}, Point{1,0,1}}}, // right
    Face{{Point{0,0,1}, Point{1,0,1}, Point{1,1,1}, Point{0,1,1}}}, // top
    Face{{Point{0,0,0}, Point{0,1,0}, Point{1,1,0}, Point{1,0,0}}}, // bottom*/
  };
  //get total Face
  for(int i = 0; i < 276; i++){
      vector<Point> new_a = {};
      for(int j = 0; j < 3; j++){
          p.x = instPoly.vertices[instPoly.facets[i][0]][j];
          p.y = instPoly.vertices[instPoly.facets[i][1]][j];
          p.z = instPoly.vertices[instPoly.facets[i][2]][j];
          new_a.push_back(p);
      }
      cube.push_back(Face{new_a});
  }
  p.x = 9;
  p.y = 9;
  p.z = 9;
  c.x = 100;
  c.y = 100;
  c.z = 100;
  vector<Point> new_a = {};
  new_a.push_back(p);
  new_a.push_back(c);
  //v = new_a;
  //v = new_b;
  //f[0].v = new_a;
  //f[1].v = new_b;
  cube.push_back(Face{new_a});
  

  //cube.push_back(f{{p}});

// Displaying the 2D vector 
  for (auto p : cube){

    //cout << p.x << endl;
	cout<< p.v[1].x << endl;
  }
  free(instPoly.vertices);
  free(instPoly.facets);

}

