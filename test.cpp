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
bool isInConvexPoly(Point const& p, std::vector<Face> const& fs) {
  for (Face const& f : fs) {
    Vector p2f = f.v[0] - p;        
    double d = p2f.dot(f.normal());
    d /= p2f.norm();                 // for numeric stability
    //cout << f.v[3].x <<"  " <<p.x << "  "<< d << endl;
    //cout << d << endl;
    constexpr double bound = -1e-15; // use 1e15 to exclude boundaries
    if (d < bound)
      return false;
  }

  return true;
}
int main(int argc, char* argv[]) {
  assert(argc == 3+1);
  char* end;
  Point p;

  p.x = std::strtod(argv[1], &end);
  p.y = std::strtod(argv[2], &end);
  p.z = std::strtod(argv[3], &end);
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

  Point init;
  //Point c;
  vector<Face> f;
  Face f2;
  vector<Face> cube{};
  //get total Face
  for(int i = 0; i < 276; i++){
      vector<Point> new_a = {};
      for(int j = 0; j < 3; j++){
          init.x = instPoly.vertices[instPoly.facets[i][j]][0];
          init.y = instPoly.vertices[instPoly.facets[i][j]][1];
          init.z = instPoly.vertices[instPoly.facets[i][j]][2];
          new_a.push_back(init);
          // cout << instPoly.facets[i][j] << "  "<<  "" <<" "<< init.x <<" "<< init.y<< " " << init.z << endl;
      }
      
      cube.push_back(Face{new_a});
  }
  /*
  for (auto p : cube){

    //cout << p.x << endl;
	cout<< p.v[0].x << p.v[0].y << p.v[0].z << endl;
  }*/
  cout << (isInConvexPoly(p, cube) ? "inside" : "outside") << endl;

  
  free(instPoly.vertices);
  free(instPoly.facets);
  return 0;
}

