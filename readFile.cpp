#include <fstream>
#include <sstream>
#include <string>
#include<iostream>

using namespace std;
struct polyhedron{
    int vertex_num;
    int facet_num;
    vector<vector<float> > vertices; 
    vector<vector<int> > facets;
};

int main()
{
  vector<float> data[3];
  vector<int> dataInt[3];
  struct polyhedron instPoly;
  std::ifstream infile("0.off");
  bool isInt=false;;
  string line;
  int i=0;

  while(getline(infile, line))
  {
    //istringstream iss(line);
    float a, b,c;
    int ia,ib,ic;
    if (!(infile >> a >> b >> c)) { 
        break; 
    } // error
	
	if(a % 1 == 0)
	{
	 isInt = true;
	 dataInt[0]= ia;
	 dataInt[1]= ib;
	 dataInt[2]= ic;
	instPoly.facets.push_back(data);
	}
       else
       {
	data[0]= a;
	data[1]= b;
	data[2]= c;
	instPoly.vertices.push_back(data);
       }
	
	
	//cout << a << b << c << endl;
    // process pair (a,b)
  }

  for(i=0;i<10;i++)
  {
	cout<< instPoly.facets[i][0]] << endl;
  }
}

