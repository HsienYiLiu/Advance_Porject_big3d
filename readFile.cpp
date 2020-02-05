#include <fstream>
#include <sstream>
#include <string>
#include<iostream>

using namespace std;


int main()
{
  float data[500][3];
  std::ifstream infile("0.off");

  string line;
  int i=0;

  while(getline(infile, line))
  {
    //istringstream iss(line);
    float a, b,c;
    if (!(infile >> a >> b >> c)) { 
        break; 
    } // error
	data[i][0]= a;
	data[i][1]= b;
	data[i][2]= c;
	i++;
	//cout << a << b << c << endl;
    // process pair (a,b)
  }

  for(i=0;i<10;i++)
  {
   cout << data[i][0] << " " << data[i][1] << " " <<  data[i][2] << endl;
  }
}

