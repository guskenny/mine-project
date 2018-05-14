// generates a random graph file for testing
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

int main (int argc, char *argv[]) {
  ofstream ofile;
  if (argc != 4)
    cout << "usage: " << argv[0] << " <filename> <width> <depth>\n";
  else
    ofile.open (argv[1]);

  int width = atoi(argv[2]);
  int depth = atoi(argv[3]);
  int n_nodes = width*depth;
  int n_arcs = (width*3 - 2)*(depth-1);

  ofile << "d " << width << " " << depth << endl;
  ofile << "p max " << n_nodes << " " << n_arcs << endl;

 
  float MIN = -3.0;
  float MAX = 6.0;
  srand(time(NULL));
 
  for (int i=0; i < n_nodes; i++){
    float profit = 0;
    while (profit == 0){
//      profit = MIN + (rand() % (int)(MAX-MIN + 1));
      profit = ((MAX-MIN)*((float)rand()/RAND_MAX))+MIN;
    }
    ofile << "s " << i << " " << profit << std::endl;
  }
 
  for (int i=width;i<n_nodes;i++){
    int k,l;
    if (i % width == 0){
      k = 0;
      l = 2;
    }
    else if (i % width == (width-1)){
      k = -1;
      l = 1;
    }
    else{
      k = -1;
      l = 2;
    }
    
    for (int j=k; j<l; j++)
      ofile << "a " << i << " " << (i-width)+j << " 1000\n";
    }

  ofile.close();
  return 0;
}
