#include <iostream>
#include <ctime>
#include "Vector.hpp"

using namespace std;

int main(void) {
  unsigned dim = 10000000, start, end;
  Vector<double> v1(dim), v2(dim), v3(dim);
  
  for(unsigned i=0;i<dim;i++) {
    v1[i] = i*(-1.0)+1.0;
    v2[i] = i;
    if(i%2==1) v2[i] *= -1.0;
  }
  
//   cout << "v1:\n" << v1 << "\n";
//   cout << "v2:\n" << v2 << "\n";
  cout << "START:\n";
  start = clock();
  v3 = v1*v2;
  end = clock();
  cout << "END:\n";
  cout << "Time for dim=" << dim << " needed: " << (double)(end-start)/CLOCKS_PER_SEC << "\n";
  return 0;
}
