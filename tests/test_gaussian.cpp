#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>

#include "Vector.hpp"
#include "Matrix2D.hpp"
#include "gaussian_elimination.hpp"


typedef double value_type;

value_type rand_number() {
  unsigned bkomma = 1000;
  unsigned neg = 1000;
  unsigned pos = 1000;
  return ((value_type)(rand()%((neg+pos)*bkomma)))/bkomma - neg;
}

int main() {
  unsigned dim = 400;
  Matrix2D<value_type> A(dim,dim),test(dim,1);
  Vector<value_type> b(dim), x(dim);
  
  for(unsigned i=0;i<dim;i++) {
    for(unsigned j=0;j<dim;j++) {
      A[i][j] = rand_number();
    }
    b[i] = rand_number();
  }
  
  x = gaussian_elimination(A,b);
  test = A*x;
  value_type sum = 0.0;
  for(unsigned i=0;i<dim;i++) {
    sum += std::abs(test[i][0]-b[i]);
  }
  std::cout << "Differenz: " << sum << "\n";
  
  return 0;
}