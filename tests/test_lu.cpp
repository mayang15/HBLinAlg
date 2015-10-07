#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>

#include "Vector.hpp"
#include "Matrix2D.hpp"
#include "lu_factors.hpp"


typedef double value_type;

value_type rand_number() {
  unsigned bkomma = 1;
  unsigned neg = 20;
  unsigned pos = 20;
  return ((value_type)(rand()%((neg+pos)*bkomma)))/bkomma - neg;
}

int main() {
  unsigned dim = 1000;
  Matrix2D<value_type> A(dim,dim),test(dim,dim);
  
  for(unsigned i=0;i<dim;i++) {
    for(unsigned j=0;j<dim;j++) {
      A[i][j] = rand_number();
    }
  }
  
  std::vector<Matrix2D<value_type> > luf;
  luf = lu_factors(A);
//   std::cout << "P:\n" << luf[0] << "L:\n" << luf[1] << "R:\n" << luf[2];
  
  test =luf[0]*luf[1]*luf[2];
  value_type sum = 0.0;
  for(unsigned i=0;i<dim;i++) {
    for(unsigned j=0;j<dim;j++) sum += std::abs(test[i][j]-A[i][j]);
  }
  std::cout << "Differenz: " << sum << "\n";
  
  return 0;
}