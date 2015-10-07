#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>

#include "Vector.hpp"
#include "Matrix2D.hpp"
#include "lin_solve_lu.hpp"


typedef double value_type;

value_type rand_number() {
  srand(clock()%1000);
  unsigned bkomma = 1000;
  unsigned neg = 1000;
  unsigned pos = 1000;
  return ((value_type)(rand()%((neg+pos)*bkomma)))/bkomma - neg;
}

int main() {
  unsigned startdim = 1000, dist = 500, count = 2, end, start;
  double time;
  for(unsigned t=0;t<count;t++) {
    unsigned dim = startdim + t*dist; 
    Matrix2D<value_type> A(dim,dim),test(dim,1);
    Vector<value_type> b(dim), x(dim);
    
    for(unsigned i=0;i<dim;i++) {
      for(unsigned j=0;j<dim;j++) {
	A[i][j] = rand_number();
      }
      b[i] = rand_number();
    }
    start = clock();
    x = lin_solve_lu(A,b);
    end = clock();
    
    test = A*x;
//     std::cout << "b=\n" << b << "test=\n" << test;
    value_type sum = 0.0, max = 0.0;
    for(unsigned i=0;i<dim;i++) {
      sum += std::abs(test[i][0]-b[i]);
      if(std::abs(test[i][0]-b[i]) > max) max = std::abs(test[i][0]-b[i]);
    }
    
    std::cout << startdim+t*dist << "	" << (double)(end-start)/CLOCKS_PER_SEC << "	" << sum << "	" << max <<"\n";
  }
  return 0;
}