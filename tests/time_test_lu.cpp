#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "Vector.hpp"
#include "Matrix2D.hpp"
#include "MatrixSP.hpp"
#include "lu_factors.hpp"


typedef long double value_type;

value_type rand_number() {
  srand(clock()%1000);
  unsigned bkomma = 1000;
  unsigned neg = 1000;
  unsigned pos = 1000;
  return ((value_type)(rand()%((neg+pos)*bkomma)))/bkomma - neg;
}

int main() {
  unsigned startdim = 100, dist = 100, count = 1, end, start;
  double time;
  for(unsigned t=0;t<count;t++) {
    unsigned dim = startdim + t*dist; 
    Matrix2D<value_type> A(dim,dim),test(dim,dim);
    std::vector<Matrix2D<value_type> > lus;
    
    for(unsigned i=0;i<dim;i++) {
      for(unsigned j=0;j<dim;j++) {
	A[i][j] = rand_number();
      }
    }
    start = clock();
    lus = lu_factors(A);
    end = clock();
    
    test = lus[0]*lus[1]*lus[2];
    value_type sum = 0.0, max = 0.0;
    for(unsigned i=0;i<dim;i++) {
      for(unsigned j=0;j<dim;j++) {
	sum += std::abs(test[i][j]-A[i][j]);
	if(std::abs(test[i][j]-A[i][j]) >max) max = std::abs(test[i][j]-A[i][j]);
      }
    }
    
    std::cout << startdim+t*dist << "	" << (double)(end-start)/CLOCKS_PER_SEC << "	" << sum << "	" << max <<"\n";
  }
  return 0;
}