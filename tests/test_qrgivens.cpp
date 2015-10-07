#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <utility>

#include "Vector.hpp"
#include "Matrix2D.hpp"
#include "qr_givens.hpp"

typedef double value_type;
unsigned seed = 5;

value_type rand_number() {
  unsigned bkomma = 10;
  unsigned neg = 1;
  unsigned pos = 1;
  value_type res = ((value_type)(rand()%((neg+pos)*bkomma)))/bkomma - neg;
  return res;
}

int main() {
  unsigned dim1 = 2000, dim2 = 2000;
  Matrix2D<value_type> A(std::max(dim1,dim2),dim2),test(std::max(dim1,dim2),dim2);
  
  for(unsigned i=0;i<dim1;i++) {
    for(unsigned j=0;j<dim2;j++) {
      A[i][j] = rand_number();
    }
  }

if(dim1 < dim2)
{
	for(unsigned i=dim1;i<std::max(dim1,dim2);i++) {
    	for(unsigned j=0;j<dim2;j++) {
    	  A[i][j] = i*A[0][j];
    	}
  	}
}
unsigned start_t, end_t;
  start_t = clock();
  QR_Givens<value_type> sol(A); 
  end_t = clock();
  // R = G2*G1*A
  // (G2*G1)T * R = A
  
  
  test = sol.get_q()*sol.get_r();
 // std::cout << "A:\n" << A << sol.get_q() << sol.get_r();
  
  value_type sum = 0.0;
  for(unsigned i=0;i<dim1;i++) {
    for(unsigned j=0;j<std::max(dim1,dim2);j++) {
      sum += std::abs(test[i][j]-A[i][j]);
    }
  }
//	std::cout << "diff mat:\n" << test-A << "\n";
  std::cout << "Differenz: " << sum << "\n";
std::cout << "Time: " << (double)(end_t-start_t)/CLOCKS_PER_SEC << "\n";
  
  return 0;
}
