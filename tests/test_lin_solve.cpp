#include "lin_solve.hpp"
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <utility>

typedef double value_type;

value_type rand_number() {
  unsigned bkomma = 10;
  unsigned neg = 1;
  unsigned pos = 1;
  return ((value_type)(rand()%((neg+pos)*bkomma)))/bkomma - neg;
}

int main() {
	unsigned dim1 = 3, dim2 = 2;
	Matrix2D<value_type> A(dim1,dim2);
	Vector<value_type> b(dim1), x(dim2);
	  
	for(unsigned i=0;i<dim1;i++) {
		for(unsigned j=0;j<dim2;j++) {
		  	A[i][j] = rand_number();
		}
		b[i] = rand_number();
	}

	x = b/A;

	std::cout << "Defekt: \n" << A*x-b << "\n";

	return 0;
}
