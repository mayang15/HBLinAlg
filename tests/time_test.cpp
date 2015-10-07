#include "Matrix2D.hpp"
#include <iostream>
#include <ctime>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;
typedef double value_type;

long double flops(unsigned dim) {
  long double n = dim;
  if(dim<=256) return (2*n*n*n-n*n)*1.0e-09;
  if(dim%2==1) return (18*((n+1)*(n+1))/4*1.0e-09 + 7*flops((dim+1)/2));
  return (18*(n*n)/4*1.0e-09 + 7*flops(dim/2));
}
	  

int main() {
    unsigned dim1, dim2, dim3, startdim, count, dist, start, end;
    count = 1;
    dist = 0;
    startdim =2048-20;
    
    double time[count], mtltime[count];
    for(unsigned k=0;k<count;k++) {
	dim1 = startdim + k*dist;
	dim2 = dim1;
	dim3 = dim1;
	
	cout << "Dim: " << dim1 << "\n";
	cout << "own start\n";
	Matrix2D<value_type> A(dim1,dim2), B(dim2,dim3), C(dim1,dim3);
	for(unsigned i=0;i<dim1;++i) {
	    for(unsigned j=0;j<dim2;++j) {
		A[i][j] = (i * 17)%89 + (j*23)%101;
	    }
	}
	
	for(unsigned i=0;i<dim2;++i) {
	    for(unsigned j=0;j<dim3;++j) {
		B[i][j] = (i * 19)%83 + (j*29)%113;
	    }
	}
	
	cout << "start measure:\n";
	start = clock();
 	C = A*B;
	end = clock();
	cout << "end measure:\n";
	time[k] = (value_type)(end - start)/(value_type)CLOCKS_PER_SEC;
	
	cout << "MTL start\n";
	mtl::dense2D<value_type> Am(dim1,dim2),Bm(dim2,dim3),Cm(dim1,dim3);
	for(unsigned i=0;i<dim1;++i) {
	    for(unsigned j=0;j<dim2;++j) {
		Am[i][j] = (i * 17)%89  + (j*23)%101;
	    }
	}
	
	for(unsigned i=0;i<dim2;++i) {
	    for(unsigned j=0;j<dim3;++j) {
		Bm[i][j] = (i * 19)%83 + (j*29)%113;
	    }
	}
	
	start = clock();
  	Cm = Am*Bm;
	end = clock();
	mtltime[k] = (double)(end - start)/(double)CLOCKS_PER_SEC;
	
    }
    cout << "Dim	Mine	MTL	GFLOPS	GTFLOPS\n";
    for(unsigned k=0;k<count;k++) {
	cout << startdim+k*dist << "	" << time[k] << "	" << mtltime[k] << "	" << flops(startdim+k*dist)/time[k] << "	" << flops(startdim+k*dist)/mtltime[k] << "\n";
    }
    return 0;
}
