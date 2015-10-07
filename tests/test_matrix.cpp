#include "Matrix2D_block.hpp"
#include <iostream>
#include <ctime>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;
typedef double value_type;

int main() {
    
    
    unsigned start1, start2, end, dim1 = 3000, dim2 = 3000, dim3 = 3000;
    
    
    start1 = clock();
    
    Matrix<value_type> A(dim1,dim2), B(dim2,dim3), C(dim1,dim3);
    
    for(unsigned i=0;i<dim1;++i) {
	for(unsigned j=0;j<dim2;++j) {
	    A[i][j] = (i * 17)%89 * 0.4 + (j*23)%101 * 0.3;
	}
    }
    
    for(unsigned i=0;i<dim2;++i) {
	for(unsigned j=0;j<dim3;++j) {
	    B[i][j] = (i * 19)%83 * 0.3 + (j*29)%113 * 0.2;
	}
    }
    
    start2 = clock();
    C = A*B;
    end = clock();
    cout << "Zeit Alloc: " << (value_type)(start2 - start1)/(value_type)CLOCKS_PER_SEC << "\n";
    cout << "Zeit MUL: " << (value_type)(end - start2)/(value_type)CLOCKS_PER_SEC << "\n";
//     cout << "A:\n" << A;
//     cout << "B:\n" << B;
//          cout << "C:\n" << C;
// 	 cout << "Csimp:\n" << Matrix<value_type>::simple_mul(A,B);
    
//     exit(1);
    start1 = clock();
    mtl::matrix::dense2D<value_type> Am(dim1,dim2),Bm(dim2,dim3),Cm(dim1,dim3);
    
    for(unsigned i=0;i<dim1;++i) {
	for(unsigned j=0;j<dim2;++j) {
	    Am[i][j] = (i * 17)%89 * 0.4 + (j*23)%101 * 0.3;
	}
    }
    
    for(unsigned i=0;i<dim2;++i) {
	for(unsigned j=0;j<dim3;++j) {
	    Bm[i][j] = (i * 19)%83 * 0.3 + (j*29)%113 * 0.2;
	}
    }
    
    start2 = clock();
    Cm = Am*Bm;
    end = clock();
    cout << "Zeit MTL Alloc: " << (value_type)(start2 - start1)/(value_type)CLOCKS_PER_SEC << "\n";
    cout << "Zeit MTL MUL: " << (value_type)(end - start2)/(value_type)CLOCKS_PER_SEC << "\n";
//         cout << "Cm:\n" << Cm;
    bool bitch = true;
    unsigned cou = 0;
//     cout << "D:\n";
    for(int i=0;i<dim1 && bitch;i++) {
	for(int j=0;j<dim3 && bitch;j++) {
	    
	    if(std::abs(Cm[i][j]-C[i][j]) > 1.0e-5) {
//  		cout << std::abs(Cm[i][j]-C[i][j]) << " ";
		cou++;
	    }
// 	    else cout << 0 << " ";
	}
// 	cout << "\n";
    }
    if(cou == 0) cout << "succer\n";
    else cout << "looser " << cou << "/" << dim1*dim3 << "\n";
    
    return 0;
}