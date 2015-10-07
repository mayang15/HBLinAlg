#ifndef MC_LIN_SOLVE_HPP
#define MC_LIN_SOLVE_HPP

#include <vector>
#include "Vector.hpp"
#include "Matrix.hpp"
#include "qr_givens.hpp"

typedef unsigned size_type;


template<typename value_type>
Vector<value_type> lin_solve(Matrix<value_type> const & A, Vector<value_type> const & bin) {
	
	//Create A=QR 
	
	if(A.dimension().first > A.dimension().second)
		return lin_solve(A.trans()*A,A.trans()*bin);
	Matrix<value_type> Ac(std::max(A.dimension().first,A.dimension().second),A.dimension().second);
	Ac = A;
	QR_Givens<value_type> qr(Ac);
	/*std::cout << "A:\n " << A << " b: " << bin << "\n\n" ;
	std::cout << "Q: " << qr.get_q() << " R: " << qr.get_r() << "\n";
	std::cout << "a=QR: " << qr.get_q() * qr.get_r() << "\n";*/

	// calc b=QT*bin
	Vector<value_type> b = qr.get_q().trans()*bin;

	//std::cout << "Qb: \n" << b << "\n\n";

	// backtrack Rx=b
	Vector<value_type> x(A.dimension().second);
		// If check if zerolines for overdefined
	size_type dim1 = A.dimension().first, dim2 = A.dimension().second;
			
	for(size_type i=dim1-1;i!=0;--i)
		if(qr.get_r().is_zero(i))
		{
			x[i] = 1.0;
			//std::cout << "found zeroline\n";
		}
		else
		{
			dim1 = i+1;
			break;
		}
	

	if(dim1 < dim2)
	{
		//=> extend R to nxn and set k xs to 1 and beging at n-k line in backtracking
		for(size_type i = dim1; i<dim2;++i)
		{
			x[i] = 1.0;
		}
	}

	// real backtracking:
	for(size_type i = dim1-1;;--i)
	{
		x[i] = b[i];
		for(size_type j = i+1;j<dim2;++j)
		{
			x[i] -= qr.get_r()[i][j]*x[j];
		}
		x[i]/= qr.get_r()[i][i];
		if(i==0)
			break;
	}
	//std::cout << "x: \n" << x << "\n";
	return x;
}


template<typename value_type>
Vector<value_type> operator/(Vector<value_type> const & bin,Matrix<value_type> const & A)
{
	return lin_solve(A,bin);
}

#endif
