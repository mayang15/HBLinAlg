#ifndef MC_MATRIX
#define MC_MATRIX 1


#include <iostream>
#include <utility>
#include <cmath>
#include "algebra.hpp"

template<typename value_type>
class Matrix {
    
public:
    typedef unsigned int size_type;
    
    Matrix() : dim1(0), dim2(0){}
    
    Matrix(size_type dim1, size_type dim2) : dim1(dim1), dim2(dim2) {
		init_alloc();
    }
    
    Matrix(size_type dim1, size_type dim2, value_type initval) : dim1(dim1), dim2(dim2) {
		init_alloc();
		init_diag(initval);
    }
    
    Matrix(Matrix<value_type> const & IN) : dim1(IN.dim1), dim2(IN.dim2) {
		init_alloc();
		init_array(IN.elems);
    }
    Matrix(Matrix<value_type> const & IN, bool copy) : dim1(IN.dim1), dim2(IN.dim2) {
		if(copy) {
			init_alloc();
			init_array(IN.elems);
		}
		else {
			elems = IN.elems;
			sub = false;
			copied = false;
		}
    }
    
    ~Matrix() {
		if(copied) {
			for(size_type i=0;i<dim1;++i) {
				delete [] elems[i];
			}
			delete [] elems;
		}
		else {
			if(sub) {
				delete [] elems;
			}
		}
    }
    
    std::pair<size_type,size_type> dimension() const {
		return std::pair<size_type,size_type>(dim1,dim2);
    }


	bool is_zero(size_type line)
	{
		for(size_type i=0;i<dim2;++i)
			if(std::abs(elems[line][i]) > eps )
				return false;

		return true;
	}
    
    void swap_cols(size_type x1, size_type x2) {
      value_type *help = elems[x1];
      elems[x1] = elems[x2];
      elems[x2] = help;
    }
    
    Matrix<value_type> trans() const {
		Matrix<value_type> res(dim2,dim1);
		for(size_type i=0;i<dim1;i++) {
			for(size_type j=0;j<dim2;j++) res[j][i] = elems[i][j];
		}
		return res;
    }
    
    value_type* operator[](size_type index) const {
		return elems[index];
    }
    
    value_type & operator()(size_type i, size_type j) const {
		return elems[i][j];
    }
    
    Matrix<value_type>& operator=(Matrix<value_type> IN) {
		size_type min1 = std::min(IN.dim1,dim1), min2 = std::min(IN.dim2,dim2);
		for(size_type i=0;i<min1;++i) {
			for(size_type j=0;j<min2;++j) elems[i][j] = IN.elems[i][j];
		}
		return *this;
    }
    friend std::ostream & operator<<(std::ostream & lhs, Matrix<value_type> const & rhs) {
		for(size_type i=0;i<rhs.dim1;i++) {
			for(size_type j=0;j<rhs.dim2;j++) {
				lhs << rhs.elems[i][j] << " ";
			}
			lhs << "\n";
		}
		lhs << "\n";
		return lhs;
    }
    
    friend bool operator==(Matrix<value_type> const & lhs, Matrix<value_type> const & rhs) {
		for(size_type i=0;i<lhs.dim1;i++)
			for(size_type j=0;j<lhs.dim2;j++) 
				if(lhs[i][j] != rhs[i][j])
					return false;
		return true;
    }
    
    Matrix<value_type> get_sub(std::pair<size_type,size_type> d1,std::pair<size_type,size_type> d2) const {
		Matrix<value_type> res(d1.second-d1.first,d2.second-d2.first, false);
		res.elems = new value_type*[res.dim1];
		for(size_type i=d1.first;i<d1.second;++i) {
			res.elems[i-d1.first] = &(elems[i][d2.first]);
		}
		return res;
    }
    
    friend Matrix<value_type> operator*(Matrix<value_type> const & lhs, Matrix<value_type> const & rhs) {
		return Matrix::strassen_mul(lhs,rhs);
    }
    
    
    friend Matrix<value_type> operator*(value_type sk, Matrix<value_type> const & rhs) {
		Matrix<value_type> res(rhs,true);
		for(size_type i=0;i<rhs.dim1;++i) {
			for(size_type j=0;j<rhs.dim2;++j) 
				res[i][j] *= sk;
		}
		return res;
    }
    friend Matrix<value_type> operator*(Matrix<value_type> const & rhs, value_type sk)
	 {
		Matrix<value_type> res(rhs,true);
		for(size_type i=0;i<rhs.dim1;++i) {
			for(size_type j=0;j<rhs.dim2;++j) 
				res[i][j] *= sk;
		}
		return res;
    }
    
    friend Matrix<value_type> operator+(Matrix<value_type> const & lhs, Matrix<value_type> const & rhs) {
		if(lhs.dim1 < rhs.dim1 || lhs.dim2 < rhs.dim2) {
			Matrix<value_type> res(rhs);
			for(size_type i=0;i<lhs.dim1;++i) {
				for(size_type j=0;j<lhs.dim2-(lhs.dim2%4);j+=4) 
				{
					res.elems[i][j] += lhs.elems[i][j];
					res.elems[i][j+1] += lhs.elems[i][j+1];
					res.elems[i][j+2] += lhs.elems[i][j+2];
					res.elems[i][j+3] += lhs.elems[i][j+3];
				}
				for(size_type j=lhs.dim2-(lhs.dim2%4);j<lhs.dim2;++j) 
				{
					res.elems[i][j] += lhs.elems[i][j];
				}
			}
			return res;
		}
		else {
			Matrix<value_type> res(rhs);
			for(size_type i=0;i<rhs.dim1;++i) {
				for(size_type j=0;j<rhs.dim2-(rhs.dim2%4);j+=4) 
				{
					res.elems[i][j] += lhs.elems[i][j];
					res.elems[i][j+1] += lhs.elems[i][j+1];
					res.elems[i][j+2] += lhs.elems[i][j+2];
					res.elems[i][j+3] += lhs.elems[i][j+3];
				}
				for(size_type j=rhs.dim2-(rhs.dim2%4);j<rhs.dim2;++j) 
				{
					res.elems[i][j] += lhs.elems[i][j];
				}
			}
			return res;
		}
    }
    friend Matrix<value_type> operator-(Matrix<value_type> const & lhs, Matrix<value_type> const & rhs) {
		if(lhs.dim1 < rhs.dim1 || lhs.dim2 < rhs.dim2) {
			Matrix<value_type> res(-1*rhs);
			for(size_type i=0;i<lhs.dim1;++i) {
				for(size_type j=0;j<lhs.dim2;++j) res.elems[i][j] += lhs.elems[i][j];
			}
			return res;
		}
		else {
			Matrix<value_type> res(lhs);
			for(size_type i=0;i<rhs.dim1;++i) {
				for(size_type j=0;j<rhs.dim2;++j) res.elems[i][j] -= rhs.elems[i][j];
			}
			return res;
		}
    }
#ifdef MC_VECTOR
    Matrix(Vector<value_type> vin) : dim1(vin.size()), dim2(1) {
      init_alloc();
      for(size_type i=0;i<dim1;i++) 
		  elems[i][0] = vin[i];
    }
    friend Matrix<value_type> operator*(Vector<value_type> const & lhs, Matrix<value_type> const & rhs) {
		Vector<value_type> res(lhs.dim1);
		for(size_type i=0;i<lhs.dim1;i++) {
			res[i] = 0.0;
			for(size_type j=0;j<lhs.dim2;j++) {
				res[i] += lhs[j]*rhs[j][i];
			}
		}
		return res;
    }
    
    friend Vector<value_type> operator*(Matrix<value_type> const & lhs, Vector<value_type> const & rhs) {
		Vector<value_type> res(lhs.dim2);
		for(size_type i=0;i<lhs.dim1;i++) {
			res[i] = 0.0;
			for(size_type j=0;j<lhs.dim2;j++) {
				res[i] += lhs[i][j]*rhs[j];
			}
		}
		return res;
    }
    void assigne(size_type index, Vector<value_type> const & in)
	{
		for(size_type i=0;i<std::min(dim2,in.size());++i)
			elems[index][i] = in[i];
	}
#endif 
    
    void assigne(size_type index, value_type * const & in)
	{
		for(size_type i=0;i<dim2;++i)
			elems[index][i] = in[i];
	}
    
    
private:
    
    value_type **elems;
    bool copied,sub;
    const size_type dim1,dim2;
#if __cplusplus > 199711L
	static constexpr value_type eps = 1.0e-10;
#else
	static const value_type eps = 1.0e-10;
#endif
    
    Matrix(size_type dim1, size_type dim2, bool copy) : dim1(dim1), dim2(dim2) {
		if(copy) {
			init_alloc();
		}
		else {
			copied = false;
			sub = true;
		}
    }
    
    Matrix(Matrix<value_type> const & IN, 
			 std::pair<size_type,size_type> xdim, 
			 std::pair<size_type,size_type> ydim, 
			 size_type dim1, 
			 size_type dim2, 
			 value_type inval) : dim1(dim1), dim2(dim2) {
		init_alloc();
		unsigned x = xdim.second-xdim.first, y = ydim.second-ydim.first;
		for(size_type i=0;i<dim1;i++) {
			for(size_type j=0;j<dim2;j++) {
				if(x > i && y > j) {
					elems[i][j] = IN[xdim.first+i][ydim.first+j];
				}
				else {
					elems[i][j] = inval;
				}
			}
		}    
    }
    
    void init_alloc() {
		copied = true;
		sub = false;
		elems = new value_type*[dim1];
		for(size_type i=0;i<dim1;++i) {
			elems[i] = new value_type[dim2];
			__builtin_assume_aligned(elems[i],16);
		}
    }
    
    void init_diag(value_type initval) {
		for(size_type i=0;i<dim1;++i) {
			for(size_type j=0;j<dim2;++j) {
			if(i==j) elems[i][i] = initval;
			else {
				elems[i][j] = 0.0;
			}
			}
		}
    }
    
    void init_array(value_type ** const IN) {
		for(size_type i=0;i<dim1;++i) {
			for(size_type j=0;j<dim2;++j) {
				elems[i][j] = IN[i][j];
			}
		}
    } 
    
	
    static Matrix<value_type> simple_mul(Matrix<value_type> const & lhs, Matrix<value_type> const & rhs) {
	//        std::cout << "SIMPLE_MUL(" << lhs.dim1 << "x" << lhs.dim2 << "," << rhs.dim1 << "x" << rhs.dim2 << ")\n";
		Matrix<value_type> res(lhs.dim1,rhs.dim2,(value_type)0.0); 
		for(size_type i=0;i<lhs.dim1;++i) {
			for(size_type k=0;k<lhs.dim2-(lhs.dim2%4);k+=4) {
				for(size_type j=0;j<rhs.dim2;j++) {
					res.elems[i][j] += lhs.elems[i][k]*rhs.elems[k][j] + lhs.elems[i][k+1]*rhs.elems[k+1][j] + lhs.elems[i][k+2]*rhs.elems[k+2][j] + lhs.elems[i][k+3]*rhs.elems[k+3][j];
				}
			}
			for(size_type k=lhs.dim2-(lhs.dim2%4);k<lhs.dim2;k++) {
				for(size_type j=0;j<rhs.dim2;j++) {
					res.elems[i][j] += lhs.elems[i][k]*rhs.elems[k][j];
				}
			}
		}
		
		return res;
    }
	 
    static Matrix<value_type> strassen_mul(Matrix<value_type> const & lhs, Matrix<value_type> const & rhs) {
		// 1x1 Matrix
		unsigned block = 256;
		if((lhs.dim1+lhs.dim2+rhs.dim1+rhs.dim2)/4 <= block || std::min(lhs.dim1,lhs.dim2)<2 || std::min(rhs.dim1,rhs.dim2)<2) return simple_mul(lhs,rhs);
		// ungerade? 
		Matrix<value_type> A(lhs,false);
		Matrix<value_type> B(rhs,false);
		std::pair<size_type,size_type> a1,a2,a_1,a_2,b1,b2,b_1,b_2;
		
		if(lhs.dim1%2==1) {
			a1 = std::make_pair(0,A.dim1/2+1);
			a2 = std::make_pair(A.dim1/2+1,A.dim1);
		}
		else {
			a1 = std::make_pair(0,A.dim1/2);
			a2 = std::make_pair(A.dim1/2,A.dim1);
		}
		if(lhs.dim2%2==1) {
			a_1 = std::make_pair(0,A.dim2/2+1);
			a_2 = std::make_pair(A.dim2/2+1,A.dim2);
		}
		else {
			a_1 = std::make_pair(0,A.dim2/2);
			a_2 = std::make_pair(A.dim2/2,A.dim2);
		}
		if(rhs.dim1%2==1) {
			b1 = std::make_pair(0,B.dim1/2+1);
			b2 = std::make_pair(B.dim1/2+1,B.dim1);
		}
		else {
			b1 = std::make_pair(0,B.dim1/2);
			b2 = std::make_pair(B.dim1/2,B.dim1);
		}
		if(rhs.dim2%2==1) {
			b_1 = std::make_pair(0,B.dim2/2+1);
			b_2 = std::make_pair(B.dim2/2+1,B.dim2);
		}
		else {
			b_1 = std::make_pair(0,B.dim2/2);
			b_2 = std::make_pair(B.dim2/2,B.dim2);
		}
		
		// aufteilung
		
		// Formeln mit rekursion
		Matrix<value_type> P1(strassen_mul(A.get_sub(a1,a_1)+A.get_sub(a2,a_2),B.get_sub(b1,b_1)+B.get_sub(b2,b_2)));
		
		Matrix<value_type> P2(strassen_mul(A.get_sub(a2,a_1)+A.get_sub(a2,a_2),B.get_sub(b1,b_1)));
				
		Matrix<value_type> P3(strassen_mul(A.get_sub(a1,a_1),B.get_sub(b1,b_2)-B.get_sub(b2,b_2)));
			
		Matrix<value_type> P4(strassen_mul(Matrix<value_type>(A,a2,a_2,a2.second-a2.first,a_2.second-a_2.first+(lhs.dim2%2),0.0),B.get_sub(b2,b_1)-B.get_sub(b1,b_1)));
		
		
		Matrix<value_type> P5(strassen_mul(A.get_sub(a1,a_1)+A.get_sub(a1,a_2),Matrix<value_type>(B,b2,b_2,b2.second-b2.first+(rhs.dim1%2),b_2.second-b_2.first,0.0)));
		
				
		Matrix<value_type> P6(strassen_mul(A.get_sub(a2,a_1)-A.get_sub(a1,a_1),B.get_sub(b1,b_1)+B.get_sub(b1,b_2)));
				
		Matrix<value_type> P7(strassen_mul(A.get_sub(a1,a_2)-A.get_sub(a2,a_2),B.get_sub(b2,b_1)+B.get_sub(b2,b_2)));
		
		Matrix<value_type> res(A.dim1,B.dim2);
		res.get_sub(a1,b_1) = (P1+P4-P5+P7);
		res.get_sub(a1,b_2) = P3+P5;
		res.get_sub(a2,b_1) = P2+P4;
		res.get_sub(a2,b_2) = P1+P3-P2+P6;
		
	// 	delete P4_1;
	// 	delete P5_2;
		return res;
    }
};


template<typename value_type>
Matrix<value_type> diag(unsigned dim, value_type val = 1.0)
{
	Matrix<value_type> res(dim,dim);
	for(unsigned i=0;i<dim;++i)
		res[i][i] = val;
	
	return res;
}


#endif //MATRIX2D-Macro
