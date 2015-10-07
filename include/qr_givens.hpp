#ifndef MC_QR_GIVENS
#define MC_QR_GIVENS

#include <utility>
#include <cmath>

#include "Matrix.hpp"

template<typename value_type>
class QR_Givens {
  
  typedef typename Matrix<value_type>::size_type size_type;
  
public:
  QR_Givens(Matrix<value_type> const & IN) : R(IN),G(Matrix<value_type>(IN.dimension().first,IN.dimension().first,(value_type)1.0)) {
    eps = 1.0e-10;
    calc();
  }
  
  void set_rotation(value_type a1, value_type a2) {
    value_type t, temp;
    if(a1>=a2) {
      t = a2/std::abs(a1);
      temp = sqrt(1+std::abs(t)*std::abs(t));
      s = t/temp;
      c = (a1/std::abs(a1))/temp;
    }
    else {
      t = a1/std::abs(a2);
      temp = sqrt(1+std::abs(t)*std::abs(t));
      s = (a2/std::abs(a2))/temp;
      c = t/temp;
    }
  }
  
  void calc() {
    value_type temp;
     for(size_type j=0;j<R.dimension().second;j++) {
       for(size_type i=R.dimension().first-1;i>j;i--) {
	if(std::abs(R[i][j])>eps) {
	  set_rotation(R[j][j],R[i][j]);
	  for(size_type k=0;k<R.dimension().second;k++) {
	    temp = R[j][k];
	    R[j][k] = c*R[j][k] + s*R[i][k];
	    R[i][k] = -s*temp + c*R[i][k];
	    temp = G[j][k];
	    G[j][k] = c*G[j][k] + s*G[i][k];
	    G[i][k] = -s*temp + c*G[i][k];
	  }
	}
 	R[i][j] = 0.0;
      }
    }
    G = G.trans(); 
  }
  
  Matrix<value_type> get_q() {
    return G;
  }
  
  Matrix<value_type> get_r() {
    return R;
  }
  
  std::pair<Matrix<value_type>, Matrix<value_type> > get() {
    return std::pair<Matrix<value_type>, Matrix<value_type> >(G,R);
  }
  
private:
  Matrix<value_type> G,R;
  value_type c,s, eps;
  
};

template<typename value_type>
std::pair<Matrix<value_type>, Matrix<value_type> > qr_givens(Matrix<value_type> const & IN) {
  QR_Givens<value_type> res(IN);
  return res.get();
}

#endif

