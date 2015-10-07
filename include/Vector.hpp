#ifndef MC_VECTOR
#define MC_VECTOR 1

#include <iostream>
#include <cmath>
#include <utility>
#include <algorithm>
#include "algebra.hpp"

template<typename value_type>
class Vector {
public:
  typedef unsigned int size_type;
  
  Vector() : dim(0) {}
  
  Vector(size_type dim) : dim(dim) {
    init_alloc();
  }
  
  Vector(size_type dim, value_type initval) : dim(dim) {
    init_alloc();
    for(unsigned i=0;i<dim;i++) elems[i] = initval;
  }
  
  Vector(Vector<value_type> const & IN) : dim(IN.dim) {
    init_alloc();
    init_array(IN.elems);
  }
  
  ~Vector() {
    delete [] elems;
  }
  
  size_type size() const {
    return dim;
  }
  
  Vector<value_type> & operator=(Vector<value_type> const & IN) {
    if(dim!=IN.dim) {
      if(dim!=0) delete[] elems;
      dim = IN.dim;
      init_alloc();
    }
    for(size_type i=0;i<dim;i++) {
      elems[i] = IN.elems[i];
    }
    return *this;
  }
#if __cplusplus > 199711L
  Vector<value_type> & operator=(std::initializer_list<value_type> const in)
  {
	std::copy_n( in.begin(), in.size(), elems );
	return *this;
  }
#else
  Vector<value_type> & operator=(std::vector<value_type> const in)
  {
	for(unsigned i=0;i<in.size();++i)
		elems[i] = in[i];
	return *this;
  }
#endif
  
  value_type & operator[](size_type index) const {
    return elems[index];
  }
  
  friend value_type operator*(Vector<value_type> const & v1, Vector<value_type> const & v2) {
    value_type sum = 0;
    for(size_type i=0;i<v1.dim;i++) sum += v1[i]*v2[i];
    return sum;
  }
  
  friend Vector<value_type> operator+(Vector<value_type> const & v1, Vector<value_type> const & v2) {
    Vector<value_type> res(v1.dim);
    for(size_type i=0;i<v1.dim;i++) res[i] = v1[i]+v2[i];
    return res;
  }

  friend Vector<value_type> operator-(Vector<value_type> const & v1, Vector<value_type> const & v2) {
    Vector<value_type> res(v1.dim);
    for(size_type i=0;i<v1.dim;i++) res[i] = v1[i]-v2[i];
    return res;
  }
  
  friend Vector<value_type> operator+(Vector<value_type> const & v1, value_type sk) {
    Vector<value_type> res(v1.dim);
    for(size_type i=0;i<v1.dim;i++) res[i] += sk;
    return res;
  }
  
  friend Vector<value_type> operator+(value_type sk, Vector<value_type> const & v1) {
    Vector<value_type> res(v1.dim);
    for(size_type i=0;i<v1.dim;i++) res[i] += sk;
    return res;
  }
  
  friend Vector<value_type> operator*(Vector<value_type> const & v, value_type sk) {
    Vector<value_type> res(v.dim);
    for(size_type i=0;i<v.dim;i++) res = sk*v[i];
    return res;
  }
  
  friend Vector<value_type> operator*(value_type sk,Vector<value_type> const & v) {
    Vector<value_type> res(v.dim);
    for(size_type i=0;i<v.dim;i++) res = sk*v[i];
    return res;
  }
  
  friend Vector<value_type> operator/(Vector<value_type> const & v,value_type sk) {
    Vector<value_type> res(v.dim);
    for(size_type i=0;i<v.dim;i++) res = v[i]/sk;
    return res;
  }
  
  friend std::ostream & operator<<(std::ostream & lhs, Vector<value_type> const & v) {
	lhs << "(";
    for(size_type i=0;i<v.dim;i++) lhs << v[i] << " ";
    return lhs << ")";
  }
  
  friend bool operator==(Vector<value_type> const & lhs, Vector<value_type> const & rhs) {
    for(size_type i=0;i<lhs.dim;i++) if(lhs[i]!=rhs[i]) return false;
    return true;
  }
  
  value_type norm() const
  {
	  value_type res = 0;
	  for(size_type i=0;i<dim;++i)
	  {
		  res += elems[i]*elems[i];
	  }
	  return sqrt(res);
  }
  
private:
  size_type dim;
  value_type *elems;
  
  void init_alloc() {
    elems = new value_type[dim];
	//__builtin_assume_aligned(elems,16);
  }
  
  void init_array(value_type *in) {
    for(size_type i=0;i<dim;i++) elems[i] = in[i];
  }
    
};

#endif // VECTOR-Macro
