#ifndef MC_PARAMETER_FUNCTION_HPP
#define MC_PARAMETER_FUNCTION_HPP

#include <vector>
#include <cmath>
#include "Vector.hpp"
#include "Matrix.hpp"

typedef unsigned size_type;


/**
 * \brief For LM-Method, takes an vector X and its defining Parameters A which will be determined be LM
 */
template<typename value_type>
class ParameterFunction 
{
public:
	
	virtual ~ParameterFunction() { }
	
	virtual Vector<value_type> fx(Vector<value_type> const & x) const = 0;
	
	virtual Vector<value_type> fxa(Vector<value_type> const & x, Vector<value_type> const a) const = 0;
	
	virtual Matrix<value_type> dfa(Vector<value_type> const & x) const = 0;
	
	virtual value_type operator[](size_type a_index) const = 0;
	
	virtual value_type & operator[](size_type a_index) = 0;
	
	virtual Vector<value_type> A() const = 0;
	
	virtual size_type size() const = 0;
	
	virtual size_type x_dim() const = 0;
	
	virtual size_type y_dim() const = 0;

private:

};


#endif
