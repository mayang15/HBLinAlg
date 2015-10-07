#include "lm_algo.hpp"
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <utility>

typedef double value_type;

template<typename value_type>
class ParameterPolynome : public ParameterFunction<value_type>
{
public:
	
	ParameterPolynome(Vector<value_type> const & in)
	{
		a = in;
	}
	
	~ParameterPolynome() { }
	
	Vector<value_type> fx(Vector<value_type> const & x) const
	{
		Vector<value_type> res(1,0.0);
		res[0] = a[0];
		for(size_type i=1;i<a.size();++i)
		{
			res[0] *= x[0];
			res[0] += a[i];
		}
		return res;
	}
	
	Matrix<value_type> dfa(Vector<value_type> const & x) const
	{
		Matrix<value_type> res(1,a.size(),0.0);
		Vector<value_type> a_cp(a.size(),0.0);
		a_cp = a;
		const value_type h = 1.0e-8;
		for(size_type i=0;i<res.dimension().second;++i)
		{
			a_cp[i] += h;
			res[0][i] = (fxa(x,a_cp)-fx(x))[0]/h;
			a_cp[i] -= h;
		}
		return res;
	}
	
	Vector<value_type> fxa(Vector<value_type> const & x, Vector<value_type> const a) const
	{
		Vector<value_type> res(1,0.0);
		res[0] = a[0];
		for(size_type i=1;i<a.size();++i)
		{
			res[0] *= x[0];
			res[0] += a[i];
		}
		return res;
	}
	
	Vector<value_type> A() const
	{
		return a;
	}
	
	value_type operator[](size_type a_index) const
	{
		return a[a_index];
	}
	
	value_type & operator[](size_type a_index)
	{
		return a[a_index];		
	}
	
	size_type size() const
	{
		return a.size();
	}
	
	size_type x_dim() const
	{
		return 1;
	}
	
	size_type y_dim() const
	{
		return 1;
	}
	
private:
	Vector<value_type> a;
	
};

value_type rand_number() {
  unsigned bkomma = 10;
  unsigned neg = 1;
  unsigned pos = 1;
  return ((value_type)(rand()%((neg+pos)*bkomma)))/bkomma - neg;
}

int main() {
	std::vector<Vector<value_type> > x(3,Vector<value_type>(1,0.0)), y(3,Vector<value_type>(1,0.0));
	
	std::cout << "here-1\n";
	x[0] = {-2}; y[0] = {0};
	x[1] = {-1}; y[1] = {-1};
	x[2] = {1}; y[2] = {0}; 
	std::cout << "here0\n";
	Vector<value_type> a(3);
	for(unsigned i=0;i<a.size();++i)
		a[i] = i;
	std::cout << "here1\n";
	ParameterPolynome<value_type> pp(a);
	std::cout << "here2\n";
	LM_Solver<value_type> lms(pp,x,y);
	std::cout << "here3\n";
	lms.solve();
	std::cout << "here4\n";
	
	std::cout << "polynom durch: \n";
	for(unsigned i=0;i<x.size();++i)
		std::cout << "(" << x[i] << "," << y[i] << ") ";
	std::cout << "\ngegeben durch:\n";
	
	std::cout << pp.A() << " ";
	std::cout << "\n";
	
	std::cout << "x1:" << x[0] << " f:=" << pp.fx(x[0]) << "\n";
	std::cout << "x2:" << x[1] << " f:=" << pp.fx(x[1]) << "\n";
	std::cout << "x3:" << x[2] << " f:=" << pp.fx(x[2]) << "\n";

	return 0;
}
