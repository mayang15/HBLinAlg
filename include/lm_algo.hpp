#ifndef MC_LM_ALGO_HPP
#define MC_LM_ALGO_HPP

#include "Vector.hpp"
#include "Matrix.hpp"
#include "lin_solve.hpp"
#include "ParameterFunction.hpp"
#include <vector>

typedef unsigned size_type;


template<typename value_type, typename value_base = double>
class LM_Solver
{
public:
	LM_Solver(ParameterFunction<value_type> & in, std::vector<Vector<value_type> > const & xs, std::vector<Vector<value_type> > const & ys)
	{
		x_samps = xs;
		y_samps = ys;
		
		pf = &in;
		beta1 = 0.3;
		beta2 = 0.7;
		
		mu0 = 1;
		
		tol = 1.0e-12;
	}
	
	void solve()
	{
		Vector<value_type> F(x_samps.size()*pf->y_dim()), sk(pf->size()*pf->y_dim()+pf->size(),0.0), 
						   b(x_samps.size()*pf->y_dim()), y(x_samps.size()*pf->y_dim(),0);
		Matrix<value_type> dF(x_samps.size()*pf->y_dim(),pf->size()), E(pf->size(),pf->size(),1.0), A(x_samps.size()*pf->y_dim()+pf->size(),pf->size());
		sk[0] = 10;
		value_base mu = mu0, eps_mu;
		for(size_type i=0;i<y.size();++i)
			for(size_type j=0;j<pf->y_dim();++j)
				y[i*pf->y_dim()+j] = y_samps[i][j];
		while(sk.norm()>tol)
		{
			funcVal(F,dF);
			A = dF;
			for(size_type i=0;i<pf->size();++i)
				A[x_samps.size()*pf->y_dim()+i][i] = mu;
			b = y-F;
			sk = b/A;
			eps_mu = testUpdate(F,dF,sk);
			
			if(eps_mu <= beta1)
			{
				mu *= 2.0;
			}
			else
			{
				if(eps_mu >= beta2)
				{
					mu *= 0.5;
				}
				else
				updatePF(sk);
			}
		}
		
	}
	
	void setBetas(value_type b1, value_type b2)
	{
		beta1 = b1;
		beta2 = b2;
	}
	
	void setTolerance(value_type t)
	{
		tol = t;
	}
	
	void setMu0(value_type mu)
	{
		mu0 = mu;
	}

private:
	std::vector<Vector<value_type> > x_samps, y_samps;
	value_base beta1, beta2, tol, mu0;
	
	ParameterFunction<value_type> * pf;
	
	void funcVal(Vector<value_type> & F, Matrix<value_type> & dF) const
	{
		Matrix<value_type> dftemp(pf->y_dim(),pf->size());
		Vector<value_type> ftemp;
		for(size_type i = 0;i<x_samps.size();++i)
		{
			ftemp = pf->fx(x_samps[i]);
			dftemp = pf->dfa(x_samps[i]);
			for(size_type j=0;j<pf->y_dim();++j) 
			{
				F[i*pf->y_dim()+j] = ftemp[j];
				dF.assigne(i*pf->y_dim()+j,dftemp[j]); 
			}
		}
	}
	
	void funcVal(Vector<value_type> & F, Vector<value_type> const & a) const
	{
		Vector<value_type> ftemp;
		for(size_type i = 0;i<x_samps.size();++i)
		{
			ftemp = pf->fxa(x_samps[i],a);
			for(size_type j=0;j<pf->y_dim();++j) 
			{
				F[i*pf->y_dim()+j] = ftemp[j];
			}
		}
	}
	
	void updatePF(Vector<value_type> const & a)
	{
		for(size_type i=0;i<pf->size();++i)
			(*pf)[i] += a[i];
	}
	
	value_type testUpdate(Vector<value_type> const & F, Matrix<value_type> dF, Vector<value_type> const & a) const
	{
		value_base normF, normF_new, normF_df, res;
		Vector<value_type> F_new(F.size());
		
		normF = F.norm();
		normF *= normF;
		
		funcVal(F_new,a+pf->A());
		normF_new =F_new.norm();
		normF_new *= normF_new;
		
		normF_df = (F+dF*a).norm();
		normF_df *= normF_df;
		
		res = (normF - normF_new)/(normF - normF_df);
		
		return res;
	}

	
};


#endif
