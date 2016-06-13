#ifndef NEWTON_METHOD_HPP_
#define NEWTON_METHOD_HPP_

#define FAILSAVE

#include <iostream>
#include <cmath>


class newton_method
{
public :
	newton_method(const double &Tolerance = 0.05, unsigned const int &MaxIterations = 20) :
		eps(Tolerance), Nmax(MaxIterations)
	{	}

	template<typename Func_f, typename Func_df>
	double solve(double x, Func_f f, Func_df df)
	{
		double omega = 1.;
		double g = 0.8;
		double dx = 0.;
		double dx_alt = f(x)/df(x);
       #if defined(FAILSAVE)
		double x0 = x;
       #endif

		int i=0;
		do
		{
			x -= omega*dx;

			dx = f(x)/df(x);
			omega = fabs(dx)>fabs(dx_alt) ?
					g*omega : (1.-g)*omega + g*1.;
			dx_alt = dx;
			++i;
		} while ( fabs(dx)>eps && i<Nmax );

       #if defined(FAILSAVE)
    	if(fabs(dx)<eps) // everthing ok
    		return x-dx;

    	if(dx!=dx)
    	{
    		std::cout << "WARNING: Newton-method went nan" << std::endl;
    		return x;
    	}

    	if(fabs(f(x0)/df(x0)) < fabs(f(x)/df(x)))
    	{
    		std::cout << "WARNING: Newton-method diverged" << std::endl;
    		return x0;
    	}

    	std::cout << "WARNING: didn't reach convergence" << std::endl;
       #endif

		return x-dx;
	}


private :
	const double eps;
	const unsigned int Nmax;
};




#endif
