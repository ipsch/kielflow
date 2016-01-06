#include "newton_method.hpp"

newton_method::newton_method(double (&function)(double), double (&derivative)(double), const double &accuracy) :
	f(function), df(derivative), eps(accuracy)
{

}

double newton_method::solve(const double &x_i) const
{
	double eps_dump[20];
    try
    {
    	double x_ip1;

    	for(int i=0; i<20; ++i)
    	{
    		double dx = f(x_i)/df(x_i);
    		x_ip1 = x_i - dx;


    		eps_dump[i] = eps;
    		if( fabs(dx) < eps  ) return x_ip1;

    	}
    	throw 1;
    }
    catch(int e)
    {
    	std::cout << "An exception occurred. Exception Nr. " << e << std::endl;
    	for(int i=0; i<10; ++i)
    		std::cout << eps_dump[i] << std::endl;
    }

	return x_i;
}
