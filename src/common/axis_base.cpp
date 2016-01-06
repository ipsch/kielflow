#include "axis_base.hpp"

axis::axis()
{
   #ifdef _MY_VERBOSE_TEDIOUS
	logger log("axis");
	log << "axis()";
   #endif
	e_i = e_x;
	type_id = -1;
	N = 0;
	l0 = 0.;
	L = 0.;
	m = 0;
}

double axis::dS(const int &index) const
{
	std::cout << "ERROR: axis is using fallback dS method\n";
	std::cout << "inherited from class axis! this shouldn't happen!\n";
	throw("ERROR: axis is using default dS method");
	return 0;
}


axis::~axis()
{
   #ifdef _MY_VERBOSE_TEDIOUS
	logger log("axis");
	log << "~axis()";
   #endif
}


/*
int axis::index_at(const double &val) const
{


}
*/

/*
double& axis::operator() (const int &i)
{
	throw("not implemented: double& axis::operator() (const int &i)");
}


double  axis::operator() (const int &i) const
{
	throw("not implemented: double  axis::operator() (const int &i) const");
}


int& axis::operator() (const double &x)
{
	throw("not implemented: int& axis::operator() (const double &x)");
}


int  axis::operator() (const double &x) const
{
	throw("not implemented: int  axis::operator() (const double &x) const");
}



int  axis::operator() (const double &x, double &lambda) const
{

	throw("axis::index_at(): something went terribly wrong");
}


int  axis::operator() (const double &x, double &lambda) const
{
	double x_eq = x;
	int shift = 0;

	if(x_eq<l0)
	{
		do
		{
			shift -= 1;
			x_eq = x_eq + L;
		} while(x_eq<l0);
	}

	if(x>=l0+L)
	{
		do
		{
			shift += 1;
			x_eq = x_eq - L;
		} while(x_eq>=l0+L);
	}

	for(int i=0; i<N-1; ++i)
	{
		if( (val_at(i)<=x_eq) && (x_eq < val_at(i+1))  )
		{
			lambda = (val_at(i+1)-x)/(val_at(i+1) - val_at(i));
			return i + N*shift;
		}

	}
	throw("axis::index_at(): something went terribly wrong");
}
*/
