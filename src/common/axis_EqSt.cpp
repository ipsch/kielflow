#include "axis_EqSt.hpp"


// ############ equidistant coordinate Axis ###################################

axis_CoEqSt::axis_CoEqSt(const double &l0_, const double &L_, const int &N_)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("axis_CoEqSt");
	my_log << "axis_CoEqSt(const double &l0_, const double &L_, const int &N_)";
   #endif
	l0 = l0_;
	L = L_;
	N = N_;
	type_id = 1;
}

axis_CoEqSt::axis_CoEqSt(const axis_CoEqSt &that)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("axis_CoEqSt");
	my_log << "axis_CoEqSt(const axis_CoEqSt &that)";
   #endif
	l0 = that.l0;
	L = that.L;
	N = that.N;
	type_id = that.type_id;
}

axis_CoEqSt::~axis_CoEqSt()
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("axis_CoEqSt");
	my_log << "~axis_CoEqSt()";
   #endif
}

axis * axis_CoEqSt::clone () const
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("axis_CoEqSt");
	my_log << "clone()";
   #endif
	return new axis_CoEqSt(*this);
}

void axis_CoEqSt::resize(const int &N_)
{
	this->N = N_;
}

double axis_CoEqSt::val_at(const int &index) const
{
	return l0 + index*(L/N);
}

double axis_CoEqSt::k_val_at(const int &index) const
{
	double pi = acos(-1.);
	return ( (index<=N/2) ? (2.*pi*index/L) : (2.*pi*(index-N)/L) );
}

int axis_CoEqSt::index_at(const double &val) const
{
	return N*(val-l0)/L;
}

int axis_CoEqSt::k_index_at(const double &val) const
{
	std::cout << "axis_CoEqSt::k_index_at(const double &val) const not implemented\n";
	throw("axis_CoEqSt::k_val_at(const int &index) const not implemented\n");
	return 0;
}


double axis_CoEqSt::S(const double &x) const
{
	return x;
}

double axis_CoEqSt::dS(const int &index) const
{
	return 1.;
}



// functors
/*
double axis_CoEqSt::operator() (const int &i) const
{
	return l0 + i*(L/N);
}

int axis_CoEqSt::operator() (const double &x) const
{
	// secure access to axis value
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

	if(x_eq>=l0+L)
	{
		do
		{
			shift += 1;
			x_eq = x_eq - L;

		} while(x_eq>=l0+L);
	}

	return floor(N*(x_eq-l0)/L) + N*shift;
}

int axis_CoEqSt::operator() (const double &x, double &lambda) const
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

	if(x_eq>=l0+L)
	{
		do
		{
			shift += 1;
			x_eq = x_eq - L;
		} while(x_eq>=l0+L);
	}

	for(int i=0; i<N; ++i)
	{
		if( (operator()(i)<=x_eq) && (x_eq < operator()(i+1))  )
		{
			lambda = (operator()(i+1) - x_eq) /
					 (operator()(i+1) - operator()(i));
			return floor(N*(x_eq-l0)/L) + N*shift;
		}

	}
	throw("axis::index_at(): something went terribly wrong");
}

*/

