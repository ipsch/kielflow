#include "axis_CoEqSt.hpp"


// ############ equidistant coordinate Axis ###################################

axis_CoEqSt::axis_CoEqSt(const double &l0_, const double &L_, const int &N_)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoEqSt");
	log << "axis_CoEqSt(const double &l0_, const double &L_, const int &N_)";
   #endif
	l0 = l0_;
	L = L_;
	N = N_;
	type_id = 1;
}

axis_CoEqSt::axis_CoEqSt(const axis_CoEqSt &that)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoEqSt");
	log << "axis_CoEqSt(const axis_CoEqSt &that)";
   #endif
	l0 = that.l0;
	L = that.L;
	N = that.N;
	type_id = that.type_id;
}

axis_CoEqSt::~axis_CoEqSt()
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoEqSt");
	log << "~axis_CoEqSt()";
   #endif
}

axis * axis_CoEqSt::clone () const
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoEqSt");
	log << "clone()";
   #endif
	return new axis_CoEqSt(*this);
}

axis * axis_CoEqSt::create_reciprocal() const
{
  #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoEqSt");
	log << "create_reciprocal()";
   #endif

	return new axis_FoEqSt(l0,L,N);
}

double axis_CoEqSt::val_at(const int &index) const
{
	return l0 + index*(L/N);
}

int axis_CoEqSt::index_at(const double &val) const
{
	return N*(val-l0)/L;
}




// functors

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

