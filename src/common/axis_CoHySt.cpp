#include "axis_CoHySt.hpp"



axis_CoHySt::axis_CoHySt(const double &l0_, const double &L_, const int &N_, const int &m_)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoHySt");
	log << "axis_CoHySt(const double &l0_, const double &L_, const int &N_, const int &m_)";
   #endif
	l0 = l0_;
	L = L_;
	N = N_;
	m = m_;
	type_id = 2;
}

axis_CoHySt::axis_CoHySt(const axis_CoHySt &that)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoHySt");
	log << "axis_CoHySt(const axis_CoHySt &that)";
   #endif
	l0 = that.l0;
	L = that.L;
	N = that.N;
	m = that.m;
	type_id = that.type_id;
}

axis_CoHySt::~axis_CoHySt()
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoHySt");
	log << "~axis_CoHySt()";
   #endif
}

axis * axis_CoHySt::clone () const
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoHySt");
	log << "clone() const";
   #endif
	return new axis_CoHySt(*this);
}

axis * axis_CoHySt::create_reciprocal() const
{
  #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoHySt");
	log << "create_reciprocal()";
   #endif

	return new axis_FoHySt(l0,L,N,m);
}

double axis_CoHySt::val_at(const int &index) const
{
	double x_val;
	x_val = l0 + L*double(index)/(N);

	int shift = 0;
	int tmp_i = index;


	if(x_val< -L/2.)
	{
		do
		{
			x_val += L;
			shift -= 1;
		} while(x_val<-L/2.);
	}

	if(x_val>= L/2.)
	{
		do
		{
			x_val -= L;
			shift += 1;
		} while(x_val>=l0+L);
	}


	 x_val =  shift*L + (L/2.)*sinh(2*(x_val)*m/L)/sinh(m);

	return x_val;
}

int axis_CoHySt::index_at(const double &val) const
{
	double x_eq = val;
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


	 x_eq = (L/(2.*m))*asinh(2.*x_eq*sinh(m)/L);
	 int r_int = nearbyint( N*(x_eq-l0)/L );

	 return r_int + N*shift;
}



double axis_CoHySt::operator() (const int &i) const
{
	double x_val;
	x_val = l0 + L*double(i)/(N);

	int shift = 0;
	int tmp_i = i;


	if(x_val< -L/2.)
	{
		do
		{
			x_val += L;
			shift -= 1;
		} while(x_val<-L/2.);
	}

	if(x_val>= L/2.)
	{
		do
		{
			x_val -= L;
			shift += 1;
		} while(x_val>=l0+L);
	}


	 x_val =  shift*L + (L/2.)*sinh(2*(x_val)*m/L)/sinh(m);

	return x_val;
}


int axis_CoHySt::operator() (const double &x) const
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


	 x_eq = (L/(2.*m))*asinh(2.*x_eq*sinh(m)/L);
	 int r_int = floor( N*(x_eq-l0)/L );

	 return r_int + N*shift;
}


int axis_CoHySt::operator() (const double &x, double &lambda) const
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
			x_eq = (L/(2.*m))*asinh(2.*x_eq*sinh(m)/L);
			return floor(N*(x_eq-l0)/L) + N*shift;
		}

	}
	throw("axis::index_at(): something went terribly wrong");
}

