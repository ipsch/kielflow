#include "axis_FoEqSt.hpp"


axis_FoEqSt::axis_FoEqSt(const double &l0_, const double &L_, const int &N_)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoEqSt");
	log << "axis_FoEqSt(const double &l0_, const double &L_, const int &N_)";
   #endif
	l0 = l0_;
	L = L_;
	N = N_;
	type_id = 1;
}

axis_FoEqSt::axis_FoEqSt(const axis_FoEqSt &that)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoEqSt");
	log << "axis_FoEqSt(const axis_FoEqSt &that)";
   #endif
	l0 = that.l0;
	L = that.L;
	N = that.N;
	type_id = that.type_id;
}

axis_FoEqSt::~axis_FoEqSt()
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoEqSt");
	log << "~axis_FoEqSt()";
   #endif
}

axis * axis_FoEqSt::clone () const
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoEqSt");
	log << "clone() const";
   #endif
	return new axis_FoEqSt(*this);
}

axis * axis_FoEqSt::create_reciprocal() const
{
  #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoEqSt");
	log << "create_reciprocal()";
   #endif

	return new axis_CoEqSt(l0,L,N);
}

double axis_FoEqSt::val_at(const int &index) const
{
	double pi = acos(-1.);
	return ( (index<=N/2) ? (2.*pi*index/L) : (2.*pi*(index-N)/L) );
}

int axis_FoEqSt::index_at(const double &val) const
{
	double pi = acos(-1.);
	return ( (val>=0.)) ? (L*val/(2.*pi)) : (N + L*val/(2.*pi));
}



double axis_FoEqSt::operator() (const int &i) const
{
	double pi = acos(-1.);
	return ( (i<=N/2) ? (2.*pi*i/L) : (2.*pi*(i-N)/L) );
}


int axis_FoEqSt::operator() (const double &x) const
{
	double pi = acos(-1.);
	return ( (x>=0.)) ? (L*x/(2.*pi)) : (N + L*x/(2.*pi));
}


int axis_FoEqSt::operator() (const double &x, double &lambda) const
{
	// ToDo : implementiere axis_FoEqSt::operator() (const double &x, double &lambda) const
	throw("axis_FoEqSt::operator() (const double &x, double &lambda) const - not implemented yet");
}






