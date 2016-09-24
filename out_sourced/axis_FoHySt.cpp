#include "axis_FoHySt.hpp"








// ############## Hyperbolic axis system ####################################


axis_FoHySt::axis_FoHySt(const double &l0_, const double &L_, const double &N_, const int &m_)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoHySt");
	log << "axis_FoHySt(const double &l0_, const double &L_, const int &N_, const int &m_)";
   #endif
	l0 = l0_;
	L = L_;
	N = N_;
	m = m_;
	type_id = 2;
}


axis_FoHySt::axis_FoHySt(const axis_FoHySt &that)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoHySt");
	log << "axis_FoHySt(const axis_FoHySt &that)";
   #endif
	l0 = that.l0;
	L = that.L;
	N = that.N;
	m = that.m;
	type_id = that.type_id;
}

axis_FoHySt::~axis_FoHySt()
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoHySt");
	log << "~axis_FoHySt()";
   #endif
}



axis * axis_FoHySt::clone () const
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoHySt");
	log << "clone() const";
   #endif
	return new axis_FoHySt(*this);
}

axis * axis_FoHySt::create_reciprocal() const
{
  #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoHySt");
	log << "create_reciprocal()";
   #endif

	return new axis_CoHySt(l0,L,N,m);
}

double axis_FoHySt::val_at(const int &index) const
{
	double pi = acos(-1.);
	return ( (index<=N/2) ? (2.*pi*index/L) : (2.*pi*(index-N)/L) );
}


int axis_FoHySt::index_at(const double &val) const
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoHySt");
	log << "index_at(const double &val) const";
	log << "ERROR";
   #endif
	std::cout << "Nicht implementiert: ";
	std::cout << "axis_FoHySt::index_at(const double &val) const \n";
	return 0;
}

double axis_FoHySt::k_val_at(const int &index) const
{
	double pi = acos(-1.);
	return ( (index<=N/2) ? (2.*pi*index/L) : (2.*pi*(index-N)/L) );
}


int axis_FoHySt::k_index_at(const double &val) const
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoHySt");
	log << "index_at(const double &val) const";
	log << "ERROR";
   #endif
	std::cout << "Nicht implementiert: ";
	std::cout << "axis_FoHySt::index_at(const double &val) const \n";
	return 0;
}

double axis_FoHySt::dS(const int &index) const
{
	double x_val = l0 + L*double(index)/double(N);
	double y_val;

	if     (x_val>= L/2.) y_val = (m/sinh(m))*cosh(2.*(x_val-L)*m/L);
	else if(x_val< -L/2.) y_val = (m/sinh(m))*cosh(2.*(x_val+L)*m/L);
	else                  y_val = (m/sinh(m))*cosh(2.*(x_val  )*m/L);

	return y_val;
}



double axis_FoHySt::operator() (const int &i) const
{
	double pi = acos(-1.);
	return ( (i<=N/2) ? (2.*pi*i/L) : (2.*pi*(i-N)/L) );
}


int axis_FoHySt::operator() (const double &x) const
{
	double pi = acos(-1.);
	return ( (x>=0.)) ? (L*x/(2.*pi)) : (N + L*x/(2.*pi));
}


int axis_FoHySt::operator() (const double &x, double &lambda) const
{
	// ToDo : implementiere axis_FoEqSt::operator() (const double &x, double &lambda) const
	throw("axis_FoEqSt::operator() (const double &x, double &lambda) const - not implemented yet");
}


