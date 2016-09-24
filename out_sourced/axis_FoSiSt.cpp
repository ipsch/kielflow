#include "axis_FoSiSt.hpp"








// ############## Hyperbolic axis system ####################################


axis_FoSiSt::axis_FoSiSt(const double &l0_, const double &L_, const double &N_, const double &m_)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("axis_FoSiSt");
	my_log << "axis_FoSiSt(const double &l0_, const double &L_, const int &N_, const int &m_)";
   #endif
	l0 = l0_;
	L = L_;
	N = N_;
	m = m_;
	type_id = 3;
}


axis_FoSiSt::axis_FoSiSt(const axis_FoSiSt &that)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("axis_FoSiSt");
	my_log << "axis_FoSiSt(const axis_FoSiSt &that)";
   #endif
	l0 = that.l0;
	L = that.L;
	N = that.N;
	m = that.m;
	type_id = that.type_id;
}

axis_FoSiSt::~axis_FoSiSt()
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("axis_FoSiSt");
	my_log << "~axis_FoSiSt()";
   #endif
}



axis * axis_FoSiSt::clone () const
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("axis_FoSiSt");
	my_log << "clone() const";
   #endif
	return new axis_FoSiSt(*this);
}

axis * axis_FoSiSt::create_reciprocal() const
{
  #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("axis_FoSiSt");
	my_log << "create_reciprocal()";
   #endif

	return new axis_CoSiSt(l0,L,N,m);
}

double axis_FoSiSt::val_at(const int &index) const
{
	double pi = acos(-1.);
	return ( (index<=N/2) ? (2.*pi*index/L) : (2.*pi*(index-N)/L) );
}

double axis_FoSiSt::k_val_at(const int &index) const
{
	double pi = acos(-1.);
	return ( (index<=N/2) ? (2.*pi*index/L) : (2.*pi*(index-N)/L) );
}


int axis_FoSiSt::index_at(const double &val) const
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("axis_FoSiSt");
	my_log << "index_at(const double &val) const";
	my_log << "ERROR";
   #endif
	std::cout << "Nicht implementiert: ";
	std::cout << "axis_FoSiSt::index_at(const double &val) const \n";
	return 0;
}

int axis_FoSiSt::k_index_at(const double &val) const
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("axis_FoSiSt");
	my_log << "index_at(const double &val) const";
	my_log << "ERROR";
   #endif
	std::cout << "Nicht implementiert: ";
	std::cout << "axis_FoSiSt::k_index_at(const double &val) const \n";
	return 0;
}



double axis_FoSiSt::operator() (const int &i) const
{
	double pi = acos(-1.);
	return ( (i<=N/2) ? (2.*pi*i/L) : (2.*pi*(i-N)/L) );
}


int axis_FoSiSt::operator() (const double &x) const
{
	double pi = acos(-1.);
	return ( (x>=0.)) ? (L*x/(2.*pi)) : (N + L*x/(2.*pi));
}


int axis_FoSiSt::operator() (const double &x, double &lambda) const
{
	// ToDo : implementiere axis_FoEqSt::operator() (const double &x, double &lambda) const
	throw("axis_FoEqSt::operator() (const double &x, double &lambda) const - not implemented yet");
}


double axis_FoSiSt::dS(const int &index) const
{
	const double pi = acos(-1.);
	double x = l0 + (L/N)*index;
	return 1 - (2.*pi*m/L)*cos(2.*pi*x/L);
}
