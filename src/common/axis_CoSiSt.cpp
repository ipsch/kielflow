#include "axis_CoSiSt.hpp"




double axis_CoSiSt::bisect(const double &a, const double &b, const double &y) const
{
	double x = 0.5*(a+b);
	double eps = 1.e-6;
	if(fabs(y-S(x)) < eps)
		return x;
	return (S(x) < y) ? bisect(x,b,y) : bisect(a,x,y);
}



axis_CoSiSt::axis_CoSiSt(const double &l0_, const double &L_, const int &N_, const double &m_)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoSiSt");
	log << "axis_CoSiSt(const double &l0_, const double &L_, const int &N_, const int &m_)";
   #endif
	l0 = l0_;
	L = L_;
	N = N_;
	m = m_;
	type_id = 3;
}

axis_CoSiSt::axis_CoSiSt(const axis_CoSiSt &that)
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoSiSt");
	log << "axis_CoSiSt(const axis_CoSiSt &that)";
   #endif
	l0 = that.l0;
	L = that.L;
	N = that.N;
	m = that.m;
	type_id = that.type_id;
}

axis_CoSiSt::~axis_CoSiSt()
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoSiSt");
	log << "~axis_CoSiSt()";
   #endif
}

axis * axis_CoSiSt::clone () const
{
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_FoSiSt");
	log << "clone() const";
   #endif
	return new axis_CoSiSt(*this);
}

axis * axis_CoSiSt::create_reciprocal() const
{
  #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoSiSt");
	log << "create_reciprocal()";
   #endif

	return new axis_FoSiSt(l0,L,N,m);
}

inline double axis_CoSiSt::S(const double &x) const
{
	const double pi = acos(-1.);
	return x - m*sin(2.*pi*x/L);
}

double axis_CoSiSt::val_at(const int &index) const
{
	double x;
	x = l0 + L*double(index)/(N);

	int shift = 0;


	if(x< -L/2.)
	{
		do
		{
			x += L;
			shift -= 1;
		} while(x<-L/2.);
	}

	if(x>= L/2.)
	{
		do
		{
			x -= L;
			shift += 1;
		} while(x>=l0+L);
	}

	x =  shift*L + S(x);
	return x;
}



int axis_CoSiSt::index_at(const double &val) const
{
	// ToDo : shift wenn val out of range
	std::cout << "WARNING: not tested yet:" << std::endl;
	std::cout << "int axis_CoSiSt::index_at(const double &val) const" << std::endl;
	double x_eq;
	x_eq = bisect(l0,l0+L,val);
	int r_int = nearbyint(N*(x_eq-l0)/L);
	return r_int;
}



double axis_CoSiSt::operator() (const int &i) const
{
	double x;
	x = l0 + L*double(i)/(N);

	int shift = 0;


	if(x< -L/2.)
	{
		do
		{
			x += L;
			shift -= 1;
		} while(x<-L/2.);
	}

	if(x>= L/2.)
	{
		do
		{
			x -= L;
			shift += 1;
		} while(x>=l0+L);
	}

	x =  shift*L + S(x);
	return x;
}


int axis_CoSiSt::operator() (const double &x) const
{
	// ToDo : shift wenn val out of range
	std::cout << "WARNING: not tested yet:" << std::endl;
	std::cout << "axis_CoSiSt::operator() (const double &x) const" << std::endl;

	double x_eq;
	x_eq = bisect(l0,l0+L,x);
	int r_int = nearbyint(N*(x_eq-l0)/L);
	return r_int;
}


int axis_CoSiSt::operator() (const double &x, double &lambda) const
{
	// ToDo : invertierung dieser Achse implementieren
	std::cout << "ERROR in Invertierung von CoSiSt nicht implementiert" << std::endl;
	throw("ERROR in Invertierung von CoSiSt nicht implementiert");

	 return 0;
}

