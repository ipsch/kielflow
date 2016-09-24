#include "axis_SiSt.hpp"




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
	/* m_ = Der Anstieg der  Skalierungsfunktion bei x=0; Werte für m_ ]0,1[*/
   #ifdef _MY_VERBOSE_MORE
	logger log("axis_CoSiSt");
	log << "axis_CoSiSt(const double &l0_, const double &L_, const int &N_, const int &m_)";
   #endif
	double pi = acos(-1.);
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
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("axis_CoSiSt");
	my_log << "clone()";
   #endif
	return new axis_CoSiSt(*this);
}


double axis_CoSiSt::val_at(const int &index) const
{

	int shift = 0;
	if(index<0)
	{
		do
		{
			shift += 1;
		} while(index+shift*N <0);
	}

	if(index>=N)
	{
		do
		{
			shift -= 1;
		} while(index+shift*N>=N-1);
	}


	double x;
	x = l0 + L*double(index+shift*N)/(N);
	x =  S(x)-shift*L;
	return x;
}


double axis_CoSiSt::k_val_at(const int &index) const
{
	double pi = acos(-1.);
	return ( (index<=N/2) ? (2.*pi*index/L) : (2.*pi*(index-N)/L) );
}

int axis_CoSiSt::index_at(const double &val) const
{
	double xtmp = val;

	while(xtmp<l0)
		xtmp += L;
	while(xtmp >= l0+L)
		xtmp -= L;

	double x_eq;
	x_eq = bisect(l0,l0+L,xtmp);
	int r_int = nearbyint(N*(x_eq-l0)/L);
	if(r_int==N)
		r_int = 0;
	return r_int;
}

int axis_CoSiSt::k_index_at(const double &val) const
{
	std::cout << "error: axis_CoSiSt::k_index_at not implemented\n";
	throw("error: axis_CoSiSt::k_index_at not implemented\n");
	return -1;
}


inline double axis_CoSiSt::S(const double &x) const
{
	const double pi = acos(-1.);
	return x - m*sin(2.*pi*x/L);
}

double axis_CoSiSt::dS(const int &index) const
{
	const double pi = acos(-1.);
	double x = l0 + (L/N)*index;
	return 1 - (2.*pi*m/L)*cos(2.*pi*x/L);
}


/*
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
	double xtmp = x;

	while(xtmp<l0)
		xtmp += L;
	while(xtmp >= l0+L)
		xtmp -= L;

	double x_eq;
	x_eq = bisect(l0,l0+L,xtmp);
	int r_int = nearbyint(N*(x_eq-l0)/L);
	if(r_int==N) // ToDo : This step is missing in all other axis types
		r_int = 0;
	return r_int;
}


int axis_CoSiSt::operator() (const double &x, double &lambda) const
{
	// ToDo : invertierung dieser Achse implementieren
	std::cout << "ERROR in Invertierung von CoSiSt nicht implementiert" << std::endl;
	throw("ERROR in Invertierung von CoSiSt nicht implementiert");

	 return 0;
}

*/
