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

