#include "init_cond.hpp"
#include "o_math.hpp"





fkt3d_staticPlsm::fkt3d_staticPlsm(const int &N, double * X, double * Y, double * f) :
	N_(N)
{
	//double L = X[N_-1] - X[N_];
	h = X[1] - X[0];

	iterations = 100000;

	theta = 10.;

	X_ = new double[N_];
	Y_ = new double[N_];
	f_ = new double[N_];

	Tt_ = 0L;

	// copy grid
	for(int i=0; i<N_; ++i)
		X_[i] = X[i];

	// copy data
	for(int i=0; i<N_; ++i)
		Y_[i] = Y[i];

	// copy density
	for(int i=0; i<N_; ++i)
		f_[i] = f[i];

}



fkt3d_staticPlsm::~fkt3d_staticPlsm()
{
	delete[] X_;
	delete[] Y_;
	delete[] f_;
}


void fkt3d_staticPlsm::get_results(double * Y) const
{
	for(int i=0; i<N_; ++i)
	{
		Y[i] = Y_[i];
	}
	return;
}


void fkt3d_staticPlsm::solve(void)
{
	Tt_ = new double[N_];


	for(int iter=0; iter<iterations; ++iter)
	{
		double sup = 0.;
		//std::cout << "pre-solver iteration-loop " << iter << std::endl;
		lin_jacobi();
		for(int n=0; n<N_; ++n)
		{
			sup = max<double>(sup, fabs(Y_[n] - Tt_[n]));
			Y_[n] = Tt_[n];
		}
		//std::cout << "sup= " << sup << std::endl;
	}



	delete[] Tt_;
}






void fkt3d_staticPlsm::lin_jacobi(void)
{

	for(int i=1; i<N_-1; ++i)
	{
		double r = X_[i];
		double A = 1. + h/r;
		double B = 1. - h/r;

		Tt_[i] = 0.5*( A*Y_[i+1] + B*Y_[i-1] + h*h*NL_(i) );
	}


	// Border right side (Y_[N_] = 0)
	double r = X_[N_-1];
	double A = 1. + h/r;
	Tt_[N_-1] = 0.5*(  A*Y_[N_-2] + h*h*NL_(N_-1) );

	// Border left side (Y_[-1] = Y_[1])
	Tt_[0]    = 0.5*(  2*Y_[1]  + h*h*NL_(0) );

	return;
}

double fkt3d_staticPlsm::NL_(const int& i)
{
	return exp(-theta*Y_[i]) - exp(Y_[i]) + f_[i];
}




void fkt3d_staticPlsm::nlin_jacobi(void)
{

	for(int i=0; i<N_; ++i)
	{
		Tt_[i] = newton(i);
	}

	return;
}

double fkt3d_staticPlsm::newton(const int &i)
{

	double A, B;

	double r = X_[i];
	double x_i = Y_[i];
	double x_ip1; //
	double rho = f_[i];

	double omega_SOR = 1.;
	double eps = 1.e-8;


	if( (i!=0) && (i!=N_-1))
	{
		A = Y_[i-1];
		B = Y_[i+1];
	}
	else
	{
		if(i==0)
			A = B = Y_[1];
		if(i==N_-1)
			A = Y_[N_-2];
			B = 0.;
	}


	// N_ewton Method
	for(int it=0; it<20; ++it)
	{
		double dx;
		dx = this->f_df(x_ip1,r,A,B,rho);
		x_ip1 -= omega_SOR*dx;
		if( fabs(dx) < eps  ) return x_ip1;
	}

	throw("newton_method: maximum number of iterations exceeded");
	return x_i;
}

double fkt3d_staticPlsm::f_df(const double &x, const double &r, const double &A, const double &B, const double &nd)
{

	// Zähler
	double f = (A -2*x + B) / (h*h);
	       //f+= 2.*(B - A)/(2.*h*r);
		   f+= (B - A)/(h*r);
	       f+= exp(-theta*x) - exp(x) + nd;
	       f+=nd;


	// Nenner
	double df = 0.;

	df -= 2./(h*h);
	df +=  -theta*exp(-theta*x) - exp(x);

	return (f/df);

}


