#ifndef INIT_COND_HPP_
#define INIT_COND_HPP_

#include <cmath>
#include "interface_1d_fkt.hpp"
#include "interface_3d_fkt.hpp"
#include <iostream>
#include <fstream>

void main_loop(void);




class fkt3d_staticPlsm
{
public :
	fkt3d_staticPlsm(const int &N, double * X, double * Y, double * f);
	~fkt3d_staticPlsm();
	void solve(void);

	void get_results(double * Y) const;
	void set_iterations(const int &max_iter){iterations = max_iter;}
	void set_theta(const double &theta_) {theta = theta_;}
private :
	fkt3d_staticPlsm();


	double h;
	double theta;
	int iterations;

	int N_;
	double * X_;
	double * Y_;
	double * f_;
	double * Tt_;


	void lin_jacobi(void);
	double NL_(const int& i); // nonlinearity

	void nlin_jacobi(void);
	double newton(const int &i);
	double f_df(const double &x, const double &r,
			const double &A, const double &B, const double &nd);


};


#endif
