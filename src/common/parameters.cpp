#include "parameters.hpp"




parameters::parameters(double M_, double tau_, double theta_, double mu_, double beta_) :
		pi(acos(-1)), M(M_), tau(tau_), theta(theta_), mu(mu_), beta(beta_)
{	}

parameters::parameters(const parameters &that) :
		pi(acos(-1)), M(that.M), tau(that.tau), theta(that.theta), mu(that.mu), beta(that.beta)
{	}


parameters& parameters::operator= (parameters const& rhs)
{
	M     = rhs.M;
	tau   = rhs.tau;
	theta = rhs.theta;
	mu    = rhs.mu;
	beta  = rhs.beta;

	return *this;
}
