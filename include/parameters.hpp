#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include<cmath> // für die Berechnung von PI


class parameters
{
public :
	parameters(double M_, double tau_, double theta_, double mu_, double beta_);
	parameters(const parameters &that);
	double pi;  // Kreiszahl [1]
	double M;   // Mach-Number
	double tau;
	double theta;
	double mu;
	double beta;

	parameters& operator= (parameters const& rhs);
};


#endif  // Ende der Parameters Klasse
