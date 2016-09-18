#ifndef AXIS_COSIST_HPP_
#define AXIS_COSIST_HPP_

#if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif

#include "axis.hpp"


// ####################################
class axis_CoSiSt : public virtual axis
{
public :
	axis_CoSiSt(const double &l0_, const double &L_, const int &N_, const double &m_=0.);
	axis_CoSiSt(const axis_CoSiSt &that);
	~axis_CoSiSt();
	axis * clone () const;

	double val_at(const int &index) const;
	int index_at(const double &val) const;
	double k_val_at(const int &index) const;
	int k_index_at(const double &val) const;

	void resize(const int &N_) {N = N_;}

	// ToDo : remove
	//virtual double operator() (const int &i) const;
	//virtual int operator() (const double &x) const;
	//virtual int operator() (const double &x, double &lambda) const;
private :
	double bisect(const double &a, const double &b, const double &y) const;
	double S(const double &x) const;
	double dS(const int &index) const;
};


#endif
