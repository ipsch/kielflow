#ifndef AXIS_COHYST_HPP_
#define AXIS_COHYST_HPP_

#if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif

#include "axis.hpp"


// ####################################
class axis_CoHySt : public virtual axis
{
public :
	axis_CoHySt(const double &l0_, const double &L_, const int &N_, const int &m_=1);
	axis_CoHySt(const axis_CoHySt &that);
	~axis_CoHySt();
	axis * clone () const;

	double val_at(const int &index) const;
	int index_at(const double &val) const;
	double k_val_at(const int &index) const;
	int k_index_at(const double &val) const;

	void resize(const int &N_);

	// ToDo : remove
	//virtual double operator() (const int &i) const;
	//virtual int operator() (const double &x) const;
	//virtual int operator() (const double &x, double &lambda) const;
};


#endif
