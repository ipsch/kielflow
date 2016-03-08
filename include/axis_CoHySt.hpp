#ifndef AXIS_COHYST_HPP_
#define AXIS_COHYST_HPP_

#if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif

#include "axis_Co.hpp"
#include "axis_FoHySt.hpp"


// ####################################
class axis_CoHySt : public virtual axis_Co
{
public :
	axis_CoHySt(const double &l0_, const double &L_, const int &N_, const int &m_=1);
	axis_CoHySt(const axis_CoHySt &that);
	~axis_CoHySt();
	double val_at(const int &index) const;
	int index_at(const double &val) const;
	axis * clone () const;
	axis * create_reciprocal() const;

	virtual double operator() (const int &i) const;
	virtual int operator() (const double &x) const;
	virtual int operator() (const double &x, double &lambda) const;
};


#endif
