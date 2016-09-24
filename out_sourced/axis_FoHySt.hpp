#ifndef AXIS_FOHYST_HPP_
#define AXIS_FOHYST_HPP_


#include "../include/axis_HySt.hpp"
#include "axis_Fo.hpp"


// ############## Hyperbolic axis system ####################################
class axis_FoHySt : public virtual axis_Fo
{
public :
	axis_FoHySt(const double &l0_, const double &L_, const double &N_, const int &m_=1);
	axis_FoHySt(const axis_FoHySt &that);
	~axis_FoHySt();
	double val_at(const int &index) const;
	int index_at(const double &val) const;
	double k_val_at(const int &index) const;
	int k_index_at(const double &val) const;
	axis * clone () const;
	axis * create_reciprocal() const;

	double dS(const int &index) const;

	virtual double operator() (const int &i) const;
	virtual int operator() (const double &x) const;
	virtual int operator() (const double &x, double &lambda) const;

};




#endif
