#ifndef AXIS_FOSIST_HPP_
#define AXIS_FOSIST_HPP_


#include "axis_Fo.hpp"
#include "axis_CoSiSt.hpp"


// ############## Hyperbolic axis system ####################################
class axis_FoSiSt : public virtual axis_Fo
{
public :
	axis_FoSiSt(const double &l0_, const double &L_, const double &N_, const double &m_=0.);
	axis_FoSiSt(const axis_FoSiSt &that);
	~axis_FoSiSt();
	double val_at(const int &index) const;
	int index_at(const double &val) const;
	axis * clone () const;
	axis * create_reciprocal() const;

	virtual double operator() (const int &i) const;
	virtual int operator() (const double &x) const;
	virtual int operator() (const double &x, double &lambda) const;

	double dS(const int &index) const;
private :
};


#endif
