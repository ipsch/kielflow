#ifndef AXIS_COSIST_HPP_
#define AXIS_COSIST_HPP_


#include "axis_Co.hpp"
#include "axis_FoSiSt.hpp"


// ####################################
class axis_CoSiSt : public virtual axis_Co
{
public :
	axis_CoSiSt(const double &l0_, const double &L_, const int &N_, const double &m_=0.);
	axis_CoSiSt(const axis_CoSiSt &that);
	~axis_CoSiSt();
	double val_at(const int &index) const;
	int index_at(const double &val) const;
	axis * clone () const;
	axis * create_reciprocal() const;

	virtual double operator() (const int &i) const;
	virtual int operator() (const double &x) const;
	virtual int operator() (const double &x, double &lambda) const;
private :
	double bisect(const double &a, const double &b, const double &y) const;
	double S(const double &x) const;
};


#endif
