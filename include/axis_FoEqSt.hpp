#ifndef AXIS_FOEQST_HPP_
#define AXIS_FOEQST_HPP_

#include "axis_Fo.hpp"
#include "axis_CoEqSt.hpp"



class axis_CoEqSt;
// ####################################
class axis_FoEqSt : public axis_Fo
{
public :
	axis_FoEqSt(const double &l0_, const double &L_, const int &N_);
	axis_FoEqSt(const axis_FoEqSt &that);
	~axis_FoEqSt();
	double val_at(const int &index) const;
	int index_at(const double &val) const;

	axis * clone () const;
	axis * create_reciprocal() const;


	virtual double operator() (const int &i) const;
	virtual int operator() (const double &x) const;
	virtual int operator() (const double &x, double &lambda) const;
};

#endif
