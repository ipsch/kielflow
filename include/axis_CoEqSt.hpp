#ifndef AXIS_COEQST_HPP_
#define AXIS_COEQST_HPP_

#include "axis_Co.hpp"
#include "axis_FoEqSt.hpp"



// ####################################
class axis_CoEqSt : public axis_Co
{
public :
	axis_CoEqSt(const double &l0_, const double &L_, const int &N_);
	axis_CoEqSt(const axis_CoEqSt &that);
	~axis_CoEqSt();
	double val_at(const int &index) const;
	int index_at(const double &val) const;

	axis * clone () const;
	axis * create_reciprocal() const;


	virtual double operator() (const int &i) const;
	virtual int operator() (const double &x) const;
	virtual int operator() (const double &x, double &lambda) const;
};


#endif
