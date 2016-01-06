#ifndef INTERFACE_1D_FKT_HPP_
#define INTERFACE_1D_FKT_HPP_

class interface_1d_fkt
{
public :
	virtual ~interface_1d_fkt() { };
	virtual double operator() (const double &r) const = 0;
};

#endif /* END OF INTERFACE_1D_FKT_HPP_ */
