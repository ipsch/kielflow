#ifndef INTERFACE_3D_FKT_HPP_
#define INTERFACE_3D_FKT_HPP_

class interface_3d_fkt
{
public :
	virtual ~interface_3d_fkt() { };
	virtual double operator() (const double &x, const double &y, const double &z) const = 0;
};

#endif /* END OF INTERFACE_1D_FKT_HPP_ */
