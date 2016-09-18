#ifndef AXIS_INTERFACE_HPP_
#define AXIS_INTERFACE_HPP_

class axis_interface
{
public :
	virtual ~axis_interface() {	};

	// unsecure access
	virtual double val_at(const int &index) const = 0;
	virtual int index_at(const double &val) const = 0;
	virtual double k_val_at(const int &index) const = 0;
	virtual int k_index_at(const double &val) const = 0;

	// secure access
	//virtual double operator() (const int &i) const = 0;
	//virtual int operator() (const double &x) const = 0;
	//virtual int operator() (const double &x, double &lambda) const = 0;




};



#endif
