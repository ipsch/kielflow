#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <cmath>

class particle
{
public :
	particle();
	particle(const double &x, const double &y, const double &z,
			 const double &q, const double &r);
	double get_x() const {return pos_x;}
	double get_y() const {return pos_y;}
	double get_z() const {return pos_z;}
	double get_q() const {return charge;}
	double get_r() const {return radius;}
	double Yukawa(const double &x, const double &y, const double &z);
	double Yukawa_periodic(const double &x, const double &y, const double &z,
			const double &Lx, const double &Ly, const double &Lz);

private :
	double pos_x;
	double pos_y;
	double pos_z;
	double charge;
	double radius;
};


#endif
