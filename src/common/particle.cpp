#include "particle.hpp"


particle::particle(const double &x, const double &y, const double &z,
		           const double &q, const double &r) :
pos_x(x), pos_y(y), pos_z(z), charge(q), radius(r)
{

}


double particle::Yukawa(const double &x, const double &y, const double &z)
{
	double pi = acos(-1.);
	double r_sub_r0 = sqrt( (x-pos_x)*(x-pos_x) + (y-pos_y)*(y-pos_y) + (z-pos_z)*(z-pos_z) );




	if(r_sub_r0 > radius)
		return -(1./(4.*pi))*charge*(exp(-r_sub_r0)/r_sub_r0);
	return -(1./(4.*pi))*charge*(exp(-radius)/radius);
}



double particle::Yukawa_periodic(const double &x, const double &y, const double &z,
		const double &Lx, const double &Ly, const double &Lz)
{
	double sum = Yukawa(x,y,z);

	// periodische Randbedingungen 1.Ordnung (6)
	sum +=Yukawa(x+Lx, y   , z   );
	sum +=Yukawa(x-Lx, y   , z   );
	sum +=Yukawa(x   , y+Ly, z   );
	sum +=Yukawa(x   , y-Ly, z   );
	sum +=Yukawa(x   , y   , z+Lz);
	sum +=Yukawa(x   , y   , z-Lz);

	// periodische Randbedingungen 2.Ordnung (12)
	sum +=Yukawa(x+Lx, y+Ly, z   );
	sum +=Yukawa(x+Lx, y-Ly, z   );
	sum +=Yukawa(x+Lx,   Ly, z+Lz);
	sum +=Yukawa(x+Lx,   Ly, z-Lz);
	sum +=Yukawa(x-Lx, y+Ly, z   );
	sum +=Yukawa(x-Lx, y-Ly, z   );
	sum +=Yukawa(x-Lx,   Ly, z+Lz);
	sum +=Yukawa(x-Lx,   Ly, z-Lz);
	sum +=Yukawa(x   , y+Ly, z+Lz);
	sum +=Yukawa(x   , y+Ly, z-Lz);
	sum +=Yukawa(x   , y-Ly, z+Lz);
	sum +=Yukawa(x   , y-Ly, z-Lz);

	// periodische Randbedingungen 3.Ordnung (8)
	sum +=Yukawa(x+Lx, y+Ly, z+Lz);
	sum +=Yukawa(x+Lx, y+Ly, z-Lz);
	sum +=Yukawa(x+Lx, y-Ly, z+Lz);
	sum +=Yukawa(x+Lx, y-Ly, z-Lz);
	sum +=Yukawa(x-Lx, y+Ly, z+Lz);
	sum +=Yukawa(x-Lx, y+Ly, z-Lz);
	sum +=Yukawa(x-Lx, y-Ly, z+Lz);
	sum +=Yukawa(x-Lx, y-Ly, z-Lz);

	return sum;

}
