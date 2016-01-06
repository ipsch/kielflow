#include "masks.hpp"









mask_box::mask_box(const double &x0, const double &y0, const double &z0, const double &Lx, const double &Ly, const double &Lz ) :
	my_x(x0), my_y(y0), my_z(z0), Lx_(Lx), Ly_(Ly), Lz_(Lz)
{

}

double mask_box::operator()(const double &x, const double &y, const double &z) const
{
	if( (x> my_x - 0.5*Lx_) && (x<my_x + 0.5*Lx_) )
		if( (y> my_y - 0.5*Ly_) && (y<my_y + 0.5*Ly_) )
			if( (z> my_z - 0.5*Lz_) && (z<my_z + 0.5*Lz_) )
				return 1.;
	return 0.;
}



mask_LetterA::mask_LetterA(const double &x0, const double &y0, const double &z0, const double &H, const double &W, const double &L) :
		my_x(x0), my_y(y0), my_z(z0), my_height(H), my_width(W), my_length(L)
{

}

double mask_LetterA::operator()(const double &x, const double &y, const double &z) const
{
	double d = 0.33333*my_width;

	if( (z > my_z + 0.5*my_length) && (z < my_z - 0.5*my_length) )
		return 0.;

	if ( (y < my_y + 0.5*my_height) && (y > my_y - 0.5*my_height) )
		if ( (x > my_x - d +  d*(y-my_y)/my_height ) && (x < my_x  +  d*(y-my_y)/my_height) )
			return 1.;

	if ( (y < my_y ) && (y > my_y - 0.5*my_height) )
		if ( (x > my_x - d*(y-my_y)/my_height ) && (x < my_x +d - d*(y-my_y)/my_height) )
			return 1.;

	return 0.;
}



mask_LetterU::mask_LetterU(const double &x0, const double &y0, const double &z0, const double &H, const double &W, const double &L) :
	my_x(x0), my_y(y0), my_z(z0), my_height(H), my_width(W), my_length(L)
{

}

double mask_LetterU::operator()(const double &x, const double &y, const double &z) const
{
	if( (z > my_z + 0.5*my_length) && (z < my_z - 0.5*my_length) )
		return 0.;

	if( (x>= my_x -0.5*my_width) && (x< my_x - 0.25*my_width) )
		if( (y< my_y + 0.5*my_height) && (y> my_y - 0.5*my_height + 0.5*my_width) )
			return 1.;

	if( (x< my_x +0.5*my_width) && (x> my_x + 0.25*my_width) )
		if( (y< my_y + 0.5*my_height) && (y> my_y - 0.5*my_height + 0.5*my_width) )
			return 1.;

	if( y <= my_y - 0.5*my_height + 0.5*my_width)
		if ( sqrt( pow(x-my_x, 2.) + pow(my_y - 0.5*my_height + 0.5*my_width - y , 2.) ) < 0.5*my_width )
			if ( sqrt( pow(x-my_x, 2.) + pow(my_y - 0.5*my_height + 0.5*my_width - y , 2.) ) > 0.25*my_width )
				return 1.;

	return 0.;
}



mask_LetterC::mask_LetterC(const double &x0, const double &y0, const double &z0, const double &H, const double &W, const double &L) :
	my_x(x0), my_y(y0), my_z(z0), my_height(H), my_width(W), my_length(L)
{

}

double mask_LetterC::operator()(const double &x, const double &y, const double &z) const
{
	double d = min<double>(0.25*my_height, 0.25*my_width);

	if( (z > my_z + 0.5*my_length) && (z < my_z - 0.5*my_length) )
		return 0.;

	if( (x>my_x) && (y< my_y + d*0.95) && (y > my_y - d*0.95) )
		return 0;


	if( sqrt( pow((x-my_x)/my_width ,2.) + pow( (y-my_y)/my_height,2.)  ) <.5)
		if( sqrt( pow((x-my_x)/my_width ,2.) + pow( (y-my_y)/my_height,2.)  ) > 0.5*(1.-d))
			return 1.;

	return 0.;

}



mask_CAU::mask_CAU(const double &x0, const double &y0, const double &z0) :
		my_x(x0), my_y(y0), my_z(z0), C(-1.,0.,0.,1.,1.,1.), A(0.,0.,0.,1.,1.,1.), U(1.,0.,0.,1.,1.,1.)
{

}

double mask_CAU::operator()(const double &x, const double &y, const double &z) const
{


	if( (C(x,y,z)==1.) || (A(x,y,z)==1.) || (U(x,y,z)==1.) )
	{
		return 1.;
	}

	return 0;


}



double set_box(const double &x, const double &y, const double &z)
{
	double my_x = 0.;
	double my_y = 0.;
	double my_z = 0.;

	double Lx_ = 4.;
	double Ly_ = 4.;
	double Lz_ = 4.;


	if( (x> my_x - 0.5*Lx_) && (x<my_x + 0.5*Lx_) )
		if( (y> my_y - 0.5*Ly_) && (y<my_y + 0.5*Ly_) )
			if( (z> my_z - 0.5*Lz_) && (z<my_z + 0.5*Lz_) )
				return 10.;

	return 0.;
}

double set_const(const double &x, const double &y, const double &z)
{
	return 1.;
}


/*
 *
 *
 class fkt3d_Yukawa : public mask
{
public :
	fkt3d_Yukawa(const double &x_pos = 0., const double &y_pos = 0., const double &z_pos = 0., \
		const double &radius = 1., const double &potential = 10000.);
	double operator()(const double &x, const double &y, const double &z) const;
private :
	const double R;
	const double Phi;
};

 *
 */



fkt3d_Poisson_periodic::fkt3d_Poisson_periodic(const double &x0, const double &y0, const double &z0, \
		const double &radius, const double &potential, \
		const double &Lx, const double &Ly, const double &Lz) :
				my_x(x0), my_y(y0), my_z(z0), R(radius), P(potential), Lx(Lx), Ly(Ly), Lz(Lz)
{
	const double pi = acos(-1.);
	Q = -4.*pi*P*R;
}

double fkt3d_Poisson_periodic::operator() (const double &x, const double &y, const double &z) const
{
	const double pi = acos(-1.);
	const double N  = .25*Q/pi;

	double sum = 0;

	double r = sqrt( (x-my_x)*(x-my_x) + (y-my_y)*(y-my_y) + (z-my_z)*(z-my_z));

	if(r < R) return P;

	// Ladung
	sum -= one_over_r( x-my_x, y-my_y, z-my_z );
	// O1 - Spiegelladung
	sum += one_over_r( x-my_x+Lx, y-my_y   , z-my_z    );
	sum += one_over_r( x-my_x-Lx, y-my_y   , z-my_z    );
	sum += one_over_r( x-my_x   , y-my_y+Ly, z-my_z    );
	sum += one_over_r( x-my_x   , y-my_y-Ly, z-my_z    );
	sum += one_over_r( x-my_x   , y-my_y   , z-my_z+Lz );
	sum += one_over_r( x-my_x   , y-my_y   , z-my_z-Lz );
	// O2 - Spiegelladung
	sum -= one_over_r( x-my_x   , y-my_y+Ly, z-my_z+Lz );
	sum -= one_over_r( x-my_x   , y-my_y+Ly, z-my_z-Lz );
	sum -= one_over_r( x-my_x   , y-my_y-Ly, z-my_z+Lz );
	sum -= one_over_r( x-my_x   , y-my_y-Ly, z-my_z-Lz );
	sum -= one_over_r( x-my_x+Lx, y-my_y   , z-my_z+Lz );
	sum -= one_over_r( x-my_x+Lx, y-my_y   , z-my_z-Lz );
	sum -= one_over_r( x-my_x-Lx, y-my_y   , z-my_z+Lz );
	sum -= one_over_r( x-my_x-Lx, y-my_y   , z-my_z-Lz );
	sum -= one_over_r( x-my_x+Lx, y-my_y+Ly, z-my_z    );
	sum -= one_over_r( x-my_x+Lx, y-my_y-Ly, z-my_z    );
	sum -= one_over_r( x-my_x-Lx, y-my_y+Lz, z-my_z    );
	sum -= one_over_r( x-my_x-Lx, y-my_y-Lz, z-my_z    );
	// O3 - Spiegelladung
	sum += one_over_r( x-my_x+Lx, y-my_y+Ly, z-my_z+Lz );
	sum += one_over_r( x-my_x+Lx, y-my_y+Ly, z-my_z-Lz );
	sum += one_over_r( x-my_x+Lx, y-my_y-Ly, z-my_z+Lz );
	sum += one_over_r( x-my_x+Lx, y-my_y-Ly, z-my_z-Lz );
	sum += one_over_r( x-my_x-Lx, y-my_y+Ly, z-my_z+Lz );
	sum += one_over_r( x-my_x-Lx, y-my_y+Ly, z-my_z-Lz );
	sum += one_over_r( x-my_x-Lx, y-my_y-Ly, z-my_z+Lz );
	sum += one_over_r( x-my_x-Lx, y-my_y-Ly, z-my_z-Lz );

	return N*sum;
}




double set_zero(const double &x, const double &y, const double &z)
{
	return 0.;
}

double set_sin(const double &x, const double &y, const double &z)
{
	double pi = acos(-1.);
	double A = -1.;
	double w = 2.*pi/10.;
	return A*sin(w*x);
	//return sin(2.*pi*x/8)*sin(2.*pi*y/10.);
}

double set_sinx(const double &x, const double &y, const double &z)
{
	double pi = acos(-1.);
	double L = 10.;
	double A = .5;
	return A*sin(2.*pi*x/L);
	//return sin(2.*pi*x/8)*sin(2.*pi*y/10.);
}

double set_siny(const double &x, const double &y, const double &z)
{
	double pi = acos(-1.);
	double L = 6;
	return -0.5*sin(2.*pi*y/10.);
	//return sin(2.*pi*x/8)*sin(2.*pi*y/10.);
}

double one_over_r(const double &x, const double &y, const double &z)
{
	return 1./sqrt(x*x + y*y + z*z);
}

double Yukawa_1d(const double &r, const double &q, const double &mu)
{
	const double pi = acos(-1.);
	const double N = .25/pi;
	return -N*q*(exp(-mu*r)/r);
}

double gaussXY(const double &x, const double &y, const double &z)
{
	return 10.*exp(-.7*(x*x+y*y));
	//return sin(2.*pi*x/8)*sin(2.*pi*y/10.);
}



double naca_profile(const double &x, const double &d)
{
	double a0 = 1.4845;
	double a1 = -0.63;
	double a2 = -1.758;
	double a3 = 1.4215;
	double a4 = -0.5075;

	double l=1.;

	return (d/l)*(a0*sqrt(x) + a1*x + a2*x*x + a3*x*x*x + a4*x*x*x*x);
}

double define_H_naca(const double &x, const double &y, const double &z)
{
	if(x<0.) return 0.;
	if(x>1.) return 0.;

	double y_profil_pos =  naca_profile(x,.3);
	double y_profil_neg = -naca_profile(x,.1);

	if( (y>y_profil_pos) || (y<y_profil_neg) )
		return 0;

	return 1;
}
