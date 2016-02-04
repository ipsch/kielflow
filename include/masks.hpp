#ifndef MASKS_HPP_
#define MASKS_HPP_


#include <cmath>
#include "interface_1d_fkt.hpp"
#include "interface_3d_fkt.hpp"
#include "o_math.hpp"

#ifdef _MY_VERBOSE
#include "logger.hpp"
#endif



// 1d functions
double naca_profile(const double &x, const double &d);
double Yukawa_1d(const double &r, const double &q = 1., const double &mu = 1.);

class const_fkt : public interface_1d_fkt
{
public :
	const_fkt(const double &const_val) : C(const_val) { };
	~const_fkt() { };
	double operator() (const double &r) {return C;}
private :
	const double C;
};

class fkt1d_theta : public interface_1d_fkt
{
public :
	fkt1d_theta() : r0(0.), A1(0.), A2(1.) { }
	fkt1d_theta(const double &x0, const double &amplitude1 = 0., const double &amplitude2 = 1.) :
		r0(x0), A1(amplitude1), A2(amplitude2) { };
	~fkt1d_theta() { };
	double operator() (const double &r) const {return (r<r0) ? A1 : A2;}
private :
	const double r0;
	const double A1;
	const double A2;
};

class smooth_rectangle : public interface_1d_fkt
{
public :
	smooth_rectangle(const double &amplitude, const double &steepness, const double &shift) :
		A(amplitude), S(steepness), r0(shift) { };
	~smooth_rectangle() { };
	double operator() (const double &r) const {return A*exp(-S* pow(fabs(r-r0),S));}
private :
	const double r0;
	const double S;
	const double A;
};

class Gauss_1d_fkt : public interface_1d_fkt
{
public :
	Gauss_1d_fkt(const double & amplitude = 1., const double &varianz = 1.) :
		A(amplitude), V(varianz)
		{ }
	~Gauss_1d_fkt() { }
	double operator()(const double &r) const {return A*exp(-V*(r*r));}
private :
	const double A;
	const double V;
};

class fkt1d_Yukawa : public interface_1d_fkt
{
public :
	fkt1d_Yukawa(const double &x_pos = 0., const double &y_pos = 0., const double &z_pos = 0., \
		const double &radius = 1., const double &potential = 10000., const double &screening_length = 1.) :
			pi(acos(-1.)), R(radius), P(potential), mu(screening_length), conducting(false)
				{Q = 4.*pi*P*R*exp(R); }
	void set_by_potential() {Q = 4.*pi*P*R*exp(R*mu); }
	void set_y_charge() { P = .25*Q*exp(-R*mu)/(R*pi); }
	void switch_conductivity(){conducting = !conducting;}
	~fkt1d_Yukawa() { };
	double operator()(const double &r) const
	{
		if(r > R)
			return (.25/pi)*Q*(exp(-r*mu)/r);
		if(conducting)
			return P;
		return 0.5*P*(3. - r*r/(R*R));
	}
private :
	bool conducting; // default = false
	const double pi;
	const double R;

	const double mu;
	double P;
	double Q;
};


class fkt1d_interpolated : public interface_1d_fkt
{
public :
	fkt1d_interpolated(double *X, double *Y, const int &N) : N_(N)
	{
		X_ = new double[N_];
		Y_ = new double[N_];

		for(int i=0; i<N_; ++i)
		{
			X_[i] = X[i];
			Y_[i] = Y[i];
		}
		for(int i=1; i<N_; ++i)
			if(!(X[i-1]<X[i]))
				throw("bad grid");
	}

	double operator()(const double &r) const
	{
		try
		{
			if( (X_[0]>r) || (r>=X_[N_-1]) )
			{
				std::cout << "out of range" << std::endl;
				std::cout << r << " !in[" << X_[0] << ":" << X_[N_-1] << "]" << std::endl;
				throw("x out of range");
			}
		}
		catch(char *e)
		{
			std::cout << e << std::endl;
		}


		int i = 0;
		while(X_[i]<r)
			++i;

		if(r==X_[i-1]) return Y_[i-1];

		double gamma = (r-X_[i-1]) / (X_[i]-X_[i-1]);

		return (1.-gamma)*Y_[i-1] + (gamma)*Y_[i];
	}


	~fkt1d_interpolated()
	{
		delete[] X_;
		delete[] Y_;
	}

private :
	fkt1d_interpolated() : N_(0), X_(0L), Y_(0L) { throw("initialization not allowed"); }
	double * X_;
	double * Y_;
	int N_;


};


class fkt3d_from_fkt1d : public interface_3d_fkt
{
public :
	fkt3d_from_fkt1d(interface_1d_fkt &profile,
			const double &x0 = 0., const double &y0 = 0., const double &z0 = 0.) :
				P(profile), my_x(x0), my_y(y0), my_z(z0) {}
	~fkt3d_from_fkt1d() { };
	double operator() (const double &x, const double &y, const double &z) const
	{
		double R = sqrt((x-my_x)*(x-my_x) + (y-my_y)*(y-my_y) + (z-my_z)*(z-my_z));
		return P(R);
	}
private :
	interface_1d_fkt &P;
	const double my_x;
	const double my_y;
	const double my_z;
};


class fkt3d_shift : public interface_3d_fkt
{
public :
	fkt3d_shift(interface_3d_fkt &fkt, const double &x, const double &y, const double &z) :
		my_fkt(fkt), x0(x), y0(y), z0(z)
{ }
	double operator() (const double &x, const double &y, const double &z) const
	{
		return my_fkt(x-x0,y-y0,z-z0);
	}
private :
	interface_3d_fkt &my_fkt;
	const double x0;
	const double y0;
	const double z0;
};


class mask_sphere : public interface_3d_fkt
{
public :
	mask_sphere(const double &x0, const double &y0, const double &z0, const double &R);
	~mask_sphere() { };
	double operator()(const double &x, const double &y, const double &z) const;
private :
	const double my_R;
};

class mask_box : public interface_3d_fkt
{
public :
	mask_box(const double &x0, const double &y0, const double &z0, const double &Lx, const double &Ly, const double &Lz);
	~mask_box() { };
	double operator() (const double &x, const double &y, const double &z) const;
private :
	const double Lx_;
	const double Ly_;
	const double Lz_;
	const double my_x;
	const double my_y;
	const double my_z;
};

class mask_LetterA : public interface_3d_fkt
{
public :
	mask_LetterA(const double &x0, const double &y0, const double &z0, const double &H = 1., const double &W = 1., const double &L = 1.);
	~mask_LetterA() { };
	double operator()(const double &x, const double &y, const double &z) const;
private :
	const double my_height;
	const double my_width;
	const double my_length;
	const double my_x;
	const double my_y;
	const double my_z;
};


class mask_LetterC : public interface_3d_fkt
{
public :
	mask_LetterC(const double &x0, const double &y0, const double &z0, const double &H = 1., const double &W = 1., const double &L = 1.);
	~mask_LetterC() { };
	double operator()(const double &x, const double &y, const double &z) const;
private :
	const double my_height;
	const double my_width;
	const double my_length;
	const double my_x;
	const double my_y;
	const double my_z;
};


class mask_LetterU : public interface_3d_fkt
{
public :
	mask_LetterU(const double &x0, const double &y0, const double &z0, const double &H = 1., const double &W = 1., const double &L = 1.);
	~mask_LetterU() { };
	double operator()(const double &x, const double &y, const double &z) const;
private :
	const double my_height;
	const double my_width;
	const double my_length;
	const double my_x;
	const double my_y;
	const double my_z;
};


class mask_CAU : public interface_3d_fkt
{
public :
	mask_CAU(const double &x0, const double &y0, const double &z0);
	~mask_CAU() { };
	double operator()(const double &x, const double &y, const double &z) const;
private :
	mask_LetterC C;
	mask_LetterA A;
	mask_LetterU U;
	const double my_x;
	const double my_y;
	const double my_z;
};



class fkt3d_const : public interface_3d_fkt
{
public :
	fkt3d_const(const double &val = 1.) : my_val(val) {	}
	~fkt3d_const() { }
	double operator()(const double &x, const double &y, const double &z) const {return my_val;}
private :
	const double my_val;
};

class fkt3d_const_of_t : public interface_3d_fkt
{
public :
	fkt3d_const_of_t(const double &val = 1.) : my_val(val) {	}
	~fkt3d_const_of_t() { }
	double operator()(const double &x, const double &y, const double &z) const {return my_val;}
	void set_val(const double val) {my_val = val;}
private :
	double my_val;
};

class fkt3d_sin : public interface_3d_fkt
{
public :
	fkt3d_sin(const double &amplitude = 1., const double &kx=1., const double &ky=0., const double &kz=0.) :
		A(amplitude), kx(kx), ky(ky), kz(kz)
		{ }
	~fkt3d_sin() { }
	double operator()(const double &x, const double &y, const double &z) const
	{return A*sin(kx*x + ky*y + kz*z);}
private :
	const double A;
	const double kx;
	const double ky;
	const double kz;
};

class fkt3d_Gauss : public interface_3d_fkt
{
public :
	fkt3d_Gauss(const double & amplitude = 1., const double &varianz_x = 1., const double &varianz_y = 1., const double &varianz_z = 1.) :
		A(amplitude), Vx(varianz_x), Vy(varianz_y), Vz(varianz_z)
		{ }
	~fkt3d_Gauss() { }
	double operator()(const double &x, const double &y, const double &z) const
	{return A*exp(-pow(x/Vx,2.) - pow(y/Vy,2.) - pow(z/Vz,2.));}
private :
	const double A;
	const double Vx;
	const double Vy;
	const double Vz;
};

class fkt3d_Poisson_periodic : public interface_3d_fkt
{
public :
	fkt3d_Poisson_periodic(const double &x_pos = 0., const double &y_pos = 0., const double &z_pos = 0., \
		const double &radius = 1., const double &potential = 10000., \
		const double &Lx = 10., const double &Ly = 10., const double &Lz = 10.);
	~fkt3d_Poisson_periodic() { };
	double operator()(const double &x, const double &y, const double &z) const;
private :
	const double R;
	const double P;
	const double Lx;
	const double Ly;
	const double Lz;
	double Q;
	const double my_x;
	const double my_y;
	const double my_z;
};




// 3d functions
double set_box(const double &x, const double &y, const double &z);
double set_const(const double &x, const double &y, const double &z);
double set_zero(const double &x, const double &y, const double &z);
double set_sin(const double &x, const double &y, const double &z);
double set_sinx(const double &x, const double &y, const double &z);
double set_siny(const double &x, const double &y, const double &z);
double gaussXY(const double &x, const double &y, const double &z);
double gaussXYZ(const double &x, const double &y, const double &z);
double define_H_naca(const double &x, const double &y, const double &z);
double one_over_r(const double &x, const double &y, const double &z);














#endif
