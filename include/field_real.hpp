#ifndef FIELD_REAL_HPP_
#define FIELD_REAL_HPP_

#include "field_base.hpp"
#include "interface_3d_fkt.hpp"

class field_real : public field
{
public :
	field_real(const grid &Omega);
	field_real(const field_real &that);
	~field_real();

	void fill(double (*fill_fkt)(const double &,const double &, const double &)) const;
	void fill(interface_3d_fkt const &rhs);
	bool IsNan();

	template<typename Func>
	void fill2(Func fkt)
	{
   	   #ifdef _MY_VERBOSE_MORE
		logger log("field_real");
		log << "fill(double (*fill_fkt)(double,double,double))";
   	   #endif
		for( int i=0; i < Nx; ++i)
		{
			double x = my_grid.x_axis->val_at(i);
			for( int j=0; j < Ny; ++j)
			{
				double y = my_grid.y_axis->val_at(j);
				for( int k=0; k < Nz; ++k)
				{
					double z = my_grid.z_axis->val_at(k);
					double fkt_val =  fkt(x,y,z);
					val[index(i,j,k)] = fkt_val;
				}
			}
		}
		return;
	}


	void add(double (*fill_fkt)(const double &,const double &, const double &)) const;
	void mulitply(const double &lambda);

	double val_at(int ix, int iy, int iz) const;

	double& operator() (const int &ix, const int &iy, const int &iz); // Subscript operators often come in pairs
	double  operator() (const int &ix, const int &iy, const int &iz) const; // Subscript operators often come in pairs
	double  operator() (const double &x, const double &y, const double &z) const;
	double& operator() (int N); // Subscript operators often come in pairs
	//double  operator() (int N) const; // Subscript operators often come in pairs

	double * val;

	//void resize(const int &NX, const int &NY, const int &NZ);
	void resize(int NX, int NY, int NZ);
	// ToDO : functoren implementieren

	void swap(field_real &that);
	field_real& operator= (const field_real &rhs);
	field_real& operator+= (const field_real &rhs);

	//field_real& operator/= (double const& lambda);
};





#endif
