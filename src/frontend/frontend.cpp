#define MOVABLE

// set (git)-repository-Version to unknown
// (if not defined in makefile during compile process)
#ifndef VERSION_STRING
#define VERSION_STRING "unknown"
#endif


// gcc standard C++ libraries (common stuff)
#include <iostream>  // standard I/O-Operations (I/O to console)
#include <fstream>   // same as iostream but for files
#include <cmath>     // basic math stuff (sqrt(), cabs(), pow(a,n) etc.)
#include <string>    // class : string (here: used for filenames only; Maybe I should avoid this library)
#include <sstream>   // class : string-stream (MORE OVERHEAD) I'am using this to interpret data read from a file


// gcc standard C libraries (common stuff)
//#include <complex.h> // simple arithmetics with complex numbers (defines symbol "I" for imaginary unit)
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <unistd.h>
#include <limits.h>
#include <omp.h>

// uncommon libraries
#include "fftw3.h"   // google it (something about discret fourier Transformations)

// own includes
#include "operations.hpp"
#include "parameters.hpp"
#include "grid.hpp"
//#include "particle.hpp"
#include "field.hpp"
#include "masks.hpp"
#include "init_cond.hpp"
#include "subdim.hpp"
#include "solver_poisson_jacobi_lin.hpp"
#include "solver_poisson_jacobi_nlin.hpp"
#include "solver_poisson_multigrid.hpp"

#include "IO.hpp"


//double omp_get_wtime(void);

#if defined (_MY_VERBOSE) || defined (_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif

// number of grid-points
// in different space directions

int Nx = 256;
int Ny = 128;
int Nz = 128;

double Lx = 20.;
double Ly = 10.;
double Lz = 10.;

double M = 0.3;
double tau = 0.1;
double theta = 30.;
double mu = 0.;
double beta = 0.0;






// main #######################################################################
int main(void)
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("frontend");
	my_log << "running version";
	my_log << VERSION_STRING;
   #endif

	double pi = acos(-1.);
	axis_CoSiSt x_axis(-0.5*Lx, Lx, Nx, Lx*(1.-0.25)/(2*pi) );
	axis_CoSiSt y_axis(-0.5*Ly, Ly, Ny, Ly*(1.-0.25)/(2*pi) );
	axis_CoSiSt z_axis(-0.5*Lz, Lz, Nz, Lz*(1.-0.25)/(2*pi) );

	/*
	axis_CoSiSt x_axis(-5.,10.,Nx,1.24);
	axis_CoSiSt y_axis(-5.,10.,Ny,1.24);
	axis_CoSiSt z_axis(-5.,10.,Nz,1.24);

	axis_CoHySt x_axis(-5.,10.,Nx,3);
	axis_CoHySt y_axis(-5.,10.,Ny,3);
	axis_CoHySt z_axis(-5.,10.,Nz,3);

	axis_CoEqSt x_axis(-5.,10.,128);
	axis_CoEqSt y_axis(-5.,10.,128);
	axis_CoEqSt z_axis(-5.,10.,128);
	*/

	grid_Co Omega(x_axis,y_axis,z_axis);
	parameters Params(M,tau,theta,mu,beta);


	subdim my_dim;
	my_dim.xpos = Nx/2.;
	my_dim.ypos = Ny/2.;
	my_dim.zpos = Nz/2.;
	my_dim.direction = 0;
	my_dim.plane = 2;



   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "create parameters ";
   #endif
	double v0 = M/tau;

   #ifdef _MY_VERBOSE
	my_log << "init fields";
   #endif
	field_real Ux(Omega), Uy(Omega), Uz(Omega);
	field_real ni(Omega), nd(Omega);
	field_real Ph(Omega);

	fkt3d_const H_zero(0.);
	Ux.fill(&set_zero);
	Uy.fill(&set_zero);
	Uz.fill(&set_zero);


	//theta_fkt dust_1d_fkt(R, Q, 0.);
	//Gauss_1d_fkt dust_1d_fkt(Q, 10000);
	//theta_fkt dust_1d_fkt(R, Q);
	//smooth_rectangle dust_1d_fkt(nd0, 56, -.8);
	//fkt3d_from_fkt1d dust_3d_fkt(dust_1d_fkt);
	// ##### DUST #####
	double Rd = 0.1183;
	double nd0 = -11.940;
	double shift = 0;
	fkt3d_Gauss dust_3d_fkt(-14.,0.15,0.15,0.15);
	//fkt3d_shift H_3d_shifted_fkt(dust_3d_fkt, shift, 0., 0.);
	nd.fill(dust_3d_fkt);
	//save_2d(nd, my_dim, "./data/nd2d.dat");
	//save_1d(nd, my_dim, "./data/nd1d.dat");


	// ##### RELAXATIONS-SOLVER #####
	fkt3d_const boundary_shape(0.); // there are no additional boundarys (except Domain borders)
	fkt3d_const boundary_value(0.); // the value of the boundarys is always zero
	double SOR = .9; // Overrelaxation parameter
	solver_poisson_jacobi_nlin NLJ(boundary_shape, boundary_value, SOR);
	NLJ.set_max_iterations(80);
	NLJ.limit_max = 1.e-5;
	NLJ.limit_sum = 1.e-7;

	// ##### MULTIGRID #####
	solver_poisson_multigrid MG(NLJ);
	{ // setup MG-cycle
		const int MG_N = 6;
		int * MG_steps_sizes = new int[MG_N]
                  {1,1,1,-1,-1,-1};
	    MG_lvl_control * cycle_shape = new MG_lvl_control[MG_N]
			      {lvl_keep, lvl_down, lvl_keep, lvl_down, lvl_up, lvl_up};
	    MG.set_level_control(MG_N, MG_steps_sizes, cycle_shape);
	    // Sketch of cycle:
	    // _
	    //  \_  /
	    //    \/
	} // configure done

	MG.solve(Ph,nd);


	//Ph.fill(Potential_shifted);
	for(int i=0; i<ni.N; ++i)
		ni.val[i] = exp(-theta*Ph.val[i]);


   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	std::cout << "save all" << std::endl;
	my_log << "save all";
   #endif
	file_create("./data/fields.h5");
	save_grid(Omega, "./data/fields.h5");
	save_time(0.);
	save_parameters(Params, "./data/fields.h5");
	save_field_real("Ux", Ux, "./data/fields.h5");
	save_field_real("Uy", Uy, "./data/fields.h5");
	save_field_real("Uz", Uz, "./data/fields.h5");
	save_field_real("ni", ni, "./data/fields.h5");
	save_field_real("Ph", Ph, "./data/fields.h5");

   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
	std::cout  << "finished - frontend terminated";

	return 0;
}







/*
{ // #### DEALAISING-FILTER TO FILE ####
	field_imag dummy(Omega);
	for(int i=0; i<dummy.N; ++i)
	{
		dummy.val[i][0] = 1.;
		dummy.val[i][1] = 1.;
	}
handle_dealiasing(dummy);
save_major_wavevectors("./diagnostics/DealaisingFilter.dat", dummy);
} // #### DEND DEALAISING-FILTER TO FILE ####
*/

/*
// #### PRESOLVER ####
int N = 40000;
double L = 30;
double * X = new double[N];
double * Y = new double[N];
double * f = new double[N];
// load previous data
std::ifstream input("./data/Phi_1d.dat");
for(int i=0; i<N; ++i)
{
	double tmp_data1;
	double tmp_data2;
	double tmp_data3;
	input >> tmp_data1;
	input >> tmp_data2;
	input >> tmp_data3;

	X[i] = (L/N)*i;
	Y[i] = 0.;
	Y[i] = tmp_data2;
	f[i] = dust_1d_fkt(X[i]);
}
input.close();

// presolver
fkt3d_staticPlsm Phi_init(N,X,Y,f);
Phi_init.set_iterations(1000);
Phi_init.set_theta(10.);
Phi_init.solve();
Phi_init.get_results(Y);

// save result from presolver for next time
std::ofstream output("./data/Phi_1d.dat", std::ios::trunc);
output << std::scientific << std::setprecision(8);
for(int i=0; i<N; ++i)
{
	output << X[i] << "\t" << Y[i] << "\t" << f[i] << std::endl;
}
output.close();
std::cout << "pre-solver done" << std::endl;

fkt1d_interpolated Potential_profile(X, Y, N);
fkt3d_from_fkt1d Potential_fkt3d(Potential_profile);
fkt3d_shift Potential_shifted(Potential_fkt3d, shift,0.,0.);

delete[] X;
delete[] Y;
delete[] f;

*/

