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
#include "field_integrate.hpp"
#include "field_interpolation.hpp"
#include "masks.hpp"
#include "subdim.hpp"
#include "solver_poisson_jacobi_lin.hpp"
#include "solver_poisson_jacobi_nlin.hpp"
#include "solver_poisson_multigrid.hpp"
#include "OP_FFT.hpp"
#include "OP_iFFT.hpp"
#include "operations.hpp"
#include "IO.hpp"


//double omp_get_wtime(void);

#if defined (_MY_VERBOSE) || defined (_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif


double pi = acos(-1.);
// number of grid-points
// in different space directions

int Nx = 512;
//int Nx = 384;
int Ny = 4;
int Nz = 4;
// Box dimensions
double Lx = 26.;
double Ly = 2.;
double Lz = 2.;

// physical parameters
double M = .0;
double tau = 0.1;
double global_theta = 30.;
double mu = 0.;
double beta = 0.0;
//double global_Q = -11498.5;
double global_Q = 0.;
double radius_a = 0.15;

void create_input_from_old_data(field_real &Ux, field_real &Uy, field_real &Uz, field_real &ni, field_real &Ph)
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("create_input_from_old_data(..)");
	my_log << "start";
   #endif

   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "load old data (in fourier space)";
   #endif
	field_imag FUx = load_field_imag("Ux", "./data/input1.h5");
	field_imag FUy = load_field_imag("Uy", "./data/input1.h5");
	field_imag FUz = load_field_imag("Uz", "./data/input1.h5");
	field_imag Fni = load_field_imag("ni", "./data/input1.h5");
	field_imag FPh = load_field_imag("Ph", "./data/input1.h5");

   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "allocate memory for old data in real space ";
   #endif
	field_real AUx(FUx.my_grid);
	field_real AUy(FUx.my_grid);
	field_real AUz(FUx.my_grid);
	field_real Ani(FUx.my_grid);
	field_real APh(FUx.my_grid);

   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "transform fourier space data to real space";
   #endif
	iFFT(FUx,AUx);
	iFFT(FUy,AUy);
	iFFT(FUz,AUz);
	iFFT(Fni,Ani);
	iFFT(FPh,APh);

   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "interpolate old data onto new grids";
   #endif
	OP_XhtoYh_lvl1(AUx,Ux,2,0.,-3.);
	OP_XhtoYh_lvl1(AUy,Uy,2,0.,-3.);
	OP_XhtoYh_lvl1(AUz,Uz,2,0.,-3.);
	OP_XhtoYh_lvl1(Ani,ni,2,1.,-3.);
	OP_XhtoYh_lvl1(APh,Ph,2,0.,-3.);

	OP_smoothing_lvl2(Ux,Ux);
	OP_smoothing_lvl2(Uy,Uy);
	OP_smoothing_lvl2(Uz,Uz);
	OP_smoothing_lvl2(ni,ni);

   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";// number of grid-points
	// in different space directions
   #endif

	return;
}

void create_input_from_MGsolver(field_real &ni, field_real &Ph)
{
	// ##### DUST #####
	field_real nd(ni.my_grid);
	//double Rd = 0.1183;
	double default_Q = global_Q; //-11498.5;
	double scale_Q = 8.1720e-06;
	// -Q/2299.7
	//fkt3d_Gauss dust_3d_fkt(Q*scale_Q,0.15,0.15,0.15);
	fkt3d_Gauss dust_3d_fkt(default_Q*scale_Q,radius_a,radius_a,radius_a);
	//fkt3d_shift H_3d_shifted_fkt(dust_3d_fkt, shift, 0., 0.);
	nd.fill(dust_3d_fkt);

	// ##### Output Dust #####
	subdim my_dim;
	my_dim.xpos = Nx/2.;
	my_dim.ypos = Ny/2.;
	my_dim.zpos = Nz/2.;
	my_dim.direction = 0;
	my_dim.plane = 2;
	save_2d(nd, my_dim, "./data/nd2d.dat");
	save_1d(nd, my_dim, "./data/nd1d.dat");

	field_integrate IntegratorTrapez(nd.my_grid, 1);
	field_integrate IntegratorSimpson(nd.my_grid, 2);
	double Qges;
	Qges=IntegratorTrapez.execute(nd);
	std::cout << "Qges= " << Qges << " (by trapez rule)\n";
	Qges=IntegratorSimpson.execute(nd);
	std::cout << "Qges= " << Qges << " (by Simpson rule)\n";

	// ##### RELAXATIONS-SOLVER #####
	fkt3d_const boundary_shape(0.); // there are no additional boundarys (except Domain borders)
	fkt3d_const boundary_value(0.); // the value of the boundarys is always zero
	double SOR = .8; // Overrelaxation parameter
	solver_poisson_jacobi_nlin NLJ(boundary_shape, boundary_value, SOR);
	NLJ.set_max_iterations(80);
	NLJ.set_limit_NormMax(1.e-5);
	NLJ.set_limit_NormSum(1.e-7);

	// ##### MULTIGRID #####
	solver_poisson_multigrid MG(NLJ);
	const int MG_N = 4;
	int * MG_steps_sizes = new int[MG_N]
                 {0,-1,-1,-1};
	MG_lvl_control * cycle_shape = new MG_lvl_control[MG_N]
			  { lvl_down3, lvl_up, lvl_up, lvl_up};
	double * MG_error = new double[MG_N]
			  {.01,
			   0.0005,
			   0.0001,
			   0.0001};
	MG.set_level_control(MG_N, MG_steps_sizes, cycle_shape, MG_error);

	delete[] MG_steps_sizes;
	delete[] cycle_shape;


	MG.solve(Ph,nd);

	for(int i=0; i<ni.N; ++i)
	{
		ni.val[i] = exp(-global_theta*Ph.val[i]);
	}




	return;
}


// main #######################################################################
int main(int argc,char **argv)
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("frontend");
	my_log << "running version";
	my_log << VERSION_STRING;
   #endif

	parameters Params(M,tau,global_theta,mu,beta);

	int switch_opt;
	int opt_indent = 0;
    opterr = 0;
	 while ((switch_opt = getopt (argc, argv, "a:hM:Q:T:t:b:")) != -1)
	    {

	    	switch (switch_opt)
	    	{

	        case 'a':
	        	radius_a = std::atof(optarg);
	        	opt_indent+=2;
	        	break;

	        case 'M':
	        	Params.M = std::atof(optarg);
	        	opt_indent+=2;
	        	break;

	        case 'Q':
	        	global_Q = std::atof(optarg);
	        	opt_indent+=2;
	        	break;

	        case 'T':
	        	global_theta = std::atof(optarg);
	        	Params.theta = global_theta;
	        	opt_indent+=2;
	        	break;

	        case 't':
	        	Params.tau = std::atof(optarg);
	        	opt_indent+=2;
	        	break;

	        case 'b':
	        	Params.beta = std::atof(optarg);
	        	opt_indent+=2;
	        	break;

	    	case 'h':
	    		// ToDo : write help fkt
	    		//help();
	    		return 0;
	    		break;

	        case '?':
	        	if (optopt == 'f')
	        	{
	        		std::cout << "Option" << optopt << " requires an argument.\n";
	        		break;
	        	}
	        	else if (isprint (optopt))
	        	{
	        		std::cout << "Unknown option" << optopt << ".\n";
	        	}
	        	else
	        	{
	        		fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
	        	}
	        	return 1;
	        	break;

	        default:
	        	abort ();
	        	break;

	    	} // END SWITCH
	    } // END WHILE




		axis_CoSiSt x_axis(-0.5*Lx, Lx, Nx, Lx*(1.-0.3)/(2*pi) );
		axis_CoSiSt y_axis(-0.5*Ly, Ly, Ny, Ly*(1.-0.3)/(2*pi) );
		axis_CoSiSt z_axis(-0.5*Lz, Lz, Nz, Lz*(1.-0.3)/(2*pi) );


	/*
	axis_CoSiSt x_axis(-0.5*Lx, Lx, Nx, Lx*(1.-0.25)/(2*pi) );
	axis_CoSiSt y_axis(-0.5*Ly, Ly, Ny, Ly*(1.-0.25)/(2*pi) );
	axis_CoSiSt z_axis(-0.5*Lz, Lz, Nz, Lz*(1.-0.25)/(2*pi) );

	axis_CoHySt x_axis(-5.,10.,Nx,3);
	axis_CoHySt y_axis(-5.,10.,Ny,3);
	axis_CoHySt z_axis(-5.,10.,Nz,3);

	axis_CoEqSt x_axis(-0.5*Lx,Lx,Nx);
	axis_CoEqSt y_axis(-0.5*Ly,Ly,Ny);
	axis_CoEqSt z_axis(-0.5*Lz,Lz,Nz);
	*/


	grid Omega(x_axis,y_axis,z_axis);



   #if defined(_MY_VERBOSE)
	my_log << "init fields";
   #endif
	field_real Ux(Omega), Uy(Omega), Uz(Omega);
	field_real ni(Omega);
	field_real Ph(Omega);

	field_imag Fni(Omega);
	field_imag FPh(Omega);

	//create_input_from_MGsolver(ni, Ph);
	//create_input_from_old_data(Ux, Uy, Uz, ni, Ph);

	OP_FFT my_FFT(Omega);
	OP_iFFT my_iFFT(Omega);

	my_FFT(ni,Fni);
	my_FFT(Ph,FPh);

	dealiasing_undesignated(Fni, [] (double k) {return exp(-2000.*pow(k,15.));});
	dealiasing_undesignated(FPh, [] (double k) {return exp(-2000.*pow(k,15.));});

	my_iFFT(Fni,ni);
	my_iFFT(FPh,Ph);

	ni.fill2([&] (double x, double y, double z) {return 1.;});
	//ni.fill2([&] (double x, double y, double z) {return 1.+exp(-(x*x)/0.1);});
	//Ux.fill2([&] (double x, double y, double z) {return 0.2*sin(2*3.14152*x/8.);});
	Ph.fill2([&] (double x, double y, double z) {return 0.;});
	Ux.fill2([&] (double x, double y, double z) {return 0.;});
	Uy.fill2([&] (double x, double y, double z) {return 0.;});
	Uz.fill2([&] (double x, double y, double z) {return 0.;});

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

