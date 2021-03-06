 //#define MOVABLE

// set (git)-repository-Version to unknown
// (if not defined in makefile during compile process)
#ifndef VERSION_STRING
#define VERSION_STRING "unknown"
#endif

#if defined(_MY_VERBOSE_LESS) || defined(_MY_VERBOSE) || \
	defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif

// gcc standard C++ libraries (common stuff)
#include <iostream>  // standard I/O-Operations (I/O to console)
#include <cmath>     // basic math stuff (sqrt(), cabs(), pow(a,n) etc.)
#include <string>    // class : string (here: used for filenames only; Maybe I should avoid this library)
#include <sstream>   // class : string-stream (MORE OVERHEAD) I'am using this to interpret data read from a file



// components
#include "parameters.hpp"
#include "grid.hpp"
#include "field_real.hpp"
#include "field_imag.hpp"
#include "IO.hpp"
#include "counter.hpp"
#include "operations.hpp"

// solver for poission equation
#include "solver_poisson_multigrid.hpp"
#include "solver_poisson_jacobi_lin.hpp"
#include "solver_poisson_jacobi_nlin.hpp"

// evaluation of right hand side
#include "rhs.hpp"
#include "rhs_standard.hpp"
//#include "rhs_lab.hpp"

// time integration
#include "euler_method.hpp"
#include "RKO4.hpp"





#include <omp.h>







// main #######################################################################

int main(int argc,char **argv)
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("main");
	my_log.set_file("./kielstream.log");
	my_log << "running version";
	my_log << VERSION_STRING;
   #endif

	int iterations= 3000;
	counter i_output(10);
	counter i_backup(50);

	double charge_Q = -11498.5;
	double radius_a = 0.15;
    bool override_M = false;
    bool override_theta = false;
    bool override_tau = false;
    bool override_beta = false;
    parameters override_Params;
    std::string IO_file = "./data/fields.h5";

	int switch_opt;
	int opt_indent = 0;
    opterr = 0;
    while ((switch_opt = getopt (argc, argv, "hM:Q:T:t:b:")) != -1)
    {
    	opt_indent+=2;
    	switch (switch_opt)
    	{
        case 'a':
        	radius_a = std::atof(optarg);
        	break;

        case 'Q':
        	charge_Q = std::atof(optarg);
        	break;

        case 'M':
        	override_M = true;
        	override_Params.M = std::atof(optarg);
        	break;

        case 'T':
        	override_theta = true;
        	override_Params.theta = std::atof(optarg);
        	break;

        case 't':
        	override_tau = true;
        	override_Params.tau = std::atof(optarg);
        	break;

        case 'b':
        	override_beta = true;
        	override_Params.beta = std::atof(optarg);
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

    for(int i=opt_indent+1; i<argc; ++i)
    	IO_file = argv[i];


	std::cout << "IO_file " << IO_file << std::endl;
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "loading input from";
	my_log << IO_file;
   #endif
	grid Omega     = load_grid(IO_file);
	double t_total    = load_time(IO_file);
	parameters Params = load_parameters(IO_file);
	field_imag FUx    = load_field_imag("Ux", IO_file);
	field_imag FUy    = load_field_imag("Uy", IO_file);
	field_imag FUz    = load_field_imag("Uz", IO_file);
	field_imag Fni    = load_field_imag("ni", IO_file);
	field_imag FPh    = load_field_imag("Ph", IO_file);


	if(override_M) Params.M = override_Params.M;
	if(override_theta) Params.theta = override_Params.theta;
	if(override_tau) Params.tau = override_Params.tau;
	if(override_beta) Params.beta = override_Params.beta;


    // ##### ION DENSITY TO LOGSCALE ######
    field_real ni_tmp(Omega);
	field_imag Fni_tmp(Omega);

	iFFT(Fni,ni_tmp);

	for(int i=0; i<ni_tmp.N; ++i)
	{
		ni_tmp.val[i] = log(ni_tmp.val[i]);
	}
	FFT(ni_tmp,Fni);
	//dealiasing_undesignated(Fni, [] (const double &k) {return k<(1./3.) ? 1. : 0.; });

	// ##### DUST #####
	double scale_Q = 8.1720e-06; // umrechnungsfaktor auf dimensionslose Einheiten
	field_real nd(Omega);
	fkt3d_Gauss dust_3d_fkt(charge_Q*scale_Q,radius_a,radius_a,radius_a);
	nd.fill(dust_3d_fkt);

	field_real Hd(Omega);
	fkt3d_Gauss dust_3d_mask(1.,radius_a,radius_a,radius_a);
	Hd.fill(dust_3d_mask);


	// ##### RELAXATIONS-SOLVER #####
	fkt3d_const boundary_shape(0.); // there are no additional boundarys (except Domain borders)
	fkt3d_const boundary_value(0.); // the potential at the boundary is zero
	double SOR = .8; // Overrelaxation parameter
	solver_poisson_jacobi_nlin NLJ(boundary_shape, boundary_value, SOR);
	NLJ.set_max_iterations(80);
	NLJ.set_limit_NormMax(.5e-5);
	NLJ.set_limit_NormSum(.5e-7);


	// ##### MULTIGRID #####
	solver_poisson_multigrid MG(NLJ);
	const int MG_N = 4;

	int * MG_steps_sizes = new int[MG_N]
              {0,-1,-1,-1};
	MG_lvl_control * cycle_shape = new MG_lvl_control[MG_N]
			  { lvl_down3, lvl_up, lvl_up, lvl_up};
	double * MG_error = new double[MG_N]
			  { 0.01, 0.0005, 0.0001, 0.0001};

	MG.set_level_control(MG_N, MG_steps_sizes, cycle_shape, MG_error);

	delete[] MG_steps_sizes;
	delete[] MG_error;
	delete[] cycle_shape;




	// ##### RHS #####
	fkt3d_barrier Barrier_fkt(nd.my_grid.x_axis->val_at(0), -3.);
	field_real Barrier(Omega);
	Barrier.fill2([] (const double &x, const double &y, const double &z) {return exp(-pow(x/0.15,2.) - pow(y/0.15,2.) - pow(z/0.15,2.));;});
    //Barrier.fill2([] (const double &x, const double &y, const double &z) {return 0.;});
	field_real Ph(FPh.my_grid);
	iFFT(FPh, Ph);
	rhs_standard rhs(Params, MG, Ph, nd, Barrier, Omega);
	//rhs_lab rhs(Params, MG, Ph, nd, Barrier);
   #if defined(TEST_RHS)
	//rhs.solve(FUx, FUy, FUz, Fni);
   #endif
	// ##### TIME-INTEGRATOR #####
	double t_delta = 0.01;
	Runge_kutta_O4 time_integrator(rhs, t_delta);
	//euler_method time_integrator(rhs, t_delta);


   //#define TEST_MULTIGRID
   #ifdef TEST_MULTIGRID
	// ##### MULTIGRID TESTING #####
	// uncomment for testing purposes only
	std::cout << "MG testing" << std::endl;
	field_real ni(*FPh.my_grid);
	field_real rho(*FPh.my_grid);
	iFFT(Fni,ni);
	for(int i=0; i<ni.N; ++i)
	{
		Ph.val[i] = 0.;
		ni.val[i] += nd.val[i];
	}
	MG.solve(Ph,ni);
   #endif

    /*
	// Remember to set #define _MONITOR_EVOLUTION in RKO4.cpp
	//iteration.set_range(1);
	field_real dummy(Omega);
	dummy.fill2([&] (double x, double y, double z) {return 0.;});
	//FFT(dummy, FUx);
	dummy.fill2([&] (double x, double y, double z) {return 0.;});
	//FFT(dummy, FUy);
	dummy.fill2([&] (double x, double y, double z) {return 0.;});
	//FFT(dummy, FUz);
	dummy.fill2([&] (double x, double y, double z) {return cos(2.*acos(-1.)/8.*x);});
	//FFT(dummy, Fni);
	*/

   // #define KILL_MAIN_LOOP
   #ifndef KILL_MAIN_LOOP
	// #### Time-Step #########################################################
	for(unsigned int step=0; step<iterations; step++, ++i_output, ++i_backup)
	{
		std::cout << "step(time) (iteration: " << step << ")" << std::endl;
       #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
		std::stringstream sstr;
		sstr << "step(time) (iteration: " << step << ")";
		my_log << sstr.str();
       #endif

		time_integrator.solve(FUx, FUy, FUz, Fni);


		// ##### SLICE OUTPUT #####
		if(!i_output.good())
		{
           #if defined (_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
			my_log << "saving slice";
           #endif

			FFT(Ph,FPh);

			// ION DENSITY TO LINEAR SCALE
			iFFT(Fni,ni_tmp);
			for(int i=0; i<ni_tmp.N; ++i)
				ni_tmp.val[i] = exp(ni_tmp.val[i]);
			FFT(ni_tmp, Fni_tmp);
			save_slice(step, Omega, t_total, Params, FUx, FUy, FUz, Fni_tmp, FPh);

			i_output.reset();
		}

		// ##### BACKUP OUTPUT #####
		if(!i_backup.good())
		{
           #if defined (_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	        my_log << "saving backup";
           #endif

			// ION DENSITY TO LINEAR SCALE
			iFFT(Fni,ni_tmp);
			for(int i=0; i<ni_tmp.N; ++i)
				ni_tmp.val[i] = exp(ni_tmp.val[i]);
			FFT(ni_tmp, Fni_tmp);

			std::stringstream sstr;
	        sstr.str(std::string());
			sstr << "./data/backup_at_" << step << ".h5";
			save_all(Omega, t_total, Params, FUx, FUy, FUz, Fni_tmp, FPh, sstr.str());

			i_backup.reset();
		}


		t_total += t_delta;
	}
   #endif // KILL_MAIN_LOOP



	// ##### ION DENSITY TO LINEAR SCALE #####
	iFFT(Fni,ni_tmp);
	for(int i=0; i<ni_tmp.N; ++i)
		ni_tmp.val[i] = exp(ni_tmp.val[i]);
	FFT(ni_tmp,Fni);

	//save_all(Omega, t_total, Params, FUx, FUy, FUz, Fni, FPh, IO_file);

	std::cout  << "finished - kielflow terminated" << std::endl;
   #if defined (_MY_VERBOSE) || defined (_MY_VERBOSE_MORE) || defined (_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return 0;
}



