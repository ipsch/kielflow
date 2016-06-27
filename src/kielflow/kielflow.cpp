 //#define MOVABLE

// set (git)-repository-Version to unknown
// (if not defined in makefile during compile process)
#ifndef VERSION_STRING
#define VERSION_STRING "unknown"
#endif

#if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif

// gcc standard C++ libraries (common stuff)
#include <iostream>  // standard I/O-Operations (I/O to console)

#include <cmath>     // basic math stuff (sqrt(), cabs(), pow(a,n) etc.)
#include <string>    // class : string (here: used for filenames only; Maybe I should avoid this library)
#include <sstream>   // class : string-stream (MORE OVERHEAD) I'am using this to interpret data read from a file

#include <vector>


// uncommon libraries
#include "fftw3.h"   // google it (something about discret fourier Transformations)

// own includes

#include "parameters.hpp"
#include "field.hpp"
#include "IO.hpp"
#include "counter.hpp"
#include "operations.hpp"


#include "particle.hpp"


#include "solver_poisson_jacobi_lin.hpp"
#include "solver_poisson_jacobi_nlin.hpp"
#include "solver_poisson_multigrid.hpp"

#include "rhs.hpp"
#include "rhs_standard.hpp"

#include "euler_method.hpp"
#include "RKO4.hpp"





#include <omp.h>
//double omp_get_wtime(void);








// main #######################################################################

int main(int argc,char **argv)
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("main");
	my_log << "running version";
	my_log << VERSION_STRING;
   #endif


    bool override_M = false;
    double default_Q = -11498.5;
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

        case 'M':
        	override_M = true;
        	override_Params.M = std::atof(optarg);
        	break;

        case 'Q':
        	default_Q = std::atof(optarg);
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

	counter iteration(5000);   // set counter for how many iterations are allowed
	counter i_output(1);
	counter i_backup(10);


   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "loading input";
   #endif
	grid_Co Omega     = load_grid(IO_file);
	double t_total    = load_time(IO_file);
	parameters Params = load_parameters(IO_file);
	field_imag FUx    = load_field_imag("Ux", IO_file);
	field_imag FUy    = load_field_imag("Uy", IO_file);
	field_imag FUz    = load_field_imag("Uz", IO_file);
	field_imag Fni    = load_field_imag("ni", IO_file);
	field_imag FPh    = load_field_imag("Ph", IO_file);
	std::vector<particle> particle_list;
	load_particles(particle_list, "./config/particles.dat");

	if(override_M) Params.M = override_Params.M;
	if(override_theta) Params.theta = override_Params.theta;
	if(override_tau) Params.tau = override_Params.tau;
	if(override_beta) Params.beta = override_Params.beta;


	// ##### DUST #####
	double scale_Q = 8.1720e-06; // umrechnungsfaktor auf dimensionslose Einheiten
	field_real nd(*FPh.my_grid);
	fkt3d_Gauss dust_3d_fkt(default_Q*scale_Q,0.15,0.15,0.15);
	nd.fill(dust_3d_fkt);

	field_real Hd(*FPh.my_grid);
	fkt3d_Gauss dust_3d_mask(1.,0.15,0.15,0.15);
	Hd.fill(dust_3d_mask);


	// ##### RELAXATIONS-SOLVER #####
	fkt3d_const boundary_shape(0.); // there are no additional boundarys (except Domain borders)
	fkt3d_const boundary_value(0.); // the potential at the boundary is zero
	double SOR = .8; // Overrelaxation parameter
	solver_poisson_jacobi_nlin NLJ(boundary_shape, boundary_value, SOR);
	NLJ.set_max_iterations(80);
	NLJ.limit_max = .5e-5; // precision in percent of amplitude
	NLJ.limit_sum = .5e-7;


	// ##### MULTIGRID #####
	solver_poisson_multigrid MG(NLJ);
	{ // setup MG-cycle
		const int MG_N = 6;
		int * MG_steps_sizes = new int[MG_N]
                  {1,1,1,-1,-1,-1};
	    MG_lvl_control * cycle_shape = new MG_lvl_control[MG_N]
				  {lvl_keep, lvl_down, lvl_keep, lvl_down, lvl_up, lvl_up};
	    double * MG_error = new double[MG_N]
				  {.01,.01,0.01,2.e-5,5.e-5,1.e-4};
	    MG.set_level_control(MG_N, MG_steps_sizes, cycle_shape, MG_error);
	    // Sketch of cycle:
	    // _
	    //  \_  /
	    //    \/
	    delete[] MG_steps_sizes;
	    delete[] MG_error;
	    delete[] cycle_shape;
	} // configure done



	// ##### RHS #####
	fkt3d_barrier Barrier_fkt(nd.my_grid->x_axis->val_at(0), -3.);
	field_real Barrier(*FPh.my_grid);
    Barrier.fill(Barrier_fkt);
	field_real Ph(*FPh.my_grid);
	iFFT(FPh, Ph);
	rhs_standard rhs(Params, MG, Ph, nd, Barrier);
   #if defined(TEST_RHS)
	//rhs.solve(FUx, FUy, FUz, Fni);
   #endif
	// ##### TIME-INTEGRATOR #####
	double t_delta = 0.01;
	Runge_kutta_O4 time_integrator(rhs, t_delta);


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



   //#define KILL_MAIN_LOOP
   #ifndef KILL_MAIN_LOOP
	// #### Time-Step #########################################################
	while(iteration.good())
	{
       #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
		std::stringstream sstr;
		sstr << "step(time) (iteration: " << (iteration.show()) << ")";
		my_log << sstr.str();
       #endif

		std::cout << "step(time) (iteration: " << (iteration.show()) << ")" << std::endl;
		time_integrator.solve(FUx, FUy, FUz, Fni);
		//rhs.solve(FUx, FUy, FUz, Fni);
		//MG.solve(Ph,ni);
		//NLJ.solve(Ph,ni);

		iteration.up();
		i_output.up();
		i_backup.up();

		// ##### SLICE OUTPUT #####
		if(!i_output.good())
		{
           #if defined (_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
			sstr << "step(time) (iteration: " << (iteration.show()) << ") saving step to file.";
			my_log << sstr.str();
			std::cout << sstr.str() << std::endl;
           #endif

			FFT(Ph,FPh);
			save_slice(iteration.show(),particle_list, Omega, t_total, Params, FUx, FUy, FUz, Fni, FPh);
			i_output.reset();

			/*
			sstr.str(std::string());
			sstr << "./data/FUx_" << (iteration.show()) << ".dat";
			save_major_wavevectors(sstr.str(), FUx);

			sstr.str(std::string());
			sstr << "./data/FUy_" << (iteration.show()) << ".dat";
			save_major_wavevectors(sstr.str(), FUy);

			sstr.str(std::string());
			sstr << "./data/FUz_" << (iteration.show()) << ".dat";
			sa#endif /* KILL_MAIN_LOOPve_major_wavevectors(sstr.str(), FUz);

			sstr.str(std::string());
			sstr << "./data/Fni_" << (iteration.show()) << ".dat";
			save_major_wavevectors(sstr.str(), Fni);
*/
		}

		// ##### BACKUP OUTPUT #####
		if(!i_backup.good())
		{
			// ToDo : also overwrite fields.h5
           #if defined (_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	        sstr << "step(time) (iteration: " << (iteration.show()) << ") saving step to file.";
	        my_log << sstr.str();
           #endif
	        save_all(particle_list, Omega, t_total, Params, FUx, FUy, FUz, Fni, FPh, "./data/fields.h5");
	        sstr.str(std::string());
			sstr << "./data/backup_at_" << (iteration.show()) << ".h5";
			save_all(particle_list, Omega, t_total, Params, FUx, FUy, FUz, Fni, FPh, sstr.str());
			i_backup.reset();
		}

		t_total += t_delta;
	}
   #endif // KILL_MAIN_LOOP

	save_all(particle_list, Omega, t_total, Params, FUx, FUy, FUz, Fni, FPh, IO_file);

	// ##### DONE #####



	std::cout  << "finished - kielflow terminated" << std::endl;
   #if defined (_MY_VERBOSE) || defined (_MY_VERBOSE_MORE) || defined (_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return 0;
}



