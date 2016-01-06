 //#define MOVABLE

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




#ifdef _MY_VERBOSE
#include "logger.hpp"
#endif





#include <omp.h>
//double omp_get_wtime(void);




double charge(const double t)
{
	double val = -0.55*0.5*t;
	if(val < -0.55)
		return -0.55;
	return val;
}





// main #######################################################################

int main(void)
{
   #ifdef _MY_VERBOSE
	logger my_log("main");
	my_log << "start";
   #endif

	counter iteration(5000);   // set counter for how many iterations are allowed
	counter i_output(1);
	counter i_backup(10);


   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "loading input";
   #endif
	grid_Co Omega     = load_grid("./data/fields.h5");
	double t_total    = load_time("./data/fields.h5");
	parameters Params = load_parameters("./data/fields.h5");
	field_imag FUx    = load_field_imag("Ux", "./data/fields.h5");
	field_imag FUy    = load_field_imag("Uy", "./data/fields.h5");
	field_imag FUz    = load_field_imag("Uz", "./data/fields.h5");
	field_imag Fni    = load_field_imag("ni", "./data/fields.h5");
	field_imag FPh    = load_field_imag("Ph", "./data/fields.h5");
	std::vector<particle> particle_list;
	load_particles(particle_list, "./config/particles.dat");



	// ##### DUST #####
	field_real nd(*FPh.my_grid);
	fkt3d_Gauss dust_3d_fkt(-14.,0.15,0.15,0.15);
	nd.fill(dust_3d_fkt);

	field_real Hd(*FPh.my_grid);
	fkt3d_Gauss dust_3d_mask(1.,0.15,0.15,0.15);
	Hd.fill(dust_3d_mask);


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


	// ##### RHS #####
	field_real Ph(*FPh.my_grid);
	iFFT(FPh, Ph);
	rhs_standard rhs(Params, MG, Ph, nd, Hd);


	// ##### TIME-INTEGRATOR #####
	double t_delta = 0.01;
	Runge_kutta_O4 time_integrator(rhs, t_delta);


	// ##### MULTIGRID TESTING #####
	// uncomment for testing purposes only
	/*
	std::cout << "MG testing" << std::endl;
	field_real ni(*FPh.my_grid);
	field_real rho(*FPh.my_grid);
	iFFT(Fni,ni);
	for(int i=0; i<ni.N; ++i)
	{
		Ph.val[i] = 0.;
		rho.val[i] = (ni.val[i]) + nd.val[i];
	}
	MG.solve(Ph,rho);
	*/





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
			save_major_wavevectors(sstr.str(), FUz);

			sstr.str(std::string());
			sstr << "./data/Fni_" << (iteration.show()) << ".dat";
			save_major_wavevectors(sstr.str(), Fni);
*/
		}

		// ##### BACKUP OUTPUT #####
		if(!i_backup.good())
		{
           #if defined (_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	        sstr << "step(time) (iteration: " << (iteration.show()) << ") saving step to file.";
	        my_log << sstr.str();
           #endif
	        sstr.str(std::string());
			sstr << "./data/backup_at_" << (iteration.show()) << ".h5";
			save_all(particle_list, Omega, t_total, Params, FUx, FUy, FUz, Fni, FPh, sstr.str());
			i_backup.reset();
		}

		save_all(particle_list, Omega, t_total, Params, FUx, FUy, FUz, Fni, FPh);

		t_total += t_delta;
	}



	//rhs.solve(FUx, FUy, FUz, Fni);

	save_all(particle_list, Omega, t_total, Params, FUx, FUy, FUz, Fni, FPh);

	// ##### DONE #####
	std::cout  << "finished - kielflow terminated" << std::endl;
   #if defined (_MY_VERBOSE) || defined (_MY_VERBOSE_MORE) || defined (_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return 0;
}



