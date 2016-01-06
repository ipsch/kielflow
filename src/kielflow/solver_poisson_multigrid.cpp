#include "solver_poisson_multigrid.hpp"



solver_poisson_multigrid::solver_poisson_multigrid()
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_multigrid::solver_poisson_multigrid()");
	my_log << "start";
   #endif

	ptr_relaxation_method = 0L;

	I_hto2h = &OP_hto2h;
	I_2htoh = &OP_2htoh_lvl0;

	my_cascades = 0;
	my_lvl = 0L;
	my_steps = 0L;

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
}


solver_poisson_multigrid::solver_poisson_multigrid(interface_relaxation_solver &method)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_multigrid::solver_poisson_multigrid(interface_relaxation_solver &method)");
	my_log << "start";
   #endif

	ptr_relaxation_method = &method;

	I_hto2h = &OP_hto2h;
	I_2htoh = &OP_2htoh_lvl0;

	my_cascades = 0;
	my_lvl = new  MG_lvl_control{lvl_keep};
	my_steps = new int{1};

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
}


solver_poisson_multigrid::~solver_poisson_multigrid()
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_multigrid::~solver_poisson_multigrid()");
	my_log << "start";
   #endif
}


void solver_poisson_multigrid::set_level_control(const int &N, int * lvl_steps, MG_lvl_control * lvl_C)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_multigrid::set_level_control(..)");
	my_log << "start";
   #endif

	delete[] my_lvl;
	delete[] my_steps;

	my_cascades = N;
	my_steps = new int[my_cascades];
	my_lvl = new MG_lvl_control[my_cascades];

	int lvl_current = 0;

	for(int cascade=0; cascade<my_cascades; cascade++)
	{ // copy cycle shape
		my_steps[cascade] = lvl_steps[cascade];
		my_lvl[cascade] = lvl_C[cascade];

		// following stuff is for to avoid input errors
		if(my_lvl[cascade]==lvl_up)
			lvl_current++;
		if(my_lvl[cascade]==lvl_down)
			lvl_current--;
		if(lvl_current>0)
		{
			std::stringstream error_msg;
			std::string error_string;
			error_msg << "error while setting multigrid cycle:\n";
			error_msg << "cycle level became greater than base zero (base level).\n";
			error_string = error_msg.str();
           #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	        my_log << error_string;
           #endif
	        throw(error_string);
		}
	}  // END copy cycle shape

   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
	return;
}


void solver_poisson_multigrid::solve(field_real &Phi_IO, field_real &rho)
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_multigrid::solve(..)");
	my_log << "start";
   #endif

	int lvl_current = 0;
	int cascade_current = 0;
	field_real Phi_n(*Phi_IO.my_grid);
	field_real rho_n(*rho.my_grid);
	rho_n = rho; // only important for the first run per call

	for(cascade_current=0; cascade_current<my_cascades; ++cascade_current)
	{
		if(my_lvl[cascade_current]==lvl_keep)
		{
			Phi_n = Phi_IO; // just copy the input
		}

		if(my_lvl[cascade_current]==lvl_down)
		{
           #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
			my_log << "moving Phi_n to coarse grid";
		   #endif
			// average potential on current mesh
			// to coarser mesh if cycle[i]=0
			I_2htoh(Phi_IO,Phi_n);
			lvl_current--;
		}

		if(my_lvl[cascade_current]==lvl_up)
		{
           #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	        my_log << "moving Phi_n to fine grid";
           #endif
			// interpolate potential on current mesh
			// to finer mesh if cycle[i]=1
			I_hto2h(Phi_IO,Phi_n);
			lvl_current++;
		}

       #if defined (_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
		std::stringstream IO_stream;
		std::string IO_string;
		IO_stream << "MG cascade " << (cascade_current+1) << "/" << my_cascades;
		IO_stream << " level " << lvl_current;
		IO_string = IO_stream.str();
		std::cout << IO_string << std::endl;
	    my_log << IO_string;
       #endif

		/*
		 * charge-density on finer mesh can never be recovered from
		 * coarse grid. Therefore the density has to be derived
		 * from the highest level of resolution down to the current level
		 */

		if(rho_n.N != Phi_n.N)
		{
		   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
			my_log << "moving rho_n to original grid";
		   #endif
			// check if current level is the highest level of resoultion
			rho_n.resize(rho.Nx,rho.Ny,rho.Nz);
			rho_n = rho;
		}

		while(rho_n.N != Phi_n.N)
		{
           #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
			my_log << "moving rho_n to coarse grid";
           #endif
			// average density on current mesh
			// to coarser mesh until
			// resolution of potential and density
			// are equal
			field_real bc(*rho_n.my_grid);
			bc = rho_n;
			I_2htoh(bc,rho_n);
		};

		// run relaxation solver on current level
		ptr_relaxation_method->set_max_iterations(my_steps[cascade_current]);
		ptr_relaxation_method->solve(Phi_n,rho_n);



		// ToDo : Fast Forward code zum überspringen von Kaskaden.
		/*
		// code zum fast forwarden: Die Rechnung auf einem gröberen Gitter kann
		// übersprungen werden, wenn:
		// - die Relaxation auf diesem Gitter bereits konvergiert ist
		// - in dieser Kaskade die Gittergröße nicht geändert wurde
		// - noch mindestens zwei Kaskaden aus zu führen sind
		// - die nächste Kaskade das Gitter verkleinern würde
		if( (my_lvl[cascade_current]==lvl_keep) && ptr_relaxation_method->converged() &&
				((cascade_current+1) < my_cascades) )
		{
			int lvl_change = 0;
			int cascade_skips = 0;

			do
			{
				cascade_skips++;
				if(my_lvl[cascade_current+cascade_skips]==lvl_up)
					lvl_change++;
				if(my_lvl[cascade_current+cascade_skips]==lvl_down)
					lvl_change--;

			} while(lvl_change!=0);

			my_lvl[cascade_current+cascade_skips] = lvl_keep;
			cascade_current = cascade_current+cascade_skips-1;

		}
		*/

		Phi_IO.resize(Phi_n.Nx, Phi_n.Ny, Phi_n.Nz);
		Phi_IO = Phi_n;
	}

	if(Phi_IO.N != rho.N)
	{
		std::cout << "MG kehrte nicht zum original Gitter zurück" << std::endl;
		throw("MG kehrte nicht zum original Gitter zurück");
	}

   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return;

}


