#include "solver_poisson_multigrid.hpp"



solver_poisson_multigrid::solver_poisson_multigrid()
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_multigrid::solver_poisson_multigrid()");
	my_log << "start";
   #endif

	ptr_relaxation_method = 0L;

	I_xto2x = &OP_xto2x;
	I_yto2y = &OP_yto2y;
	I_zto2z = &OP_zto2z;
	I_hto2h = &OP_hto2h;

	I_2xtox = &OP_2xtox;
	I_2ytoy = &OP_2ytoy;
	I_2ztoz = &OP_2ztoz;
	I_2htoh = &OP_2htoh_lvl0;

	my_cascades = 0;
	my_lvl = 0L;
	my_steps = 0L;
	my_tolerance = 0L;

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

	I_xto2x = &OP_xto2x;
	I_yto2y = &OP_yto2y;
	I_zto2z = &OP_zto2z;
	I_hto2h = &OP_hto2h;

	I_2xtox = &OP_2xtox;
	I_2ytoy = &OP_2ytoy;
	I_2ztoz = &OP_2ztoz;
	I_2htoh = &OP_2htoh_lvl0;

	my_cascades = 1;
	my_lvl = new MG_lvl_control[1]{lvl_keep};
	my_steps = new int[1] {1};
	my_tolerance = new double[1] {1};

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
}


solver_poisson_multigrid::~solver_poisson_multigrid()
{
	delete[] my_lvl;
	delete[] my_steps;
	delete[] my_tolerance;

   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_multigrid::~solver_poisson_multigrid()");
	my_log << "~solver_poisson_multigrid()";
   #endif
}


void solver_poisson_multigrid::set_level_control(const int &N, int * lvl_steps, MG_lvl_control * lvl_C,
		double * lvl_tol)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_multigrid::set_level_control(..)");
	my_log << "start";
   #endif

	delete[] my_lvl;
	delete[] my_steps;
	delete[] my_tolerance;

	my_cascades = N;
	my_steps = new int[my_cascades];
	my_lvl = new MG_lvl_control[my_cascades];
	my_tolerance = new double[my_cascades];

	int lvl_current = 0;

	for(int cascade=0; cascade<my_cascades; cascade++)
	{ // copy cycle shape
		my_steps[cascade] = lvl_steps[cascade];
		my_lvl[cascade] = lvl_C[cascade];
		my_tolerance[cascade] = lvl_tol[cascade];

		// following stuff is for to avoid input errors
		if(my_lvl[cascade]==lvl_up_x)
			lvl_current++;
		if(my_lvl[cascade]==lvl_up_y)
			lvl_current++;
		if(my_lvl[cascade]==lvl_up_z)
			lvl_current++;
		if(my_lvl[cascade]==lvl_up)
			lvl_current++;
		if(my_lvl[cascade]==lvl_down_x)
			lvl_current--;
		if(my_lvl[cascade]==lvl_down_y)
			lvl_current--;
		if(my_lvl[cascade]==lvl_down_z)
			lvl_current--;
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


int solver_poisson_multigrid::refine_mesh(const MG_lvl_control &step, const field_real &in, field_real &out)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_multigrid::refine_mesh(..)");
	my_log << "start";
   #endif
	switch(step)

	{
	// actions that refine the grid in all three space directions
    case lvl_keep :
    	out = in;        // just copy the input
       #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
        my_log << "lvl_keep";
       #endif
        return 0;

	case lvl_up :
		I_hto2h(in,out);
       #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
        my_log << "lvl_up";
       #endif
		return 1;

	case lvl_down :
		I_2htoh(in,out); // average potential on current mesh
       #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
        my_log << "lvl_down";
       #endif
		return -1;       // to coarser mesh if cycle[i]=0

	// actions that apply only to one direction in space
	case lvl_up_x :
		I_xto2x(in,out);
       #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
        my_log << "lvl_up_x";
       #endif
		return 1;

	case lvl_up_y :
		I_yto2y(in,out);
       #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
        my_log << "lvl_up_y";
       #endif
		return 1;

	case lvl_up_z :
		I_zto2z(in,out);
       #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
        my_log << "lvl_up_z";
       #endif
		 return 1;

	case lvl_down_x :
		I_2xtox(in,out);
       #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
        my_log << "lvl_down_x";
       #endif
		return -1;

	case lvl_down_y :
		I_2ytoy(in,out);
       #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
        my_log << "lvl_down_y";
       #endif
		return -1;

	case lvl_down_z :
		I_2ztoz(in,out);
       #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
        my_log << "lvl_down_z";
       #endif
		return -1;

	default :
		throw("no instruction to handle MG");
		break;
	} // end of switch

	return 0;
}

void solver_poisson_multigrid::solve(field_real &Phi_IO, field_real &rho)
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_multigrid::solve(..)");
	my_log << "start";
   #endif

	int lvl_current = 0;
	int cascade_current = 0;
	field_real Phi_n(Phi_IO.my_grid);
	field_real rho_n(rho.my_grid);
	rho_n = rho; // only important for the first run per call

   #if defined (_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
    std::stringstream IO_stream;
    std::string IO_string;
    IO_stream << "MG cascade " << (cascade_current+1) << "/" << my_cascades;
    IO_stream << " level " << lvl_current;
    IO_string = IO_stream.str();
    std::cout << IO_string << std::endl;
    my_log << IO_string;
#endif

	for(cascade_current=0; cascade_current<my_cascades; ++cascade_current)
	{
       #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
		std::stringstream msg;
		msg << "cascade: " << cascade_current;
	    my_log << msg.str();
       #endif
		// refine mesh for the Potential according to given plan saved in my_lvl
		refine_mesh(my_lvl[cascade_current],Phi_IO,Phi_n);

		// important note :
		// charge density on a mesh can never be recovered from coarser grid.
		// Therefore the density has to be derived
		// from the highest level of resolution down to the current level

		// check if current grid of Phi is equal to grid of charge density
		// otherwise return to finest grid for charge density
		if( (rho_n.Nx != Phi_n.Nx) || (rho_n.Ny != Phi_n.Ny) || (rho_n.Nz != Phi_n.Nz) )
		{
		   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
			my_log << "moving rho_n to original grid";
		   #endif
			rho_n.resize(rho.Nx,rho.Ny,rho.Nz);
			rho_n = rho;
		}

		while( (rho_n.Nx != Phi_n.Nx) &&
			   (rho_n.Ny != Phi_n.Ny) &&
			   (rho_n.Nz != Phi_n.Nz) )
		{
           #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	        my_log << "reduce rho_n in all directions";
           #endif
			field_real bc(rho_n.my_grid);
			bc = rho_n;
			I_2htoh(bc,rho_n);
		}

		while( (rho_n.Nx != Phi_n.Nx) )
		{
           #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
            my_log << "reduce rho_n in x-directions";
           #endif
			field_real bc(rho_n.my_grid);
			bc = rho_n;
			I_2xtox(bc,rho_n);
		}

		while( (rho_n.Ny != Phi_n.Ny) )
		{
           #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
            my_log << "reduce rho_n in x-directions";
           #endif
			field_real bc(rho_n.my_grid);
			bc = rho_n;
			I_2ytoy(bc,rho_n);
		}

		while( (rho_n.Nz != Phi_n.Nz) )
		{
           #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
            my_log << "reduce rho_n in z-directions";
           #endif
			field_real bc(rho_n.my_grid);
			bc = rho_n;
			I_2ztoz(bc,rho_n);
		}

		// run relaxation solver on current level
		ptr_relaxation_method->set_max_iterations(my_steps[cascade_current]);
		ptr_relaxation_method->set_tolerance(my_tolerance[cascade_current]);
		ptr_relaxation_method->solve(Phi_n,rho_n);

		/*
		 * ToDo : why didn't it work properly
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

			//my_lvl[cascade_current+cascade_skips] = lvl_keep;
			//cascade_current = cascade_current+cascade_skips-1;

			//my_lvl[cascade_current+cascade_skips] = lvl_keep;
           #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	        my_log << "skipping cascades:";
	        my_log << cascade_skips;
           #endif
			cascade_current = cascade_current+cascade_skips;
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


