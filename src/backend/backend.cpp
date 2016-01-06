// gcc standard C++ libraries (common stuff)
#include <iostream>  // standard I/O-Operations (I/O to console)
#include <fstream>   // same as iostream but for files
#include <cmath>     // basic math stuff (sqrt(), cabs(), pow(a,n) etc.)
#include <string>    // class : string (here: used for filenames only; Maybe I should avoid this library)
#include <sstream>   // class : string-stream (MORE OVERHEAD) I'am using this to interpret data read from a file
#include <iomanip>   // sets format of output-stream (scientific & fixed) (is this really necessary??)

#include "parameters.hpp"
#include "grid.hpp"
//#include "particle.hpp"
#include "plot.hpp"
#include "IO.hpp"
#include "o_string.hpp"

#include "fftw3.h"













void help(void)
{
	std::cout << "backend to kielflow v0.1" << std::endl;
}



int main(int argc, char *argv[])
{
   #ifdef _MY_VERBOSE
	logger log("backend");
	log << "start";
   #endif


	std::string file_name_h5 = "./data/fields.h5";
	int slice_flag = 0;


	int switch_opt;

    opterr = 0;

    while ((switch_opt = getopt (argc, argv, "hf:s")) != -1)
    	switch (switch_opt)
    	{

    	case 'h':
        help();
        return 0;
          break;

        case 'f':
        file_name_h5 = optarg;
          break;

        case 's':
        slice_flag = 1;
          break;

        case '?':
        if (optopt == 'f')
        {
        	std::cout << "Option" << optopt << " requires an argument.\n";
        }
        else if (isprint (optopt))
        {
        	std::cout << "Unknown option" << optopt << ".\n";
        }
        else
        {
        	fprintf (stderr,
        			"Unknown option character `\\x%x'.\n",
                     optopt);
        }
        return 1;
          break;

        default:
        abort ();
          break;

    	} // END SWITCH





   #ifdef _MY_VERBOSE
	log << "load domain";
   #endif
	grid_Co Omega = load_grid(file_name_h5);


	std::cout << "achsen" << std::endl;
	std::cout << Omega.x_axis->m << std::endl;
	std::cout << Omega.y_axis->m << std::endl;
	std::cout << Omega.z_axis->m << std::endl;


	for(int i=0; i<Omega.Nx; ++i)
		std::cout << Omega.x_axis->val_at(i) << std::endl;

   #ifdef _MY_VERBOSE
	log << "load parameters";
   #endif
	parameters Params = load_parameters(file_name_h5);


	// ToDo : Particle missing here
   #ifdef _MY_VERBOSE
	log << "load particles";
   #endif
	//particle::load();






	if(slice_flag==true)
	{
	   #ifdef _MY_VERBOSE
		log << "slice_flag==true";
	   #endif

		double * Ux = new double [Omega.Nx*Omega.Ny];
		double * Uy = new double [Omega.Nx*Omega.Ny];
		double * Uz = new double [Omega.Nx*Omega.Ny];
		double * ni = new double [Omega.Nx*Omega.Ny];
		double * Ph = new double [Omega.Nx*Omega.Ny];



		load_slice("Ux", Ux,file_name_h5);
		load_slice("Uy", Uy,file_name_h5);
		load_slice("Uz", Uz,file_name_h5);
		load_slice("ni", ni,file_name_h5);
		load_slice("Ph", Ph,file_name_h5);

		std::string file_name_dat = ChangeFileExtension(file_name_h5,".dat");
		save_frame(Omega, Ux, Uy, Uz, ni, Ph, file_name_dat);
		save_density_pl4 (Omega, Params, file_name_dat);
		save_velocity_pl4(Omega, Params, file_name_dat);

		std::cout << "done.";
       #ifdef _MY_VERBOSE
		log << "done";
       #endif
		return 0;
	}


	// Load fields (real-space) from file
   #ifdef _MY_VERBOSE
	log << "load fields (real-space)";
   #endif
	field_real Ux = load_field_real(file_name_h5,"Ux");
	field_real Uy = load_field_real(file_name_h5,"Uy");
	field_real Uz = load_field_real(file_name_h5,"Uz");
	field_real ni = load_field_real(file_name_h5,"ni");
	field_real Ph = load_field_real(file_name_h5,"Ph");



	plot_2d::default_plane = 2;       // 0=x-fixed; 1=y fixed; 2 = z fixed
	plot_2d::default_direction = 0;   // along x-axis (0), y-axis(1), z-axis(2)
	plot_2d::default_xpos = Omega.x_axis->index_at(0.0);
	plot_2d::default_ypos = Omega.y_axis->index_at(0.0);
	plot_2d::default_zpos = Omega.z_axis->index_at(0.0);


	subdim save_opt;

	save_opt.default_plane = plot_2d::default_plane;
	save_opt.default_direction = plot_2d::default_direction;
	save_opt.default_xpos = plot_2d::default_xpos;
	save_opt.default_ypos = plot_2d::default_ypos;
	save_opt.default_zpos = plot_2d::default_zpos;

	#ifdef _MY_VERBOSE
	log << "saving 2D plot data";
   #endif
	save_2d(Ux, Uy, Uz, ni, Ph, save_opt, "./data/splot_data.dat");

   #ifdef _MY_VERBOSE
	log << "saving 1D plot data";
   #endif
	save_1d(Ux, Uy, Uz, ni, Ph, save_opt, "./data/plot_data.dat");





	density_plot PL_density(ni,Params);
	PL_density.auto_set(ni);
	PL_density.save();


	potential_plot PL_potential(Ph,Params);
	PL_potential.auto_set(Ph);
	PL_potential.save();


	velocity_plot PL_velocity(Ux,Uy,Uz,Params);
	PL_velocity.save();



   #ifdef _MY_VERBOSE
	log << "done. backend terminated";
   #endif
	std::cout << "backend terminated" << std::endl;

	return 0;
}

