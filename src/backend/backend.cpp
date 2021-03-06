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
#include "IO.hpp"
#include "o_string.hpp"
#include "subdim.hpp"


#include "fftw3.h"
#include "operations.hpp"
#include "OP_FFT.hpp"

#include "field_integrate.hpp"


/*

 #define POISSON
 #define ADVECTION
 #define DIFFUSION
 #define E_INT
 #define E_EXT
 #define DISSIPATION
 //#define PENALIZATION_U
 #define CONTINUITY
 //#define DEALAISING_MORE
 #define SPECTRAL_VISCOSITY
#define COMBINED
#define SAVE_TOTAL

*/




void help(void)
{
	std::cout << "backend to kielflow v0.1" << std::endl;
	std::cout << "o | output : path to output." << std::endl;
	std::cout << "    filename will be appended with data type" << std::endl;
}



int main(int argc, char *argv[])
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("backend");
	my_log << "start";
   #endif


	std::string file_input = "./data/fields.h5";
	std::string file_output_1d = "./data/plot_data.dat";
	std::string file_output_2d = "./data/splot_data.dat";
	std::string file_output_F1d = "./data/plot_Fdata.dat";
	std::string file_output_F2d = "./data/splot_Fdata.dat";
	bool slice_flag = false;
	bool use_target_dir = false;
	bool fourier_flag = false;
	std::string target_dir = "./data/";

	int N_files = 1;
	std::string * files_in = new std::string[1] {"./data/fields.h5"};

	int switch_opt;
	int opt_indent;
    opterr = 0;

    while ((switch_opt = getopt (argc, argv, "Fhsd:")) != -1)
    {

    	switch (switch_opt)
    	{

    	case 'F': {
    		fourier_flag = true;
    		opt_indent+=1;
    		break;
    	}

    	case 'h': {
    		help();
    		return 0;
    		break;
    	}

    	case 'd' : {
    		use_target_dir = true;
    		target_dir = optarg;
    		opt_indent+=2;
    		break;
    	}

    	case 's': {
    		slice_flag = true;
    		opt_indent+=1;
    		break;
    	}

        case '?': {
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
        		fprintf (stderr,
        			"Unknown option character `\\x%x'.\n",
                     optopt);
        	}
        	return 1;
        	break;
        }

        default: {
        	abort ();
        	break;
        }

    	} // END SWITCH
	} // END WHILE

    if(opt_indent+1<argc)
    {
    	delete[] files_in;
    	N_files = argc-(opt_indent+1);
    	files_in = new std::string[N_files];
        for(int i=0; i<N_files; ++i)
        	files_in[i] = argv[opt_indent+1+i];
    }

    for(int i=0; i<N_files; ++i)
    {

    	file_input = files_in[i];

    	std::string t_dir_;
    	std::string t_file_;
    	t_dir_ = ExtractDirectory(file_input);
    	t_file_ = ExtractFilename(file_input);
		t_file_ = ChangeFileExtension(t_file_, ".dat");

		if(use_target_dir)
		{
			file_output_1d = target_dir + "plot_data " + t_file_;
			file_output_2d = target_dir + "splot_data " + t_file_;
			file_output_F1d = target_dir + "plot_Fdata " + t_file_;
			file_output_F2d = target_dir + "splot_Fdata " + t_file_;
		}
		else
		{
			file_output_1d = t_dir_ + "plot_data " + t_file_;
			file_output_2d = t_dir_ + "splot_data " + t_file_;
			file_output_F1d = t_dir_ + "plot_Fdata " + t_file_;
			file_output_F2d = t_dir_ + "splot_Fdata " + t_file_;
		}

	    // status msg to console
	    std::cout << "backend: input: " << file_input << std::endl;
	    std::cout << "backend: output: " << file_output_1d << std::endl;
	    std::cout << "                 " << file_output_2d << std::endl;

	    // load parameters and domain
		grid Omega = load_grid(file_input);
		parameters Params = load_parameters(file_input);
		//particle::load();

		// set backend to process sliced data (2d instead of 3d data-files)
		if(slice_flag==true)
		{
			std::cout << "create data from sliced h5-file" << std::endl;
		   #if defined (_MY_VERBOSE_MORE) || defined (_MY_VERBOSE_TEDIOUS)
			my_log << "slice_flag==true";
		   #endif

			double * Ux = new double [Omega.Nx*Omega.Ny];
			double * Uy = new double [Omega.Nx*Omega.Ny];
			double * Uz = new double [Omega.Nx*Omega.Ny];
			double * ni = new double [Omega.Nx*Omega.Ny];
			double * Ph = new double [Omega.Nx*Omega.Ny];

			load_slice("Ux", Ux,file_input);
			load_slice("Uy", Uy,file_input);
			load_slice("Uz", Uz,file_input);
			load_slice("ni", ni,file_input);
			load_slice("Ph", Ph,file_input);

			save_frame_1d(Omega, Ux, Uy, Uz, ni, Ph, file_output_1d);
			save_frame_2d(Omega, Ux, Uy, Uz, ni, Ph, file_output_2d);

			delete[] Ux;
			delete[] Uy;
			delete[] Uz;
			delete[] ni;
			delete[] Ph;
		} // end processing 2d (slice) data

		if(slice_flag==false)
		{
			std::cout << "create data from full h5-file" << std::endl;
		   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
			my_log << "load fields (real-space)";
		   #endif

			// ToDo : probably memory leak in following block
			field_real Ux = load_field_real(file_input,"Ux");
			field_real Uy = load_field_real(file_input,"Uy");
			field_real Uz = load_field_real(file_input,"Uz");
			field_real ni = load_field_real(file_input,"ni");
			field_real Ph = load_field_real(file_input,"Ph");

			subdim save_opt;

			save_opt.plane = 2;
			save_opt.direction = 0;
			save_opt.xpos = Omega.x_axis->index_at(0.0);
			save_opt.ypos = Omega.y_axis->index_at(0.0);
			save_opt.zpos = Omega.z_axis->index_at(0.0);
			save_opt.Nx = Omega.x_axis->N;
			save_opt.Ny = Omega.y_axis->N;
			save_opt.Nz = Omega.z_axis->N;
			save_opt.set();

		   #if defined(_MY_VERBOSE_MORE) || defined (_MY_VERBOSE_TEDIOUS)
			my_log << "saving 2D plot data";
		   #endif
			save_2d(Ux, Uy, Uz, ni, Ph, save_opt, file_output_2d);
			save_1d(Ux, Uy, Uz, ni, Ph, save_opt, file_output_1d);

			if(fourier_flag)
			{
				save_opt.zpos = 0;
				std::cout << "create additional Fourier data" << std::endl;
				OP_FFT my_FTT(Omega);
				field_imag FUx(Omega);
				field_imag FUy(Omega);
				field_imag FUz(Omega);
				field_imag Fni(Omega);
				field_imag FPh(Omega);

				FFT(Ux,FUx);
				FFT(Uy,FUy);
				FFT(Uz,FUz);
				FFT(ni,Fni);
				FFT(Ph,FPh);

				save_F2d(FUx, FUy, FUz, Fni, FPh, save_opt, file_output_F2d);
				save_F1d(FUx, FUy, FUz, Fni, FPh, save_opt, file_output_F1d);

			}
		} // End processing 3d data

    } // END loop over file list



   #if defined(ADVECTION) || defined(E_EXT) || defined(E_INT) || defined(DISSIPATION)
	field_imag FUx(*Ux.my_grid);
	field_imag FUy(*Ux.my_grid);
	field_imag FUz(*Ux.my_grid);
	field_imag Fni(*Ux.my_grid);
	field_imag FPh(*Ux.my_grid);
	field_imag Buffer_FUx(*Ux.my_grid);
	field_imag Buffer_FUy(*Ux.my_grid);
	field_imag Buffer_FUz(*Ux.my_grid);
	field_imag Buffer_Fni(*Ux.my_grid);
	field_imag Buffer_FPh(*Ux.my_grid);

	field_imag FUx_total(*Ux.my_grid);
	field_imag FUy_total(*Uy.my_grid);
	field_imag FUz_total(*Uz.my_grid);

	field_real out(*Ux.my_grid);
	field_real Ux_total(*Ux.my_grid);
	field_real Uy_total(*Ux.my_grid);
	field_real Uz_total(*Ux.my_grid);

	FFT(Ux,FUx);
	FFT(Uy,FUy);
	FFT(Uz,FUz);
	FFT(ni,Fni);
	FFT(Ph,FPh);

	for(int i=0; i<Ux_total.N; ++i)
	{
		Ux_total.val[i] = 0.;
		Uy_total.val[i] = 0.;
		Uz_total.val[i] = 0.;
	}

	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		FUx_total.val[i][0] = 0.;
		FUx_total.val[i][1] = 0.;
		FUy_total.val[i][0] = 0.;
		FUy_total.val[i][1] = 0.;
		FUz_total.val[i][0] = 0.;
		FUz_total.val[i][1] = 0.;
	}
   #endif


   #if defined(ADVECTION) // ##########################
	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		Buffer_FUx.val[i][0] = 0.;
		Buffer_FUx.val[i][1] = 0.;
		Buffer_FUy.val[i][0] = 0.;
		Buffer_FUy.val[i][1] = 0.;
		Buffer_FUz.val[i][0] = 0.;
		Buffer_FUz.val[i][1] = 0.;
	}
	static effect_advection advection;
	advection.execute(FUx, FUy, FUz, FUx, Buffer_FUx);
	advection.execute(FUx, FUy, FUz, FUy, Buffer_FUy);
	advection.execute(FUx, FUy, FUz, FUz, Buffer_FUz);

	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		FUx_total.val[i][0] += Buffer_FUx.val[i][0];
		FUx_total.val[i][1] += Buffer_FUx.val[i][1];
		FUy_total.val[i][0] += Buffer_FUy.val[i][0];
		FUy_total.val[i][1] += Buffer_FUy.val[i][1];
		FUz_total.val[i][0] += Buffer_FUz.val[i][0];
		FUz_total.val[i][1] += Buffer_FUz.val[i][1];
	}

	iFFT(Buffer_FUx,out);
	save_2d(out, save_opt, "./data/splot_ADVECTION_Ux.dat");
	save_1d(out, save_opt, "./data/plot_ADVECTION_Ux.dat");
	for(int i=0; i<Ux_total.N; ++i)
		Ux_total.val[i] += out.val[i];

	iFFT(Buffer_FUy,out);
	save_2d(out, save_opt, "./data/splot_ADVECTION_Uy.dat");
	save_1d(out, save_opt, "./data/plot_ADVECTION_Uy.dat");
	for(int i=0; i<Ux_total.N; ++i)
		Uy_total.val[i] += out.val[i];

	iFFT(Buffer_FUz,out);
	save_2d(out, save_opt, "./data/splot_ADVECTION_Uz.dat");
	save_1d(out, save_opt, "./data/plot_ADVECTION_Uz.dat");
	for(int i=0; i<Ux_total.N; ++i)
		Uz_total.val[i] += out.val[i];
   #endif // END ADVECTION


   #if defined(E_EXT) // ##########################
	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		Buffer_FUx.val[i][0] = 0.;
		Buffer_FUx.val[i][1] = 0.;
		Buffer_FUy.val[i][0] = 0.;
		Buffer_FUy.val[i][1] = 0.;
		Buffer_FUz.val[i][0] = 0.;
		Buffer_FUz.val[i][1] = 0.;
	}
	static effect_translation translation_Ex(Params.M, e_x);
	// translation_Ex is also applied to ion density
	translation_Ex.execute(FUx, Buffer_FUx);
	translation_Ex.execute(FUy, Buffer_FUy);
	translation_Ex.execute(FUz, Buffer_FUz);

	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		FUx_total.val[i][0] += Buffer_FUx.val[i][0];
		FUx_total.val[i][1] += Buffer_FUx.val[i][1];
		FUy_total.val[i][0] += Buffer_FUy.val[i][0];
		FUy_total.val[i][1] += Buffer_FUy.val[i][1];
		FUz_total.val[i][0] += Buffer_FUz.val[i][0];
		FUz_total.val[i][1] += Buffer_FUz.val[i][1];
	}

	iFFT(Buffer_FUx,out);
	save_2d(out, save_opt, "./data/splot_E_EXT_Ux.dat");
	save_1d(out, save_opt, "./data/plot_E_EXT_Ux.dat");
	for(int i=0; i<Ux_total.N; ++i)
		Ux_total.val[i] += out.val[i];

	iFFT(Buffer_FUy,out);
	save_2d(out, save_opt, "./data/splot_E_EXT_Uy.dat");
	save_1d(out, save_opt, "./data/plot_E_EXT_Uy.dat");
	for(int i=0; i<Ux_total.N; ++i)
		Uy_total.val[i] += out.val[i];

	iFFT(Buffer_FUz,out);
	save_2d(out, save_opt, "./data/splot_E_EXT_Uz.dat");
	save_1d(out, save_opt, "./data/plot_E_EXT_Uz.dat");
	for(int i=0; i<Ux_total.N; ++i)
		Uz_total.val[i] += out.val[i];
   #endif


   #if defined(DIFFUSION) // ##########################
	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		Buffer_FUx.val[i][0] = 0.;
		Buffer_FUx.val[i][1] = 0.;
		Buffer_FUy.val[i][0] = 0.;
		Buffer_FUy.val[i][1] = 0.;
		Buffer_FUz.val[i][0] = 0.;
		Buffer_FUz.val[i][1] = 0.;
	}
	static effect_force_E DivP;

	for(int i=0; i<ni.N; ++i)
	{
		out.val[i] = log(ni.val[i])/Params.theta;
		if(fabs(ni.val[i]) < 1.e-10)
			out.val[i] = 0.;
	}
	field_imag FBuffer(*ni.my_grid);
	FFT(out,Buffer_Fni);
	DivP.execute(Buffer_Fni, Buffer_FUx, Buffer_FUy, Buffer_FUz);

	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		FUx_total.val[i][0] += Buffer_FUx.val[i][0];
		FUx_total.val[i][1] += Buffer_FUx.val[i][1];
		FUy_total.val[i][0] += Buffer_FUy.val[i][0];
		FUy_total.val[i][1] += Buffer_FUy.val[i][1];
		FUz_total.val[i][0] += Buffer_FUz.val[i][0];
		FUz_total.val[i][1] += Buffer_FUz.val[i][1];
	}

	iFFT(Buffer_FUx,out);
	save_2d(out, save_opt, "./data/splot_DIFFUSION_Ux.dat");
	save_1d(out, save_opt, "./data/plot_DIFFUSION_Ux.dat");
	//for(int i=0; i<Ux_total.N; ++i)
	//	Ux_total.val[i] += out.val[i];

	iFFT(Buffer_FUy,out);
	save_2d(out, save_opt, "./data/splot_DIFFUSION_Uy.dat");
	save_1d(out, save_opt, "./data/plot_DIFFUSION_Uy.dat");
	for(int i=0; i<Uy_total.N; ++i)
		Uy_total.val[i] += out.val[i];

	iFFT(Buffer_FUz,out);
	save_2d(out, save_opt, "./data/splot_DIFFUSION_Uz.dat");
	save_1d(out, save_opt, "./data/plot_DIFFUSION_Uz.dat");
	for(int i=0; i<Uz_total.N; ++i)
		Uz_total.val[i] += out.val[i];
   #endif // END DIFFUSION


   #if defined(E_INT) // ##########################
	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		Buffer_FUx.val[i][0] = 0.;
		Buffer_FUx.val[i][1] = 0.;
		Buffer_FUy.val[i][0] = 0.;
		Buffer_FUy.val[i][1] = 0.;
		Buffer_FUz.val[i][0] = 0.;
		Buffer_FUz.val[i][1] = 0.;
	}
	static effect_force_E DivE;
	for(int i=0; i<ni.N; ++i)
	{
		out.val[i] += Ph.val[i];
	}
	FFT(out,Buffer_Fni);
	DivE.execute(Buffer_Fni, Buffer_FUx, Buffer_FUy, Buffer_FUz);


	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		FUx_total.val[i][0] += Buffer_FUx.val[i][0];
		FUx_total.val[i][1] += Buffer_FUx.val[i][1];
		FUy_total.val[i][0] += Buffer_FUy.val[i][0];
		FUy_total.val[i][1] += Buffer_FUy.val[i][1];
		FUz_total.val[i][0] += Buffer_FUz.val[i][0];
		FUz_total.val[i][1] += Buffer_FUz.val[i][1];
	}

	iFFT(Buffer_FUx,out);
	save_2d(out, save_opt, "./data/splot_E_INT_Ux.dat");
	save_1d(out, save_opt, "./data/plot_E_INT_Ux.dat");
	//for(int i=0; i<Ux_total.N; ++i)
	//	Ux_total.val[i] += out.val[i];

	iFFT(Buffer_FUy,out);
	save_2d(out, save_opt, "./data/splot_E_INT_Uy.dat");
	save_1d(out, save_opt, "./data/plot_E_INT_Uy.dat");
	for(int i=0; i<Uy_total.N; ++i)
		Uy_total.val[i] += out.val[i];

	iFFT(Buffer_FUz,out);
	save_2d(out, save_opt, "./data/splot_E_INT_Uz.dat");
	save_1d(out, save_opt, "./data/plot_E_INT_Uz.dat");
	for(int i=0; i<Uz_total.N; ++i)
		Uz_total.val[i] += out.val[i];
   #endif // END E_INT



#if defined(COMBINED) // ##########################
	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		Buffer_FUx.val[i][0] = 0.;
		Buffer_FUx.val[i][1] = 0.;
		Buffer_FUy.val[i][0] = 0.;
		Buffer_FUy.val[i][1] = 0.;
		Buffer_FUz.val[i][0] = 0.;
		Buffer_FUz.val[i][1] = 0.;
	}
	static effect_force_E DivPP;

	for(int i=0; i<ni.N; ++i)
	{
		out.val[i] = log(ni.val[i])/Params.theta + Ph.val[i];;
		if(fabs(ni.val[i]) < 1.e-10)
			out.val[i] = 0.;
	}

	save_2d(out, save_opt, "./data/splot_COMBINED_Div.dat");
	save_1d(out, save_opt, "./data/plot_COMBINED_Div.dat");

	field_imag FFBuffer(*ni.my_grid);
	FFT(out,Buffer_Fni);
	DivPP.execute(Buffer_Fni, Buffer_FUx, Buffer_FUy, Buffer_FUz);

	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		FUx_total.val[i][0] += Buffer_FUx.val[i][0];
		FUx_total.val[i][1] += Buffer_FUx.val[i][1];
		FUy_total.val[i][0] += Buffer_FUy.val[i][0];
		FUy_total.val[i][1] += Buffer_FUy.val[i][1];
		FUz_total.val[i][0] += Buffer_FUz.val[i][0];
		FUz_total.val[i][1] += Buffer_FUz.val[i][1];
	}

	iFFT(Buffer_FUx,out);
	save_2d(out, save_opt, "./data/splot_COMBINED_Ux.dat");
	save_1d(out, save_opt, "./data/plot_COMBINED_Ux.dat");
	for(int i=0; i<Ux_total.N; ++i)
		Ux_total.val[i] += out.val[i];

	iFFT(Buffer_FUy,out);
	save_2d(out, save_opt, "./data/splot_COMBINED_Uy.dat");
	save_1d(out, save_opt, "./data/plot_COMBINED_Uy.dat");
	for(int i=0; i<Uy_total.N; ++i)
		Uy_total.val[i] += out.val[i];

	iFFT(Buffer_FUz,out);
	save_2d(out, save_opt, "./data/splot_COMBINED_Uz.dat");
	save_1d(out, save_opt, "./data/plot_COMBINED_Uz.dat");
	for(int i=0; i<Uz_total.N; ++i)
		Uz_total.val[i] += out.val[i];
   #endif // END COMBINED


   #if defined(DISSIPATION) // ##########################
	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		Buffer_FUx.val[i][0] = 0.;
		Buffer_FUx.val[i][1] = 0.;
		Buffer_FUy.val[i][0] = 0.;
		Buffer_FUy.val[i][1] = 0.;
		Buffer_FUz.val[i][0] = 0.;
		Buffer_FUz.val[i][1] = 0.;
	}
	static effect_force_linear dissipation(-Params.tau);
	dissipation.execute(FUx, Buffer_FUx);
	dissipation.execute(FUy, Buffer_FUy);
	dissipation.execute(FUz, Buffer_FUz);

	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		FUx_total.val[i][0] += Buffer_FUx.val[i][0];
		FUx_total.val[i][1] += Buffer_FUx.val[i][1];
		FUy_total.val[i][0] += Buffer_FUy.val[i][0];
		FUy_total.val[i][1] += Buffer_FUy.val[i][1];
		FUz_total.val[i][0] += Buffer_FUz.val[i][0];
		FUz_total.val[i][1] += Buffer_FUz.val[i][1];
	}

	iFFT(Buffer_FUx,out);
	save_2d(out, save_opt, "./data/splot_DISSIPATION_Ux.dat");
	save_1d(out, save_opt, "./data/plot_DISSIPATION_Ux.dat");
	for(int i=0; i<Ux_total.N; ++i)
		Ux_total.val[i] += out.val[i];

	iFFT(Buffer_FUy,out);
	save_2d(out, save_opt, "./data/splot_DISSIPATION_Uy.dat");
	save_1d(out, save_opt, "./data/plot_DISSIPATION_Uy.dat");
	for(int i=0; i<Uy_total.N; ++i)
		Uy_total.val[i] += out.val[i];

	iFFT(Buffer_FUz,out);
	save_2d(out, save_opt, "./data/splot_DISSIPATION_Uz.dat");
	save_1d(out, save_opt, "./data/plot_DISSIPATION_Uz.dat");
	for(int i=0; i<Uz_total.N; ++i)
		Uz_total.val[i] += out.val[i];
   #endif // END DISSIPATION


   #if defined(SPECTRAL_VISCOSITY)
	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		Buffer_FUx.val[i][0] = 0.;
		Buffer_FUx.val[i][1] = 0.;
		Buffer_FUy.val[i][0] = 0.;
		Buffer_FUy.val[i][1] = 0.;
		Buffer_FUz.val[i][0] = 0.;
		Buffer_FUz.val[i][1] = 0.;
	}
	static effect_spectral_viscosity VS(0.66,*FUx.my_grid);
	VS.execute(FUx,Buffer_FUx);
	VS.execute(FUy,Buffer_FUy);
	VS.execute(FUz,Buffer_FUz);

	iFFT(Buffer_FUx,out);
	save_2d(out, save_opt, "./data/splot_SPECTRAL_VISCOSITY_Ux.dat");
	save_1d(out, save_opt, "./data/plot_SPECTRAL_VISCOSITY_Ux.dat");
	for(int i=0; i<Ux_total.N; ++i)
		Ux_total.val[i] += out.val[i];

	iFFT(Buffer_FUy,out);
	save_2d(out, save_opt, "./data/splot_SPECTRAL_VISCOSITY_Uy.dat");
	save_1d(out, save_opt, "./data/plot_SPECTRAL_VISCOSITY_Uy.dat");
	for(int i=0; i<Uy_total.N; ++i)
		Uy_total.val[i] += out.val[i];

	iFFT(Buffer_FUz,out);
	save_2d(out, save_opt, "./data/splot_SPECTRAL_VISCOSITY_Uz.dat");
	save_1d(out, save_opt, "./data/plot_SPECTRAL_VISCOSITY_Uz.dat");
	for(int i=0; i<Uz_total.N; ++i)
		Uz_total.val[i] += out.val[i];
   #endif

#if defined(SAVE_TOTAL)
	save_2d(Ux_total, save_opt, "./data/splot_TOTAL_Ux.dat");
	save_1d(Ux_total, save_opt, "./data/plot_TOTAL_Ux.dat");

	save_2d(Uy_total, save_opt, "./data/splot_TOTAL_Uy.dat");
	save_1d(Uy_total, save_opt, "./data/plot_TOTAL_Uy.dat");

	save_2d(Uz_total, save_opt, "./data/splot_TOTAL_Uz.dat");
	save_1d(Uz_total, save_opt, "./data/plot_TOTAL_Uz.dat");
#endif

	/*
	field_integrate IntO(*Ux.my_grid);
	double V_ges = 0;
	for(int i=0; i<Ux.N; ++i)
	{
		V_ges += IntO.V[i];
	}
		std::cout << IntO.execute(Ux);
	std::cout << "V_ges" << V_ges << "\n";
	*/

   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
	std::cout << "done" << std::endl;

	return 0;
}

