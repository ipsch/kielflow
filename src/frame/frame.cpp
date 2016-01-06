// gcc standard C++ libraries (common stuff)
#include <iostream>  // standard I/O-Operations (I/O to console)
#include <fstream>   // same as iostream but for files
#include <cmath>     // basic math stuff (sqrt(), cabs(), pow(a,n) etc.)
#include <string>    // class : string (here: used for filenames only; Maybe I should avoid this library)
#include <sstream>   // class : string-stream (MORE OVERHEAD) I'am using this to interpret data read from a file
#include <iomanip>   // sets format of output-stream (scientific & fixed) (is this really necessary??)

#include "parameters.hpp"
#include "domain.hpp"
#include "particle.hpp"
#include "plot.hpp"
#include "IO.hpp"

#include "fftw3.h"
























void save_pl4(std::string file_name)
{

	std::stringstream s_stream;
	double label_pos_x, label_pos_y;

	s_stream << "# File for ploting electric-potential via plot4" << std::endl;
	s_stream << "#PLOT -l" << std::endl;
	s_stream << "#PLOT -s 6" << std::endl;
	s_stream << "#PLOT -w " << 16*(domain::Lx/domain::Ly) <<" -h 16" << std::endl;
	s_stream << std::endl;
	s_stream << "########## quick Options #######################################################" << std::endl;
	s_stream << "cb_scale = 1." << std::endl;
	s_stream << std::endl;

	s_stream << "########## SPLOT Options #######################################################" << std::endl;
	s_stream << "set palette model RGB defined ( 0\"black\", 1\"blue\", 2\"cyan\",3\"green\",4\"yellow\", 5\"red\",6\"black\" )" << std::endl;
	s_stream << "unset surface" << std::endl;
	s_stream << "set pm3d map" << std::endl;
	s_stream << "set pm3d implicit at b" << std::endl;
		s_stream << "set contour base" << std::endl;
		s_stream << "set cntrparam bspline" << std::endl;
		s_stream << "set cntrparam levels auto 3";
		s_stream << std::endl;

		s_stream << "########## Canvas ##############################################################" << std::endl;
		s_stream << "set lmargin 1" << std::endl;
		s_stream << "set tmargin 1" << std::endl;
		s_stream << "set rmargin 3" << std::endl;
		s_stream << "set bmargin 2" << std::endl;
		s_stream << std::endl;

		s_stream << "########## Format and Scaling ##################################################" << std::endl;
		s_stream << "# set format x ''" << std::endl;
		s_stream << "# set format y ''" << std::endl;
		s_stream << "# set format z ''" << std::endl;
		s_stream << "set format cb \"%.1e\"" << std::endl;
		s_stream << "# set logscale x" << std::endl;
		s_stream << "# set logscale y" << std::endl;
		s_stream << "# set logscale z" << std::endl;
		s_stream << std::endl;

		s_stream << "########## Ranges ##############################################################" << std::endl;
		s_stream << "set xrange ["<< domain::x0 << ":" << domain::x0+domain::Lx << "]" << std::endl;
		s_stream << "set yrange ["<< domain::y0 << ":" << domain::y0+domain::Ly << "]" << std::endl;
		s_stream << "set cbrange [cb_scale*0.00845726:cb_scale*-0.00845726]" << std::endl;
		s_stream << std::endl;

		s_stream << "########## Labels #############################################################" << std::endl;
		s_stream << "set xlabel  'x'" << std::endl;
		s_stream << "set ylabel  'y' offset -2,0" << std::endl;
		s_stream << "set cblabel 'cb'" << std::endl;
		s_stream << "set title 'title'" << std::endl;
		s_stream << std::endl;

		s_stream << "########## Legend #############################################################" << std::endl;
		label_pos_x = domain::x0 + domain::Lx - 0.10*domain::Lx;
		label_pos_y = domain::y0 + domain::Ly - 0.10*domain::Ly;
		s_stream << "set key at " << label_pos_x << "," << label_pos_y   << std::endl;
		s_stream << "set key title '";
		s_stream << "$M="       << parameters::M     << "$ \\\t";
		s_stream << "$\\theta=" << parameters::theta << "$ \\\t";
		s_stream << "$\\tau="   << parameters::tau   << "$";
		s_stream << "'" << std::endl;
		s_stream << std::endl;

		s_stream << "########## Labels #############################################################" << std::endl;
		// label elektrisches Feld
		label_pos_x = domain::x0 + 0.1*domain::Lx;
		label_pos_y = domain::y0 + domain::Ly - 0.10*domain::Ly;
		s_stream << "set label '[b]{$\\stackrel{\\vec{E}_0}{\\longrightarrow}$ }' at ";
		s_stream << label_pos_x <<","<< label_pos_y <<" front" << std::endl;
		// label Schnitt-Ebene
		label_pos_x = domain::x0 + 0.05*domain::Lx;
		label_pos_y = domain::y0 + 0.2*domain::Ly;
		s_stream << "set label '[b]{Schnitt in z-Ebene @ $z=0$ }' at ";
		s_stream << label_pos_y <<","<< label_pos_x <<" front" << std::endl;
		s_stream << std::endl;

		s_stream << "########## Style Definitions ###################################################" << std::endl;
		s_stream << "set style line 1 lt 1 lw 2 lc rgb \"#aaaaaa\"" << std::endl;
		s_stream << "set style line 2 lt 1 lw 2 lc rgb \"#aaaaaa\"" << std::endl;
		s_stream << "set style line 3 lt 1 lw 2 lc rgb \"#aaaaaa\"" << std::endl;
		s_stream << "set style line 4 lt 1 lw 2 lc rgb \"#aaaaaa\"" << std::endl;
		s_stream << "set style increment userstyles" << std::endl;
		s_stream << std::endl;

		s_stream << "########## Plot-Command ########################################################" << std::endl;
		s_stream << "splot \"" + file_name + ".dat\" using ($1):($2):($7) notitle" << std::endl;


	std::string pl4_content = s_stream.str();

	std::string output = "./data/"+ file_name + ".pl4";

	std::ofstream output_stream(output.c_str(), std::ofstream::trunc);



	output_stream << pl4_content;
	output_stream.close();

	return;

}



int main(int argc, char *argv[])
{
	std::string current_exec_name = argv[0]; // Name of the current exec program
	std::string first_arge;
	std::vector<std::string> all_args;

	std::string raw_name;
	std::string file_name_h5;
	std::string file_name_dat;

	if (argc > 1)
	{
		first_arge = argv[1];
	    all_args.assign(argv + 1, argv + argc);


	    raw_name = all_args[0];
	    file_name_h5  = "./data/"+all_args[0]+".h5";
	    file_name_dat = "./data/"+all_args[0]+".dat";


	}
	else
	{
		raw_name = "fields";
		file_name_h5  = "./data/fields.h5";
		file_name_dat = "./data/fields.dat";

	}



	load_domain(file_name_h5);

	domain::Nx = 128;
	domain::Ny = 128;
	domain::Nz = 128;
	domain::N_real_ges = domain::Nx*domain::Ny*domain::Nz;
	domain::N_imag_ges = domain::Nx*domain::Ny*(domain::Nz/2 +1);
	domain::Lx = 10;
	domain::Ly = 10;
	domain::Lz = 10;
	domain::x0 = -5;
	domain::y0 = -5;
	domain::z0 = -5;
	domain::dt = 0.01;

	load_parameters(file_name_h5);


	particle::load();




	load_slice(Ux,"Ux",file_name_h5);
	load_slice(Uy,"Uy",file_name_h5);
	load_slice(Uz,"Uz",file_name_h5);
	load_slice(ni,"ni",file_name_h5);
	load_slice(Ph,"Ph",file_name_h5);






	return 0;
}
