#include "plot_2d_velocity.hpp"



velocity_plot::velocity_plot(const field_real &xFdata, const field_real &yFdata, const field_real &zFdata, const parameters &Pdata) :
plot_2d(xFdata, Pdata)
{
	this->file = "splot_velocity.pl4";
}



std::string velocity_plot::script_velocity(void)
{
	std::stringstream s_stream;

	s_stream << "# File for ploting Ion-Flux via plot4" << std::endl;

	s_stream << "#PLOT -l" << std::endl;
	s_stream << "#PLOT -s 6" << std::endl;
	s_stream << "#PLOT -w 15 -h 10" << std::endl;

	s_stream << "########## SPLOT Options #######################################################" << std::endl;
	s_stream << "scale = 1." << std::endl;


	s_stream << "#set palette model RGB defined ( 0\"black\", 1\"blue\", 2\"cyan\",3\"green\",4\"yellow\", 5\"red\",8\"purple\", 9\"white\" ) " << std::endl;
	s_stream << "#unset surface" << std::endl;
	s_stream << "#set pm3d map" << std::endl;
	s_stream << "#set pm3d implicit at b" << std::endl;

	s_stream << "#set contour base" << std::endl;
	s_stream << "#set cntrparam bspline" << std::endl;
	s_stream << "#set cntrparam levels discrete -1,0,1" << std::endl;

	s_stream << "########## Canvas ############################################################## " << std::endl;
	s_stream << "set lmargin 12" << std::endl;
	s_stream << "set tmargin 4" << std::endl;
	s_stream << "set rmargin 2" << std::endl;
	s_stream << "set bmargin 4" << std::endl;

	s_stream << "########## Format and Scaling ##################################################" << std::endl;
	s_stream << "# set format x ''" << std::endl;
	s_stream << "# set format y ''" << std::endl;
	s_stream << "# set format z ''" << std::endl;
	s_stream << "#set format cb \"%.1e\"" << std::endl;
	s_stream << "# set logscale x" << std::endl;
	s_stream << "# set logscale y" << std::endl;
	s_stream << "# set logscale z" << std::endl;



	s_stream << "########## Ranges ##############################################################" << std::endl;
	s_stream << "#set xrange [-2:12]" << std::endl;
	s_stream << "#set yrange [0:6]" << std::endl;
	s_stream << "#set cbrange [2.5e-4:-5e-5]" << std::endl;

	s_stream << "########## Labels #############################################################" << std::endl;
	s_stream << "set xlabel  '" << this->  x_label << "'" << std::endl;
	s_stream << "set ylabel  '" << this->  y_label << "' offset -2,0" << std::endl;
	s_stream << "set cblabel '" << this-> cb_label << "'" << std::endl;
	s_stream << "set title '"<< this->title << "'" << std::endl;
	s_stream << std::endl;

	s_stream << "########## Style Definitions ###################################################" << std::endl;
	s_stream << "set style arrow 1 head filled size screen 0.03,5 ls 1" << std::endl;
	s_stream << "set style line 1 lt 1 lw 6 lc rgb \"black\"" << std::endl;
	s_stream << "set style line 2 lt 1 lw 1.5 lc rgb \"yellow\"" << std::endl;
	s_stream << "set style line 3 lt 1 lw 1.5 lc rgb \"#006400\"" << std::endl;
	s_stream << "set style line 4 lt 1 lw 1.5 lc rgb \"blue\"" << std::endl;
	s_stream << "set style increment userstyles" << std::endl;

	s_stream << "########## Plot-Command ########################################################" << std::endl;
	s_stream << "plot \"splot_data.dat\" using " + this->gnuplot_plane + ":" + this->gnuplot_v_vec + " with vectors arrowstyle 1 notitle" << std::endl;


	return s_stream.str();
}
