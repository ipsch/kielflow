#include "plot_2d_density.hpp"




density_plot::density_plot(const field_real &Fdata, const parameters &Pdata) :
plot_2d(Fdata, Pdata)
{
	this->file = "splot_density.pl4";
}


void density_plot::auto_set(const field_real &field)
{
	diagnostics field_stats;

	field_stats.analyze(field);
	int index = 0;


	// autoscaling color-bar axis
	cb_max = field_stats.get_supremum();
	cb_min = field_stats.get_infinum();

	std::cout << "cb_max: " << cb_max << std::endl;
	std::cout << "cb_min: " << cb_min << std::endl;
	std::cout << std::endl;

	std::cout << "supremum@: (";
	std::cout << field_stats.get_i_supremum() << ",";
	std::cout << field_stats.get_j_supremum() << ",";
	std::cout << field_stats.get_k_supremum() << ")" << std::endl;

	index = field.my_grid->index_at(field_stats.i_supremum,field_stats.j_supremum,field_stats.k_supremum);
	std::cout << "supremum = " << field.val[index] << std::endl;
	std::cout << std::endl;


	std::cout << "infinum@: (";
	std::cout << field_stats.get_i_infinum() << ",";
	std::cout << field_stats.get_j_infinum() << ",";
	std::cout << field_stats.get_k_infinum() << ")" << std::endl;

	index = field.my_grid->index_at(field_stats.i_infinum,field_stats.j_infinum,field_stats.k_infinum);
	std::cout << "infinum = " << field.val[index] << std::endl;
	std::cout << std::endl;


	cb_max = max(  2.*fabs(cb_min),cb_max);
	cb_min = min(-0.5*fabs(cb_max),cb_min);

	return;
}


std::string density_plot::script_density(void)
{
	std::stringstream s_stream;

	s_stream << "# File for ploting Ion-Densities via plot4" << std::endl;
	s_stream << std::endl;

	s_stream << "#PLOT -l" << std::endl;
	s_stream << "#PLOT -s 6" << std::endl;
	s_stream << "#PLOT -w 15 -h 10" << std::endl;
	s_stream << std::endl;

	s_stream << "########## quick Options #######################################################" << std::endl;
	s_stream << "cb_scale = 1." << std::endl;
	s_stream << std::endl;

    s_stream << "########## SPLOT Options #######################################################" << std::endl;
    s_stream << "set palette model RGB defined ( 0\"black\", 1\"blue\", 2\"cyan\",3\"green\",4\"yellow\", 5\"red\",8\"purple\", 9\"white\" )" << std::endl;
    s_stream << "unset surface" << std::endl;
    s_stream << "set pm3d map" << std::endl;
    s_stream << "set pm3d implicit at b" << std::endl;
    s_stream << "# set contour base" << std::endl;
    s_stream << "# set cntrparam bspline" << std::endl;
    s_stream << "# set cntrparam levels discrete -1,0,1" << std::endl;
    s_stream << std::endl;

    s_stream << "########## Canvas ##############################################################" << std::endl;
    s_stream << "set lmargin 2" << std::endl;
    s_stream << "set tmargin 2" << std::endl;
    s_stream << "set rmargin 2" << std::endl;
    s_stream << "set bmargin 2" << std::endl;
    s_stream << std::endl;

    s_stream << "########## Format and Scaling ##################################################" << std::endl;
    s_stream << "# set format x ''" << std::endl;
    s_stream << "# set format y ''" << std::endl;
    s_stream << "# set format z ''" << std::endl;
    s_stream << "# set format cb \"%.1e\"" << std::endl;
    s_stream << "# set logscale x" << std::endl;
    s_stream << "# set logscale y" << std::endl;
    s_stream << "# set logscale z" << std::endl;
    s_stream << std::endl;

    s_stream << "########## Ranges ##############################################################" << std::endl;
    s_stream << "#set xrange [-2:12]" << std::endl;
    s_stream << "#set yrange [0:6]" << std::endl;
    s_stream << "set cbrange [cb_scale*" << cb_min << ":cb_scale*" << cb_max << "]" << std::endl;
    s_stream << std::endl;

	s_stream << "########## Labels #############################################################" << std::endl;
	s_stream << "set xlabel  '" << this->  x_label << "'" << std::endl;
	s_stream << "set ylabel  '" << this->  y_label << "' offset -2,0" << std::endl;
	s_stream << "set cblabel '" << this-> cb_label << "'" << std::endl;
	s_stream << "set title '"<< this->title << "'" << std::endl;
	s_stream << std::endl;

    s_stream << "########## Style Definitions ###################################################" << std::endl;
    s_stream << "set style line 1 lt 1 lw 6 lc rgb \"black\" " << std::endl;
    s_stream << "set style line 2 lt 1 lw 1.5 lc rgb \"yellow\" " << std::endl;
    s_stream << "set style line 3 lt 1 lw 1.5 lc rgb \"#006400\" " << std::endl;
    s_stream << "set style line 4 lt 1 lw 1.5 lc rgb \"blue\" " << std::endl;
    s_stream << "set style increment userstyles" << std::endl;
    s_stream << std::endl;

    s_stream << "########## Plot-Command ########################################################" << std::endl;
    s_stream << "splot \"splot_data.dat\" using " + this->gnuplot_plane + ":7 notitle" << std::endl;
	return s_stream.str();
}



