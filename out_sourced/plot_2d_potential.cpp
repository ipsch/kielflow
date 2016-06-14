#include "plot_2d_potential.hpp"




potential_plot::potential_plot(const field_real &Fdata, const parameters &Pdata) :
plot_2d(Fdata, Pdata)
{
	this->file = "splot_potential.pl4";
	this->title = "Potential $\\Phi$";
	this->cb_label = "$\\Phi$ in $\\frac{q_i^2}{m_i c_{si}^2}$";
}

void potential_plot::auto_set(const field_real &field)
{
	diagnostics field_stats;

	field_stats.analyze(field);

	// autoscaling color-bar axis
	cb_max = field_stats.get_supremum();
	cb_min = field_stats.get_infinum();


	cb_max = max(fabs(cb_min),cb_max);
	cb_min =  -cb_max;

	return;
}

std::string potential_plot::script_potential(void)
{
	std::stringstream s_stream;
	double label_pos_x, label_pos_y;

	s_stream << "# File for ploting electric-potential via plot4" << std::endl;
	s_stream << std::endl;

	s_stream << "#PLOT -l" << std::endl;
	s_stream << "#PLOT -s 6" << std::endl;
	s_stream << "#PLOT -w " << 16*(my_field.my_grid->x_axis->L/my_field.my_grid->y_axis->L) <<" -h 16" << std::endl;
	s_stream << std::endl;

	s_stream << "########## quick Options #######################################################" << std::endl;
	s_stream << "if (!exists(\"filename\")) filename='splot_data.dat'" << std::endl;
	s_stream << "cb_scale = 1." << std::endl;
	s_stream << std::endl;

	s_stream << "########## SPLOT Options #######################################################" << std::endl;
	s_stream << "set palette model RGB defined ( 0\"black\", 1\"blue\", 2\"cyan\",3\"green\",4\"yellow\", 5\"red\",6\"black\" )" << std::endl;
	s_stream << "unset surface" << std::endl;
	s_stream << "set pm3d map" << std::endl;
	s_stream << "set pm3d implicit at b" << std::endl;
	s_stream << "#set contour base" << std::endl;
	s_stream << "#set cntrparam bspline" << std::endl;
	s_stream << "#set cntrparam levels discrete -1,0,1" << std::endl;
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


	double x_a = my_field.my_grid->x_axis->l0;
	double x_b = my_field.my_grid->x_axis->l0 + my_field.my_grid->x_axis->L;
	s_stream << "########## Ranges ##############################################################" << std::endl;
	s_stream << "set xrange ["<< x_a << ":" << x_b << "]" << std::endl;
	s_stream << "#set yrange [0:6]" << std::endl;
	s_stream << "set cbrange [cb_scale*" << cb_max << ":cb_scale*" << cb_min << "]" << std::endl;
	s_stream << std::endl;

	s_stream << "########## Labels #############################################################" << std::endl;
	s_stream << "set xlabel  '" << this->  x_label << "'" << std::endl;
	s_stream << "set ylabel  '" << this->  y_label << "' offset -2,0" << std::endl;
	s_stream << "set cblabel '" << this-> cb_label << "'" << std::endl;
	s_stream << "set title '"<< this->title << "'" << std::endl;
	s_stream << std::endl;

	s_stream << "########## Legend #############################################################" << std::endl;
	label_pos_x = my_field.my_grid->x_axis->l0 + 0.9*my_field.my_grid->x_axis->L;
	label_pos_y = my_field.my_grid->y_axis->l0 + 0.9*my_field.my_grid->y_axis->L;
	s_stream << "set key at " << label_pos_x << "," << label_pos_y   << std::endl;
	s_stream << "set key title '";
	s_stream << "$M="       << my_parameters.M << "$ \\\t";
	s_stream << "$\\theta=" << my_parameters.theta << "$ \\\t";
	s_stream << "$\\tau="   << my_parameters.tau   << "$";
	s_stream << "'" << std::endl;
	s_stream << std::endl;

	s_stream << "########## Labels #############################################################" << std::endl;
	// label elektrisches Feld
	label_pos_x = my_field.my_grid->x_axis->l0 + 0.1*my_field.my_grid->x_axis->L;
	label_pos_y = my_field.my_grid->y_axis->l0 + 0.9*my_field.my_grid->y_axis->L;
	s_stream << "set label '[b]{$\\stackrel{\\vec{E}_0}{\\longrightarrow}$ }' at ";
	s_stream << label_pos_x <<","<< label_pos_y <<" front" << std::endl;
	// label Schnitt-Ebene
	label_pos_x = my_field.my_grid->x_axis->l0 + 0.05*my_field.my_grid->x_axis->L;
	label_pos_y = my_field.my_grid->y_axis->l0 + 0.2*my_field.my_grid->y_axis->L;
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
	s_stream << "splot filename using " + this->gnuplot_plane + ":($8) notitle" << std::endl;


	return s_stream.str();
}


