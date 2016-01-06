
#include "plot_2d.hpp"


int plot_2d::default_plane = 2;
int plot_2d::default_direction = 0;
int plot_2d::default_axis = 0;
int plot_2d::default_xpos = 0;
int plot_2d::default_ypos = 0;
int plot_2d::default_zpos = 0;

//int plot::dim2_plane_int = 2;
//std::string plot::dim2_plane_str = "($1):($2)";



plot_2d::plot_2d(const field_real &Fdata, const parameters &Pdata) :
my_parameters(Pdata), my_field(Fdata)
{
   #ifdef _MY_VERBOSE
	logger log("plot");
	log << "plot(const field_real &Fdata, const parameters &Pdata)";
   #endif

	this->path = "./data/";
	this->file = "splot_dummy.pl4";
	this->title = "";
	this->cb_min = 0;
	this->cb_max = 0;

	try
	{
		switch(plot_2d::default_plane)
		{
		case  0:
			this->dim2_plane_int = 0;
			this->gnuplot_plane = "($2):($3)";
			this->gnuplot_v_vec = "(($5)*scale):(($6)*scale)";
			this->x_label = "y in $\\lambda_{De}$";
			this->y_label = "z in $\\lambda_{De}$";
		break;
		case  1:
			this->dim2_plane_int = 1;
			this->gnuplot_plane = "($1):($3)";
			this->gnuplot_v_vec = "(($4)*scale):(($6)*scale)";
			this->x_label = "x in $\\lambda_{De}$";
			this->y_label = "z in $\\lambda_{De}$";
		break;
		case  2:
			this->dim2_plane_int = 2;
			this->gnuplot_plane = "($1):($2)";
			this->gnuplot_v_vec = "(($4)*scale):(($5)*scale)";
			this->x_label = "x in $\\lambda_{De}$";
			this->y_label = "y in $\\lambda_{De}$";
		break;
		default:
			throw -1;
		break;
		}
	} // END Try
	catch (int e)
	{
		std::cout << "An exception occurred: Tow-dimensional plane for output undefined!\n";
	}

}


plot_2d::~plot_2d()
{
   #ifdef _MY_VERBOSE
	logger log("plot");
	log << "~plot()";
   #endif
}

void plot_2d::save_(std::string script)
{
	std::string file_name = path+file;
	std::ofstream output_stream(file_name.c_str(), std::ofstream::trunc);
	output_stream << script;
	output_stream.close();
	return;
}







void save_density_pl4(const grid_Co &Omega, const parameters &Params, std::string file_name)
{

	std::stringstream s_stream;
	double label_pos_x, label_pos_y;

	s_stream << "# File for ploting ion-density via plot4" << std::endl;
	s_stream << "#PLOT -l" << std::endl;
	s_stream << "#PLOT -s 6" << std::endl;
	s_stream << "#PLOT -w " << 16*(Omega.x_axis->L/Omega.y_axis->L) <<" -h 16" << std::endl;
	s_stream << std::endl;
	s_stream << "########## quick Options #######################################################" << std::endl;
	s_stream << "cb_scale = 1." << std::endl;
	s_stream << std::endl;

	s_stream << "########## SPLOT Options #######################################################" << std::endl;
	s_stream << "set palette model RGB defined ( 0\"black\", 1\"blue\", 2\"cyan\",3\"green\",4\"yellow\", 5\"red\",6\"black\" )" << std::endl;
	s_stream << "unset surface" << std::endl;
	s_stream << "set pm3d map" << std::endl;
	s_stream << "set pm3d implicit at b" << std::endl;
	s_stream << "#set contour base" << std::endl;
	s_stream << "#set cntrparam bspline" << std::endl;
	s_stream << "#set cntrparam levels auto 3";
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
	s_stream << "set xrange ["<< Omega.x_axis->val_at(0) << ":" << Omega.x_axis->val_at(Omega.Nx -1) << "]" << std::endl;
	s_stream << "set yrange ["<< Omega.y_axis->val_at(0) << ":" << Omega.y_axis->val_at(Omega.Nx -1) << "]" << std::endl;
	s_stream << "set cbrange [cb_scale*0.00845726:cb_scale*-0.00845726]" << std::endl;
	s_stream << std::endl;

	s_stream << "########## Labels #############################################################" << std::endl;
	s_stream << "set xlabel  'x'" << std::endl;
	s_stream << "set ylabel  'y' offset -2,0" << std::endl;
	s_stream << "set cblabel 'cb'" << std::endl;
	s_stream << "set title 'density'" << std::endl;
	s_stream << std::endl;

	s_stream << "########## Legend #############################################################" << std::endl;
	label_pos_x = Omega.x_axis->val_at(0) + Omega.x_axis->L - 0.10*Omega.x_axis->L;
	label_pos_y = Omega.y_axis->val_at(0) + Omega.y_axis->L - 0.10*Omega.y_axis->L;
	s_stream << "set key at " << label_pos_x << "," << label_pos_y   << std::endl;
	s_stream << "set key title '";
	s_stream << "$M="       << Params.M     << "$ \\\t";
	s_stream << "$\\theta=" << Params.theta << "$ \\\t";
	s_stream << "$\\tau="   << Params.tau   << "$";
	s_stream << "'" << std::endl;
	s_stream << std::endl;

	s_stream << "########## Labels #############################################################" << std::endl;
	// label elektrisches Feld
	label_pos_x = Omega.x_axis->val_at(0) + 0.1*Omega.x_axis->L;
	label_pos_y = Omega.y_axis->val_at(0) + Omega.y_axis->L - 0.10*Omega.y_axis->L;
	s_stream << "set label '[b]{$\\stackrel{\\vec{E}_0}{\\longrightarrow}$ }' at ";
	s_stream << label_pos_x <<","<< label_pos_y <<" front" << std::endl;
	// label Schnitt-Ebene
	label_pos_x = Omega.x_axis->val_at(0) + 0.05*Omega.x_axis->L;
	label_pos_y = Omega.y_axis->val_at(0) + 0.05*Omega.y_axis->L;
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
	s_stream << "splot \"" + ExtractFilename(file_name) + "\" using ($1):($2):($6) notitle" << std::endl;

	std::string pl4_content = s_stream.str();

	file_name = ChangeFileExtension(file_name,".pl4");
	FindAndReplace(file_name,".pl4","_d.pl4");

	std::ofstream output_stream(file_name.c_str(), std::ofstream::trunc);
	output_stream << pl4_content;
	output_stream.close();

	return;

}




void save_velocity_pl4(const grid_Co &Omega, const parameters &Params, std::string file_name)
{


	std::stringstream s_stream;
	double label_pos_x, label_pos_y;

	s_stream << "# File for ploting a velocity field via plot4" << std::endl;
	s_stream << "#PLOT -l" << std::endl;
	s_stream << "#PLOT -s 6" << std::endl;
	s_stream << "#PLOT -w " << 16*(Omega.x_axis->L/Omega.y_axis->L) <<" -h 16" << std::endl;
	s_stream << std::endl;
	s_stream << "########## quick Options #######################################################" << std::endl;
	s_stream << "vec_scale = 1." << std::endl;
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
	s_stream << "#set format cb \"%.1e\"" << std::endl;
	s_stream << "# set logscale x" << std::endl;
	s_stream << "# set logscale y" << std::endl;
	s_stream << "# set logscale z" << std::endl;
	s_stream << std::endl;

	s_stream << "########## Ranges ##############################################################" << std::endl;
	s_stream << "set xrange ["<< Omega.x_axis->val_at(0) << ":" << Omega.x_axis->val_at(Omega.Nx -1) << "]" << std::endl;
	s_stream << "set yrange ["<< Omega.y_axis->val_at(0) << ":" << Omega.y_axis->val_at(Omega.Nx -1) << "]" << std::endl;
	s_stream << std::endl;

	s_stream << "########## Labels #############################################################" << std::endl;
	s_stream << "set xlabel  'x'" << std::endl;
	s_stream << "set ylabel  'y' offset -2,0" << std::endl;
	s_stream << "set cblabel 'cb'" << std::endl;
	s_stream << "set title 'density'" << std::endl;
	s_stream << std::endl;
	s_stream << "########## Legend #############################################################" << std::endl;
	label_pos_x = Omega.x_axis->val_at(0) + Omega.x_axis->L - 0.10*Omega.x_axis->L;
	label_pos_y = Omega.y_axis->val_at(0) + Omega.y_axis->L - 0.10*Omega.y_axis->L;
	s_stream << "set key at " << label_pos_x << "," << label_pos_y   << std::endl;
	s_stream << "set key title '";
	s_stream << "$M="       << Params.M     << "$ \\\t";
	s_stream << "$\\theta=" << Params.theta << "$ \\\t";
	s_stream << "$\\tau="   << Params.tau   << "$";
	s_stream << "'" << std::endl;
	s_stream << std::endl;

	s_stream << "########## Labels #############################################################" << std::endl;
	// label elektrisches Feld
	label_pos_x = Omega.x_axis->val_at(0) + 0.1*Omega.x_axis->L;
	label_pos_y = Omega.y_axis->val_at(0) + Omega.y_axis->L - 0.10*Omega.y_axis->L;
	s_stream << "set label '[b]{$\\stackrel{\\vec{E}_0}{\\longrightarrow}$ }' at ";
	s_stream << label_pos_x <<","<< label_pos_y <<" front" << std::endl;
	// label Schnitt-Ebene
	label_pos_x = Omega.x_axis->val_at(0) + 0.05*Omega.x_axis->L;
	label_pos_y = Omega.y_axis->val_at(0) + 0.05*Omega.y_axis->L;
	s_stream << "set label '[b]{Schnitt in z-Ebene @ $z=0$ }' at ";
	s_stream << label_pos_y <<","<< label_pos_x <<" front" << std::endl;
	s_stream << std::endl;

	s_stream << "########## Style Definitions ###################################################" << std::endl;
	s_stream << "set style line 1 lt 1 lw 2 lc rgb \"#aaaaaa\"" << std::endl;
	s_stream << "set style line 2 lt 1 lw 2 lc rgb \"#aaaaaa\"" << std::endl;
	s_stream << "set style line 3 lt 1 lw 2 lc rgb \"#aaaaaa\"" << std::endl;
	s_stream << "set style line 4 lt 1 lw 2 lc rgb \"#aaaaaa\"" << std::endl;
	s_stream << "set style increment userstyles" << std::endl;
	s_stream << "set style arrow 1 head filled size screen 0.03,5 ls 1" << std::endl;
	s_stream << std::endl;

	s_stream << "########## Plot-Command ########################################################" << std::endl;
	s_stream << "plot \"" + ExtractFilename(file_name) + "\" using ($1):($2):(($3)*vec_scale):(($4)*vec_scale) every 4:4 with vectors arrowstyle 1 notitle" << std::endl;


	std::string pl4_content = s_stream.str();

	file_name = ChangeFileExtension(file_name,".pl4");
	FindAndReplace(file_name,".pl4","_v.pl4");

	std::ofstream output_stream(file_name.c_str(), std::ofstream::trunc);
	output_stream << pl4_content;
	output_stream.close();



	return;

}





