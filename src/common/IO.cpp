#include "IO.hpp"



inline bool file_exist (const std::string& path)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "file_exist (const std::string& name)";
	std::string status_text = "file: " + path;
	my_log << status_text;
   #endif

    std::ifstream file_to_check(path.c_str());
    if (file_to_check.good())
    {
    	file_to_check.close();
        return true;
    }
    else
    {
    	file_to_check.close();
        return false;
    }
}

void file_create(const std::string &path)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "file_create(const std::string &path)";
	std::string status_text = "file: " + path;
	my_log << status_text;
   #endif

	echelon::file hdf5_file(path, echelon::file::create_mode::truncate);

	echelon::group gr_parameters = hdf5_file.create_group("parameters");
	echelon::group gr_domain     = hdf5_file.create_group("domain");
	echelon::group gr_data       = hdf5_file.create_group("data");

	return;
}



void save_axis(const axis * const A, const std::string &C, const std::string &path)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "save_axis(const axis * const A, const std::string &C, const std::string &path)";
	std::string status_text = "path: " + path;
	my_log << status_text;
   #endif

    if (!file_exist(path))
    	file_create(path);

	echelon::file hdf5_file(path, echelon::file::open_mode::read_write);
	echelon::group gr_domain = hdf5_file["domain"];

	echelon::scalar_dataset ds_type = gr_domain.create_scalar_dataset("type_"+C, A->type_id);
	echelon::scalar_dataset ds_m    = gr_domain.create_scalar_dataset("m_"+C   , A->m);
	echelon::scalar_dataset ds_l0   = gr_domain.create_scalar_dataset("l0_"+C  , A->l0);
	echelon::scalar_dataset ds_N    = gr_domain.create_scalar_dataset("N_"+C   , A->N);
	echelon::scalar_dataset ds_L    = gr_domain.create_scalar_dataset("L_"+C   , A->L);


   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return;
}

void load_axis(const std::string &C, const std::string &path, axis_Co * &A)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "load_axis(...)";
	std::string status_text = "path: " + path;
	my_log << status_text;
    #endif

	echelon::file hdf5_file(path, echelon::file::open_mode::read_only);
	echelon::group gr_domain = hdf5_file["domain"];

	int type, N;
	double l0, L, m;

    echelon::scalar_dataset ds_type = gr_domain["type_"+C];
    type <<= ds_type;
    echelon::scalar_dataset ds_m   = gr_domain["m_"+C];
    m <<= ds_m;
    echelon::scalar_dataset ds_l0   = gr_domain["l0_"+C];
    l0 <<= ds_l0;
    echelon::scalar_dataset ds_N    = gr_domain["N_"+C];
    N <<= ds_N;
    echelon::scalar_dataset ds_L    = gr_domain["L_"+C];
    L <<= ds_L;

    switch(type)
    {
    case 1:
    	A = new axis_CoEqSt(l0,L,N);
    	break;
    case 2:
    	A = new axis_CoHySt(l0,L,N,m);
    	break;
    case 3:
    	A = new axis_CoSiSt(l0,L,N,double(m));
    	break;
    }

   #if defined(_MY_VERBOSE_TEDIOUS)
    my_log << "done";
   #endif

	return;
}



void save_grid(const grid_Co &Target, const std::string &path)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "save_grid(const std::string &path, const grid &Target)";
	std::string status_text = "path: " + path;
	my_log << status_text;
   #endif

	save_axis(Target.x_axis,"x",path);
	save_axis(Target.y_axis,"y",path);
	save_axis(Target.z_axis,"z",path);

	// Lx / Ly / Lz #########################################################
	//echelon::scalar_dataset ds_dt = gr_domain.create_scalar_dataset<double>("dt", Target.dt);
	//echelon::scalar_dataset ds_t  = gr_domain.create_scalar_dataset<double>("t" , Target.t);

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return;
}

grid_Co load_grid(const std::string &path)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "load_grid(const std::string &path)";
	std::string status_text = "path: " + path;
	my_log << status_text;
   #endif


	axis_Co * x_axis;
	axis_Co * y_axis;
	axis_Co * z_axis;

	load_axis("x", path, x_axis);
	load_axis("y", path, y_axis);
	load_axis("z", path, z_axis);

    grid_Co Omega(*x_axis, *y_axis , *z_axis);

   #if defined(_MY_VERBOSE_TEDIOUS)
	std::stringstream sstr;
	my_log << "domain(Nx,Ny,Nz,x0,y0,z0,Lx,Ly,Lz,dt):";
	sstr << "(" << Omega.x_axis->N;
	sstr << "," << Omega.y_axis->N;
	sstr << "," << Omega.z_axis->N;
	sstr << "," << Omega.x_axis->l0;
	sstr << "," << Omega.y_axis->l0;
	sstr << "," << Omega.z_axis->l0;
	sstr << "," << Omega.x_axis->L;
	sstr << "," << Omega.y_axis->L;
	sstr << "," << Omega.z_axis->L;
	//sstr << "," << Omega.dt;
	sstr << ")";
	my_log << sstr.str();
	my_log << "done";
   #endif

	delete x_axis;
	delete y_axis;
	delete z_axis;

	return Omega;
}



void save_time(const double &Time, const std::string &path)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "save_time(const double &TIME, const std::string &path)";
	std::string status_text = "path: " + path;
	my_log << status_text;
   #endif

	echelon::file hdf5_file(path, echelon::file::open_mode::read_write);
	echelon::group gr_parameters = hdf5_file["parameters"];

	echelon::scalar_dataset ds_T     = hdf5_file.create_scalar_dataset("time", Time);

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return;
}

double load_time(const std::string &path)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO.CPP: load_parameters");
	std::string status_text = "start - path: " + path;
	my_log << status_text;
   #endif

	double Time;

	echelon::file hdf5_file(path, echelon::file::open_mode::read_only);
	echelon::scalar_dataset ds_T = hdf5_file["time"];
	Time  <<= ds_T;


   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return Time;
}



void save_parameters(const parameters &Target, const std::string &path)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "save_parameters(const std::string &path, const parameters &Target)";
	std::string status_text = "path: " + path;
	my_log << status_text;
   #endif

	echelon::file hdf5_file(path, echelon::file::open_mode::read_write);
	echelon::group gr_parameters = hdf5_file["parameters"];

	echelon::scalar_dataset ds_M     = gr_parameters.create_scalar_dataset("mach number M",           Target.M);
	echelon::scalar_dataset ds_tau   = gr_parameters.create_scalar_dataset("damping ratio tau",       Target.tau);
	echelon::scalar_dataset ds_theta = gr_parameters.create_scalar_dataset("temperature ratio theta", Target.theta);
	echelon::scalar_dataset ds_mu    = gr_parameters.create_scalar_dataset("velocity diffusion mu",   Target.mu);
	echelon::scalar_dataset ds_beta  = gr_parameters.create_scalar_dataset("magnetization beta",      Target.beta);

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return;
}

parameters load_parameters(const std::string &path)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO.CPP: load_parameters");
	std::string status_text = "start - path: " + path;
	my_log << status_text;
   #endif

	double M, tau, theta, mu, beta;

	echelon::file hdf5_file(path, echelon::file::open_mode::read_only);
	echelon::group gr_parameters = hdf5_file["parameters"];

	double var_buffer = 0.;

	echelon::scalar_dataset ds_M = gr_parameters["mach number M"];
	M  <<= ds_M;

    echelon::scalar_dataset ds_tau      = gr_parameters["damping ratio tau"];
    tau <<= ds_tau;

    echelon::scalar_dataset ds_theta    = gr_parameters["temperature ratio theta"];
	theta <<= ds_theta;

    echelon::scalar_dataset ds_mu       = gr_parameters["velocity diffusion mu"];
    mu <<= ds_mu;

    echelon::scalar_dataset ds_beta     = gr_parameters["magnetization beta"];
    beta <<= ds_beta;

    parameters Params(M,tau,theta,mu,beta);

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return Params;
}




void save_field_real(const std::string &data_set, const field_real &Target, const std::string &path)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "save_field_real(const std::string &path, const std::string &data_set, const field_real &Target)";
	std::string status_text = "path: " + path + " data_set: " + data_set;
	my_log << status_text;
   #endif

	long unsigned int Nx_ = static_cast<int>(Target.Nx);
	long unsigned int Ny_ = static_cast<int>(Target.Ny);
	long unsigned int Nz_ = static_cast<int>(Target.Nz);


	echelon::multi_array<double> conversion_array({Nx_, Ny_, Nz_}, 0.);
	for( int i=0; i<Target.Nx; ++i)
		for( int j=0; j<Target.Ny; ++j)
			for( int k=0; k<Target.Nz; ++k)
				conversion_array(i,j,k) = Target.val[Target.my_grid->index_at(i,j,k)];

	echelon::file hdf5_file(path, echelon::file::open_mode::read_write);
	echelon::group data = hdf5_file["data"];
	echelon::dataset ds = data.create_dataset<double>(data_set, {Nx_, Ny_, Nz_});
    ds <<= conversion_array;

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done.";
   #endif

	return;
}

field_real load_field_real(const std::string &path, const std::string &data_set)
//void load_field_real(double * ptr_data, const std::string &data_set, const std::string &path)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "load_field_real(const std::string &path, const std::string &data_set)";
	std::string status_text = "path: " + path + " data_set " + data_set;
	my_log << status_text;
   #endif

	grid_Co Omega = load_grid(path);

	field_real Target = field_real(Omega);
	long unsigned int Nx_ = static_cast<int>(Target.Nx);
	long unsigned int Ny_ = static_cast<int>(Target.Ny);
	long unsigned int Nz_ = static_cast<int>(Target.Nz);

	echelon::file hdf5_file(path, echelon::file::open_mode::read_only);
	echelon::group data = hdf5_file["data"];
	echelon::dataset ds = data[data_set];

	echelon::multi_array<double> conversion_array({Nx_ ,Ny_ ,Nz_ }, 0.);
    echelon::auto_reshape(conversion_array) <<= ds;

	for( int i=0; i<Target.Nx; ++i)
		for( int j=0; j<Target.Ny; ++j)
			for( int k=0; k<Target.Nz; ++k)
				Target.val[Target.my_grid->index_at(i,j,k)] = conversion_array(i,j,k);

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
	return Target;
}




void save_field_imag(const std::string &data_set, const field_imag &Target, const std::string &path)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "save_field_imag(const std::string &path, const std::string &data_set, const field_imag &Target, const grid &Omega)";
	std::string status_text = "path: " + data_set + " data_set " + path;
	my_log << status_text;
   #endif

	field_real tmp = field_real(*Target.my_grid);

	iFFT(Target,tmp);

	long unsigned int Nx_ = static_cast<int>(tmp.Nx);
	long unsigned int Ny_ = static_cast<int>(tmp.Ny);
	long unsigned int Nz_ = static_cast<int>(tmp.Nz);

	echelon::multi_array<double> conversion_array({Nx_, Ny_, Nz_}, 0.);
	for( int i=0; i<tmp.Nx; ++i)
		for( int j=0; j<tmp.Ny; ++j)
			for( int k=0; k<tmp.Nz; ++k)
				conversion_array(i,j,k) = tmp.val[tmp.my_grid->index_at(i,j,k)];


	echelon::file hdf5_file(path, echelon::file::open_mode::read_write);
	echelon::group data = hdf5_file["data"];
	echelon::dataset ds = data.create_dataset<double>(data_set, {Nx_, Ny_, Nz_});
    ds <<= conversion_array;

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done.";
   #endif

	return;
}

field_imag load_field_imag(const std::string &data_set, const std::string &path)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "load_field_imag(const std::string &path, const std::string &data_set)";
	std::string status_text = "path: " + data_set + " data_set " + path;
	my_log << status_text;
   #endif

	grid_Co Omega = load_grid(path);
	field_real RF = load_field_real(path, data_set);
	field_imag IF(Omega);
	FFT(RF,IF);

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done.";
   #endif

	return IF;
}




void load_particles(std::vector<particle> &particle_list, const std::string &path)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "load_particle";
   #endif



	std::ifstream input_stream(path.c_str(), std::ofstream::in);
	std::string input_line;

	if(!input_stream)
	{ // file couldn't be opened
		std::cout << "   load(): Error, file wasn't found or couldn't be opened. Abort!\n";
		return;
	}

	double x,y,z,q,R;

	getline(input_stream,input_line); // read header & discard

	if(!input_stream.eof()) do
	{
		getline(input_stream,input_line); // read next line
		std::stringstream sstr(input_line);
		sstr >> x; // x-coordinate)
		sstr >> y; // y-coordinate
		sstr >> z; // z-coordinate
		sstr >> q; // q (charge)
		sstr >> R;

		particle_list.push_back(particle(x,y,z,q,R));

	} while(!input_stream.eof()) ;

	input_stream.close();

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
}

void save_particles(std::vector<particle> &particle_list, const std::string &path)
{
   #ifdef _MY_VERBOSE
	logger my_log("particle::save");
	my_log << "save_particles(const std::string &path, std::vector<particle> &particle_list)";
   #endif

	std::ofstream output_stream(path.c_str(), std::ofstream::trunc);

	// write header
	output_stream << "x\ty\tz\tq\tR";

	// write data
	for (std::vector<particle>::iterator it=particle_list.begin(); it != particle_list.end(); ++it)
	{
       #ifdef _MY_VERBOSE
		my_log << "saving a particle";
       #endif

		output_stream << "\n";
		output_stream << (*it).get_x() << "\t";
		output_stream << (*it).get_y() << "\t";
		output_stream << (*it).get_z() << "\t";
		output_stream << (*it).get_q() << "\t";
		output_stream << (*it).get_r();
	}

	output_stream.close();

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return;
}




void save_1d(const field_real &Ux, const field_real &Uy, const field_real  &Uz, const field_real  &ni, const field_real &Ph,
		     const subdim & sdim, const std::string &path)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("save_1d");
	my_log << "start";
   #endif

	 int first_fixed;
	 int second_fixed;

	 int i_fast;
	 int N_fast;

	 int *i, *j, *k;

	int direction = sdim.default_direction;

	switch( direction )
	{
	case  0:
		first_fixed = sdim.default_ypos;
		second_fixed = sdim.default_zpos;
		i = &i_fast;
		j = &first_fixed;
		k = &second_fixed;
	    N_fast = Ux.Nx;
		break;
	case 1:
		first_fixed = sdim.default_xpos;
		second_fixed = sdim.default_zpos;
		i = &first_fixed;
		j = &i_fast;
		k = &second_fixed;
	    N_fast = Ux.Ny;
	    break;
	case 2:
		first_fixed = sdim.default_xpos;
		second_fixed = sdim.default_ypos;
		i = &first_fixed;
		j = &second_fixed;
	    k = &i_fast;
	    N_fast = Ux.Nz;
	    break;
	} // END of switch

	std::ofstream output_stream(path.c_str(), std::ofstream::trunc);

	// write header to file (3-lines of header)
	output_stream << "output from PFluidDy\n";
	output_stream << "Domain: ";
	output_stream << "[" << Ux.my_grid->x_axis->val_at(0) << "," <<  Ux.my_grid->x_axis->val_at(Ux.Nx) << "]x";
	output_stream << "[" << Ux.my_grid->y_axis->val_at(0) << "," <<  Ux.my_grid->y_axis->val_at(Ux.Nx) << "]x";
	output_stream << "[" << Ux.my_grid->z_axis->val_at(0) << "," <<  Ux.my_grid->z_axis->val_at(Ux.Nx) << "]\t";
	output_stream << "Nx/Ny/Nz ";
	output_stream << Ux.Nx << "/";
	output_stream << Ux.Ny << "/";
	output_stream << Uy.Nz << "\n";
	output_stream << "x\ty\tz\tUx\tUy\tUz\tni\n\n";

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "writing to file";
   #endif

	// write data to file
	for(i_fast=0; i_fast<N_fast; ++i_fast)
	{
		//std::cout << "(" << *i << "," << *j << "," << *k << ")" << std::endl;
		//i_fast << " " << domain::index_3d_re(*i,*j,*k) << std::endl;
		output_stream << std::scientific << std::setprecision(3);
		output_stream << Ux.my_grid->x_axis->val_at(*i) << "\t";
		output_stream << Ux.my_grid->y_axis->val_at(*j) << "\t";
		output_stream << Ux.my_grid->z_axis->val_at(*k) << "\t";
		output_stream << std::scientific << std::setprecision(6);
		int index = Ux.my_grid->index_at(*i,*j,*k);
		output_stream << Ux.val[index] << "\t";
		output_stream << Uy.val[index] << "\t";
		output_stream << Uz.val[index] << "\t";
		output_stream << ni.val[index] << "\t";
		output_stream << Ph.val[index] << "\n";
		//output_stream << std::endl;
	}

	output_stream.close();

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done.";
   #endif

	return;
}

void save_1d(const field_real &Ux, const subdim & sdim, const std::string &path)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("save_1d");
	my_log << "start";
   #endif

	 int first_fixed;
	 int second_fixed;

	 int i_fast;
	 int N_fast;

	 int *i, *j, *k;

	int direction = sdim.default_direction;

	switch( direction )
	{
	case  0:
		first_fixed = sdim.default_ypos;
		second_fixed = sdim.default_zpos;
		i = &i_fast;
		j = &first_fixed;
		k = &second_fixed;
	    N_fast = Ux.Nx;
		break;
	case 1:
		first_fixed = sdim.default_xpos;
		second_fixed = sdim.default_zpos;
		i = &first_fixed;
		j = &i_fast;
		k = &second_fixed;
	    N_fast = Ux.Ny;
	    break;
	case 2:
		first_fixed = sdim.default_xpos;
		second_fixed = sdim.default_ypos;
		i = &first_fixed;
		j = &second_fixed;
	    k = &i_fast;
	    N_fast = Ux.Nz;
	    break;
	} // END of switch



	std::ofstream output_stream(path.c_str(), std::ofstream::trunc);


	// write header to file (3-lines of header)
	output_stream << "output from PFluidDy\n";
	output_stream << "Domain: ";
	output_stream << "[" << Ux.my_grid->x_axis->val_at(0) << "," <<  Ux.my_grid->x_axis->val_at(Ux.Nx) << "]x";
	output_stream << "[" << Ux.my_grid->y_axis->val_at(0) << "," <<  Ux.my_grid->y_axis->val_at(Ux.Nx) << "]x";
	output_stream << "[" << Ux.my_grid->z_axis->val_at(0) << "," <<  Ux.my_grid->z_axis->val_at(Ux.Nx) << "]\t";
	output_stream << "Nx/Ny/Nz ";
	output_stream << Ux.Nx << "/";
	output_stream << Ux.Ny << "/";
	output_stream << Ux.Nz << "\n";
	output_stream << "x\ty\tz\tUx\tUy\tUz\tni\n\n";

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "writing to file";
   #endif

	// write data to file
	for(i_fast=0; i_fast<N_fast; ++i_fast)
	{
		//std::cout << "(" << *i << "," << *j << "," << *k << ")" << std::endl;
		//i_fast << " " << domain::index_3d_re(*i,*j,*k) << std::endl;
		output_stream << std::scientific << std::setprecision(3);
		output_stream << Ux.my_grid->x_axis->val_at(*i) << "\t";
		output_stream << Ux.my_grid->y_axis->val_at(*j) << "\t";
		output_stream << Ux.my_grid->z_axis->val_at(*k) << "\t";
		output_stream << std::scientific << std::setprecision(6);
		int index = Ux.my_grid->index_at(*i,*j,*k);
		output_stream << Ux.val[index] << "\n";
		//output_stream << std::endl;
	}

	output_stream.close();

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done.";
   #endif

	return;
}

void save_2d(const field_real &Ux, const subdim & sdim,  const std::string &path)
{
	int i_fixed;
	int i_slow;
	int i_fast;

	int N_slow;
	int N_fast;

	int *i, *j, *k;

	int direction = sdim.default_plane;
	    // 0 : x fixiert => out: y-z-ebene
		// 1 : y fixiert => out: x-z-ebene
	    // 2 : z fixiert => out: x-y-ebene

	switch( direction )
	{
	case 0:
		i_fixed = sdim.default_xpos;
		i = &i_fixed;
		j = &i_slow;
	    k = &i_fast;
	    N_slow = Ux.Ny;
	    N_fast = Ux.Nz;
	    break;
	case 1:
		i_fixed = sdim.default_ypos;
		i = &i_slow;
		j = &i_fixed;
		k = &i_fast;
	    N_slow = Ux.Nx;
	    N_fast = Ux.Nz;
	    break;
	case  2:
		i_fixed = sdim.default_zpos;
		i = &i_slow;
		j = &i_fast;
		k = &i_fixed;
	    N_slow = Ux.Nx;
	    N_fast = Ux.Ny;
		break;
	} // END of switch


	std::ofstream output_stream(path.c_str(), std::ofstream::trunc);


	// write header to file (2-lines of header)
	output_stream << "Domain: ";
	output_stream << "[" << Ux.my_grid->x_axis->val_at(0) << "," << Ux.my_grid->x_axis->val_at(Ux.Nx)<< "]x";
	output_stream << "[" << Ux.my_grid->y_axis->val_at(0) << "," << Ux.my_grid->y_axis->val_at(Ux.Nx)<< "]x";
	output_stream << "[" << Ux.my_grid->z_axis->val_at(0) << "," << Ux.my_grid->z_axis->val_at(Ux.Nx)<< "]x";
	output_stream << "Nx/Ny/Nz ";
	output_stream << Ux.Nx << "/";
	output_stream << Ux.Ny << "/";
	output_stream << Ux.Nz << "\n";
	output_stream << "x\ty\tz\tUx\tUy\tUz\tni\n\n";


   #if defined(MY_VERBOSE_MORE) || defined(MY_VERBOSE_TEDIOUS)
	std::cout << "   save(): saving data to " << path;
	std::cout << " .. this may take a while\n";
   #endif

	// write data to file
	for(i_slow=0; i_slow<N_slow; ++i_slow)
	{
		for(i_fast=0; i_fast<N_fast; ++i_fast)
		{
			output_stream << std::scientific << std::setprecision(3);
			output_stream << Ux.my_grid->x_axis->val_at(*i) << "\t";
			output_stream << Ux.my_grid->y_axis->val_at(*j) << "\t";
			output_stream << Ux.my_grid->z_axis->val_at(*k) << "\t";
			output_stream << std::scientific << std::setprecision(6);
			 int index = Ux.my_grid->index_at(*i,*j,*k);
			output_stream << Ux.val[index] << "\n";
		}
		output_stream << std::endl;
	}
	output_stream.close();

	return;
}

void save_2d(const field_real &Ux, const field_real &Uy, const field_real  &Uz, const field_real  &ni, const field_real &Ph, const subdim & sdim,
		     const std::string &path)
{
	int i_fixed;
	int i_slow;
	int i_fast;

	int N_slow;
	int N_fast;

	int *i, *j, *k;

	int direction = sdim.default_plane;
	    // 0 : x fixiert => out: y-z-ebene
		// 1 : y fixiert => out: x-z-ebene
	    // 2 : z fixiert => out: x-y-ebene

	switch( direction )
	{
	case 0:
		i_fixed = sdim.default_xpos;
		i = &i_fixed;
		j = &i_slow;
	    k = &i_fast;
	    N_slow = Ux.Ny;
	    N_fast = Ux.Nz;
	    break;
	case 1:
		i_fixed = sdim.default_ypos;
		i = &i_slow;
		j = &i_fixed;
		k = &i_fast;
	    N_slow = Ux.Nx;
	    N_fast = Ux.Nz;
	    break;
	case  2:
		i_fixed = sdim.default_zpos;
		i = &i_slow;
		j = &i_fast;
		k = &i_fixed;
	    N_slow = Ux.Nx;
	    N_fast = Ux.Ny;
		break;
	} // END of switch


	std::ofstream output_stream(path.c_str(), std::ofstream::trunc);

	// write header to file (3-lines of header)
	output_stream << "output from PFluidDy\n";
	output_stream << "Domain: ";
	output_stream << "[" << Ux.my_grid->x_axis->val_at(0) << "," << Ux.my_grid->x_axis->val_at(Ux.Nx)<< "]x";
	output_stream << "[" << Ux.my_grid->y_axis->val_at(0) << "," << Ux.my_grid->y_axis->val_at(Ux.Nx)<< "]x";
	output_stream << "[" << Ux.my_grid->z_axis->val_at(0) << "," << Ux.my_grid->z_axis->val_at(Ux.Nx)<< "]x";
	output_stream << "Nx/Ny/Nz ";
	output_stream << Ux.Nx << "/";
	output_stream << Ux.Ny << "/";
	output_stream << Uy.Nz << "\n";
	output_stream << "x\ty\tz\tUx\tUy\tUz\tni\n\n";

   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	std::cout << "   save(): saving data to ";
	std::cout << path;
	std::cout << " .. this may take a while\n";
   #endif

	// write data to file
	for(i_slow=0; i_slow<N_slow; ++i_slow)
	{
		for(i_fast=0; i_fast<N_fast; ++i_fast)
		{
			output_stream << std::scientific << std::setprecision(3);
			output_stream << Ux.my_grid->x_axis->val_at(*i) << "\t";
			output_stream << Ux.my_grid->y_axis->val_at(*j) << "\t";
			output_stream << Ux.my_grid->z_axis->val_at(*k) << "\t";
			output_stream << std::scientific << std::setprecision(6);
			 int index = Ux.my_grid->index_at(*i,*j,*k);
			output_stream << Ux.val[index] << "\t";
			output_stream << Uy.val[index] << "\t";
			output_stream << Uz.val[index] << "\t";
			output_stream << ni.val[index] << "\t";
			output_stream << Ph.val[index] << "\n";
			//output_stream << std::endl;
		}
		output_stream << std::endl;
	}
	output_stream.close();

	return;
}

void save_3d(const std::string &file_name , const field_real &Ux, const field_real &Uy, const field_real  &Uz, const field_real  &ni, const field_real &Ph)
{
	std::string path = "./data/"+file_name;
	std::ofstream output_stream(path.c_str(), std::ofstream::trunc);

	// write header to file (3-lines of header)
	output_stream << "output from PFluidDy\n";
	output_stream << "Domain: ";
	output_stream << "Nx= " << Ux.my_grid->x_axis->N <<" ";
	output_stream << "Ny= " << Ux.my_grid->y_axis->N <<" ";
	output_stream << "Nz= " << Ux.my_grid->z_axis->N <<" ";
	output_stream << "x0= " << Ux.my_grid->x_axis->l0 <<" ";
	output_stream << "y0= " << Ux.my_grid->y_axis->l0 <<" ";
	output_stream << "z0= " << Ux.my_grid->z_axis->l0 <<" ";
	output_stream << "Lx= " << Ux.my_grid->x_axis->L <<" ";
	output_stream << "Ly= " << Ux.my_grid->y_axis->L <<" ";
	output_stream << "Lz= " << Ux.my_grid->z_axis->L <<"\n";
	output_stream << "x\ty\tz\tUx\tUy\tUz\tni\n\n";

	std::cout << "   save(): saving data to ";
	std::cout << file_name;
	std::cout << " .. this may take a while\n";
	// write data to file
	for( int i=0; i<Ux.Nx; ++i)
	{
		for( int j=0; j<Ux.Ny; ++j)
		{
			for( int k=0; k<Ux.Nz; ++k)
			{
				output_stream << std::scientific << std::setprecision(3);
				output_stream << Ux.my_grid->x_axis->val_at(i) << "\t";
				output_stream << Ux.my_grid->y_axis->val_at(j) << "\t";
				output_stream << Ux.my_grid->z_axis->val_at(k) << "\t";
				output_stream << std::scientific << std::setprecision(6);
				 int index = Ux.my_grid->index_at(i,j,k);
				output_stream << Ux.val[index] << "\t";
				output_stream << Uy.val[index] << "\t";
				output_stream << Uz.val[index] << "\t";
				output_stream << ni.val[index] << "\t";
				output_stream << Ph.val[index] << "\n";
			}
			output_stream << std::endl;
		}
		//output_stream << std::endl;
	}

	return;
}

void save_3d(const field_real &XX, const std::string &path)
{
	std::ofstream output_stream(path.c_str(), std::ofstream::trunc);

	std::cout << "   save(): saving data to ";
	std::cout << path;
	std::cout << " .. this may take a while\n";
	// write data to file
	for( int i=0; i<XX.Nx; ++i)
	{
		for( int j=0; j<XX.Ny; ++j)
		{
			for( int k=0; k<XX.Nz; ++k)
			{
				output_stream << std::scientific << std::setprecision(3);
				output_stream << XX.my_grid->x_axis->val_at(i) << "\t";
				output_stream << XX.my_grid->y_axis->val_at(j) << "\t";
				output_stream << XX.my_grid->z_axis->val_at(k) << "\t";
				output_stream << std::scientific << std::setprecision(6);
				int index = XX.my_grid->index_at(i,j,k);
				output_stream << XX.val[index] << "\n";
			}
			output_stream << std::endl;
		}
		//output_stream << std::endl;
	}

	return;
}




void save_slice(const std::string &data_set, const field_imag &Target, const std::string &path)
{
   #if defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	std::string status_text = "start - path: " + path + "data_set: " + data_set;
	my_log << status_text;
   #endif


	field_real tmp_3d_out(*Target.my_grid);
	iFFT(Target,tmp_3d_out);



	long unsigned int Nx_ = static_cast<int>(tmp_3d_out.Nx);
	long unsigned int Ny_ = static_cast<int>(tmp_3d_out.Ny);
	long unsigned int Nz_ = static_cast<int>(tmp_3d_out.Nz);

/*
	int Nx_ = tmp_3d_out.Nx;
	int Ny_ = tmp_3d_out.Ny;
	int Nz_ = tmp_3d_out.Nz;
*/

	echelon::multi_array<double> conversion_array({Nx_, Ny_}, 0.);

	int k = tmp_3d_out.my_grid->z_axis->index_at(0.);
	for(int i=0; i<tmp_3d_out.Nx; ++i)
		for(int j=0; j<tmp_3d_out.Ny; ++j)
		{
			int index = tmp_3d_out.my_grid->index_at(i,j,k);
			conversion_array(i,j) = tmp_3d_out.val[index];

		}

	echelon::file hdf5_file(path, echelon::file::open_mode::read_write);
	echelon::group data = hdf5_file["data"];
	echelon::dataset ds = data.create_dataset<double>(data_set, {Nx_, Ny_});

	ds <<= conversion_array;

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return;
}

void load_slice(const std::string &data_set, double * ptr_data, const std::string &path)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("IO");
	my_log << "load_slice(const std::string &path, const std::string &data_set)" ;
   #endif

	grid_Co Omega = load_grid(path);

	/*
	long unsigned int Nx_ = static_cast<int>(Omega.Nx);
	long unsigned int Ny_ = static_cast<int>(Omega.Ny);
	long unsigned int Nz_ = static_cast<int>(Omega.Nz);
*/

	long unsigned int Nx_ = (long unsigned int) Omega.x_axis->N;
	long unsigned int Ny_ = (long unsigned int) Omega.y_axis->N;
	long unsigned int Nz_ = (long unsigned int) Omega.z_axis->N;

	echelon::file hdf5_file(path, echelon::file::open_mode::read_only);
	echelon::group data = hdf5_file["data"];
	echelon::dataset ds = data[data_set];

	echelon::multi_array<double> conversion_array({Nx_,Ny_}, 0.);
    echelon::auto_reshape(conversion_array) <<= ds;




   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "allocate temporary real field and copy data from echelon::multi_array";
   #endif
	int k = Omega.z_axis->index_at(0.);
	for(int i=0; i<Omega.Nx; ++i)
		for(int j=0; j<Omega.Ny; ++j)
		{

			int index = j + i*Ny_;
			ptr_data[index] = conversion_array(i,j);

           #if defined(_MY_VERBOSE_TEDIOUS)
			std::stringstream sstr;
			sstr << i << "\t" << j << "\t" << conversion_array(i,j);
			std::string lineout = sstr.str();
			my_log << lineout;
           #endif
		}


   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
	return;
}


std::string I4(int number){
  std::stringstream s_stream;
  s_stream.width(4);
  s_stream.fill('0');
  s_stream.setf(std::ios::right);
  s_stream << number << std::flush;
  return s_stream.str();}







void save_all(std::vector<particle> &particle_list,
		      const grid_Co &Omega,
			  const double &Time,
		      const parameters &Params,
		      const field_imag &FUx, const field_imag &FUy, const field_imag &FUz,
		      const field_imag &Fni, const field_imag &FPh,
		      const std::string &path
		)
{

	file_create(path);
	save_grid(Omega, path);
	save_time(Time, path);
	save_parameters(Params, path);
	save_field_imag("Ux", FUx, path);
	save_field_imag("Uy", FUy, path);
	save_field_imag("Uz", FUz, path);
	save_field_imag("ni", Fni, path);
	save_field_imag("Ph", FPh, path);
	return;
}


void save_slice(const int &step,
		  std::vector<particle> &particle_list,
	      const grid_Co &Omega,
		  const double &Time,
	      const parameters &Params,
	      const field_imag &FUx, const field_imag &FUy, const field_imag &FUz,
	      const field_imag &Fni, const field_imag &FPh
	)
{
	std::string path = "./data/slice_"+I4(step)+".h5";

	save_particles(particle_list, "./config/particles.dat");
	file_create(path);
	save_grid(Omega, path);
	save_time(Time, path);
	save_parameters(Params, path);
	save_slice("Ux", FUx, path);
	save_slice("Uy", FUy, path);
	save_slice("Uz", FUz, path);
	save_slice("ni", Fni, path);
	save_slice("Ph", FPh, path);

	return;
}


void save_frame(const grid_Co &Omega,
		        double* Ux, double* Uy, double* Uz, double* ni, double* Ph,
			    const std::string &path = "./data/frame.dat"
			    )
{
	std::ofstream output_stream(path.c_str(), std::ofstream::trunc);

	// write data to file
	for(int i=0; i<Omega.Nx; ++i)
	{
		for(int j=0; j<Omega.Ny; ++j)
		{
			output_stream << std::scientific << std::setprecision(3);
			output_stream << Omega.x_axis->val_at(i) << "\t";
			output_stream << Omega.y_axis->val_at(j) << "\t";
			output_stream << std::scientific << std::setprecision(6);
			int index = j + Omega.Ny*i;
			output_stream << Ux[index] << "\t";
			output_stream << Uy[index] << "\t";
			output_stream << Uz[index] << "\t";
			output_stream << ni[index] << "\t";
			output_stream << Ph[index] << "\n";

		}
		output_stream << std::endl;
	}
	output_stream.close();
	std::cout << "   save(): success\n";
	return;
}


void save_major_wavevectors(const std::string &filename, field_imag &data)
{
	std::ofstream output_stream(filename.c_str(), std::ofstream::trunc);

	// Axis Labels
	output_stream << "ik\tRe_Fx\tIm_Fx";
	output_stream << "jk\tRe_Fy\tIm_Fy";
	output_stream << "kk\tRe_Fz\tIm_Fz\n\n";

	std::stringstream sstr;
	std::string output_line;

	double NMax = max<int>(data.Nx, data.Ny);
	       NMax = max<int>(NMax   , data.Nz);

	for(int iter=0; iter<NMax; ++iter)
	{
		int index;

		sstr << iter << "\t";


		// write info about major kx into stream
		if(iter<data.Nx)
		{
			index = data.my_grid->index_at(iter,0,0);
			sstr << data.val[index][0] << "\t" << data.val[index][1];
		}
		else
		{
			sstr << " \t ";
		}
		sstr << "\t";

		// write info about major ky into stream
		if(iter<data.Ny)
		{
			index = data.my_grid->index_at(0,iter,0);
			sstr << data.val[index][0] << "\t" << data.val[index][1];
		}
		else
		{
			sstr << "- -";
		}
		sstr << "\t";

		// write info about major kz into stream
		if(iter<data.Nz)
		{
			index = data.my_grid->index_at(0,0,iter);
			sstr << data.val[index][0] << "\t" << data.val[index][1];
		}
		else
		{
			sstr << "- -";
		}
		sstr << "\n";

		// write stream to file (via string)
		output_line = sstr.str();
		sstr.str(std::string()); // clear stream
		output_stream << output_line;
	}
	output_stream.close();

	return;
}

/*




void load_field_imag(field_imag &field, const std::string &data_set, const std::string &path)
{
   #ifdef _MY_VERBOSE
	logger my_log("IO.CPP: load_field_imag");
	std::string status_text = "start - path: " + path + "data_set: " + data_set;
	my_log << status_text;
   #endif

   #ifdef _MY_VERBOSE
	my_log << "open .h5 file and load data to echelon::multi_array";
   #endif
	echelon::file hdf5_file(path, echelon::file::open_mode::read_only);
	echelon::group data = hdf5_file["data"];
	echelon::dataset ds = data[data_set];

	 int Nx = field.my_domain->x_axis->N;
	 int Ny = field.my_domain->y_axis->N;
	 int Nz = field.my_domain->z_axis->N;

	echelon::multi_array<double> conversion_array({Nx,Ny,Nz}, 0.);
    echelon::auto_reshape(conversion_array) <<= ds;

   #ifdef _MY_VERBOSE
	my_log << "allocate temporary real field and copy data from echelon::multi_array";
   #endif
	double * ptr_data = (double*) fftw_malloc(sizeof(double) * (Nx*Ny*Nz));
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
				ptr_data[field.my_domain->index_real(i,j,k)] = conversion_array(i,j,k);

   #ifdef _MY_VERBOSE
	my_log << "Tranform temporary real field to fourier-space and clean up";
   #endif
	operation TraFo;
	TraFo.FFT(ptr_data,field.val);
	fftw_free(ptr_data);

   #ifdef _MY_VERBOSE
	my_log << "done";
   #endif
	return;
}

*/





