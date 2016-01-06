#include "solver_poisson_jacobi_nlin.hpp"

int solver_poisson_jacobi_nlin::iterations_total = 0;




static const int width = 20;

void show_percent(int i) {
     int dashes = (width * i)/100;


     std::cout << ":" << std::endl;

     std::cout << '[' << std::left << std::setw(width) << std::string(dashes, '#') << ']' << std::setw(3) << i << "%";
     return;
}








solver_poisson_jacobi_nlin::solver_poisson_jacobi_nlin(interface_3d_fkt &boundary, interface_3d_fkt &val_boundary, const double &w) :
	H(boundary), val_H(val_boundary), eps(0.0001), omega_SOR(w)
{
	my_logfile = "./diagnostics/norm_max.log";
	std::ofstream output_stream(my_logfile, std::ofstream::trunc);

	//limit_max = 2.3e-05;
	//limit_sum = 3.0e-07;
	limit_max = 1.45e-06;
	limit_sum = 7.9e-08;
	norm_max = 0.;
	norm_sum = 0.;

	iteration = 0;
	invocations = 0;
	max_iterations = 50;

	hxp = 0;
	hxm = 0;
	hxmm = 0.;
	hyp = 0;
	hym = 0;
	hymm = 0;
	hzp = 0;
	hzm = 0;
	hzmm = 0;
	Hx = Hy = Hz = 0.;

	output_stream << "s" << "\t";
	output_stream << "i" << "\t";
	output_stream << "NM" << "\t";
	output_stream << "NSt" << "\t";
	output_stream << "NSi" << "\t";
	output_stream << "NSj" << "\t";
	output_stream << "NSk" << "\n";

	output_stream.close();
}





void solver_poisson_jacobi_nlin::solve(field_real &Phi_IO, field_real &rho)
{
   #if defined(_MY_VERBOSE) ||  defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_jacobi_nlin::solve(..)");
	my_log << "start";
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	std::cout << "NLJ Nmax " << max_iterations << " ";
	std::cout << "Norm(M/S) " << limit_max;
	std::cout << "/" << limit_sum << "\n";
   #endif
   #endif





	field_real Phi_n(*Phi_IO.my_grid);

	// Abbruch Bedingungen
	bool ESC_max, ESC_sum, ESC_auto, ESC_iter;
	ESC_auto = (max_iterations < 0);

	iteration=0;
	invocations++;

	// set boundaries anew (might have changed due FFT/iFFT)
	for(int i=0; i<Phi_IO.Nx; ++i)
	{
		double x = Phi_IO.my_grid->x_axis->val_at(i);
		for(int j=0; j<Phi_IO.Ny; ++j)
		{
			double y = Phi_IO.my_grid->y_axis->val_at(j);
			for(int k=0; k<Phi_IO.Nz; ++k)
			{
				double z = Phi_IO.my_grid->z_axis->val_at(k);

				if(H(x,y,z)==1.)
				{
					Phi_IO(i,j,k) = val_H(x,y,z);
				}

				if((i==0) || (j==0) || (k==0))
				{
					Phi_IO(i,j,k) = 0.;
				}
			} // end loop over k
		} // end loop over j
	} // end loop over i

    do // this is the core of the mainloop
	{
    	iteration++;
		iteration_loop(Phi_IO, Phi_n, rho);
		check_Norms(Phi_IO, Phi_n);
		Phi_IO = Phi_n;

		ESC_max = (norm_max < limit_max);
		ESC_sum = (norm_sum < limit_sum);
		ESC_iter = (iteration<max_iterations);
	    //ESC_auto = (max_iterations < 0);  // was set before main-loop (just reminder here)

       #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_TEDIOUS)
		std::cout << invocations << " " << iteration << "\t";
		std::cout << norm_max << "/" << limit_max << "\t";
		std::cout << norm_sum << "/" << limit_sum << "\t";
		std::cout << "(" << ESC_max << ";" << ESC_sum << ";";
		std::cout << ESC_iter << ";" << ESC_auto << ")" << std::endl;
       #endif

		if(ESC_max && ESC_sum)
		{
		   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
			my_log << "done";
		   #endif
			return;
		}

	} while (ESC_iter || ESC_auto);

   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return;
}





double solver_poisson_jacobi_nlin::get_HXX(const axis * const A, const int &i, double &hp, double &hm, double &hmm) const
{
	hmm= A->val_at(i) - A->val_at(i-2);
	hm = A->val_at(i) - A->val_at(i-1);
	hp = A->val_at(i+1) - A->val_at(i);



	//if(i==A->N-1) hp = A->val_at(1) - A->val_at(0);
	//if(i==0)      hm = hp;

	return hp*(hp+hm)*hm;
}



double solver_poisson_jacobi_nlin::get_PG(const field_real &in, int i, int j, int k) const
// returns PHI or boundary value depending on position (i,j,k)
{
	if(i==in.Nx)
		i=0;
	if(i==-1)
		i = in.Nx-1;

	if(j==in.Ny)
		j=0;
	if(j==-1)
		j = in.Ny-1;

	if(k==in.Nz)
		k=0;
	if(k==-1)
		k = in.Nz-1;

	return in(i,j,k);
}





void solver_poisson_jacobi_nlin::iteration_loop(const field_real &in, field_real &out, const field_real &rho)
{
   #if defined(MY_VERBOSE_MORE) || defined(MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_relaxation");
	my_log << "iteration_loop";
   #endif

	double x, y, z;

	for(int i=0; i < out.Nx; ++i)
	{
		x = in.my_grid->x_axis->val_at(i);
		Hx = get_HXX(in.my_grid->x_axis, i, hxp, hxm, hxmm);
		for(int j=0; j<out.Ny; ++j)
		{
			y = in.my_grid->y_axis->val_at(j);
			Hy = get_HXX(in.my_grid->y_axis, j, hyp, hym, hymm);
			for(int k=0; k<out.Nz; ++k)
			{
				z = in.my_grid->z_axis->val_at(k);
				Hz = get_HXX(in.my_grid->z_axis, k, hzp, hzm, hzmm);
				int index = out.my_grid->index_at(i,j,k);

				if( (H(x,y,z)==0.) && (i!=0) && (j!=0) && (k!=0) )
				{ // if not within boundary layer : calculate new potential
					//std::cout << i << "\t" << j << "\t" << k <<  std::endl;
					double rho_ijk = rho(i,j,k);
					out.val[index] = newton(i,j,k,in,rho_ijk);



				} // END if (mask)
				else
				{ // if within boundary layer : keep boundary value
					out.val[index] = in.val[index];
					//std::cout << x << "\t" << y << "\t" << z << "\t" << in.val[index] << std::endl;
				}

			} // enf loop over k
		} // end loop over j
	} // end loop over i

/*
	subdim my_dim;
	my_dim.default_xpos = out.Nx/2;
	my_dim.default_ypos = out.Ny/2;
	my_dim.default_zpos = out.Nz/2;
	my_dim.default_direction = 0;
	my_dim.default_plane = 2;

	std::string filename = "./data/nlj_" + ConvertToString<int>(iterations_total) + ".dat";
	save_2d(out,out,out,rho_ges,out,my_dim,filename);
*/
	iterations_total++;

   #if defined(MY_VERBOSE_MORE) || defined(MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif

	return;
}




void solver_poisson_jacobi_nlin::check_Norms(const field_real &field_new, const field_real &field_old)
{
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_jacobi_nlin::check_Norms(..)");
	my_log << "start";
   #endif

	norm_max = 0.;
	norm_sum = 0.;

	for(int i=0; i<field_new.N; ++i)
	{
		double delta = fabs(field_new.val[i] - field_old.val[i] );
		norm_max = max<double>(norm_max,delta);
		norm_sum += delta;
	}

	norm_sum = norm_sum/(field_new.N);

   #if defined(_MY_VERBOSE) || defined(_MY_BERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	std::stringstream status;
	status << invocations << "\t";
	status << iteration << "\t";
	status << norm_max << "\t";
	status << norm_sum << "\n";

	std::ofstream output_stream(my_logfile, std::ofstream::app);
	output_stream << status.str();
	output_stream.close();
   #if defined(_MY_BERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
   #if defined(_MY_VERBOSE_TIDEOUS)
	std::cout << status.str() << std::endl;
	my_log << "done";
   #endif
   #endif
   #endif

	return;
}


double solver_poisson_jacobi_nlin::newton(const int &i, const int j, const int k,\
		const field_real &Phi, const double &rho_ijk) const
{

	double x_ip1, x_i;
	x_ip1 = x_i = Phi(i,j,k);
	const int max_iter = 30;
	double dx_history[max_iter];

	for(int iter=0; iter<max_iter; ++iter)
	{
		double dx = f_df(x_ip1,i,j,k,Phi,rho_ijk);
		x_ip1 -= omega_SOR*dx;
		dx_history[iter] = dx;
		if( fabs(dx) < eps  ) return x_ip1;
	}

   #if  defined(_MY_VERBOSE_LESS) || defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	logger log("solver_poisson_jacobi_nlin::newton(..)");
	std::stringstream out_stream;
	std::string out_line;
	out_stream << "WARNING@ (" << i << "," << j << "," << k << ")";
	out_line = out_stream.str();
	log << out_line;
	log << "Phi: ";
	log << get_PG(Phi,i,j,k);
	log << get_PG(Phi,i+1,j,k);
	log << get_PG(Phi,i-1,j,k);
	log << get_PG(Phi,i,j+1,k);
	log << get_PG(Phi,i,j-1,k);
	log << get_PG(Phi,i,j,k+1);
	log << get_PG(Phi,i,j,k-1);
	log << "hxmm/hxm/hxp/hymm/hym/hyp/hzmm/hzm/hzp: ";
	log << hxmm;
	log << hxm;
	log << hxp;
	log << hymm;
	log << hym;
	log << hyp;
	log << hzmm;
	log << hzm;
	log << hzp;
	log << "rho_ijk: ";
	log << rho_ijk;
	log << "history newton : ";
	for(int iter=0; iter<max_iter; ++iter)
	{
		log << dx_history[iter];
	}
	log << ":-(";
   #endif

	if(dx_history[0]>dx_history[max_iter])
		return x_ip1;

	log << "error was fatal";

	throw("newton_method: maximum number of iterations exceeded");
	return x_i;
}


double solver_poisson_jacobi_nlin::f_df(const double &x,
		const int &i, const int j, const int k,
		const field_real &Phi, const double &rho_ijk) const
{
	// Zähler



	double Ax = 2.*(hxp-hxm)/(hxmm*(hxmm+hxp)*(hxmm-hxm));
	double Bx = 2.*(hxmm-hxp)/(hxm*(hxmm-hxm)*(hxm+hxp));
	double Cx = 2.*(hxmm+hxm)/(hxp*(hxm+hxp)*(hxmm+hxp));

	double Ay = 2.*(hyp-hym)/(hymm*(hymm+hyp)*(hymm-hym));
	double By = 2.*(hymm-hyp)/(hym*(hymm-hym)*(hym+hyp));
	double Cy = 2.*(hymm+hym)/(hyp*(hym+hyp)*(hymm+hyp));

	double Az = 2.*(hzp-hzm)/(hzmm*(hzmm+hzp)*(hzmm-hzm));
	double Bz = 2.*(hzmm-hzp)/(hzm*(hzmm-hzm)*(hzm+hzp));
	double Cz = 2.*(hzmm+hzm)/(hzp*(hzm+hzp)*(hzmm+hzp));


	double f = 0.;

	f += Ax*get_PG(Phi,i-2,j,k);
	f += Bx*get_PG(Phi,i-1,j,k);
	f += -(Ax+Bx+Cx)*x;
	f += Cx*get_PG(Phi,i+1,j,k);

	f += Ay*get_PG(Phi,i,j-2,k);
	f += By*get_PG(Phi,i,j-1,k);
	f += -(Ay+By+Cy)*x;
	f += Cy*get_PG(Phi,i,j+1,k);

	f += Az*get_PG(Phi,i,j,k-2);
	f += Bz*get_PG(Phi,i,j,k-1);
	f += -(Az+Bz+Cz)*x;
	f += Cz*get_PG(Phi,i,j,k+1);

	f += rho_ijk - exp(x) + exp(-10.*x);

	// Nenner
	double df = 0.;
	df += -(Ax+Bx+Cx);
	df += -(Ay+By+Cy);
	df += -(Az+Bz+Cz);
	df += -exp(x); // -10.*exp(-10.*x);

	return (f/df);
}

