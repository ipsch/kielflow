#include "solver_poisson_jacobi_nlin.hpp"

int solver_poisson_jacobi_nlin::iterations_total = 0;



solver_poisson_jacobi_nlin::solver_poisson_jacobi_nlin(interface_3d_fkt &boundary, interface_3d_fkt &val_boundary, const double &w) :
	H(boundary), val_H(val_boundary), eps(0.0001), omega_NEWTON(w), omega_SOR(1.)
{
	//limit_max = 2.3e-05;
	//limit_sum = 3.0e-07;
	norm_sum_old = -1.;
	norm_max_old = -1.;
	limit_max = 1.45e-06;
	limit_sum = 7.9e-08;
	norm_max = 0.;
	norm_sum = 0.;
	supremum = 0.;
	infinum = 0.;
	iteration = 0;
	invocation = 0;
	max_iterations = 50;
	converged_=false;
	use_boundary_=false;

	HX = 0L;
	HY = 0L;
	HZ = 0L;

	// create logfile + Header
	my_logfile = "./diagnostics/norm_max.log";
	std::ofstream output_stream(my_logfile, std::ofstream::trunc);
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

	invocation++;
	iteration=0;
	converged_=false;
	omega_SOR = 1.;
	norm_sum_old = -1.;
	field_real Phi_n(*Phi_IO.my_grid);
	H_create(*Phi_IO.my_grid);

	if(max_iterations<0) max_iterations = 1000;

	while (!converged_&&(iteration<max_iterations))
	{
    	iteration++;
    	iteration_loop(Phi_IO, Phi_n, rho);

		Phi_IO = Phi_n;
		save_evolution(Phi_IO, rho);
	}

    H_delete();

   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
	return;
}



void solver_poisson_jacobi_nlin::H_create(const grid_Co &Omega)
{
	HX = new double[Omega.Nx+1];
	HY = new double[Omega.Ny+1];
	HZ = new double[Omega.Nz+1];

	for(int i=0; i<Omega.Nx+1;++i)
		HX[i] = Omega.x_axis->val_at(i) - Omega.x_axis->val_at(i-1);
	for(int j=0; j<Omega.Ny+1;++j)
		HY[j] = Omega.y_axis->val_at(j) - Omega.y_axis->val_at(j-1);
	for(int k=0; k<Omega.Nz+1;++k)
		HZ[k] = Omega.z_axis->val_at(k) - Omega.z_axis->val_at(k-1);

	return;
}

void solver_poisson_jacobi_nlin::H_delete()
{
	delete[] HX;
	delete[] HY;
	delete[] HZ;
}




void solver_poisson_jacobi_nlin::iteration_loop(const field_real &in, field_real &out, const field_real &rho)
{
   #if defined(MY_VERBOSE_MORE) || defined(MY_VERBOSE_TEDIOUS)
	logger my_log("solver_poisson_relaxation");
	my_log << "iteration_loop";
   #endif

	iterations_total++;
	norm_max = 0.;
	norm_sum = 0.;
	supremum = in.val[0];
	infinum = in.val[0];

	newton_method NL_solver(0.05,20);

	for(int i=0; i < out.Nx; ++i)
		for(int j=0; j<out.Ny; ++j)
			for(int k=0; k<out.Nz; ++k)
			{
				double A=0.;
				double D=0.;

				// Diagonal elemente
				D += 2. / (HX[i+1]*HX[i]);
				D += 2. / (HY[j+1]*HY[j]);
				D += 2. / (HZ[k+1]*HZ[k]);

				// Nebendiagonale
				A += in.val_at(i+1,j,k) / (HX[i+1]*(HX[i+1]+HX[i]));
				A += in.val_at(i-1,j,k) / (HX[i]*(HX[i+1]+HX[i]));
				A += in.val_at(i,j+1,k) / (HY[j+1]*(HY[j+1]+HY[j]));
				A += in.val_at(i,j-1,k) / (HY[j]*(HY[j+1]+HY[j]));
				A += in.val_at(i,j,k+1) / (HZ[k+1]*(HZ[k+1]+HZ[k]));
				A += in.val_at(i,j,k-1) / (HZ[k]*(HZ[k+1]+HZ[k]));
				A = 2.*A;
				A += rho(i,j,k);

				double P_old = in(i,j,k);
				double P_new = NL_solver.solve(P_old,
						[&] (double x) {return A -D*x - exp(x);},
						[&] (double x) {return -D -exp(x);});

				double delta = fabs(P_new - P_old);
				//out(i,j,k) = P_new;
				out(i,j,k) = (1.-omega_SOR)*P_old + omega_SOR*P_new;

				supremum = max<double>(supremum,P_new);
				infinum = min<double>(infinum,P_new);
				norm_max = max<double>(norm_max,delta);
				norm_sum += delta;
			}

	// B. A. Carré :
	// "The Determination of the Optimum Accelerating Factor for Successive Over-relaxation"
	if(iteration>1) //
	{

		//if( (norm_sum<norm_sum_old) )
		if( (norm_sum<norm_sum_old) && (norm_max<norm_max_old) )
		{
			double omega_old = omega_SOR;
			double omega_new;
			double lambda = norm_sum/norm_sum_old;
			omega_new = 2./(1.+sqrt(1.-pow(lambda+omega_SOR-1.,2.)/(lambda*pow(omega_SOR,2.))));
			omega_SOR = .5*omega_old + 0.5*omega_new;
		}
		else
		{
			omega_SOR =.5;
		}
	}

	norm_sum_old = norm_sum;
	norm_max_old = norm_max;

	if(use_boundary_)
		for(int i=0; i < out.Nx; ++i)
			for(int j=0; j<out.Ny; ++j)
				for(int k=0; k<out.Nz; ++k)
				{ // begin loop over k
					double x = in.my_grid->x_axis->val_at(i);
					double y = in.my_grid->y_axis->val_at(j);
					double z = in.my_grid->z_axis->val_at(k);
					if( H(x,y,z)==0. )
						out(i,j,k) = val_H(x,y,z);
					double delta = fabs(out(i,j,k) - in(i,j,k));
					double P_new = val_H(x,y,z);

					supremum = max<double>(supremum,P_new);
					infinum = min<double>(infinum,P_new);
					norm_max = max<double>(norm_max,delta);
					norm_sum += delta;
				}




	norm_sum = norm_sum/(out.N);
	if(norm_max/(supremum-infinum) <= limit_max)
		converged_=true;
	bool ESC_sum = (norm_sum < limit_sum);



   #if defined(MY_VERBOSE_MORE) || defined(MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
	return;
}







void solver_poisson_jacobi_nlin::save_evolution(const field_real &Phi, const field_real &rho) const
{

   // #define EVOLUTION_LOGGING
   #if defined(EVOLUTION_LOGGING)
	subdim my_dim;
	my_dim.xpos = Phi.Nx/2;
	my_dim.ypos = Phi.Ny/2;
	my_dim.zpos = Phi.Nz/2;
	my_dim.direction = 0;
	my_dim.plane = 2;

	std::string filenameOne = "./data/nlj_rho_" + ConvertToString<int>(iterations_total) + ".dat";
	std::string filenameTwo = "./data/nlj_phi_" + ConvertToString<int>(iterations_total) + ".dat";
	save_2d(rho,my_dim,filenameOne);
	save_2d(Phi,my_dim,filenameTwo);
   #endif


	std::ofstream output_stream(my_logfile, std::ofstream::app);
	output_stream << invocation << "\t";
	output_stream << iteration << "\t";
	output_stream << norm_max << "\t";
	output_stream << norm_sum << "\t";
	output_stream << supremum-infinum << "\t";
	output_stream << norm_max/(supremum-infinum) << "\t";
	output_stream << limit_max << "\t";
	output_stream << omega_SOR << std::endl;
	output_stream.close();

	return;
}




