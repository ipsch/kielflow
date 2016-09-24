#include "solver_poisson_jacobi_lin.hpp"



solver_poisson_jacobi_lin::solver_poisson_jacobi_lin(interface_3d_fkt &boundary, interface_3d_fkt &val_boundary) :
	H(boundary), val_H(val_boundary)
{
	norm_maximum   =0.;
	norm_sum_total =0.;
	norm_sum_i     =0.;
	norm_sum_j     =0.;
	norm_sum_k     =0.;

	invocations = 0;
	iteration= 0;
	max_iterations = 1000;

	hxp = 0.;
	hxm = 0.;
	hyp = 0.;
	hym = 0.;
	hzp = 0.;
	hzm = 0.;

#ifdef _MY_VERBOSE
	std::ofstream output_stream("./diagnostics/relaxation.log", std::ofstream::trunc);


	output_stream << "s" << "\t";
	output_stream << "i" << "\t";
	output_stream << "NM" << "\t";
	output_stream << "NSt" << "\t";
	output_stream << "NSi" << "\t";
	output_stream << "NSj" << "\t";
	output_stream << "NSk" << "\n";

	output_stream.close();
#endif
}


void solver_poisson_jacobi_lin::solve(field_real &Phi_IO, field_real &rho)
{
   //#ifdef _MY_VERBOSE
	//logger log("linear_jacobi::solve");
	//log << "start";
   //#endif

	iteration=0;
	invocations++;

	field_real Phi_n(Phi_IO.my_grid);

	// set boundaries anew (might have changed due FFT/iFFT)
	for(int i=0; i<Phi_IO.Nx; ++i)
	{
		double x = Phi_IO.my_grid.x_axis->val_at(i);
		for(int j=0; j<Phi_IO.Ny; ++j)
		{
			double y = Phi_IO.my_grid.y_axis->val_at(j);
			for(int k=0; k<Phi_IO.Nz; ++k)
			{
				double z = Phi_IO.my_grid.z_axis->val_at(k);
				if(H(x,y,z)==1.)
				{
					Phi_IO(i,j,k) = val_H(x,y,z);
				}
			} // end loop over k
		} // end loop over j
	} // end loop over i

	do
	{
		iteration++;

		iteration_loop(Phi_IO, Phi_n, rho);
		// ToDo :
		// compare Phi_old vs. Phi_new
		// give epsilon if converged
       #ifdef _MY_VERBOSE
		check_convergence(Phi_IO, Phi_n);
       #endif

		Phi_IO = Phi_n;

	} while (iteration<max_iterations);

  // #ifdef _MY_VERBOSE
//	log << "done";
   //#endif

	return;
}








double solver_poisson_jacobi_lin::get_HXX(const axis * const A, const int &i, double &hp, double &hm) const
{
	hp = A->val_at(i+1) - A->val_at(i);
	hm = A->val_at(i) - A->val_at(i-1);

	if(i==A->N-1) hp = A->val_at(1) - A->val_at(0);
	if(i==0)     hm = hp;

	return hp*(hp+hm)*hm;
}








double solver_poisson_jacobi_lin::get_PG(const field_real &in, const int &i, const int &j, const int &k) const
// returns PHI or boundary value depending on position (i,j,k)
{



	double x = in.my_grid.x_axis->val_at(i);
	double y = in.my_grid.y_axis->val_at(j);
	double z = in.my_grid.z_axis->val_at(k);


	if(i==0) return 0.;
	if(j==0) return 0.;
	if(k==0) return 0.;
	if(i==in.Nx) return 0.;
	if(j==in.Ny) return 0.;
	if(k==in.Nz) return 0.;




	if( (i==-1) || (j==-1) || (k==-1) )
	{
		int iT = i;
		int jT = j;
		int kT = k;

		if(i==-1) iT = in.Nx-1;
		if(j==-1) jT = in.Ny-1;
		if(k==-1) kT = in.Nz-1;

		return in(iT,jT,kT);

	}


	return in(i,j,k);


}






void solver_poisson_jacobi_lin::iteration_loop(const field_real &in, field_real &out, const field_real &rho_ges)
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE)
	logger log("linear_jacobi::iteration_loop");
	log << iteration;
   #endif

	for(int i=0; i < out.Nx; ++i)
	{
		double Hx = get_HXX(in.my_grid.x_axis, i, hxp, hxm);
		double x = in.my_grid.x_axis->val_at(i);
		for(int j=0; j<out.Ny; ++j)
		{
			double Hy = get_HXX(in.my_grid.y_axis, j, hyp, hym);
			double y = in.my_grid.y_axis->val_at(j);
			for(int k=0; k<out.Nz; ++k)
			{
				double Hz = get_HXX(in.my_grid.z_axis, k, hzp, hzm);
				double z = in.my_grid.z_axis->val_at(k);

				int index = out.index(i,j,k);

				if( (H(x,y,z)==0.) || (i==0) || (j==0) || (k==0) )
				{ // if not within boundary layer : calculate new potential

					out.val[index] = 0.;
					out.val[index]+= get_PG(in,i+1,j,k) / (hxp*(hxp+hxm));
					out.val[index]+= get_PG(in,i-1,j,k) / (hxm*(hxp+hxm));
					out.val[index]+= get_PG(in,i,j+1,k) / (hyp*(hyp+hym));
					out.val[index]+= get_PG(in,i,j-1,k) / (hym*(hyp+hym));
					out.val[index]+= get_PG(in,i,j,k+1) / (hzp*(hzp+hzm));
					out.val[index]+= get_PG(in,i,j,k-1) / (hzm*(hzp+hzm));

					out.val[index]+= 2*(rho_ges(i,j,k) + 1. - exp(in.val[index] ));



					double lhs = 1./(hxp*hxm) + 1./(hyp*hym) + 1./(hzp*hzm);
					out.val[index] = out.val[index]/lhs;

				} // END if (mask)
				else
				{ // if within boundary layer : keep boundary value
					out.val[index] = in.val[index];
				}

			} // enf loop over k
		} // end loop over j
	} // end loop over i

	return;
}


//ToDo : Still got to do this
#ifdef _MY_VERBOSE
void solver_poisson_jacobi_lin::check_convergence(const field_real &field_new, const field_real &field_old)
{
	double norm_maximum = 0.;
	double norm_sum_total = 0.;
	double norm_sum_i = 0.;
	double norm_sum_j = 0.;
	double norm_sum_k = 0.;

	/*
	double YZ[field_new.Ny*field_new.Nz];
	double XY[field_new.Nx*field_new.Ny];
	double XZ[field_new.Nx*field_new.Nz];
*/

	double * YZ = new double[field_new.Ny*field_new.Nz];
	double * XY = new double[field_new.Nx*field_new.Ny];
	double * XZ = new double[field_new.Nx*field_new.Nz];

	for(int i=0; i< field_new.Ny*field_new.Nz; ++i)
		YZ[i] = 0.;

	for(int i=0; i< field_new.Nx*field_new.Ny; ++i)
		XY[i] = 0.;

	for(int i=0; i< field_new.Nx*field_new.Nz; ++i)
		XZ[i] = 0.;

	for(int i=0; i<field_new.Nx; ++i)
		for(int j=0; j<field_new.Ny; ++j)
			for(int k=0; k<field_new.Nz; ++k)
			{
				double delta = fabs(field_new(i,j,k) - field_old(i,j,k) );
				norm_maximum = max<double>(norm_maximum,delta);
				norm_sum_total += delta;

				YZ[k + j*field_new.Nz] += delta/field_new.Nx;
				XY[j + i*field_new.Ny] += delta/field_new.Nz;
				XZ[k + i*field_new.Nz] += delta/field_new.Ny;
			}

	for(int i=0; i< field_new.Ny*field_new.Nz; ++i)
		norm_sum_i = max<double>(norm_sum_i,YZ[i]);

	for(int i=0; i< field_new.Nx*field_new.Ny; ++i)
		norm_sum_j = max<double>(norm_sum_j,XY[i]);

	for(int i=0; i< field_new.Nx*field_new.Nz; ++i)
		norm_sum_k = max<double>(norm_sum_k,XZ[i]);


	std::stringstream status;

	status << invocations << "\t";
	status << iteration << "\t";
	status << norm_maximum << "\t";
	status << norm_sum_total << "\t";
	status << norm_sum_i << "\t";
	status << norm_sum_j << "\t";
	status << norm_sum_k << "\n";

	std::ofstream output_stream("./diagnostics/relaxation.log", std::ofstream::app);
	output_stream << status.str();
	output_stream.close();

	delete[] YZ;
	delete[] XY;
	delete[] XZ;

	return;
}

#endif


