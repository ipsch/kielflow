#include "rhs_standard.hpp"





rhs_standard::rhs_standard(parameters &Params, interface_poisson_solver &Solver,
		field_real &potential, field_real &density_dust, field_real &solid_mask, const grid &domain) :
my_params(Params),
my_poisson_solver(Solver),
Phi(potential),
my_solids(solid_mask),
d_dx(domain, e_x), d_dy(domain, e_y), d_dz(domain, e_z),
my_FFT(domain), my_iFFT(domain),
// Buffers
nd(density_dust),
FBuffer_Ux(domain), FBuffer_Uy(domain), FBuffer_Uz(domain), FBuffer_ni(domain),
Ux(domain), Uy(domain), Uz(domain),
dFUx_dx(domain), dFUy_dx(domain), dFUz_dx(domain),
dFUx_dy(domain), dFUy_dy(domain), dFUz_dy(domain),
dFUx_dz(domain), dFUy_dz(domain), dFUz_dz(domain),
my_dealiasing(domain, [] (const double &k) {return k<(2./3.) ? 1. : 0.;}),
FBuffer_1st(domain),
FBuffer_2nd(domain),
Buffer_1st(domain),
Buffer_2nd(domain),
N(domain.Nx*domain.Ny*(domain.Nz/2+1))



{
	std::cout << "done pre init\n";
	bool SV_warning_upper = false;
	bool SV_warning_lower = true;


	double N_real = domain.Nx*domain.Ny*domain.Nz;

	penalization_barrier = (double*) fftw_malloc(sizeof(double) * N_real);
	for(int i=0; i<domain.Nx; i++)
		for(int j=0; j<domain.Ny; j++)
			for(int k=0; k<domain.Nz; k++)
			{
				int ijk = (k + (domain.Nz)*(j + i*domain.Ny));
				double x=domain.x_axis->val_at(i);
				double L=domain.x_axis->L;
				const double pi=acos(-1);
				penalization_barrier[ijk] = (x<-0.5*L/3.) ? pow(cos(2*(1.5)*pi/L*x),2.) : 0.;
			}

	field_SV = (double*) fftw_malloc(sizeof(double) * N);
	double kx_max = domain.x_axis->k_val_at(domain.Nx/2);
	double ky_max = domain.y_axis->k_val_at(domain.Ny/2);
	double kz_max = domain.z_axis->k_val_at(domain.Nz/2);
	double kpow2 = kx_max*kx_max + ky_max*ky_max + kz_max*kz_max;
	//double eps = .9/kpow2;
	double eps = .5;
	for(int i=0; i<domain.Nx; ++i)
		for(int j=0; j<domain.Ny; ++j)
			for(int k=0; k<domain.Nz/2+1; ++k)
			{

				int ijk = (k + (domain.Nz/2+1)*(j + i*domain.Ny));
				double kx = domain.x_axis->k_val_at(i);
				double ky = domain.y_axis->k_val_at(j);
				double kz = domain.z_axis->k_val_at(k);
				kpow2 = kx*kx + ky*ky + kz*kz;
				double knorm = sqrt(pow(kx/kx_max,2.) + pow(ky/ky_max,2.) + pow(kz/kz_max,2.));
				double filter = 1. - exp(-8.*pow(sqrt( knorm ),8.));
				//double filter = exp(-36.*pow(sqrt( knorm ),26.));
				field_SV[ijk] = eps*filter;

				if(fabs(field_SV[ijk])>1.)
					SV_warning_upper=true;
				if(fabs(field_SV[ijk])>my_params.tau)
					SV_warning_lower=false;
			}

	if(SV_warning_upper)
		std::cout << "WARNING: Spectral viscosity might be too strong !!\n";
	if(SV_warning_lower)
		std::cout << "WARNING: Spectral viscosity might be too weak !!\n";


}



void rhs_standard::solve(const double &t, field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni)
// The Right-Hand-Side of the PDE is evaluated.
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
    logger my_log("rhs_standard::solve(..)");
	my_log << "start";
   #endif


	// ########################################################################
	// ################## POISSON EQUATION ####################################
   #if defined(_RHS_POISSON) && defined(_RHS_E_INT)
	for(int i=0; i<FBuffer_1st.N; ++i)
	{
		FBuffer_1st.val[i][0] = Fni.val[i][0];
		FBuffer_1st.val[i][1] = Fni.val[i][1];
	}
	my_dealiasing(FBuffer_1st);
	my_iFFT(FBuffer_1st,Buffer_1st);
	for(int ijk=0; ijk<Buffer_1st.N; ++ijk)
	{
		Buffer_1st.val[ijk] = exp(Buffer_1st.val[ijk])+nd.val[ijk];
	}
	my_poisson_solver.solve(Phi,Buffer_1st);
   #endif

	// Linear Terms first:
    // ###########################################################################
	// ############ Dissipation & spectral viscosity #############################
	// ###########################################################################
   #if defined(_RHS_DISSIPATION) && defined(_RHS_SVISCOSITY) // ##################
	for(int ijk=0; ijk<FBuffer_Ux.N; ++ijk)
	{
		FBuffer_Ux.val[ijk][0] = (my_params.tau+field_SV[ijk])*FUx.val[ijk][0];
		FBuffer_Ux.val[ijk][1] = (my_params.tau+field_SV[ijk])*FUx.val[ijk][1];
	}
	for(int ijk=0; ijk<FBuffer_Uy.N; ++ijk)
	{
		FBuffer_Uy.val[ijk][0] = (my_params.tau+field_SV[ijk])*FUy.val[ijk][0];
		FBuffer_Uy.val[ijk][1] = (my_params.tau+field_SV[ijk])*FUy.val[ijk][1];
	}
	for(int ijk=0; ijk<FBuffer_Uz.N; ++ijk)
	{
		FBuffer_Uz.val[ijk][0] = (my_params.tau+field_SV[ijk])*FUz.val[ijk][0];
		FBuffer_Uz.val[ijk][1] = (my_params.tau+field_SV[ijk])*FUz.val[ijk][1];
	}
	for(int ijk=0; ijk<FBuffer_ni.N; ++ijk)
	{
		FBuffer_ni.val[ijk][0] = field_SV[ijk]*Fni.val[ijk][0];
		FBuffer_ni.val[ijk][1] = field_SV[ijk]*Fni.val[ijk][1];
	}
   #endif

#if defined(_RHS_DISSIPATION) && !defined(_RHS_SVISCOSITY) // ##############
	for(int ijk=0; ijk<FUx.N; ++ijk)
	{
		FBuffer_Ux.val[ijk][0] = my_params.tau*FUx.val[ijk][0];
		FBuffer_Ux.val[ijk][1] = my_params.tau*FUx.val[ijk][1];
	}
	for(int ijk=0; ijk<FUy.N; ++ijk)
	{
		FBuffer_Uy.val[ijk][0] = my_params.tau*FUy.val[ijk][0];
		FBuffer_Uy.val[ijk][1] = my_params.tau*FUy.val[ijk][1];
	}
	for(int ijk=0; ijk<FUz.N; ++ijk)
	{
		FBuffer_Uz.val[ijk][0] = my_params.tau*FUz.val[ijk][0];
		FBuffer_Uz.val[ijk][1] = my_params.tau*FUz.val[ijk][1];
	}
	for(int ijk=0; ijk<FBuffer_ni.N; ++ijk)
	{
		FBuffer_ni.val[ijk][0] = 0.;
		FBuffer_ni.val[ijk][1] = 0.;
	}
#endif

#if !defined(_RHS_DISSIPATION) && defined(_RHS_SVISCOSITY) // ##################
	for(int ijk=0; ijk<FBuffer_Ux.N; ++ijk)
	{
		FBuffer_Ux.val[ijk][0] = field_SV[ijk]*FUx.val[ijk][0];
		FBuffer_Ux.val[ijk][1] = field_SV[ijk]*FUx.val[ijk][1];
	}
	for(int ijk=0; ijk<FBuffer_Uy.N; ++ijk)
	{
			FBuffer_Uy.val[ijk][0] = field_SV[ijk]*FUy.val[ijk][0];
			FBuffer_Uy.val[ijk][1] = field_SV[ijk]*FUy.val[ijk][1];
	}
	for(int ijk=0; ijk<FBuffer_Uz.N; ++ijk)
	{
		FBuffer_Uz.val[ijk][0] = field_SV[ijk]*FUz.val[ijk][0];
		FBuffer_Uz.val[ijk][1] = field_SV[ijk]*FUz.val[ijk][1];
	}
	for(int ijk=0; ijk<FBuffer_ni.N; ++ijk)
	{
		FBuffer_ni.val[ijk][0] = field_SV[ijk]*Fni.val[ijk][0];
		FBuffer_ni.val[ijk][1] = field_SV[ijk]*Fni.val[ijk][1];
	}
#endif

#if !defined(_RHS_DISSIPATION) && !defined(_RHS_SVISCOSITY) // #################
	for(int ijk=0; ijk<FBuffer_Ux.N; ++ijk)
	{
		FBuffer_Ux.val[ijk][0] = 0.;
		FBuffer_Ux.val[ijk][1] = 0.;
	}
	for(int ijk=0; ijk<FBuffer_Uy.N; ++ijk)
	{
		FBuffer_Uy.val[ijk][0] = 0.;
		FBuffer_Uy.val[ijk][1] = 0.;
	}
	for(int ijk=0; ijk<FBuffer_Uz.N; ++ijk)
	{
		FBuffer_Uz.val[ijk][0] = 0.;
		FBuffer_Uz.val[ijk][1] = 0.;
	}
	for(int ijk=0; ijk<FBuffer_ni.N; ++ijk)
	{
		FBuffer_ni.val[ijk][0] = 0.;
		FBuffer_ni.val[ijk][1] = 0.;
	}
   #endif
    // ############# END DISSIPATION ##########################################



	// another linear term (note linear because ni rescaled by exp(ni)->ni )
	// ########################################################################
	// ################### Internal E-forces ##################################
	// ########################################################################
   #if defined(_RHS_DIFFUSION) && defined(_RHS_E_INT) // ######################
	my_FFT(Phi,FBuffer_1st);
	dealiasing_undesignated(FBuffer_1st, [] (double k) {return exp(-2000.*pow(k,15.));});
	for(int ijk=0; ijk<FBuffer_1st.N; ijk++)
	{
		FBuffer_1st.val[ijk][0] += Fni.val[ijk][0]/my_params.theta;
		FBuffer_1st.val[ijk][1] += Fni.val[ijk][1]/my_params.theta;
	}
   #endif

   #if defined(_RHS_DIFFUSION) && !defined(_RHS_E_INT) // #####################
	for(int i=0; i<FBuffer_1st.N; ++i)
	{
		FBuffer_1st.val[i][0] = Fni.val[i][0]/my_params.theta;
		FBuffer_1st.val[i][1] = Fni.val[i][1]/my_params.theta;
	}
   #endif

   #if !defined(_RHS_DIFFUSION) && defined(_RHS_E_INT) // #####################
	my_FFT(Phi,FBuffer_1st);
   #endif

   #if defined(_RHS_DIFFUSION) || defined(_RHS_E_INT) // #####################
	//my_dealiasing(FBuffer_1st);
	d_dx(FBuffer_1st,FBuffer_2nd);
	FBuffer_Ux += FBuffer_2nd;
	d_dy(FBuffer_1st,FBuffer_2nd);
	FBuffer_Uy += FBuffer_2nd;
	d_dz(FBuffer_1st,FBuffer_2nd);
	FBuffer_Uz += FBuffer_2nd;
   #endif


	// dealiasing of input prior to evaluation of nonlinear terms

	my_dealiasing(FUx);
	my_dealiasing(FUy);
	my_dealiasing(FUz);
	my_dealiasing(Fni);

	my_iFFT(FUx,Ux);
	my_iFFT(FUy,Uy);
	my_iFFT(FUz,Uz);

	// all derivatives needed for RHS evaluation
	d_dx(FUx,dFUx_dx);
	d_dy(FUx,dFUx_dy);
	d_dz(FUx,dFUx_dz);

	d_dx(FUy,dFUy_dx);
	d_dy(FUy,dFUy_dy);
	d_dz(FUy,dFUy_dz);

	d_dx(FUz,dFUz_dx);
	d_dy(FUz,dFUz_dy);
	d_dz(FUz,dFUz_dz);




	// ########################################################################
   #if defined(_RHS_E_EXT) // #################################################
	for(int i=0; i<FBuffer_Ux.N; ++i)
	{
		FBuffer_Ux.val[i][0] += my_params.M*dFUx_dx.val[i][0];
		FBuffer_Ux.val[i][1] += my_params.M*dFUx_dx.val[i][1];
	}
	for(int i=0; i<FBuffer_Uy.N; ++i)
	{
		FBuffer_Uy.val[i][0] += my_params.M*dFUy_dx.val[i][0];
		FBuffer_Uy.val[i][1] += my_params.M*dFUy_dx.val[i][1];
	}
	for(int i=0; i<FBuffer_Uz.N; ++i)
	{
		FBuffer_Uz.val[i][0] += my_params.M*dFUz_dx.val[i][0];
		FBuffer_Uz.val[i][1] += my_params.M*dFUz_dx.val[i][1];
	}
   #endif


	// ########################################################################
	// ################## CONTINUITY EQUATION #################################
	// ########################################################################


	// ########################################################################
   #if defined(_RHS_E_EXT) // #################################################
	d_dx(Fni,FBuffer_1st);
	for(int i=0; i<FBuffer_ni.N; ++i)
	{
		FBuffer_ni.val[i][0] += my_params.M*FBuffer_1st.val[i][0];
		FBuffer_ni.val[i][1] += my_params.M*FBuffer_1st.val[i][1];
	}
   #endif


	// ########################################################################
   #if defined(_RHS_CONTINUITY)
	// LaTeX: n = div(\mathbf{u}) + [..]
	FBuffer_ni += dFUx_dx;
	FBuffer_ni += dFUy_dy;
	FBuffer_ni += dFUz_dz;

	// LaTeX: \dot{n} = \left(\frac{\partial n}{\partial x}\right) u_x
	//d_dx(Fni,FBuffer_1st); // see above
	my_dealiasing(FBuffer_1st);
	my_iFFT(FBuffer_1st,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)				 // Ux*(dn/dx)
		Buffer_2nd.val[i] = Buffer_1st.val[i]*Ux.val[i];

	// LaTeX: \dot{n} = \left(\frac{\partial n}{\partial y}\right) u_y
	d_dy(Fni,FBuffer_1st);
	my_dealiasing(FBuffer_1st);
	my_iFFT(FBuffer_1st,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i] += Buffer_1st.val[i]*Uy.val[i];

	// LaTeX: \dot{n} = \left(\frac{\partial n}{\partial z}\right) u_z
	d_dz(Fni,FBuffer_1st);
	my_dealiasing(FBuffer_1st);
	my_iFFT(FBuffer_1st,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i] += Buffer_1st.val[i]*Uz.val[i];

	my_FFT(Buffer_2nd,FBuffer_1st);
	FBuffer_ni += FBuffer_1st;

   #endif



/*
#if defined(_RHS_CONTINUITY)
	my_iFFT(Fni,Buffer_1st);
	for(int i=0; i<Buffer_1st.N; ++i)
		Buffer_1st.val[i] = exp(Buffer_1st.val[i]);

	for(int i=0; i<Buffer_2nd.N; ++i)				 // Ux*(dn/dx)
		Buffer_2nd.val[i] = Buffer_1st.val[i]*(M+Ux.val[i]);
	my_FFT(Buffer_2nd,FBuffer_1st);
	d_dx(FBuffer_1st,FBuffer_1st);
	FBuffer_ni = FBuffer_1st;


	for(int i=0; i<Buffer_2nd.N; ++i)				 // Ux*(dn/dx)
		Buffer_2nd.val[i] = Buffer_1st.val[i]*Uy.val[i];
	my_FFT(Buffer_2nd,FBuffer_1st);
	d_dy(FBuffer_1st,FBuffer_1st);
	FBuffer_ni += FBuffer_1st;

	for(int i=0; i<Buffer_2nd.N; ++i)				 // Ux*(dn/dx)
		Buffer_2nd.val[i] = Buffer_1st.val[i]*Uz.val[i];
	my_FFT(Buffer_2nd,FBuffer_1st);
	d_dz(FBuffer_1st,FBuffer_1st);
	FBuffer_ni += FBuffer_1st;


	my_iFFT(FBuffer_ni,Buffer_2nd);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i] = Buffer_2nd.val[i]/Buffer_1st.val[i];
	my_FFT(Buffer_2nd,FBuffer_ni);

#endif
*/





	// ########################################################################
   #if defined(_RHS_ADVECTION) // #############################################

	my_dealiasing(dFUx_dx);
	my_dealiasing(dFUx_dy);
	my_dealiasing(dFUx_dz);

	my_dealiasing(dFUy_dx);
	my_dealiasing(dFUy_dy);
	my_dealiasing(dFUy_dz);

	my_dealiasing(dFUz_dx);
	my_dealiasing(dFUz_dy);
	my_dealiasing(dFUz_dz);


	my_iFFT(dFUx_dx,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i]  = Ux.val[i]*Buffer_1st.val[i];
	my_iFFT(dFUx_dy,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i] += Uy.val[i]*Buffer_1st.val[i];
	my_iFFT(dFUx_dz,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i] += Uz.val[i]*Buffer_1st.val[i];
	my_FFT(Buffer_2nd, FBuffer_1st);
	FBuffer_Ux += FBuffer_1st;

	my_iFFT(dFUy_dx,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i]  = Ux.val[i]*Buffer_1st.val[i];
	my_iFFT(dFUy_dy,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i] += Uy.val[i]*Buffer_1st.val[i];
	my_iFFT(dFUy_dz,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i] += Uz.val[i]*Buffer_1st.val[i];
	my_FFT(Buffer_2nd, FBuffer_1st);
	FBuffer_Uy += FBuffer_1st;

	my_iFFT(dFUz_dx,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i]  = Ux.val[i]*Buffer_1st.val[i];
	my_iFFT(dFUz_dy,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i] += Uy.val[i]*Buffer_1st.val[i];
	my_iFFT(dFUz_dz,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i] += Uz.val[i]*Buffer_1st.val[i];
	my_FFT(Buffer_2nd, FBuffer_1st);
	FBuffer_Uz += FBuffer_1st;

   #endif





	// ########################################################################
   #if defined(_RHS_PENALIZATION_U) // ########################################
	static effect_penalization UX_penalization(-.5, my_solids, 0.);
	static effect_penalization ni_penalization(-.5, my_solids, 1.);
	UX_penalization.execute(FUx,FBuffer_Ux);
	UX_penalization.execute(FUy,FBuffer_Uy);
	UX_penalization.execute(FUz,FBuffer_Uz);
	ni_penalization.execute(Fni,FBuffer_ni);
   #endif



	// ########################################################################
	// PENALIZATION ###########################################################
   #if defined(_RHS_PENALIZATION)
    // penalization is calculated for two regions:
	// fist in the proximity to the dust grain (my_soldis) and second
	// far away in upstream direction (penalization_arrier) to recover
	// the undisturbed solution from the disturbed solution after
	// fields (note: periodic boundaries) reenter the simulation box.
	// For the 4 treated fields the process far away from the dust grain is
	// same (exponentially damp field value to zero where penalization_barrier
	// is non zero).
	// The process in proximity to the grain is slightly different as
	// discribed below:
	// for dFUx_dt: damp field to -machnumber where my_solids is non zero
	// (because we are in a moving reference frame)
	// for dFUy_dt and dFUz_dt: damp field to zero where my_solids is non zero.
	// for dFni_dt: reduce field by my_solid multiplied by some constant value
	// (note: ni (ion density) is in logarithmic scaling this means:
	// an almost infinite negative value corresponds to an actual vacuum.

	double eps1 = 2.;
	double eps2 = 1.0;

	iFFT(FUx, Buffer_1st);
	iFFT(FBuffer_Ux,Buffer_2nd);
	double sup = 0.;
	for(int ijk=0; ijk<Buffer_1st.N; ijk++)
	{
		double tmp;
		tmp = eps1*my_solids.val[ijk]*(Buffer_1st.val[ijk]+my_params.M*my_solids.val[ijk]);
		tmp+= (1.-my_solids.val[ijk])*Buffer_2nd.val[ijk];
		tmp+= (penalization_barrier[ijk])*Buffer_1st.val[ijk];
		Buffer_1st.val[ijk] = tmp;
	}
	FFT(Buffer_1st,FBuffer_Ux);
	my_dealiasing(FBuffer_Ux);


	iFFT(FUy, Buffer_1st);
	iFFT(FBuffer_Uy,Buffer_2nd);
	for(int ijk=0; ijk<Buffer_1st.N; ijk++)
	{
		double tmp;
		tmp = eps1*my_solids.val[ijk]*(Buffer_1st.val[ijk]+0.);
		tmp+= (1.-my_solids.val[ijk])*Buffer_2nd.val[ijk];
		tmp+= (penalization_barrier[ijk])*Buffer_1st.val[ijk];
		Buffer_1st.val[ijk] = tmp;
	}
	FFT(Buffer_1st,FBuffer_Uy);
	my_dealiasing(FBuffer_Uy);

	iFFT(FUz, Buffer_1st);
	iFFT(FBuffer_Uz,Buffer_2nd);
	for(int ijk=0; ijk<Buffer_1st.N; ijk++)
	{
		double tmp;
		tmp = eps1*my_solids.val[ijk]*(Buffer_1st.val[ijk]+0.);
		tmp+= (1.-my_solids.val[ijk])*Buffer_2nd.val[ijk];
		tmp+= (penalization_barrier[ijk])*Buffer_1st.val[ijk];
		Buffer_1st.val[ijk] = tmp;
	}
	FFT(Buffer_1st,FBuffer_Uz);
	my_dealiasing(FBuffer_Uz);

	iFFT(Fni, Buffer_1st);
	iFFT(FBuffer_ni,Buffer_2nd);
	for(int ijk=0; ijk<Buffer_1st.N; ijk++)
	{
		double tmp;
		tmp = eps1*my_solids.val[ijk]*(Buffer_1st.val[ijk]+10.);
		tmp+= (1.-my_solids.val[ijk])*Buffer_2nd.val[ijk];
		tmp+= (penalization_barrier[ijk])*Buffer_1st.val[ijk];
		Buffer_1st.val[ijk] = tmp;
	}
	FFT(Buffer_1st,FBuffer_ni);
	my_dealiasing(FBuffer_ni);
#endif


	// ########################################################################
	// RHS-VALUE BACK TO I/O ##################################################
	for(int ijk=0; ijk<FUx.N; ++ijk)
	{
		FUx.val[ijk][0] = -FBuffer_Ux.val[ijk][0];
		FUx.val[ijk][1] = -FBuffer_Ux.val[ijk][1];
	}

	for(int ijk=0; ijk<FUy.N; ++ijk)
	{
		FUy.val[ijk][0] = -FBuffer_Uy.val[ijk][0];
		FUy.val[ijk][1] = -FBuffer_Uy.val[ijk][1];
	}

	for(int ijk=0; ijk<FUz.N; ++ijk)
	{
		FUz.val[ijk][0] = -FBuffer_Uz.val[ijk][0];
		FUz.val[ijk][1] = -FBuffer_Uz.val[ijk][1];
	}

	for(int ijk=0; ijk<Fni.N; ++ijk)
	{
		Fni.val[ijk][0] = -FBuffer_ni.val[ijk][0];
		Fni.val[ijk][1] = -FBuffer_ni.val[ijk][1];
	}

   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
	return;
}
