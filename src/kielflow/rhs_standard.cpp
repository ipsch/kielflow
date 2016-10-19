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
Buffer_FUx(domain), Buffer_FUy(domain), Buffer_FUz(domain), Buffer_Fni(domain),
Ux(domain), Uy(domain), Uz(domain),
dFUx_dx(domain), dFUy_dx(domain), dFUz_dx(domain),
dFUx_dy(domain), dFUy_dy(domain), dFUz_dy(domain),
dFUx_dz(domain), dFUy_dz(domain), dFUz_dz(domain),
dealiasing(domain),
FBuffer_1st(domain),
FBuffer_2nd(domain),
Buffer_1st(domain),
Buffer_2nd(domain),
N(domain.Nx*domain.Ny*(domain.Nz/2+1))



{
	std::cout << "done pre init\n";
	bool SV_warning_upper = false;
	bool SV_warning_lower = true;


	field_SV = (double*) fftw_malloc(sizeof(double) * N);
	double kx_max = domain.x_axis->k_val_at(domain.Nx/2);
	double ky_max = domain.y_axis->k_val_at(domain.Ny/2);
	double kz_max = domain.z_axis->k_val_at(domain.Nz/2);
	double kpow2 = kx_max*kx_max + ky_max*ky_max + kz_max*kz_max;
	//double eps = .9/kpow2;
	double eps = 0.03;

	std::cout << "epsilon in spectral viscosity set to!!\n";
	std::cout << "eps= " << eps << "\n";

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
				double filter = 1. - exp(-36.*pow(sqrt( knorm ),26.));
				//double filter = exp(-36.*pow(sqrt( knorm ),26.));
				field_SV[ijk] = eps*kpow2*filter;

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
	dealiasing_theta(FBuffer_1st, 0.333);
	iFFT(FBuffer_1st,Buffer_1st);

	for(int ijk=0; ijk<Buffer_1st.N; ++ijk)
	{
		Buffer_1st.val[ijk] = exp(Buffer_1st.val[ijk])+nd.val[ijk];
	}
	my_poisson_solver.solve(Phi,Buffer_1st);
   #endif





	// First derive all all linear terms where no dealiasing is needed


   #if defined(_RHS_DISSIPATION) && defined(_RHS_SVISCOSITY) // ##################
	for(int ijk=0; ijk<Buffer_FUx.N; ++ijk)
	{
		Buffer_FUx.val[ijk][0] = (my_params.tau+field_SV[ijk])*FUx.val[ijk][0];
		Buffer_FUx.val[ijk][1] = (my_params.tau+field_SV[ijk])*FUx.val[ijk][1];
	}
	for(int ijk=0; ijk<Buffer_FUy.N; ++ijk)
	{
		Buffer_FUy.val[ijk][0] = (my_params.tau+field_SV[ijk])*FUy.val[ijk][0];
		Buffer_FUy.val[ijk][1] = (my_params.tau+field_SV[ijk])*FUy.val[ijk][1];
	}
	for(int ijk=0; ijk<Buffer_FUz.N; ++ijk)
	{
		Buffer_FUz.val[ijk][0] = (my_params.tau+field_SV[ijk])*FUz.val[ijk][0];
		Buffer_FUz.val[ijk][1] = (my_params.tau+field_SV[ijk])*FUz.val[ijk][1];
	}
   #endif

#if defined(_RHS_DISSIPATION) && !defined(_RHS_SVISCOSITY) // ##############
	for(int ijk=0; ijk<FUx.N; ++ijk)
	{
		Buffer_FUx.val[ijk][0] = my_params.tau*FUx.val[ijk][0];
		Buffer_FUx.val[ijk][1] = my_params.tau*FUx.val[ijk][1];
	}
	for(int ijk=0; ijk<FUy.N; ++ijk)
	{
		Buffer_FUy.val[ijk][0] = my_params.tau*FUy.val[ijk][0];
		Buffer_FUy.val[ijk][1] = my_params.tau*FUy.val[ijk][1];
	}
	for(int ijk=0; ijk<FUz.N; ++ijk)
	{
		Buffer_FUz.val[ijk][0] = my_params.tau*FUz.val[ijk][0];
		Buffer_FUz.val[ijk][1] = my_params.tau*FUz.val[ijk][1];
	}
#endif

#if !defined(_RHS_DISSIPATION) && defined(_RHS_SVISCOSITY) // ##################
	for(int ijk=0; ijk<Buffer_FUx.N; ++ijk)
	{
		Buffer_FUx.val[ijk][0] = field_SV[ijk]*FUx.val[ijk][0];
		Buffer_FUx.val[ijk][1] = field_SV[ijk]*FUx.val[ijk][1];
	}


	for(int ijk=0; ijk<Buffer_FUy.N; ++ijk)
	{
			Buffer_FUy.val[ijk][0] = field_SV[ijk]*FUy.val[ijk][0];
			Buffer_FUy.val[ijk][1] = field_SV[ijk]*FUy.val[ijk][1];
	}
	for(int ijk=0; ijk<Buffer_FUz.N; ++ijk)
	{
		Buffer_FUz.val[ijk][0] = field_SV[ijk]*FUz.val[ijk][0];
		Buffer_FUz.val[ijk][1] = field_SV[ijk]*FUz.val[ijk][1];
	}
#endif

#if !defined(_RHS_DISSIPATION) && !defined(_RHS_SVISCOSITY) // #################
	for(int ijk=0; ijk<Buffer_FUx.N; ++ijk)
	{
		Buffer_FUx.val[ijk][0] = 0.;
		Buffer_FUx.val[ijk][1] = 0.;
	}
	for(int ijk=0; ijk<Buffer_FUy.N; ++ijk)
	{
		Buffer_FUy.val[ijk][0] = 0.;
		Buffer_FUy.val[ijk][1] = 0.;
	}
	for(int ijk=0; ijk<Buffer_FUz.N; ++ijk)
	{
		Buffer_FUz.val[ijk][0] = 0.;
		Buffer_FUz.val[ijk][1] = 0.;
	}
   #endif
    // ############# END DISSIPATION ##########################################

	// dealiasing of input prior to evaluation of nonlinear terms

	//dealiasing(FUx);
	//dealiasing(FUy);
	//dealiasing(FUz);
	//dealiasing(Fni);

	dealiasing_theta(FUx, 0.333);
	dealiasing_theta(FUy, 0.333);
	dealiasing_theta(FUz, 0.333);
	dealiasing_theta(Fni, 0.333);

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
	for(int i=0; i<Buffer_FUx.N; ++i)
	{
		Buffer_FUx.val[i][0] += my_params.M*dFUx_dx.val[i][0];
		Buffer_FUx.val[i][1] += my_params.M*dFUx_dx.val[i][1];
	}
	for(int i=0; i<Buffer_FUy.N; ++i)
	{
		Buffer_FUy.val[i][0] += my_params.M*dFUy_dx.val[i][0];
		Buffer_FUy.val[i][1] += my_params.M*dFUy_dx.val[i][1];
	}
	for(int i=0; i<Buffer_FUz.N; ++i)
	{
		Buffer_FUz.val[i][0] += my_params.M*dFUz_dx.val[i][0];
		Buffer_FUz.val[i][1] += my_params.M*dFUz_dx.val[i][1];
	}
   #endif

	// ########################################################################
   #if defined(_RHS_DIFFUSION) && defined(_RHS_E_INT) // ######################
	my_FFT(Phi,FBuffer_1st);
	dealiasing_theta(FBuffer_1st, 0.333333); //
	for(int ijk=0; ijk<FBuffer_1st.N; ++ijk)
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

	//dealiasing(FBuffer_1st);
	d_dx(FBuffer_1st,FBuffer_2nd);
	Buffer_FUx += FBuffer_2nd;
	d_dy(FBuffer_1st,FBuffer_2nd);
	Buffer_FUy += FBuffer_2nd;
	d_dz(FBuffer_1st,FBuffer_2nd);
	Buffer_FUz += FBuffer_2nd;
   #endif





	// ########################################################################
	// ################## CONTINUITY EQUATION #################################
	// ########################################################################



   #if defined(_RHS_E_EXT) || defined(_RHS_CONTINUITY) // shared by both
	d_dx(Fni,FBuffer_1st);
   #endif

	// ########################################################################
   #if defined(_RHS_E_EXT) // #################################################
	for(int i=0; i<Buffer_Fni.N; ++i)
	{
		Buffer_Fni.val[i][0] = my_params.M*FBuffer_1st.val[i][0];
		Buffer_Fni.val[i][1] = my_params.M*FBuffer_1st.val[i][1];
	}
   #else
	for(int i=0; i<Buffer_Fni.N; ++i)
	{
		Buffer_Fni.val[i][0] = 0.;
		Buffer_Fni.val[i][1] = 0.;
	}
   #endif


	// ########################################################################
   #if defined(_RHS_CONTINUITY)
	// LaTeX: n = div(\mathbf{u}) + [..]
	Buffer_Fni += dFUx_dx;
	Buffer_Fni += dFUy_dy;
	Buffer_Fni += dFUz_dz;

	// LaTeX: \dot{n} = \left(\frac{\partial n}{\partial x}\right) u_x
	//d_dx(Fni,FBuffer_1st); // see above
	dealiasing(FBuffer_1st);
	my_iFFT(FBuffer_1st,Buffer_1st);                    // (dn/dx)
	for(int i=0; i<Buffer_2nd.N; ++i)				 // Ux*(dn/dx)
		Buffer_2nd.val[i] = Buffer_1st.val[i]*Ux.val[i];

	// LaTeX: \dot{n} = \left(\frac{\partial n}{\partial y}\right) u_y
	d_dy(Fni,FBuffer_1st);
	dealiasing(FBuffer_1st);
	my_iFFT(FBuffer_1st,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i] += Buffer_1st.val[i]*Uy.val[i];

	// LaTeX: \dot{n} = \left(\frac{\partial n}{\partial z}\right) u_z
	d_dz(Fni,FBuffer_1st);
	dealiasing(FBuffer_1st);
	my_iFFT(FBuffer_1st,Buffer_1st);
	for(int i=0; i<Buffer_2nd.N; ++i)
		Buffer_2nd.val[i] += Buffer_1st.val[i]*Uz.val[i];

	my_FFT(Buffer_2nd,FBuffer_1st);
	Buffer_Fni += FBuffer_1st;

   #endif

	// ########################################################################
   #if defined(_RHS_ADVECTION) // #############################################

	dealiasing(dFUx_dx);
	dealiasing(dFUx_dy);
	dealiasing(dFUx_dz);

	dealiasing(dFUy_dx);
	dealiasing(dFUy_dy);
	dealiasing(dFUy_dz);

	dealiasing(dFUz_dx);
	dealiasing(dFUz_dy);
	dealiasing(dFUz_dz);

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
	Buffer_FUx += FBuffer_1st;


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
	Buffer_FUy += FBuffer_1st;


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
	Buffer_FUz += FBuffer_1st;

   #endif





	// ########################################################################
   #if defined(_RHS_PENALIZATION_U) // ########################################
	static effect_penalization UX_penalization(-.5, my_solids, 0.);
	static effect_penalization ni_penalization(-.5, my_solids, 1.);
	UX_penalization.execute(FUx,Buffer_FUx);
	UX_penalization.execute(FUy,Buffer_FUy);
	UX_penalization.execute(FUz,Buffer_FUz);
	ni_penalization.execute(Fni,Buffer_Fni);
   #endif

	// ########################################################################
	// RHS-VALUE BACK TO I/O ##################################################
	for(int ijk=0; ijk<FUx.N; ++ijk)
	{
		FUx.val[ijk][0] = -Buffer_FUx.val[ijk][0];
		FUx.val[ijk][1] = -Buffer_FUx.val[ijk][1];
	}

	for(int ijk=0; ijk<FUy.N; ++ijk)
	{
		FUy.val[ijk][0] = -Buffer_FUy.val[ijk][0];
		FUy.val[ijk][1] = -Buffer_FUy.val[ijk][1];
	}

	for(int ijk=0; ijk<FUz.N; ++ijk)
	{
		FUz.val[ijk][0] = -Buffer_FUz.val[ijk][0];
		FUz.val[ijk][1] = -Buffer_FUz.val[ijk][1];
	}

	dealiasing(Buffer_Fni);
	for(int ijk=0; ijk<Fni.N; ++ijk)
	{
		Fni.val[ijk][0] = -Buffer_Fni.val[ijk][0];
		Fni.val[ijk][1] = -Buffer_Fni.val[ijk][1];
	}

   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
	return;
}
