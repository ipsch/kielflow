#include "rhs_standard.hpp"

   #define POISSON
   #define ADVECTION
   #define DIFFUSION
   #define E_INT
   #define E_EXT
   #define DISSIPATION
   //#define PENALIZATION_U
   #define CONTINUITY
   //#define DEALAISING_MORE
   #define SPECTRAL_VISCOSITY


// ############ implementation of RHS_STANDARD #########################################



rhs_standard::rhs_standard(parameters &Params, interface_poisson_solver &Solver,
		field_real &potential, field_real &density_dust, field_real &solid_mask) :
my_params(Params),
my_poisson_solver(Solver),
Phi(potential),
nd(density_dust),
my_solids(solid_mask)
{
   #ifdef _MY_VERBOSE
	logger log("solver_basic_rhs");
	log << "solve(...)";
   #endif
}

void rhs_standard::solve(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni)
// The Right-Hand-Side of the PDE is evaluated.
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
    logger my_log("rhs_standard::solve(..)");
	my_log << "start";
   #endif

	static int N = FUx.N;
	static field_imag Buffer_FUx(*FUx.my_grid);
	static field_imag Buffer_FUy(*FUx.my_grid);
	static field_imag Buffer_FUz(*FUx.my_grid);
	static field_imag Buffer_Fni(*FUx.my_grid);
	static field_imag FPh(*FUx.my_grid);
	static field_real ni(*FUx.my_grid);



   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "RHS(): clear Buffer";
   #endif
	for(int i=0; i<N; ++i)
	{
		Buffer_FUx.val[i][0] = 0.;
		Buffer_FUx.val[i][1] = 0.;
		Buffer_FUy.val[i][0] = 0.;
		Buffer_FUy.val[i][1] = 0.;
		Buffer_FUz.val[i][0] = 0.;
		Buffer_FUz.val[i][1] = 0.;
		Buffer_Fni.val[i][0] = 0.;
		Buffer_Fni.val[i][1] = 0.;
	}


	// ################## POISSON EQUATION ####################################
   #if defined(POISSON)
	iFFT(Fni,ni);
	for(int i=0; i<ni.N; ++i)
		ni.val[i] = ni.val[i]+nd.val[i];
	my_poisson_solver.solve(Phi,ni);
	FFT(Phi,FPh);
   #endif


	// ################## IMPULS EQUATION #####################################
   #if defined(DEALAISING_MORE)
	handle_dealiasing(FUx);
	handle_dealiasing(FUy);
	handle_dealiasing(FUz);
	handle_dealiasing(Fni);
   #endif

   #if defined(ADVECTION) // ##########################
	static effect_advection advection;
	advection.execute(FUx, FUy, FUz, FUx, Buffer_FUx);
	advection.execute(FUx, FUy, FUz, FUy, Buffer_FUy);
	advection.execute(FUx, FUy, FUz, FUz, Buffer_FUz);
   #endif

   #if defined(E_EXT) // ##########################
	static effect_translation translation_Ex(my_params.M, e_x);
	// translation_Ex is also applied to ion density
	translation_Ex.execute(FUx, Buffer_FUx);
	translation_Ex.execute(FUy, Buffer_FUy);
	translation_Ex.execute(FUz, Buffer_FUz);
   #endif

   #if defined(DIFFUSION) || defined(E_INT) // ##########################
	static effect_force_E Div;
   #if defined(DIFFUSION)
	iFFT(Fni,ni);
	for(int i=0; i<ni.N; ++i)
	{
		ni.val[i] = log(ni.val[i])/my_params.theta;
		if(fabs(ni.val[i]) < 1.e-10)
			ni.val[i] = 0.;
	}
   #else
	for(int i=0; i<ni.N; ++i)
	{
		ni.val[i] = 0;
	}
   #endif
   #if defined(E_INT)
	for(int i=0; i<ni.N; ++i)
	{
		ni.val[i] += Phi.val[i];
	}
   #endif
	field_imag FBuffer(*ni.my_grid);
	FFT(ni,FBuffer);
	Div.execute(FBuffer, Buffer_FUx, Buffer_FUy, Buffer_FUz);
   #endif

   #if defined(DISSIPATION) // ##########################
	static effect_force_linear dissipation(-my_params.tau);
	dissipation.execute(FUx, Buffer_FUx);
	dissipation.execute(FUy, Buffer_FUy);
	dissipation.execute(FUz, Buffer_FUz);
   #endif

	// ToDo: implement viscosity if needed
   #if defined(VISCOSITY)
	static effect_diffusion_u viscosity(my_params.theta);
	// not implemented yet
   #endif

   #if defined(SPECTRAL_VISCOSITY)
	static effect_spectral_viscosity VS(0.66,FUx.my_grid->x_axis->L,FUx.Nx);
	VS.execute(FUx,Buffer_FUx);
	VS.execute(FUy,Buffer_FUy);
	VS.execute(FUz,Buffer_FUz);
   #endif

   #if defined(PENALIZATION_U)
	static effect_penalization Ux_penalization(0.1, my_solids, -.5);
	static effect_penalization Ur_penalization(0.01, my_solids,  0. );
	Ux_penalization.execute(FUx,Buffer_FUx);
	Ur_penalization.execute(FUy,Buffer_FUy);
	Ur_penalization.execute(FUz,Buffer_FUz);
   #endif



	// ################## CONTINUITY EQUATION #################################
   #if defined(CONTINUITY)
	static effect_continuityEQ continuityEQ;
	continuityEQ.execute(FUx, FUy, FUz, Fni, Buffer_Fni);
   #endif

   #if defined(E_EXT)
	translation_Ex.execute(Fni, Buffer_Fni);
   #endif

	// ToDo : penalization n
	// This has to be the last step in evaluation of the rhs
	//solid.penalize_density(Fni,Fni);

   #if defined(DEALAISING_MORE)
	handle_dealiasing(Buffer_FUx);
	handle_dealiasing(Buffer_FUy);
	handle_dealiasing(Buffer_FUz);
	handle_dealiasing(Buffer_Fni);
   #endif

   #if defined(_MY_VERBOSE_TEDIOUS)
	my_log << "copy Buffer";
   #endif
	for(int i=0; i<N; ++i)
	{
		FUx.val[i][0] = Buffer_FUx.val[i][0];
		FUx.val[i][1] = Buffer_FUx.val[i][1];

		FUy.val[i][0] = Buffer_FUy.val[i][0];
		FUy.val[i][1] = Buffer_FUy.val[i][1];

		FUz.val[i][0] = Buffer_FUz.val[i][0];
		FUz.val[i][1] = Buffer_FUz.val[i][1];

		Fni.val[i][0] = Buffer_Fni.val[i][0];
		Fni.val[i][1] = Buffer_Fni.val[i][1];
	}


   #if defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
	my_log << "done";
   #endif
	return;
}
