#include "rhs_naca.hpp"

solver_rhs_naca::solver_rhs_naca(const parameters &Params, solver_poisson &Solver) :
solver_rhs(Params, Solver)
{
   #ifdef _MY_VERBOSE
	logger log("solver_rhs_experimental");
	log << "solver_rhs_experimental(const parameters &Params, solver_poisson &Solver)";
   #endif
}






void solver_rhs_naca::solve(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni, field_imag &FPh)
// The Right-Hand-Side of the PDE is evaluated.
{
   #ifdef _MY_VERBOSE
    logger log("solver_rhs_experimental");
	log << "solve(...)";
   #endif

	// fields in fourier-space
	static int N = FUx.N;
	static field_imag Buffer_FUx(*FUx.my_grid);
	static field_imag Buffer_FUy(*FUx.my_grid);
	static field_imag Buffer_FUz(*FUx.my_grid);
	static field_imag Buffer_Fni(*FUx.my_grid);
	static field_imag Buffer_FPh(*FUx.my_grid);


	static field_real maskH(*FUx.my_grid);
	maskH.fill(&define_H_naca);

	static effect_penalization solid(0.01,maskH, 0.0);




	static double v_stream = my_params.M;
	static double phi_stream = 2*3.14152*(0./360.);

	static effect_force_const driver_x(v_stream*cos(phi_stream));
	//static effect_force_const driver_y(v_stream*sin(phi_stream));
	static effect_diffusion diffusion_u(1.);

	static effect_force_linear dissipation(-my_params.tau);


   #ifdef _MY_VERBOSE
	log << "RHS(): clear Buffer";
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




	// ################## IMPULS EQUATION #####################################
	evaluate_advection(FUx, FUy, FUz, FUx, Buffer_FUx);
	evaluate_advection(FUx, FUy, FUz, FUy, Buffer_FUy);
	evaluate_advection(FUx, FUy, FUz, FUz, Buffer_FUz);

	solid.execute(FUx,Buffer_FUx);
	solid.execute(FUy,Buffer_FUy);
	solid.execute(FUz,Buffer_FUz);

	driver_x.execute(FUx,Buffer_FUx);
	//driver_y.execute(FUy,Buffer_FUy);


	//diffusion_u.execute(FUx,Buffer_FUx);
	//diffusion_u.execute(FUy,Buffer_FUy);
	//diffusion_u.execute(FUz,Buffer_FUz);

	evaluate_diffusion_u(Fni, Buffer_FUx, Buffer_FUy, Buffer_FUz);
   #ifdef _DIFFUSION
	evaluate_diffusion_mu(FUx, Buffer_FUx);
	evaluate_diffusion_mu(FUy, Buffer_FUy);
	evaluate_diffusion_mu(FUz, Buffer_FUz);
   #endif
   #ifdef B_FIELD_TRUE
	evaluate_force_B(          FUy, FUz,                    Buffer_FUy, Buffer_FUz           );
   #endif



	dissipation.execute(FUx, Buffer_FUx);
	dissipation.execute(FUy, Buffer_FUy);
	dissipation.execute(FUz, Buffer_FUz);


	// ################## CONTINUITY EQUATION #################################
	evaluate_continuityEQ(FUx, FUy, FUz, Fni, Buffer_Fni);
	evaluate_diffusion_n( FUx, FUy, FUz, Buffer_Fni);





	// This has to be the last step in evaluation of the rhs
	//solid.penalize_density(Fni,Fni);


   #ifdef _MY_VERBOSE
	log << "RHS(): copy Buffer";
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


	return;
}
