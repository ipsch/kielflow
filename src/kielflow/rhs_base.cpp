#include "rhs_base.hpp"

solver_rhs::solver_rhs(const parameters &Params, solver_poisson &Solver) :
		my_params(Params), my_poisson_solver(&Solver)
{

}


void solver_rhs::evaluate_advection(const field_imag &FUx, const field_imag &FUy, const field_imag &FUz, const field_imag FXX,
		field_imag &return_Buffer)
{
   #ifdef _MY_VERBOSE
	logger log("RHS");
	log << "evaluate_advection(...)";
   #endif

	// input  F[Ux], F[Uy], F[Uz] and F[XX]  (XX denotes a wildcard that can either be Ux, Uy or Uz)
	// output F[ Ux*(d/dx)XX + Uy*(d/dy)XX + Uz*(d/dz)XX ]

	// step1: berechne hier F^{-1}[F[Ux]], F^{-1}[F[Uy]] und F^{-1}[F[Uz]]
	// dies ist identisch mit
	// Ux, Uy und Uz
	field_real Ux(FUx.my_grid);
	field_real Uy(FUy.my_grid);
	field_real Uz(FUz.my_grid);


	static OP_partial_derivative d_dx(FUx,e_x);
	static OP_partial_derivative d_dy(FUx,e_y);
	static OP_partial_derivative d_dz(FUx,e_z);

	iFFT(FUx,Ux);
	iFFT(FUy,Uy);
	iFFT(FUz,Uz);


	// step2: berechne hier F^{-1}[îkx*F[Ux]], F^{-1}[îky*F[Ux]] und F^{-1}[îkz*F[Ux]]
	// dies ist identisch mit
	// (d/dx) Ux, (d/dy) Ux und (d/dz) Ux

	field_imag FBuffer(FXX.my_grid);


	d_dx(FXX,FBuffer);
	field_real dx_XX(FBuffer.my_grid);
	iFFT(FBuffer,dx_XX);


	d_dy(FXX,FBuffer);
	field_real dy_XX(FBuffer.my_grid);
	iFFT(FBuffer,dy_XX);

	d_dz(FXX,FBuffer);
	field_real dz_XX(FBuffer.my_grid);
	iFFT(FBuffer,dz_XX);

	// step3: berechne hier F[ Ux*(d/dx)Ux + Uy*(d/dy)Ux + Uz*(d/dz)Ux ]

	field_real Buffer(FXX.my_grid);


	for(int i=0; i<Ux.N; ++i)
		Buffer.val[i] = Ux.val[i]*dx_XX.val[i] + Uy.val[i]*dy_XX.val[i] + Uz.val[i]*dz_XX.val[i];


	FFT(Buffer, FBuffer);


	// Add result to RHS-Buffer
	for(int i=0; i<FUx.N; ++i)
	{
		return_Buffer.val[i][0] += -FBuffer.val[i][0];
		return_Buffer.val[i][1] += -FBuffer.val[i][1];
	}

   #ifdef _MY_VERBOSE
	log << "done";
   #endif

	return;
}




void solver_rhs::evaluate_dissipation(const field_imag &FXX,
		field_imag &Buffer_FXX)
{
   #ifdef _MY_VERBOSE
	logger log("RHS");
	log <<"evaluate_disspation(...)";
   #endif

	for(int i=0; i<FXX.N; ++i)
	{
		Buffer_FXX.val[i][0] -= my_params.tau*FXX.val[i][0];
		Buffer_FXX.val[i][1] -= my_params.tau*FXX.val[i][1];
	}

   #ifdef _MY_VERBOSE
	log <<"done";
   #endif

	return;
}

void solver_rhs::evaluate_force_E(const field_imag &FPhi,
		field_imag &Buffer_FUx, field_imag &Buffer_FUy, field_imag &Buffer_FUz)
// Evaluate the stato-electrical force upon the velocity field of the ions.
// F = -(q/m)*grad(Phi)
// The factor q/m is already included by equations due to dimensionless quantities
// evaluation of grad(Phi) is done in fourier-space by multiplying fourier coeffcients by (ik) (pay attention to sign an re<->im beeing crossed)
{
   #ifdef _MY_VERBOSE
	logger log("RHS");
	log <<"evaluate_force_E(...)";
   #endif
	for(int i=0; i<FPhi.Nx; ++i)
		for(int j=0; j<FPhi.Ny; ++j)
			for(int k=0; k<FPhi.Nz; ++k)
			{
				int index = FPhi.index(i,j,k);
				double FPhi_re = FPhi.val[index][0];
				double FPhi_im = FPhi.val[index][1];


				Buffer_FUx.val[index][0] += FPhi.my_grid.x_axis->k_val_at(i)*FPhi_im;
				Buffer_FUx.val[index][1] -= FPhi.my_grid.x_axis->k_val_at(i)*FPhi_re;

				Buffer_FUy.val[index][0] += FPhi.my_grid.y_axis->k_val_at(j)*FPhi_im;
				Buffer_FUy.val[index][1] -= FPhi.my_grid.y_axis->k_val_at(j)*FPhi_re;

				Buffer_FUz.val[index][0] += FPhi.my_grid.z_axis->k_val_at(k)*FPhi_im;
				Buffer_FUz.val[index][1] -= FPhi.my_grid.z_axis->k_val_at(k)*FPhi_re;
			}
   #ifdef _MY_VERBOSE
	log <<"done";
   #endif
	return;
}

void solver_rhs::evaluate_force_B(const field_imag &FUy, const field_imag &FUz,
		field_imag &Buffer_FUy, field_imag &Buffer_FUz)
{
	for(int i=0; i<Buffer_FUy.N; ++i)
	{
		Buffer_FUy.val[i][0] += my_params.beta*FUz.val[i][0];
		Buffer_FUy.val[i][1] += my_params.beta*FUz.val[i][1];
	}
	for(int i=0; i<Buffer_FUz.N; ++i)
	{
		Buffer_FUz.val[i][0] -= my_params.beta*FUy.val[i][0];
		Buffer_FUz.val[i][1] -= my_params.beta*FUy.val[i][1];
	}
	return;
}

void solver_rhs::evaluate_diffusion_mu(const field_imag &FXX, field_imag &Buffer_FXX)
{
	for(int i=0; i<FXX.Nx; ++i)
		for(int j=0; j<FXX.Ny; ++j)
			for(int k=0; k<FXX.Nz; ++k)
			{
				double kx = FXX.my_grid.x_axis->k_val_at(i);
				double ky = FXX.my_grid.y_axis->k_val_at(j);
				double kz = FXX.my_grid.z_axis->k_val_at(k);

				double pow_k = kx*kx + ky*ky + kz*kz;
				int index = Buffer_FXX.index(i,j,k);
				Buffer_FXX.val[index][0] += my_params.mu* pow_k * FXX.val[index][0];
				Buffer_FXX.val[index][1] += my_params.mu *pow_k * FXX.val[index][1];

			}

	return;
}

void solver_rhs::evaluate_diffusion_u(const field_imag &Fni,
		field_imag &Buffer_FUx, field_imag &Buffer_FUy, field_imag &Buffer_FUz)
{
   #ifdef _MY_VERBOSE
	logger log("RHS");
	log <<"evaluate_diffusion_u(...)";
   #endif

	OP_partial_derivative d_dx(Buffer_FUx,e_x);
	OP_partial_derivative d_dy(Buffer_FUy,e_y);
	OP_partial_derivative d_dz(Buffer_FUz,e_z);



	field_imag FBuffer(Fni.my_grid);
	field_real Buffer(Fni.my_grid);

	field_real ni(Fni.my_grid);
	iFFT(Fni,ni);

	// x-Komponente
	d_dx(Fni,FBuffer);
	iFFT(FBuffer, Buffer);

	for(int i=0; i<ni.N; ++i)
		Buffer.val[i] = Buffer.val[i]/(my_params.theta*(1.+ni.val[i]));
	FFT(Buffer, FBuffer);

	for(int i=0; i<Buffer_FUx.N; i++)
	{
		Buffer_FUx.val[i][0] += -FBuffer.val[i][0];
		Buffer_FUx.val[i][1] += -FBuffer.val[i][1];
	}

	// y-Komponente
	d_dy(Fni,FBuffer);
	iFFT(FBuffer, Buffer);

	for(int i=0; i<ni.N; ++i)
		Buffer.val[i] = Buffer.val[i]/(my_params.theta*(1.+ni.val[i]));
	FFT(Buffer, FBuffer);

	for(int i=0; i<Buffer_FUy.N; i++)
	{
		Buffer_FUy.val[i][0] += -FBuffer.val[i][0];
		Buffer_FUy.val[i][1] += -FBuffer.val[i][1];
	}


	// z-Komponente
	d_dz(Fni,FBuffer);
	iFFT(FBuffer, Buffer);

	for(int i=0; i<ni.N; ++i)
		Buffer.val[i] = Buffer.val[i]/(my_params.theta*(1.+ni.val[i]));
	FFT(Buffer, FBuffer);

	for(int i=0; i<Buffer_FUz.N; i++)
	{
		Buffer_FUz.val[i][0] += -FBuffer.val[i][0];
		Buffer_FUz.val[i][1] += -FBuffer.val[i][1];
	}


   #ifdef _MY_VERBOSE
	log <<"done";
   #endif
	return;
}

void solver_rhs::evaluate_diffusion_n(const field_imag &FUx, const field_imag &FUy, const field_imag &FUz,
		field_imag &Buffer_Fni)
	// - grad(u)
{
	int N = FUx.N;

	OP_partial_derivative d_dx(FUx,e_x);
	OP_partial_derivative d_dy(FUy,e_y);
	OP_partial_derivative d_dz(FUz,e_z);

	field_imag FBuffer_x(FUx.my_grid);
	field_imag FBuffer_y(FUy.my_grid);
	field_imag FBuffer_z(FUz.my_grid);

	d_dx(FUx,FBuffer_x);
	d_dy(FUy,FBuffer_y);
	d_dz(FUz,FBuffer_z);

	for(int i=0; i<N;++i)
	{
		Buffer_Fni.val[i][0] += -(FBuffer_x.val[i][0] + FBuffer_y.val[i][0] + FBuffer_z.val[i][0]);
		Buffer_Fni.val[i][1] += -(FBuffer_x.val[i][1] + FBuffer_y.val[i][1] + FBuffer_z.val[i][1]);
	}

	return;
}

void solver_rhs::evaluate_continuityEQ(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni,
		field_imag &Buffer_Fni)
     // - grad(n*u)

// This function evaluates the continuity equation. First a copy of Ux,Uy,Uz and ni is created in
// real-space. In real-space the products, Ux*ni, Uy*ni, Uz*ni are evaluated.
// There after the products are transformed back to fourier-space and are multiplied by ik_x btw. ik_y, ik_z (which
// accounts for the derivative in space). At last the products are added to the output-buffer.
{
   #ifdef _MY_VERBOSE
	logger log("RHS");
	log <<"evaluate_continuityEQ(...)";
   #endif



	// Temporäre größen im Ortsraum (iFFT impliziert)
	field_real ni(Fni.my_grid);
	iFFT(Fni,ni);
	field_real Ux(FUx.my_grid);
	iFFT(FUx,Ux);
	field_real Uy(FUy.my_grid);
	iFFT(FUy,Uy);
	field_real Uz(FUz.my_grid);
	iFFT(FUz,Uz);

	// evaluate non-linearity in realspace FUx \cdot ni (and so on)
	for(int i=0; i<Ux.N; ++i)
	{
		Ux.val[i] = Ux.val[i]*ni.val[i];
		Uy.val[i] = Uy.val[i]*ni.val[i];
		Uz.val[i] = Uz.val[i]*ni.val[i];
	}

	// backtransformation
	field_imag FBuffer_x(Ux.my_grid);
	FFT(Ux,FBuffer_x);
	field_imag FBuffer_y(Uy.my_grid);
	FFT(Uy,FBuffer_y);
	field_imag FBuffer_z(Uz.my_grid);
	FFT(Uz,FBuffer_z);

	// Ableitung bilden

	static OP_partial_derivative d_dx(FUx,e_x);
	static OP_partial_derivative d_dy(FUx,e_y);
	static OP_partial_derivative d_dz(FUx,e_z);

	d_dx(FBuffer_x,FBuffer_x);
	d_dy(FBuffer_y,FBuffer_y);
	d_dz(FBuffer_z,FBuffer_z);

	for(int i=0; i<Buffer_Fni.N; ++i)
	{
		Buffer_Fni.val[i][0] += -(FBuffer_x.val[i][0] + FBuffer_y.val[i][0] + FBuffer_z.val[i][0]);
		Buffer_Fni.val[i][1] += -(FBuffer_x.val[i][1] + FBuffer_y.val[i][1] + FBuffer_z.val[i][1]);
	}


   #ifdef _MY_VERBOSE
	log << "done";
   #endif
	return;
}
