#include "RKO4.hpp"





Runge_kutta_O4::Runge_kutta_O4(interface_rhs &rhs, double dt) :
	my_rhs(rhs), t_(0.), dt_(dt)
{
	my_logfile = "./diagnostics/CFL_RKO4.log";
	std::ofstream diagnose_stream(my_logfile, std::ios::trunc);
	diagnose_stream.close();
}


void Runge_kutta_O4::solve(field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni)
{
   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_TEDIOUS)
	logger log("RK-method");
   #endif

	int N = FUx.N;

	//dt_ = cfl_all(dt_,FUx,FUy,FUz);
	//std::cout <<"t+dt = " << t_ << "+" << dt_ << std::endl;


	//ToDo : was ist sinnvolle Zeit diskretisierung ?

	// Speicher für Runge-Kutta-Koeffizienten allokieren (hier ist FniCHTS optimiert)
   #if defined(_MY_TEDIOUS)
	log << "allocating fields (fourierspace) (for RK-coefficients)";
   #endif
	static field_imag RK0_FUx(*FUx.my_grid), RK0_FUy(*FUy.my_grid), RK0_FUz(*FUz.my_grid), RK0_Fni(*Fni.my_grid);
	static field_imag RK1_FUx(*FUx.my_grid), RK1_FUy(*FUy.my_grid), RK1_FUz(*FUz.my_grid), RK1_Fni(*Fni.my_grid);
	static field_imag RK2_FUx(*FUx.my_grid), RK2_FUy(*FUy.my_grid), RK2_FUz(*FUz.my_grid), RK2_Fni(*Fni.my_grid);
	static field_imag RK3_FUx(*FUx.my_grid), RK3_FUy(*FUy.my_grid), RK3_FUz(*FUz.my_grid), RK3_Fni(*Fni.my_grid);
	static field_imag RK4_FUx(*FUx.my_grid), RK4_FUy(*FUy.my_grid), RK4_FUz(*FUz.my_grid), RK4_Fni(*Fni.my_grid);




   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_TEDIOUS)
	log << "coefficient k0";
   #endif
	for(int i=0; i<FUx.N; ++i)
	{
		RK0_FUx.val[i][0] = RK1_FUx.val[i][0] = FUx.val[i][0];
		RK0_FUx.val[i][1] = RK1_FUx.val[i][1] = FUx.val[i][1];

		RK0_FUy.val[i][0] = RK1_FUy.val[i][0] = FUy.val[i][0];
		RK0_FUy.val[i][1] = RK1_FUy.val[i][1] = FUy.val[i][1];

		RK0_FUz.val[i][0] = RK1_FUz.val[i][0] = FUz.val[i][0];
		RK0_FUz.val[i][1] = RK1_FUz.val[i][1] = FUz.val[i][1];

		RK0_Fni.val[i][0] = RK1_Fni.val[i][0] = Fni.val[i][0];
		RK0_Fni.val[i][1] = RK1_Fni.val[i][1] = Fni.val[i][1];
	}

   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_TEDIOUS)
	log << "evaluating 1.coefficient";
   #endif
	my_rhs.solve(RK1_FUx, RK1_FUy, RK1_FUz, RK1_Fni); // 1. Koeff


	{//ToDo : Testweise diagnose

	double tmp_val;
	std::ofstream diagnose_stream(my_logfile, std::ios::app);
	field_real Buffer(*FUx.my_grid);

	iFFT(RK1_FUx,Buffer);
	tmp_val = supremum(Buffer);
	diagnose_stream << tmp_val << "\t";

	iFFT(RK1_FUy,Buffer);
	tmp_val = supremum(Buffer);
	diagnose_stream << tmp_val << "\t";

	iFFT(RK1_FUz,Buffer);
	tmp_val = supremum(Buffer);
	diagnose_stream << tmp_val << "\t";

	iFFT(RK1_Fni,Buffer);
	tmp_val = supremum(Buffer);
	diagnose_stream << tmp_val << "\n";

	}


	for(int i=0; i<FUx.N; ++i)
	{
		RK2_FUx.val[i][0] = RK0_FUx.val[i][0] + 0.5*dt_*RK1_FUx.val[i][0];
		RK2_FUx.val[i][1] = RK0_FUx.val[i][1] + 0.5*dt_*RK1_FUx.val[i][1];

		RK2_FUy.val[i][0] = RK0_FUy.val[i][0] + 0.5*dt_*RK1_FUy.val[i][0];
		RK2_FUy.val[i][1] = RK0_FUy.val[i][1] + 0.5*dt_*RK1_FUy.val[i][1];

		RK2_FUz.val[i][0] = RK0_FUz.val[i][0] + 0.5*dt_*RK1_FUz.val[i][0];
		RK2_FUz.val[i][1] = RK0_FUz.val[i][1] + 0.5*dt_*RK1_FUz.val[i][1];

		RK2_Fni.val[i][0] = RK0_Fni.val[i][0] + 0.5*dt_*RK1_Fni.val[i][0];
		RK2_Fni.val[i][1] = RK0_Fni.val[i][1] + 0.5*dt_*RK1_Fni.val[i][1];
	}

   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_TEDIOUS)
	log << "evaluating 2.coefficient";
   #endif
	my_rhs.solve(RK2_FUx, RK2_FUy, RK2_FUz, RK2_Fni); // 2. Koeff

	for(int i=0; i<FUx.N; ++i)
	{
		RK3_FUx.val[i][0] = RK0_FUx.val[i][0] + 0.5*dt_*RK2_FUx.val[i][0];
		RK3_FUx.val[i][1] = RK0_FUx.val[i][1] + 0.5*dt_*RK2_FUx.val[i][1];

		RK3_FUy.val[i][0] = RK0_FUy.val[i][0] + 0.5*dt_*RK2_FUy.val[i][0];
		RK3_FUy.val[i][1] = RK0_FUy.val[i][1] + 0.5*dt_*RK2_FUy.val[i][1];

		RK3_FUz.val[i][0] = RK0_FUz.val[i][0] + 0.5*dt_*RK2_FUz.val[i][0];
		RK3_FUz.val[i][1] = RK0_FUz.val[i][1] + 0.5*dt_*RK2_FUz.val[i][1];

		RK3_Fni.val[i][0] = RK0_Fni.val[i][0] + 0.5*dt_*RK2_Fni.val[i][0];
		RK3_Fni.val[i][1] = RK0_Fni.val[i][1] + 0.5*dt_*RK2_Fni.val[i][1];
	}

   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_TEDIOUS)
	log << "evaluating 3.coefficient";
   #endif
	my_rhs.solve(RK3_FUx, RK3_FUy, RK3_FUz, RK3_Fni); // 3. Koeff

	for(int i=0; i<FUx.N; ++i)
	{
		RK4_FUx.val[i][0] = RK0_FUx.val[i][0] + 1.0*dt_*RK3_FUx.val[i][0];
		RK4_FUx.val[i][1] = RK0_FUx.val[i][1] + 1.0*dt_*RK3_FUx.val[i][1];

		RK4_FUy.val[i][0] = RK0_FUy.val[i][0] + 1.0*dt_*RK3_FUy.val[i][0];
		RK4_FUy.val[i][1] = RK0_FUy.val[i][1] + 1.0*dt_*RK3_FUy.val[i][1];

		RK4_FUz.val[i][0] = RK0_FUz.val[i][0] + 1.0*dt_*RK3_FUz.val[i][0];
		RK4_FUz.val[i][1] = RK0_FUz.val[i][1] + 1.0*dt_*RK3_FUz.val[i][1];

		RK4_Fni.val[i][0] = RK0_Fni.val[i][0] + 1.0*dt_*RK3_Fni.val[i][0];
		RK4_Fni.val[i][1] = RK0_Fni.val[i][1] + 1.0*dt_*RK3_Fni.val[i][1];
	}

   #if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_TEDIOUS)
	log << "evaluating 4.coefficient";
   #endif
	my_rhs.solve(RK4_FUx, RK4_FUy, RK4_FUz, RK4_Fni); // 4. Koeff

	// Alle Rungekutta-Koeffizienten zu iteriertem Schritt zusammen bauen
   #if defined(_MY_VERBOSE_MORE) || defined(_MY_TEDIOUS)
	log << "evaluating sum";
   #endif
	for(int i=0; i<FUx.N; ++i)
	{
		FUx.val[i][0] = RK0_FUx.val[i][0] + (dt_/6.)*(1.*RK1_FUx.val[i][0] + 2.*RK2_FUx.val[i][0] + 2.*RK3_FUx.val[i][0] + 1.*RK4_FUx.val[i][0]);
		FUx.val[i][1] = RK0_FUx.val[i][1] + (dt_/6.)*(1.*RK1_FUx.val[i][1] + 2.*RK2_FUx.val[i][1] + 2.*RK3_FUx.val[i][1] + 1.*RK4_FUx.val[i][1]);

		FUy.val[i][0] = RK0_FUy.val[i][0] + (dt_/6.)*(1.*RK1_FUy.val[i][0] + 2.*RK2_FUy.val[i][0] + 2.*RK3_FUy.val[i][0] + 1.*RK4_FUy.val[i][0]);
		FUy.val[i][1] = RK0_FUy.val[i][1] + (dt_/6.)*(1.*RK1_FUy.val[i][1] + 2.*RK2_FUy.val[i][1] + 2.*RK3_FUy.val[i][1] + 1.*RK4_FUy.val[i][1]);

		FUz.val[i][0] = RK0_FUz.val[i][0] + (dt_/6.)*(1.*RK1_FUz.val[i][0] + 2.*RK2_FUz.val[i][0] + 2.*RK3_FUz.val[i][0] + 1.*RK4_FUz.val[i][0]);
		FUz.val[i][1] = RK0_FUz.val[i][1] + (dt_/6.)*(1.*RK1_FUz.val[i][1] + 2.*RK2_FUz.val[i][1] + 2.*RK3_FUz.val[i][1] + 1.*RK4_FUz.val[i][1]);

		Fni.val[i][0] = RK0_Fni.val[i][0] + (dt_/6.)*(1.*RK1_Fni.val[i][0] + 2.*RK2_Fni.val[i][0] + 2.*RK3_Fni.val[i][0] + 1.*RK4_Fni.val[i][0]);
		Fni.val[i][1] = RK0_Fni.val[i][1] + (dt_/6.)*(1.*RK1_Fni.val[i][1] + 2.*RK2_Fni.val[i][1] + 2.*RK3_Fni.val[i][1] + 1.*RK4_Fni.val[i][1]);
	}

	t_ += dt_;



   #if defined(_MY_TEDIOUS)
	log << "finished";
   #endif
	return;
}



