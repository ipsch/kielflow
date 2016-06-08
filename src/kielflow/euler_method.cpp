#include "euler_method.hpp"




euler_method::euler_method(interface_rhs &rhs, double dt) :
my_rhs(rhs), t_(0.), dt_(dt)
{

}


void euler_method::solve(const double &t, field_imag &FUx, field_imag &FUy, field_imag &FUz, field_imag &Fni)
{
   #ifdef _MY_VERBOSE
	logger log("RK-method");
   #endif



	int N = FUx.N;

	//ToDo : Schrittweitensteuerung/-konstrolle
	// Es fehlt eine Schrittweiten-steuerung oder
	// zumindeest Schrittweitenkonstrolle

	dt_ = cfl_all(dt_,FUx,FUy,FUz);
	std::cout <<"t+dt = " << t_ << "+" << dt_ << std::endl;

	// Speicher für Runge-Kutta-Koeffizienten allokieren (hier ist FniCHTS optimiert)
   #ifdef _MY_VERBOSE
	log << "allocating fields (fourierspace) (for RK-coefficients)";
   #endif
	static field_imag RK0_FUx(*FUx.my_grid);
	static field_imag RK0_FUy(*FUy.my_grid);
	static field_imag RK0_FUz(*FUz.my_grid);
	static field_imag RK0_Fni(*Fni.my_grid);


   #ifdef _MY_VERBOSE
	log << "coefficient k0";
   #endif
	for(int i=0; i<FUx.N; ++i)
	{
		RK0_FUx.val[i][0] = FUx.val[i][0];
		RK0_FUx.val[i][1] = FUx.val[i][1];

		RK0_FUy.val[i][0] = FUy.val[i][0];
		RK0_FUy.val[i][1] = FUy.val[i][1];

		RK0_FUz.val[i][0] = FUz.val[i][0];
		RK0_FUz.val[i][1] = FUz.val[i][1];

		RK0_Fni.val[i][0] = Fni.val[i][0];
		RK0_Fni.val[i][1] = Fni.val[i][1];
	}


   #ifdef _MY_VERBOSE
	log << "evaluating 1.coefficient";
   #endif
	my_rhs.solve(t_, RK0_FUx, RK0_FUy, RK0_FUz, RK0_Fni); // 1. Koeff


	// Alle Rungekutta-Koeffizienten zu iteriertem Schritt zusammen bauen
   #ifdef _MY_VERBOSE
	log << "evaluating sum";
   #endif
	for(int i=0; i<FUx.N; ++i)
	{
		FUx.val[i][0] = FUx.val[i][0] + dt_*RK0_FUx.val[i][0];
		FUx.val[i][1] = FUx.val[i][1] + dt_*RK0_FUx.val[i][1];

		FUy.val[i][0] = FUy.val[i][0] + dt_*RK0_FUy.val[i][0];
		FUy.val[i][1] = FUy.val[i][1] + dt_*RK0_FUy.val[i][1];

		FUz.val[i][0] = FUz.val[i][0] + dt_*RK0_FUz.val[i][0];
		FUz.val[i][1] = FUz.val[i][1] + dt_*RK0_FUz.val[i][1];

		Fni.val[i][0] = Fni.val[i][0] + dt_*RK0_Fni.val[i][0];
		Fni.val[i][1] = Fni.val[i][1] + dt_*RK0_Fni.val[i][1];
	}

	t_ += dt_;


   #ifdef _MY_VERBOSE
	log << "finished";
   #endif
	return;
}



