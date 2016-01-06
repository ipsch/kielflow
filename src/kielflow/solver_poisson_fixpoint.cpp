






/*

bool fixed_point_converged(fftw_complex *FPh_old , fftw_complex *FPh)
{
   #ifdef _MY_VERBOSE
	logger log("fixed_point_converged");
   #endif
	double eps = 1.e-8;
   #ifdef _MY_VERBOSE
	log << "           fixed_point_converged(): check if method has converged... ";
   #endif
	double min, max, sum;
	sum = min = max = cabs(FPh_old[0] - FPh[0]);
	for(int i=1; i<domain::N_imag_ges; ++i)
	{
		sum += cabs(FPh_old[i] - FPh[i]);
		if(min>cabs(FPh_old[i] - FPh[i]))
			min = cabs(FPh_old[i] - FPh[i]);
		if(max< cabs(FPh_old[i] - FPh[i]))
			max = cabs(FPh_old[i] - FPh[i]);
	}


	if(sum<=eps)
	{
       #ifdef _MY_VERBOSE
		log << "           fixed_point_converged(): converged!";
       #endif
		return true;
	}
   #ifdef _MY_VERBOSE
	std::stringstream sstr;
	log << "           fixed_point_converged(): not converged";
	sstr << "           inf= " << min << " sup= " << max << " norm= " << sum << "";
	log << sstr.str();
   #endif
	return false;
}

void fixed_point_iteration(fftw_complex* FPh, fftw_complex* Fni)
// control-loop for the fixed-point-iteration method
{
   #ifdef _MY_VERBOSE
	logger log("fixed_point_iteration");
   #endif
   #ifdef _MY_VERBOSE
	log <<"         fixed_point_iteration(): start";
   #endif
	bool converged = false;  // true if fixed-point-iteration method has converged
	bool exceeded  = false;  // true if fixed-point-iteration method took too long
	bool error = false;      // true if nan
	counter iteration(40);   // set counter for how many iterations are allowed
	fftw_complex *FPh_old = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (domain::N_imag_ges)); // Backup of Potential to compare old und new Phi


	while( (!converged) && (!exceeded) && (!error))
	{
       #ifdef _MY_VERBOSE
		std::stringstream sstr;
		sstr << "fixed_point_iteration(): step(" << iteration.show() << ")";
		log << sstr.str();
       #endif
		for(int i=0; i<domain::N_imag_ges; ++i)            // create a backup of Ph
			FPh_old[i] = FPh[i];
		fixed_point_step(FPh, Fni);                      // derive new Solution for Ph
		exceeded  = iteration.up();                      // check if exceeded
		converged = fixed_point_converged(FPh_old, FPh); // check if converged (compare Ph and Ph_old)
		// ToDo : error = diagnostics_nan(FPh);
	}
	// clean up
	fftw_free(FPh_old);

	if (exceeded)
	{
       #ifdef _MY_VERBOSE
		log << "fixed_point_iteration(): Warning!!";
		log << "maximum iterations exceeded. returning without proper result!!";
		log << "fixed_point_iteration(): Warning!!";
		log << "maximum iterations exceeded. returning without proper result!!";
		return;
       #endif
	}
	else if (error)
	{
       #ifdef _MY_VERBOSE
		log <<"         fixed_point_iteration(): Error!!";
		log <<"                   Result ist nan!!";
       #endif
		return;
	}

   #ifdef _MY_VERBOSE
	log <<"         fixed_point_iteration(): method converged -> returning to RHS";
   #endif

	return;
}

*/
