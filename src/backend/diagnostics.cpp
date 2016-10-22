#include"diagnostics.hpp"




/*

double diagnostics_sup(fftw_complex *DEST)
{
	double sup = cabs(DEST[0]);
	for(int i=1; i<domain::N_imag_ges; ++i)
	    if (sup<cabs(DEST[i])) sup=cabs(DEST[i]);
	return sup;
}

double diagnostics_inf(fftw_complex *DEST)
{
	double inf = cabs(DEST[0]);
	for(int i=1; i<domain::N_imag_ges; ++i)
	    if (inf>cabs(DEST[i])) inf=cabs(DEST[i]);
	return inf;
}

double diagnostics_sum(fftw_complex *DEST)
{
	double sum = 0.0;
	for(int i=0; i<domain::N_imag_ges; ++i)
		sum += cabs(DEST[i]);
	return sum;
}

*/

diagnostics::diagnostics() :
			supremum(0.), i_supremum(0), j_supremum(0), k_supremum(0),
			infinum(0.),  i_infinum(0),  j_infinum(0),  k_infinum(0)
{

	std::cout << "create diagnosis" << std::endl;

}

/*
void diagnostics::set_ijk(int index, int &i, int &j, int &k)
{
	// wichtig: hier integer division!!
	i =  index/(domain::Ny*domain::Nz);
	j = (index - i*(domain::Ny*domain::Nz)) / domain::Nz;
	k =  index - i*(domain::Ny*domain::Nz) - j*domain::Nz;
	return;
}
*/

void diagnostics::analyze(const field_real &that)
{
	infinum = supremum = that.val[0];
	for(int i=0; i<that.N;++i)
	{
		if(supremum != max(that.val[i],supremum))
		{
			supremum = max(that.val[i],supremum);
			//that.my_grid->ijk_at(i,i_supremum,j_supremum, k_supremum);
		}

		if(infinum  != min(that.val[i],infinum))
		{
			infinum  = min(that.val[i],infinum);
			//that.my_grid->ijk_at(i,i_supremum,j_supremum, k_supremum);
		}

	}
	std::cout << supremum << std::endl;
	std::cout << infinum << std::endl;
	return;
}
