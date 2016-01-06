#include "field_base.hpp"

// ToDo : field_imag Operatoren implementieren



// parent class for all types of fields (coordinate-space & fourier-space)

field::field()
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field");
	log << "field(const grid &Omega)";
   #endif
}


field::~field()
{
   #ifdef _MY_VERBOSE_MORE
	logger log("field");
	log << "~field()";
   #endif


}
