#ifndef FIELD_INTERPOLATION_HPP_
#define FIELD_INTERPOLATION_HPP_

#include "field_real.hpp"
#include "field_imag.hpp"

#if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#include <sstream>
#include <string>
#endif


// interpolations Operatoren
void OP_xto2x(const field_real &in, field_real &out);
void OP_yto2y(const field_real &in, field_real &out);
void OP_zto2z(const field_real &in, field_real &out);
void OP_hto2h(const field_real &in, field_real &out);

//
void OP_2xtox(const field_real &in, field_real &out);
void OP_2ytoy(const field_real &in, field_real &out);
void OP_2ztoz(const field_real &in, field_real &out);
void OP_2htoh_lvl0(const field_real &in, field_real &out);
// ToDo : lvl1 ist eigentlich falsch, reparieren
void OP_2htoh_lvl1(const field_real &in, field_real &out);
void OP_XhtoYh_lvl1(const field_real &in, field_real &out,
		            const int type = 0, const double &cvalue = 0.);
void OP_smoothing(const field_real &in, field_real &out);




#endif
