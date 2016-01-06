#ifndef FIELD_INTERPOLATION_HPP_
#define FIELD_INTERPOLATION_HPP_

#include "field_real.hpp"
#include "field_imag.hpp"

void OP_hto2h(const field_real &in, field_real &out);
void OP_2htoh_lvl0(const field_real &in, field_real &out);
// ToDo : lvl1 ist eigentlich falsch, reparieren
void OP_2htoh_lvl1(const field_real &in, field_real &out);
void OP_smoothing(const field_real &in, field_real &out);


#endif
