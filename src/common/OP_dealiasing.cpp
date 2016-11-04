#include "OP_dealiasing.hpp"



void OP_dealiasing::operator() (field_imag &target) const
{
	for(int ijk=0; ijk<N; ijk++)
	{
		target.val[ijk][0]*=filter[ijk];
		target.val[ijk][1]*=filter[ijk];
	}
	return;
}

OP_dealiasing::~OP_dealiasing()
{
	delete[] filter;
}
