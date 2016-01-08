#include "subdim.hpp"



void subdim::set(void)
{

	    // 0 : x fixiert => out: y-z-ebene
		// 1 : y fixiert => out: x-z-ebene
	    // 2 : z fixiert => out: x-y-ebene

	switch( plane )
	{
	case 0:

		i = &xpos;
		j = &i_slow;
	    k = &i_fast;
	    N_slow = Ny;
	    N_fast = Nz;
	    break;
	case 1:
		i = &i_slow;
		j = &ypos;
		k = &i_fast;
	    N_slow = Nx;
	    N_fast = Nz;
	    break;
	case  2:
		i = &i_slow;
		j = &i_fast;
		k = &zpos;
	    N_slow = Nx;
	    N_fast = Ny;
		break;
	} // END of switch


	return;
}
