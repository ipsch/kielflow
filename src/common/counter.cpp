#include "counter.hpp"


counter::counter(int s_max) :
step(0),
step_max(s_max)
{	}

bool counter::up(void)
{
	step++;
	if (step > step_max)
		return true; // exceeded
	return false;
}
