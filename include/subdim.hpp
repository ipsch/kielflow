#ifndef SUBDIM_HPP_
#define SUBDIM_HPP_

class subdim
{

public:
	int plane;
	double direction;
	int xpos;
	int ypos;
	int zpos;


	int * i;
	int * j;
	int * k;
	int Nx;
	int Ny;
	int Nz;

	int i_fixed;
	int i_slow;
	int i_fast;

	int N_slow;
	int N_fast;


	int index(const int &i, const int &j, const int &k) const;
	void set(void);
private :
};

#endif /* end of SUBDIM_HPP_ */
