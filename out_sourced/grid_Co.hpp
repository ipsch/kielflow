#ifndef GRID_CO_HPP_
#define GRID_CO_HPP_

#include "grid_base.hpp"
#include "grid_Fo.hpp"

class grid_Fo;

class grid_Co : public grid
{
public :

	grid_Co(const axis_Co &x, const axis_Co &y, const axis_Co &z);
	grid_Co(const grid_Co &that);
	grid_Co(const grid_Fo &that);
	~grid_Co();

	grid_Co * clone() const;
	grid_Fo * reziprocal() const;

	int index_at(const int &i, const int &j, const int &k) const;
	void ijk_at(const int &index, int &i, int &j, int &k) const;

	int index_real(const int &i, const int &j, const int &k) const;
	int index_imag(const int &i, const int &j, const int &k) const;

	//friend bool operator==(const grid_Co &c1, const grid_Co &c2);
	void set_resolution(const int &N_x, const int &N_y, const int &N_z);

};





// ToDo : Operatoren grid_Co
//bool operator==(const grid_Co &c1, const grid_Co &c2);


#endif
