#ifndef GRID_FO_HPP_
#define GRID_FO_HPP_

#if defined(_MY_VERBOSE) || defined(_MY_VERBOSE_MORE) || defined(_MY_VERBOSE_TEDIOUS)
#include "logger.hpp"
#endif

#include "grid_base.hpp"
#include "grid_Co.hpp"

class grid_Co;

class grid_Fo : public grid
{
public :

	grid_Fo(const axis_Fo &x, const axis_Fo &y, const axis_Fo &z);
	grid_Fo(const grid_Fo &that);
	grid_Fo(const grid_Co &that);
	~grid_Fo();

	grid_Fo * clone() const;
	grid_Co * reziprocal() const;

	int index_at(const int &i, const int &j, const int &k) const;
	void ijk_at(const int &index, int &i, int &j, int &k) const;

	//friend bool operator==(const grid_Fo &c1, const grid_Fo &c2);

};


// ToDo : Operatoren grid_Fo
//bool operator==(const grid_Fo &c1, const grid_Fo &c2);



#endif
