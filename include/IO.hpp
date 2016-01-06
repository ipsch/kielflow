#ifndef IO_HPP
#define IO_HPP

// set (git)-repository-Version to unknown
// (if not defined in makefile during compile process)
#ifndef VERSION_STRING
#define VERSION_STRING "unknown"
#endif

#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include "fftw3.h"
#include <iostream>

#include <vector>

#include "parameters.hpp"
#include "field.hpp"
#include "grid.hpp"
#include "axis.hpp"
#include "operations.hpp"
#include "particle.hpp"

#include <ctime>



#include <echelon/echelon.hpp>
#include <omp.h>
#include <iomanip>   // sets format of output-stream (scientific & fixed) (is this really necessary??)



#ifdef _MY_VERBOSE
#include "logger.hpp"
#endif





class subdim
{

public:
	double default_plane;
	double default_direction;
	double default_xpos;
	double default_ypos;
	double default_zpos;
};

void save_1d(const field_real &Ux, const field_real &Uy, const field_real  &Uz, const field_real  &ni, const field_real &Ph,
		     const subdim & sdim, const std::string &path = "./data/plot_data.dat");
void save_1d(const field_real &Ux,const subdim & sdim, const std::string &path = "./data/plot_data.dat");


void save_2d(const field_real &Ux, const subdim & sdim,
		     const std::string &path = "./data/splot_data.dat");
void save_2d(const field_real &Ux, const field_real &Uy, const field_real  &Uz, const field_real  &ni, const field_real &Ph,
		     const subdim & sdim, const std::string &path = "./data/splot_data.dat");
void save_3d(const std::string &file_name , const field_real &Ux, const field_real &Uy, const field_real  &Uz, const field_real  &ni, const field_real &Ph);
void save_3d(const field_real &XX, const std::string &path);

void file_create(const std::string &path);





void       save_grid      (const grid_Co &Target,
		                   const std::string &path = "./data/fields.h5");
grid_Co    load_grid      (const std::string &path = "./data/fields.h5");



void       save_time(const double &Time, const std::string &path = "./data/fields.h5");
double     load_time(const std::string &path = "./data/fields.h5");



void       save_parameters(const parameters &Target,
		                   const std::string &path = "./data/fields.h5");
parameters load_parameters(const std::string &path = "./data/fields.h5");



void       save_field_real(const std::string &data_set, const field_real &Target,
		                   const std::string &path = "./data/fields.h5");
field_real load_field_real(const std::string &data_set,
		                   const std::string &path = "./data/fields.h5");



void       save_field_imag(const std::string &data_set, const field_imag &Target,
		                   const std::string &path = "./data/fields.h5");
field_imag load_field_imag(const std::string &data_set,
		                   const std::string &path = "./data/fields.h5");

void save_slice(const std::string &data_set, const field_imag &Target, const std::string &path);
void load_slice(const std::string &data_set, double * ptr_data,const std::string &path);


void load_particles(std::vector<particle> &particle_list,
		            const std::string &path = "./config/particles.dat");
void save_particles(std::vector<particle> &particle_list,
		            const std::string &path = "./config/particles.dat");


void save_all(std::vector<particle> &particle_list,
		      const grid_Co &Omega,
			  const double &Time,
		      const parameters &Params,
		      const field_imag &FUx, const field_imag &FUy, const field_imag &FUz,
		      const field_imag &Fni, const field_imag &FPh,
		      const std::string &path  = "./data/fields.h5");


void save_slice(const int &step,
		  std::vector<particle> &particle_list,
	      const grid_Co &Omega,
		  const double &Time,
	      const parameters &Params,
	      const field_imag &FUx, const field_imag &FUy, const field_imag &FUz,
	      const field_imag &Fni, const field_imag &FPh);



void save_frame(const grid_Co &Omega, double* Ux, double* Uy, double* Uz, double* ni, double* Ph, const std::string &path);

std::string I4(int number);


void save_major_wavevectors(const std::string &filename, field_imag &data);


#endif








