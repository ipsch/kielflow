#include "logger.hpp"







int logger::level = 0;
double logger::start = omp_get_wtime();
std::string logger::file = "./foo.log";
std::ofstream logger::log_stream(logger::file.c_str(), std::ofstream::trunc);



logger::logger(std::string new_name) :
	name(new_name)
{
	level++;
}

logger::~logger()
{
	level--;
}

void logger::log_msg(std::string msg)
{
	log_stream << msg << std::endl;
	log_stream.flush();
	return;
}



void logger::save(double * DEST)
{



	return;
}


/*
void logger::save(fftw_complex * DEST)
{


	return;
}
*/

