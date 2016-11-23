#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>   // sets format of output-stream (scientific & fixed) (is this really necessary??)

#include <omp.h>
//double omp_get_wtime(void);


class logger
{
public :
	logger(std::string new_name);
	~logger();
	static std::ofstream log_stream;
	static int level;
	static double start;
	void log_msg(std::string msg);
	void save(double * DEST);
	//void save(fftw_complex * DEST);
	static void set_file(const std::string &name_logfile);

	void operator<<(const std::string& msg)
	{
		log_stream.width(8);
		log_stream << std::fixed << std::setprecision(1) << std::right;
		double now = omp_get_wtime();
		log_stream << now - start << "s ";
		for(int i=level-1; i!=0; --i)
			log_stream << "\t";
		log_stream << "(" << name << "):" << msg << std::endl;
		log_stream.flush();
	    return;
	}

	void operator<<(const double &number)
	{
		log_stream.width(8);
		log_stream << std::fixed << std::setprecision(1) << std::right;
		double now = omp_get_wtime();
		log_stream << now - start << "s ";
		for(int i=level-1; i!=0; --i)
			log_stream << "\t";
		log_stream << std::scientific;
		log_stream << std::setprecision(6);
		log_stream << "(" << name << "):" << number;
		if(number!= number)
			log_stream << "nan! ";
		log_stream << std::endl;
		log_stream.flush();
	    return;
	}

    template<typename T>
	void operator<<(const T& msg)
	{
		log_stream.width(8);
		log_stream << std::fixed << std::setprecision(1) << std::right;
		double now = omp_get_wtime();
		log_stream << now - start << "s ";
		for(int i=level-1; i!=0; --i)
			log_stream << "\t";
		log_stream << "(" << name << "):" << msg << std::endl;
		log_stream.flush();
	    return;
	}

private:
    static std::string file;
	std::string name;



} ;








#endif
