#ifndef O_STRING_HPP_
#define O_STRING_HPP_


#include <string>
#include <sstream>

std::string FindAndReplace(std::string &s, std::string toReplace, std::string replaceWith);
std::string ExtractDirectory( const std::string& path, char delimiter = '/' );
std::string ExtractFilename( const std::string& path, char delimiter = '/' );
std::string ExtractFileExtension( const std::string& path, char delimiter = '/' );
std::string ChangeFileExtension( const std::string& path, const std::string& ext, char delimiter = '/' );

template <class T> std::string ConvertToString(const T &in)
{
	std::stringstream s_stream;
	s_stream << in;
	return s_stream.str();
}



#endif
