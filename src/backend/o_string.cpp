#include "o_string.hpp"



std::string FindAndReplace(std::string &s, std::string toReplace, std::string replaceWith)
{
    return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
}

std::string ExtractDirectory( const std::string& path, char delimiter)
  //
  // Returns everything, including the trailing path separator, except the filename
  // part of the path.
  //
  // "/foo/bar/baz.txt" --> "/foo/bar/"
  {
  return path.substr( 0, path.find_last_of( delimiter ) + 1 );
  }

std::string ExtractFilename( const std::string& path, char delimiter)
  //
  // Returns only the filename part of the path.
  //
  // "/foo/bar/baz.txt" --> "baz.txt"
  {
  return path.substr( path.find_last_of( delimiter ) + 1 );
  }

std::string ExtractFileExtension( const std::string& path, char delimiter)
  //
  // Returns the file's extension, if any. The period is considered part
  // of the extension.
  //
  // "/foo/bar/baz.txt" --> ".txt"
  // "/foo/bar/baz" --> ""
{
	std::string filename = ExtractFilename( path, delimiter );
	std::string::size_type n = filename.find_last_of( '.' );
	if (n != std::string::npos)
		return filename.substr( n );
	return std::string();
	}

std::string ChangeFileExtension( const std::string& path, const std::string& ext, char delimiter)
  //
  // Modifies the filename's extension. The period is considered part
  // of the extension.
  //
  // "/foo/bar/baz.txt", ".dat" --> "/foo/bar/baz.dat"
  // "/foo/bar/baz.txt", "" --> "/foo/bar/baz"
  // "/foo/bar/baz", ".txt" --> "/foo/bar/baz.txt"
  //
  {
  std::string filename = ExtractFilename( path, delimiter );
  return ExtractDirectory( path, delimiter )
       + filename.substr( 0, filename.find_last_of( '.' ) )
       + ext;
  }
