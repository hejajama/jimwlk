// Setup.cpp is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.
#include "Setup.h"
#include <fstream>
#include <stdexcept>

//**************************************************************************
// Setup class.


//**************************************************************************
// Parameter I/O

//reads a string
string Setup::StringFind(const char* file_name, const char *st)
{
  ifstream input(file_name);
  string key;
  string value;
  
  static int flag = 0;
  
  if(flag == 0)
    {
      if(!IsFile(file_name))
	{
	  throw runtime_error(string("The input file named ") + file_name + " is absent.");
	}/* if !IsFile */
      flag = 1;
    }/* if flag == 0 */
  
  while (input >> key)
    {
      if (key == "EndOfFile")
        {
          break;
        }

      if (!(input >> value))
        {
          break;
        }

      if (key == st)
	{
	  return value;
	}/* if right, return */
    }/* while */

  throw runtime_error(string(st) + " not found in " + file_name + ".");
 }/* StringFind */

// reads a double using stringfind:
double Setup::DFind(const char *file_name, const char *st)
{
  return stod(StringFind(file_name, st));
}/* DFind */

// reads an integer using stringfind:
int Setup::IFind(const char *file_name, const char *st)
{
  double f;
  f = DFind(file_name, st);
  
  return (int) (f + 0.5);
}/* IFind */

int Setup::IsFile(const char *file_name)
{
  ifstream temp(file_name);
  return temp.good() ? 1 : 0;
}/* IsFile */
