// Setup.cpp is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.
#include "Setup.h"

//**************************************************************************
// Setup class.


//**************************************************************************
// Parameter I/O

//reads a string
char *Setup::StringFind(const char* file_name, const char *st)
{
  // char* s = new char[80];
  // char* x = new char[80];
  char* s = new char[180];
  char* x = new char[180];
  memset(s,1,179);
 
  FILE *input, *tmp_file;
  int ind, check;
  
  static int flag = 0;
  
  if(flag == 0)
    {
      if(!IsFile(file_name))
	{
	  cerr << "The input file named " << file_name << " is absent. Exiting." << endl;
	  exit(1);
	}/* if !IsFile */
      flag = 1;
    }/* if flag == 0 */
  
  input = fopen(file_name,"r");
  
  //x = char_malloc(80);
  //s = char_malloc(80);
  
  check=fscanf(input, "%s", s);
  ind = 0;
  while(strcmp(s, "EndOfFile") != 0)
    {
      check=fscanf(input, "%s", x);
      if(strcmp(s, st) == 0)
	{
	  ind++;
	  fclose(input);
	  delete[] s;
	  return x;
	}/* if right, return */
        //delete[] s;
	//s = char_malloc(80);
      check=fscanf(input, "%s", s);
    }/* while */
  
  fclose(input);
  
  if(ind == 0)
    {
      cerr << st << " not found in " << file_name << ". Exiting." << endl;
      delete[] x;
      delete[] s;
     exit(1);
    }
 }/* StringFind */

// reads a double using stringfind:
double Setup::DFind(const char *file_name, const char *st)
{
  char *s;
  s = new char[180];
  memset(s,1,179);
 
  double x;
  
  s = StringFind(file_name, st);
  
  sscanf(s, "%lf", &x);

  delete[] s;
  return x;
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
  static int isf;
  static int ind = 0;
  char st[180];
  FILE *temp;
  
  if( (temp = fopen(file_name,"r")) == NULL) return 0;
  else 
    {
      fclose(temp);
      return 1;
    }
}/* IsFile */
