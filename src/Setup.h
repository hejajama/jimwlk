// Setup.h is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.

#ifndef Setup_H
#define Setup_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;

class Setup {

public:

  // Constructor.
  Setup() {}

  char *StringFind( const char* file_name,  const char* st);
  double DFind( const char* file_name,  const char* st);
  int IFind( const char* file_name,  const char* st);
  int IsFile(const char *file_name);

};

#endif // Setup_H
