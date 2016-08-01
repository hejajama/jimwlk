// Infrared.h is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.

#ifndef Infrared_H
#define Infrared_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "Lattice.h"
#include "Parameters.h"
#include "Matrix.h"
#include "Random.h"
#include "Group.h"
#include "FFT.h"

using namespace std;

class Infrared {

 private:
  FFT *fft;
  Matrix** A;
  double As[1024];

 public:
  
  // Constructor.
  Infrared(const int nn[]) 
    {
      fft = new FFT(nn);
    };
  
  ~Infrared() { delete fft; };
  
  void regulate(Lattice *lat, Group *group, Parameters *param, Random *random, int ids);
};

#endif // Infrared_H
