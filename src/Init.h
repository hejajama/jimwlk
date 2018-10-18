// Init.h is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.

#ifndef Init_H
#define Init_H

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

class Init {

 private:
  FFT *fft;
  Matrix** A;
  double As[1024];
 
 public:
  
  // Constructor.
  Init(const int nn[]) 
    {
      fft = new FFT(nn);
    };
  
  ~Init() { delete fft; };
  
  void initU(Lattice *lat, Group *group, Parameters *param, Random *random);
  void initU1(Lattice *lat, Group *group, Parameters *param, Random *random);
  void initU2(Lattice *lat, Group *group, Parameters *param, Random *random);
  void initU3(Lattice *lat, Group *group, Parameters *param, Random *random);
  void initU4(Lattice *lat, Group *group, Parameters *param, Random *random);
  void initFromData(Lattice *lat, Group *group, Parameters *param, Random *random);
  void initFromBinaryData(Lattice *lat, Group *group, Parameters *param, Random *random);
  void initU2old(Lattice *lat, Group *group, Parameters *param, Random *random);
  void computeAx(Lattice *lat, Group *group, Parameters *param, Random *random);
  void computeAx2(Lattice *lat, Group *group, Parameters *param, Random *random);
  double SpatialDistribution(int x,int y,double R, int N);
  

};

#endif // Init_H
