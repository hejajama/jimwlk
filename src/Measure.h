// Measure.h is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.

#ifndef Measure_H
#define Measure_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <gsl/gsl_sf_bessel.h>  // include gsl for the Bessel functions needed to get F_2 and F_L

#include "Lattice.h"
#include "Parameters.h"
#include "Matrix.h"
#include "Group.h"
#include "FFT.h"

using namespace std;

class Measure {

 private:
  Matrix **tempC;
  int size;
  int Nc;
  double ds;
  double g2mu;
  FFT *fft;
  double Pi;

 public:
  
  // Constructor.
  Measure(Parameters *param, const int nn[]) 
    {
      fft = new FFT(nn);
      size = param->getSize();
      Nc = param->getNc();
      g2mu = param->getg2mu();
      Pi = param->PI;
      ds = param->getDs();
      tempC = new Matrix*[size];
      for(int i=0; i<size; i++)
	{
	  tempC[i] = new Matrix(Nc,0.);
	}
    };
  
  ~Measure() 
    {
      for(int i=0; i<size; i++)
	{
	  delete tempC[i];
	}
      
      delete [] tempC;
    };
  
  void storeU(Parameters *param, Lattice *lat, int ids); // store the current U (when called) into the array element Uy[ids] for later unequal rap. study
  void dipoleOperator(Parameters *param, Lattice *lat, int ids); // lattice content and evolution step is passed to the measuring function  
  void fourPointFunction(Parameters *param, Lattice *lat, int ids); // measure factorization
  void sixPointFunctionSquare(Parameters *param, Lattice *lat, int ids); // measure factorization
  void sixPointFunctionLine(Parameters *param, Lattice *lat, int ids); // measure factorization
  void fourPointFunctionSquare(Parameters *param, Lattice *lat, int ids); // measure factorization
  void fourPointFunctionLine(Parameters *param, Lattice *lat, int ids); // measure factorization
  void fourPointFunctionLineUnequalY(Parameters *param, Lattice *lat, int ids); // measure factorization
  void twoPointFunctionInK(Parameters *param, Lattice *lat, int ids); // lattice content and evolution step is passed to the measuring function  
  void fourPointFunctionInK(Parameters *param, Lattice *lat, int ids); // measure the unequal rapidity correlations in momentum space 
  void twoPointFunctions(Parameters *param, Lattice *lat, int ids); // mainly what Dionysis had asked for
  void analyze(Parameters *param); // routine to analyze data read from disk  
  void analyzeS6(Parameters *param); // routine to analyze data read from disk  
  void computeF2(double ***C, double ds, int steps, Parameters *param); // routine to compute F_2 from N=(<Tr(1-U^\dag U)>/Nc)   
  void TwoDCorrelator(Parameters *param, Lattice *lat, int ids);
  double PsiSq(double z, double QSq, double r, int selector, Parameters *param); // This is |Psi_{T,L}|^2, selector=0 is T, selector=1 is L (transverse and longitudinal)
  double Chi(double (Measure::*f)(double, double, double, int, Parameters*), double QSq, double r, int selector, Parameters *param); // integral over z of Psi

};

#endif // Measure_H
