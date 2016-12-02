// Parameters.h is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.

#ifndef Parameters_H
#define Parameters_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;

class Parameters 
{
 private:
  // switches:
  int initMethod;

  int Nc;      // number of colors (SU(Nc))
  int size;    // the length of the lattice (make it 2^n, with n integer)
  int runningCoupling; // switch to decide if alpha_s should run (0 constant alpha_s, 1 running coupling)
  double R;    
  int seed;    // random seed that's added to the current time to generate the full seed
  double ds;   // 'time' step
  int Ny;      // longitudinal 'resolution' (see Lappi, Eur. Phys. J. C55,285)
  double g2mu; // g^2 mu
  double Qs;   // Q_s, to be dynamically determined
  int steps;   // number of rapidity steps
  int measureSteps; // number of steps in interval between measurements
  int mode;    // mode: (1) run the evolution, (2) analysis with files from disk 
  double mu0;  // cutoff to avoid the Landau pole in the 1-loop running coupling expression
  double LambdaQCD; // LambdaQCD in units of g^2mu
  double g; // coupling g needed in the initU3 where g^2mu does not scale out
  double kappa4Factor; // factor that multiplies the ratio of kappa4/(g^2 mu^2)^3
  double m;     // mass term to regulate infrared divergence (of order Lambda_QCD)
  double L;     // length in fm of the lattice - needed when including mass regulator - new scale
  string input_wline; // filename containing Wilson lines used as an input - note that size and NC should match!
  string output_dir;  // Directory where to save intermediate Wilson lines

 public:

  // constructor:
  Parameters() 
    {
      PI = 3.141592654;
      hbarc = 0.1973269631;
    }

   double PI;
  double hbarc;   // functions to access the private variables:
  void setSeed(int x) {seed=x;}
  int getSeed() {return seed;}
  void setNc(int x) {Nc=x;}
  int getNc() {return Nc;}
  void setNy(int x) {Ny=x;}
  int getNy() {return Ny;}
  void setSize(int x) {size=x;}
  int getSize() {return size;}
  void setR(double x) {R=x;}
  double getR() {return R;}
  void setDs(double x) {ds=x;}
  double getDs() {return ds;}
  void setg2mu(double x) {g2mu=x;}
  double getg2mu() {return g2mu;}
  void setQs(double x) {Qs=x;}
  double getQs() {return Qs;}
  void setSteps(int x) {steps=x;};
  int getSteps() {return steps;}
  void setMeasureSteps(int x) {measureSteps=x;};
  int getMeasureSteps() {return measureSteps;}
  void setMode(int x) {mode=x;};
  int getMode() {return mode;}
  void setRunningCoupling(int x) {runningCoupling=x;};
  int getRunningCoupling() {return runningCoupling;}
  void setMu0(double x) {mu0=x;}
  double getMu0() {return mu0;}
  void setg(double x) {g=x;}
  double getg() {return g;}
  void setLambdaQCD(double x) {LambdaQCD=x;}
  double getLambdaQCD() {return LambdaQCD;}
  void setkappa4Factor(double x) {kappa4Factor=x;}
  double getkappa4Factor() {return kappa4Factor;}
  void setm(double x) {m=x;}
  double getm() {return m;}
  void setL(double x) {L=x;}
  double getL() {return L;}
  void setInputWline(string f){ input_wline = f; }
  string getInputWline() { return input_wline; }
  void setOutputDir(string f){ output_dir = f; }
  string getOutputDir() { return output_dir; }
  
  // switches:
  void setInitMethod(int x) {initMethod=x;}
  int getInitMethod() {return initMethod;}

};

#endif // Parameters_H
