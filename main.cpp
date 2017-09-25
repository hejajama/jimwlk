//#include <stdio.h>
#include <string>
#include <cmath>
#include <iostream>
#include <complex>
#include <fstream>
#include <time.h>
#include <vector>
#include <sstream>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>

#include "Setup.h"
#include "Init.h"
#include "Infrared.h"
#include "Random.h"
#include "FFT.h"
#include "Parameters.h"
#include "Matrix.h"
#include "Lattice.h"
#include "Measure.h"

#define _SECURE_SCL 0
#define _HAS_ITERATOR_DEBUGGING 0
using namespace std;

int readInput(Setup *setup, Parameters *param, int argc, char *argv[]);

void split(const string &s, char delim, vector<string> &elems) {
  stringstream ss(s);
  string item;
  while (getline(ss, item, delim)) {
    elems.push_back(item);
  }
}


vector<string> split(const string &s, char delim) {
  vector<string> elems;
  split(s, delim, elems);
  return elems;
}

// main program
int main(int argc, char *argv[])
{
  // welcome
  cout << endl;
  cout << "JIMWLK solver v0.5" << endl;
  cout << endl;
  // initialize helper class objects
  Setup *setup;
  setup = new Setup();
  Random *random;
  random = new Random();
  Parameters *param;
  param = new Parameters();
  
  // read parameters from file
  readInput(setup, param, argc, argv);
  
  int cells = param->getSize()*param->getSize();
  int Nc2m1 = param->getNc()*param->getNc()-1; // N_c^2-1
  int nn[2];
  int pos;
  double x,y;
  double ds = param->getDs();
  nn[0]=param->getSize();
  nn[1]=param->getSize();
  
  //initialize FFT routines
  FFT *fft;
  fft = new FFT(nn);
  //initialize init object
  Init *init;
  init = new Init(nn);
  //initialize infrared regulator object
  Infrared *infrared;
  infrared = new Infrared(nn);
  //initialize measure object
  Measure *measure;
  measure = new Measure(param,nn);
  
  int mode = param->getMode();
  
  if(mode==2)
  {
    //      measure->analyzeS6(param);
    measure->analyze(param);
    exit(1);
  }
  else if(mode>2)
  {
    cerr << "Only mode=1 (evolution) and mode=2 (analysis) supported. Mode=" << mode << " is unsupported. Exiting." << endl;
    exit(1);
  }
  // if mode=1 continue:
  
  // initialize group
  Group *group;
  group = new Group(param->getNc());
  
  //initialize random generator
  long long rnum;
  rnum=time(0)//;+param->getSeed();
  cout << "Random seed =" << rnum << " made from time "
  << rnum-param->getSeed() << " and argument "
  << param->getSeed() << endl;
  random->init_genrand64(rnum);
  
  // allocate lattice
  Lattice *lat;
  lat = new Lattice(param, param->getNc(), param->getSize());
  cout << " done." << endl;
  
  // initialize U-fields on the lattice
  init->initU(lat, group, param, random);
  cout << " done." << endl;
  
  
  //  infrared->regulate(lat, group, param, random, 0);
  
  for(int i=0; i<nn[0]*nn[1]; i++)
  {
    lat->cells[i]->setUi(lat->cells[i]->getU());
  }
  
  // Save initial state Wilson lines
  //lat->PrintWilsonLines("wline_evolution/step_0");
  
  // allocate memory
  complex<double> ** xi;
  xi = new complex<double>*[param->getSize()*param->getSize()];
  
  complex<double> ** xi2;
  xi2 = new complex<double>*[param->getSize()*param->getSize()];
  
  for(int i=0; i<param->getSize()*param->getSize(); i++)
  {
    xi[i] = new complex<double>[Nc2m1*2];
    xi2[i] = new complex<double>[Nc2m1*2];
  }
  
  // prepare Fourier transforms of K and S
  vector<complex<double> > ** K;
  K = new vector<complex<double> >*[param->getSize()*param->getSize()];
  
  for(int i=0; i<param->getSize()*param->getSize(); i++)
  {
    K[i] = new vector<complex<double> >;
  }
  
  vector<complex<double> > ** S;
  S = new vector<complex<double> >*[param->getSize()*param->getSize()];
  
  for(int i=0; i<param->getSize()*param->getSize(); i++)
  {
    S[i] = new vector<complex<double> >;
  }
  
  // C(K,xi^a), vector in a (color)
  complex<double>  ** CKxi;
  CKxi = new complex<double> *[param->getSize()*param->getSize()];
  
  for(int i=0; i<param->getSize()*param->getSize(); i++)
  {
    CKxi[i] = new complex<double>[Nc2m1];
  }
  
  
  // C(K,U^{ab} xi^a), vector in a (color)
  // H.M. disable
  /*
  complex<double> ** CKUxi;
  CKUxi = new complex<double> *[param->getSize()*param->getSize()];
  
  for(int i=0; i<param->getSize()*param->getSize(); i++)
  {
    CKUxi[i] = new complex<double>[Nc2m1];
  }
  
  Matrix ** UA;
  UA = new Matrix *[param->getSize()*param->getSize()];
  
  for(int i=0; i<param->getSize()*param->getSize(); i++)
  {
    UA[i] = new Matrix(Nc2m1,0);
  }
  
  Matrix ** UA2;
  UA2 = new Matrix *[param->getSize()*param->getSize()];
  
  for(int i=0; i<param->getSize()*param->getSize(); i++)
  {
    UA2[i] = new Matrix(Nc2m1,0);
  }
  */
  
  // H.M. matrix V xsi V^\dagger
  Matrix ** VxsiVx;
  VxsiVx = new Matrix *[param->getSize()*param->getSize()];
  Matrix ** VxsiVy;
  VxsiVy = new Matrix *[param->getSize()*param->getSize()];
  for(int i=0; i<param->getSize()*param->getSize(); i++)
  {
    VxsiVx[i] = new Matrix(param->getNc(),0);
    VxsiVy[i] = new Matrix(param->getNc(),0);
  }
  
  
  double m = param->getm();
  double Pi = param->PI;
  double r2, Lambda2;
  int Nf = 3;
  double alphas;
  Lambda2 = param->getLambdaQCD()*param->getLambdaQCD();///(2.*param->PI*2.*param->PI);
  cout << "in this version Lambda_QCD is the coordinate space Lambda!" << endl;
  cout << "Lambda^2=" << Lambda2 << endl;
  double mu0 = param->getMu0();
  for (int i=0; i<nn[0]; i++)
  {
    for (int j=0; j<nn[1]; j++)
    {
      pos = i*nn[1]+j;
      
      x = lat->cells[pos]->getX();
      y = lat->cells[pos]->getY();
      
      
      
      
      r2 = x*x+y*y;
      // K is a 2d vector (x and y)
      if (r2==0)
      {
        K[pos]->push_back(0.);
        K[pos]->push_back(0.);
        // S is a 1d vector
        S[pos]->push_back(0.);
      }
      else
      {
        double mass_regulator = 1.0;  // if m suppresses long distance tails,
        // K is multiplied by this, which is m*r*K_1(m*r)
        if (m>0)
        {
          // Regulate
          // Note that we use here r2 which is not in lattice units
          // because m is in GeV
          double length = param->getL();
          double phys_x = x/nn[0]*length; //in fm
          double phys_y = y/nn[1]*length;
          double fmgev = 5.068;
          double mr_physical = sqrt(phys_x*phys_x + phys_y*phys_y)*m*fmgev;
          // r is in fm, m is in GeV, multiply by 5!
          
          // Lattice units
          // Here x is [-N/2, N/2]
          double lat_x = sin(M_PI*x/nn[0])/(M_PI);
          double lat_y = sin(M_PI*y/nn[1])/(M_PI);
          // lat_x and lat_y are now in [-1/2,1/2] as x/nn[0] is in [-N/2, N/2]
          double lat_r = sqrt(lat_x*lat_x + lat_y*lat_y) * nn[0];
          // lat_r now tells how many lattice units the distance is
          double a = length / nn[0];
          double lat_m = m * a * fmgev;
          double bessel_argument = lat_m * lat_r;
          
          
          gsl_sf_result bes;
          gsl_set_error_handler_off(); // so that when we have too large
          // argument, we can handle error ourselves and set regulator to 0
          int status = gsl_sf_bessel_K1_e(bessel_argument, &bes);
          if (status)
          {
            // Too large argument, so x*BesselK[1,x]=0
            mass_regulator = 0.0;
          }
          else
            mass_regulator = bessel_argument * bes.val;

        }
        if (param->getRunningCoupling() == 0)
        {
          // 	  K[pos]->push_back(x/r2);
          // 		  K[pos]->push_back(y/r2);`
          // 		  // S is a 1d vector
          // 		  S[pos]->push_back(1./r2);
          x/=nn[0];
          y/=nn[1];
          // discretization without singularities
          double tmpk1 =cos(Pi*y)*(sin(2.*Pi*x)/(2.*Pi))/((pow( sin(Pi*x)/Pi ,2.) + pow( sin(Pi*y)/Pi ,2.)))/nn[0];
          double tmpk2 = cos(Pi*x)*(sin(2.*Pi*y)/(2.*Pi))/((pow( sin(Pi*x)/Pi ,2.) + pow( sin(Pi*y)/Pi ,2.)))/nn[0];
          
          
          // Regulate long distance tails, does nothing if m=0
          tmpk1 *= mass_regulator;
          tmpk2 *= mass_regulator;
        
          
          K[pos]->push_back(tmpk1);
          K[pos]->push_back(tmpk2);
          // S is a 1d vector
          S[pos]->push_back((pow( cos(Pi*y) ,2.)*pow( sin(2.*Pi*x)/(2.*Pi) ,2.)+pow( cos(Pi*x) ,2.)*pow( sin(2.*Pi*y)/(2.*Pi) ,2.))
                            /pow( (pow( sin(Pi*x)/Pi ,2.) + pow( sin(Pi*y)/Pi ,2.)) ,2.)/nn[0]/nn[0] * mass_regulator * mass_regulator);
          
        }
        else if (param->getRunningCoupling() == 1)
        {
          double c=0.2;
          double length = param->getL();
          double phys_x = x/nn[0]*length; //in fm
          double phys_y = y/nn[1]*length;
          double fmgev = 5.068;
          double phys_r2 = phys_x*phys_x + phys_y*phys_y;
          
          // Alphas in physical units! Lambda2 is lambda_QCD^2 in GeV
          
          alphas = 4.*param->PI
          /((11.0*param->getNc()-2.0*Nf)/3.*log(pow((pow(mu0*mu0/Lambda2,1./c)+pow(4./(phys_r2*Lambda2*fmgev*fmgev),1./c)),c)));
          
          // K[pos]->push_back(sqrt(alphas)*x/r2);
          // 		  K[pos]->push_back(sqrt(alphas)*y/r2);
          // 		  // S is a 1d vector
          // 		  S[pos]->push_back(alphas/r2);
          
          // discretization without singularities
          x/=nn[0];
          y/=nn[1];
          K[pos]->push_back(sqrt(alphas)*(cos(Pi*y)*(sin(2.*Pi*x)/(2.*Pi))/((pow( sin(Pi*x)/Pi ,2.) + pow( sin(Pi*y)/Pi ,2.))))/nn[0] * mass_regulator);
          K[pos]->push_back(sqrt(alphas)*(cos(Pi*x)*(sin(2.*Pi*y)/(2.*Pi))/((pow( sin(Pi*x)/Pi ,2.) + pow( sin(Pi*y)/Pi ,2.))))/nn[0] * mass_regulator);
          // S is a 1d vector
          S[pos]->push_back(alphas*(pow( cos(Pi*y) ,2.)*pow( sin(2.*Pi*x)/(2.*Pi) ,2.)+pow( cos(Pi*x) ,4.)*pow( sin(2.*Pi*y)/(2.*Pi) ,2.))
                            /pow( (pow( sin(Pi*x)/Pi ,2.) + pow( sin(Pi*y)/Pi ,2.)) ,2.)/nn[0]/nn[0]);
          
          
          // 		  if(abs(x)<0.1 && abs(y)<0.1)
          // 		    {
          // 		      cout << "x=" << x << ", y=" << y << endl;
          // 		      cout << "S1=" << (pow( cos(Pi*y) ,4.)*pow( sin(2.*Pi*x)/(2.*Pi) ,2.)+pow( cos(Pi*x) ,4.)*pow( sin(2.*Pi*y)/(2.*Pi) ,2.))
          // 			/pow( (pow( sin(Pi*x)/Pi ,2.) + pow( sin(Pi*y)/Pi ,2.)) ,2.)/nn[0]/nn[0]
          // 			   << "S2=" << (pow( cos(Pi*y) ,2.)*pow( sin(2.*Pi*x)/(2.*Pi) ,2.)+pow( cos(Pi*x) ,2.)*pow( sin(2.*Pi*y)/(2.*Pi) ,2.))
          // 			/pow( (pow( sin(Pi*x)/Pi ,2.) + pow( sin(Pi*y)/Pi ,2.)) ,2.)/nn[0]/nn[0] << ", std S=" << 1./r2 << endl;
          // 		      cout << "Kx1=" << (cos(Pi*y)*cos(Pi*y)*(sin(2.*Pi*x)/(2.*Pi))/((pow( sin(Pi*x)/Pi ,2.) + pow( sin(Pi*y)/Pi ,2.))))/nn[0]
          // 			   << "Kx2=" << (cos(Pi*y)*(sin(2.*Pi*x)/(2.*Pi))/((pow( sin(Pi*x)/Pi ,2.) + pow( sin(Pi*y)/Pi ,2.))))/nn[0]
          // 			   << ", std Kx=" << x/r2*nn[0] << endl;
          // 		      cout << "Ky1=" << (cos(Pi*x)*cos(Pi*x)*(sin(2.*Pi*y)/(2.*Pi))/((pow( sin(Pi*x)/Pi ,2.) + pow( sin(Pi*y)/Pi ,2.))))/nn[0]
          // 			   << "Ky2=" << (cos(Pi*x)*(sin(2.*Pi*y)/(2.*Pi))/((pow( sin(Pi*x)/Pi ,2.) + pow( sin(Pi*y)/Pi ,2.))))/nn[0]
          // 			   << ", std Ky=" << y/r2*nn[0] << endl;
          // 		}
        }
      }
    }
  }
  
		
  cout << "alphas_0=" << 4.*param->PI/((11*param->getNc()-2*Nf)/3.*log(mu0*mu0/Lambda2)) << endl;
		
  fft->fftnVector(K,K,nn,2,1);
  fft->fftnVector(S,S,nn,2,1);
  // now K and S contain the Discrete Fourier transforms of \vec{r}_i/r^2 and 1/r^2 respectively
  
  int Nc = param->getNc();
  int steps = param->getSteps();
  Matrix temp2(Nc,0.);
  Matrix Udag(Nc,0.);
  Matrix UAD(Nc2m1,0);
  Matrix UAt(Nc2m1,0);
  
  fstream foutrho("rhorho.dat",ios::out);
  foutrho.close();
  
  fstream foutf("fourPointDifference.dat",ios::out);
  foutf << " " << endl;
  foutf.close();
  
  fstream foutf4("sixPointSquare.dat",ios::out);
  foutf4 << " " << endl;
  foutf4.close();
  
  fstream foutL6("sixPointLine.dat",ios::out);
  foutL6 << " " << endl;
  foutL6.close();
  
  fstream foutf2("fourPointDifference2D.dat",ios::out);
  foutf2 << " " << endl;
  foutf2.close();
  
  fstream foutL("fourPointLine.dat",ios::out);
  foutL << " " << endl;
  foutL.close();
  
  fstream foutLU("fourPointLineUnequalY.dat",ios::out);
  foutLU << " " << endl;
  foutLU.close();
  
  fstream foutS("fourPointSquare.dat",ios::out);
  foutS << " " << endl;
  foutS.close();
  
  fstream fouth("history.dat",ios::out);
  fouth << " " << endl;
  fouth.close();
  
  fstream fouti("k-corr.dat",ios::out);
  fouti << " " << endl;
  fouti.close();
  
  char outname[25];
  sprintf(outname, "k-corr-unequal-%0.5f",0*param->getMeasureSteps()*5*ds*Pi*Pi);
  fstream fouty0(outname,ios::out);
  fouty0 << " " << endl;
  fouty0.close();
  sprintf(outname, "k-corr-unequal-%0.5f",1*param->getMeasureSteps()*5*ds*Pi*Pi);
  fstream fouty1(outname,ios::out);
  fouty1 << " " << endl;
  fouty1.close();
  sprintf(outname, "k-corr-unequal-%0.5f",2*param->getMeasureSteps()*5*ds*Pi*Pi);
  fstream fouty2(outname,ios::out);
  fouty2 << " " << endl;
  fouty2.close();
  sprintf(outname, "k-corr-unequal-%0.5f",3*param->getMeasureSteps()*5*ds*Pi*Pi);
  fstream fouty3(outname,ios::out);
  fouty3 << " " << endl;
  fouty3.close();
  sprintf(outname, "k-corr-unequal-%0.5f",4*param->getMeasureSteps()*5*ds*Pi*Pi);
  fstream fouty4(outname,ios::out);
  fouty4 << " " << endl;
  fouty4.close();
  sprintf(outname, "k-corr-unequal-%0.5f",5*param->getMeasureSteps()*5*ds*Pi*Pi);
  fstream fouty5(outname,ios::out);
  fouty5 << " " << endl;
  fouty5.close();
  sprintf(outname, "k-corr-unequal-%0.5f",6*param->getMeasureSteps()*5*ds*Pi*Pi);
  fstream fouty6(outname,ios::out);
  fouty6 << " " << endl;
  fouty6.close();
  sprintf(outname, "k-corr-unequal-%0.5f",7*param->getMeasureSteps()*5*ds*Pi*Pi);
  fstream fouty7(outname,ios::out);
  fouty7 << " " << endl;
  fouty7.close();
  sprintf(outname, "k-corr-unequal-%0.5f",8*param->getMeasureSteps()*5*ds*Pi*Pi);
  fstream fouty8(outname,ios::out);
  fouty8 << " " << endl;
  fouty8.close();
  sprintf(outname, "k-corr-unequal-%0.5f",9*param->getMeasureSteps()*5*ds*Pi*Pi);
  fstream fouty9(outname,ios::out);
  fouty9 << " " << endl;
  fouty9.close();
  
  fstream fout2d("2dcorr.dat",ios::out);
  fout2d << " " << endl;
  fout2d.close();
  
  fstream fouttr("traceU.dat",ios::out);
  fouttr << " " << endl;
  fouttr.close();
  
  fstream foutia("corr2.dat",ios::out);
  foutia << "> " << param->getSize()/2 << " " << floor(static_cast<double>(steps)/static_cast<double>(param->getMeasureSteps()))+1
	 << " " << param->getDs()*param->getMeasureSteps()<< endl;
  // using size/2 bins to bin the distance between cells; printing result in the beginning (the +1) and then every 'measureSteps' steps
  foutia.close();
  
  fstream fout3("corr-rQs.dat",ios::out);
  fout3 << "> " << param->getSize()/2 << " " << floor(static_cast<double>(steps)/static_cast<double>(param->getMeasureSteps()))+1
	 << " " << param->getDs()*param->getMeasureSteps()<< endl;
  // using size/2 bins to bin the distance between cells; printing result in the beginning (the +1) and then every 'measureSteps' steps
  fout3.close();
  
  fstream foutq("Qs-coord.dat",ios::out);
  foutq << " " << endl;
  foutq.close();
  
  fstream foutql("QsQL.dat",ios::out);
  foutql << " " << endl;
  foutql.close();
  
  fstream foutqs("QsQS.dat",ios::out);
  foutqs << " " << endl;
  foutqs.close();
  
  fstream foutd("dionysis.dat",ios::out);
  foutd << " " << endl;
  foutd.close();
  
  complex<double> temp;
  complex<double> I(0.,1.);
  //  complex<double> omega[Nc2m1];
  Matrix M(param->getNc(),0.);
  
  // measure quantities
  //cout << "measuring ... " << endl;
  //measure->dipoleOperator(param,lat,0);
  //measure->TwoDCorrelator(param, lat, 0);
  //measure->fourPointFunction(param, lat, 0);
  //cout << "Measure: twoPointFunctionInK..." << endl;
  //measure->twoPointFunctionInK(param,lat,0);
  
  //cout << "storing U in Uy[" << 0 << "] ... " << endl;
  //measure->storeU(param, lat, 0);
  
  //cout << "Measure: fourPointFunctionInK..." << endl;
  //measure->fourPointFunctionInK(param,lat,0);
  
  
  //cout << "Measure: fourPointFunctionLine" << endl;
  //measure->fourPointFunctionLine(param, lat, 0); //Q
  //   cout << "Measure: fourPointFunctionLineUnequalY" << endl;
  //   measure->fourPointFunctionLineUnequalY(param, lat, 0); //Q
  // cout << "Measure: fourPointFunctionSquare" << endl;
  // measure->fourPointFunctionSquare(param, lat, 0); //Q
  // cout << "Measure: sixPointFunctionLine" << endl;
  // measure->sixPointFunctionLine(param, lat, 0); //S_6
  //  cout << "Measure: sixPointFunctionSquare" << endl;
  //  measure->sixPointFunctionSquare(param, lat, 0); //S_6
  //   cout << "Measure: Dionysis' stuff" << endl;
  //   measure->twoPointFunctions(param, lat, 0);
  cout << "done." << endl;
  
  
  // ---------------------------------------------------------------------------------------------------------------------------------------
  // begin evolution
  // here the loop over steps s begins
  cout << "Beginning evolution ..." << endl;
  
  for(int ids=1; ids<=steps; ids++)
  {
    //cout << "step " << ids << endl;
    for (int i=0; i<cells; i++)
    {
      // H.M. disable
      //lat->cells[i]->computeAdjointU();
      //*UA2[i] = lat->cells[i]->getUA();
      
      // generate random Gaussian noise in every cell for Nc^2-1 color components and 2 spatial components x and y
      for (int n=0; n<Nc2m1*2; n++)
      {
        xi2[i][n]=complex<double>(random->Gauss(),0.); // real xi
      }
    }
    //cout << "xi set" << endl;
    
    // FFT of xi
    fft->fftnArray(xi2,xi,nn,2,1,2*Nc2m1);
    // the local xi now contains the Fourier transform of xi, while the original xi is stored in the array xi2
    //cout << "xi fft done" << endl;
    
    // now compute C(K_i,xi_i^a) == F^{-1}(F(K_i)F(xi_i^a)) = F^{-1}(F(K_x)F(xi_x^a)+F(K_y)F(xi_y^a))
    for (int i=0; i<cells; i++)
    {
      for (int n=0; n<Nc2m1; n++)
      {
        CKxi[i][n]=((*K[i]).at(0)*xi[i][n]+(*K[i]).at(1)*xi[i][Nc2m1+n]); // product of x components + product of y components
      }
      
      // now compute the product \tilde(U)^{ab} xi_i^b
      // replace local xi by \tilde(U)^{ab} xi_i^b (do this to save memory)
      //H.M. disable
      /*
      for (int a=0; a<Nc2m1; a++)
      {
        xi[i][0*Nc2m1+a] = 0.;
        xi[i][1*Nc2m1+a] = 0.;
        for (int b=0; b<Nc2m1; b++)
        {
          // do the product UA (matrix) times xi (vector)
          xi[i][0*Nc2m1+a] += (*UA2[i]).get(a,b)*xi2[i][0*Nc2m1+b];
          xi[i][1*Nc2m1+a] += (*UA2[i]).get(a,b)*xi2[i][1*Nc2m1+b];
        }
      }
       */
    }
    //cout << "CKxi done" << endl;
    
    
    fft->fftnArray(CKxi,CKxi,nn,2,-1,Nc2m1);
    //cout << "CKxi fft done" << endl;
    
    // now CKxi contains C(K_i,xi_i^a) - it is a vector with a components
    // checked: is real.
    
    // checked: \tilde(U)^{ab} xi_i^b (in xi) is real.
    //cout << "Uxi set" << endl;
    
    // FFT of UA*xi
    // H.M. disable
    //fft->fftnArray(xi,xi,nn,2,1,2*Nc2m1); // remember xi contains UA times xi now
    //cout << "Uxi fft done" << endl;
    
    // now compute C(K_i,U^{ab} xi_i^a) == F^{-1}(F(K_i)F(U^{ab} xi_i^a)) = F^{-1}(F(K_x)F(U^{ab} xi_x^a)+(F(K_y)F(U^{ab} xi_y^a))
    
    // H.M. disable
    /*
    for (int i=0; i<cells; i++)
    {
      for (int n=0; n<Nc2m1; n++)
      {
        CKUxi[i][n]=((*K[i])[0]*xi[i][n]+(*K[i])[1]*xi[i][Nc2m1+n]); // product of x components + product of y components
      }
    }
    //cout << "CKUxi set" << endl;
    
    fft->fftnArray(CKUxi,CKUxi,nn,2,-1,Nc2m1);
     */
    //cout << "CKUxi fft done" << endl;
    
    // now CKUxi contains C(K_i,U^{ab} xi_i^a) - it is a vector with a components
    // checked: CKUxi is real
    
    // now take \tilde{U}^\dag times CKUxi
    // replace local xi by \tilde(U^\dag)^{ac} \tilde(U)^{cb} xi_i^b (do this to save memory)
    
    // H.M.
    // Disable: keep CKxi[i][a] as  C(K_i,xi_i^a)
    /*
    for (int i=0; i<cells; i++)
    {
      UAD = *UA2[i];
      UAD.conjg();
      for (int a=0; a<Nc2m1; a++)
      {
        xi[i][a] = 0.;
        for (int b=0; b<Nc2m1; b++)
        {
          // do the product UA (matrix) times xi (vector)
          xi[i][a] += UAD.get(a,b)*CKUxi[i][b];
        }
        // now xi contains the vector (in color) U^dag times CKUxi
        // replace CKxi by CKxi - UCKUxi ( the first [ ] that multiplies sqrt(\delta s) )
        CKxi[i][a] -= xi[i][a];
      }
    }
     */
    //cout << "Udag CKUxi set" << endl;
    // C(K_i, V xi_i^a V^\dagger) == F^{-1}(F(K_i)F(xi_i^a)) = F^{-1}(F(K_x)F(xi_x^a)+F(K_y)F(xi_y^a))
    
    
    
    
    //cout << "omega1 done" << endl;
    
    // ------------------------------------------------------------------------------
    // now do the second part in \omega^a
    
    // compute Fourier transform of \tilde{U}^{ab}
    
    // compute Fourier transform of UA
    // H.M. disable
    //fft->fftn(UA2,UA,nn,2,1);
    
    //cout << "UA fft done" << endl;
    
    // multiply by S
    // H.M. disable
    /*
    for (int i=0; i<cells; i++)
    {
      *UA[i] = (*S[i])[0]*(*UA[i]);
    }
    
    // FFT back
    fft->fftn(UA,UA,nn,2,-1);
     */
    // now UA contains the convolution of S with UA
    
    //cout << "SUA fft done" << endl;
    
    // take matrix product with U^\dag
    // H.M. disable
    /*
    for (int i=0; i<cells; i++)
    {
      // do the product UAD (matrix) times SUA (matrix)
      UAD = *UA2[i];
      UAD.conjg();
      *UA[i] = UAD*(*UA[i]);
      // now UA contains the product U^\dag C(S,U)
      
      // now take the product of 0.5*i*(\tilde{t}^a)^{bc} ( U^\dag C(S,U) )^{cb}
      // (Hadamard product, same as Tr((\tilde{t}^a)^{bc} ( U^\dag C(S,U) )^{cd}) )
      
      for (int a=0; a<Nc2m1; a++)
      {
        temp=0.;
        for (int b=0; b<Nc2m1; b++)
        {
          for (int c=0; c<Nc2m1; c++)
          {
            if(group->getTA(a).get(b,c)!=0.)
              temp += group->getTA(a).get(b,c)*(*UA[i]).get(c,b);
          }
        }
        temp*=0.5*I; // 0.5*i is in Eq.(31) but missing in Eq.(34) ... must be there
        // use xi as storage
        xi[i][a]=temp;
      }
    }
     */
    
    
    
    ///// New code for evolution without adjoint
    // H.M.
    // H.M. matrix V xsi V^\dagger
    for (int i=0; i<cells; i++)
    {
      *VxsiVx[i] = Matrix(param->getNc(),0);
      *VxsiVy[i] = Matrix(param->getNc(),0);
      //Matrix tmpx(param->getNc(),0);
      //Matrix tmpy(param->getNc(),0);
      for (int a=0; a<Nc2m1; a++)
      {
        *VxsiVx[i] = *VxsiVx[i]  + xi2[i][a]* lat->cells[i]->getU()
            *group->getT(a) * lat->cells[i]->getU().conjg();
       *VxsiVy[i] = *VxsiVy[i] + xi2[i][Nc2m1+a]* lat->cells[i]->getU()
          *group->getT(a) * lat->cells[i]->getU().conjg();
      }
    }
    
    // FFT V xi V, save output to same array
    fft->fftn(VxsiVx,VxsiVx,nn,2,1);
    fft->fftn(VxsiVy,VxsiVy,nn,2,1);
    
    // Compute in Fourier space F(VxsiV) F(K)
    // Note now *K[i].at(0) is x component of FT of K, and at(1) y comp
    // Save to VxsiV to save memory
    for (int i=0; i<cells; i++)
    {
      // Calculatedot prodcut
      *VxsiVx[i] = (*K[i])[0] * (*VxsiVx[i]) +  (*K[i])[1] * (*VxsiVy[i]) ;
    }
    
    // FFT back
    fft->fftn(VxsiVx,VxsiVx,nn,2,-1); // Contains now dot product
    //fft->fftn(VxsiVy,VxsiVy,nn,2,-1);
    
    


    
    // Evolve matrix
    for (int i=0; i<cells; i++)
    {
      
      Matrix left(param->getNc(),0);
      left = -I * sqrt(ds) * (*VxsiVx[i]);
      Matrix right(param->getNc(),0);
      
      
      for (int a=0; a<Nc2m1; a++)
      {
        // CKxi is approximately real, get more stable evolution by taking the real part?
        right = right + real(CKxi[i][a]) * group->getT(a);
      }
      right = I * sqrt(ds) * right;
      
      lat->cells[i]->setU( left.expm() * lat->cells[i]->getU() * right.expm() );
      
    }
    
    
    
    // Compute C(
    
    
    //cout << "omega2 done" << endl;
    
    // now xi holds the second [], i.e. 0.5 i Tr[t^a U^\dag C(S,U)]
    // checked: is real.
    
    // now compute omega^a and evolve the U field in every cell
    // H.M. disable
    /*
    for (int i=0; i<cells; i++)
    {
      for (int a=0; a<param->getNc()*param->getNc(); a++)
      {
        M.set(a,0.);
      }
      for (int a=0; a<Nc2m1; a++)
      {
        //	      cout << sqrt(ds)*CKxi[i][a] - ds*xi[i][a] << endl;
        M += real( sqrt(ds)*CKxi[i][a] - ds*xi[i][a])*group->getT(a);
      }
      M *= I;
      lat->cells[i]->setU(lat->cells[i]->getU()*M.expm());
    }
    */
    
    // here impose infrared regulator on newly updated U-fields
    //      infrared->regulate(lat, group, param, random, ids);
    
    //cout << " done with step" << endl;
    
    // Save intermediate step Wlines
    //if (ids%10==0)
    //{
    //  stringstream fname; fname << "wline_evolution/step_" << ids;
    //  lat->PrintWilsonLines(fname.str());
   // }
    
    // measure quantities
    if ((ids)%param->getMeasureSteps()==0)
    {
      //cout << "measuring ... " << endl;
      stringstream fname;
      vector<string> inputname =split(param->getInputWline(),'/');
      
      fname << param->getOutputDir() << "/" << inputname[inputname.size()-1] << "_steps_" << ids;
      lat->PrintWilsonLines(fname.str());
      //measure->dipoleOperator(param,lat,ids);
      //measure->TwoDCorrelator(param, lat, ids);
      //measure->fourPointFunction(param, lat, ids);
      //cout << "Measure: twoPointFunctionInK..." << endl;
      //measure->twoPointFunctionInK(param,lat,ids);
      /*
      if (ids%(param->getMeasureSteps()*5)==0 && ids/(param->getMeasureSteps()*5)<10)
      {
        cout << "storing U in Uy[" << ids/(param->getMeasureSteps()*5) << "] ... " << endl;
        measure->storeU(param, lat, ids/(param->getMeasureSteps()*5));
        //cout << "Uy=" << lat->cells[25]->getUy(1) << ", U=" << lat->cells[25]->getU()<< endl; //check
      }*/
      //  cout << "Measure: fourPointFunctionInK..." << endl;
      // measure->fourPointFunctionInK(param,lat,ids);
      
      //  cout << "Measure: fourPointFunctionLine" << endl;
      // measure->fourPointFunctionLine(param, lat, ids); //Q
      //  	  cout << "Measure: fourPointFunctionLineUnequalY" << endl;
      //  	  measure->fourPointFunctionLineUnequalY(param, lat, ids); //Q
      //  cout << "Measure: fourPointFunctionSquare" << endl;
      //  measure->fourPointFunctionSquare(param, lat, ids); //Q
      // cout << "Measure: sixPointFunctionLine" << endl;
      // measure->sixPointFunctionLine(param, lat, ids); //S_6
      // cout << "Measure: sixPointFunctionSquare" << endl;
      // measure->sixPointFunctionSquare(param, lat, ids); //S_6
      // 	  cout << "Measure: Dionysis' stuff" << endl;
      // 	  measure->twoPointFunctions(param, lat, 0);
      cout << "done." << endl;
    }
  }
  
  
  // -----------------------------------------------------------------
  
  // -----------------------------------------------------------------
  
  // finalize
  for(int i=0; i<param->getSize()*param->getSize(); i++)
  {
    delete xi[i];
    delete xi2[i];
    delete K[i];
    delete S[i];
    delete CKxi[i];
    //delete CKUxi[i];
    //delete UA[i];
    //delete UA2[i];
    delete VxsiVx[i];
    delete VxsiVy[i];
  }
  
  delete [] xi;
  delete [] xi2;
  delete [] K;
  delete [] S;
  delete [] CKxi;
  //delete [] CKUxi;
  //delete [] UA;
  //delete [] UA2;
  delete [] VxsiVx;
  delete []VxsiVy;
  
  delete measure;
  delete group;
  delete init;
  delete setup;
  delete random;
  delete fft;
  delete param;
  delete lat;
}/* main */


int readInput(Setup *setup, Parameters *param, int argc, char *argv[])
{
  // the first given argument is taken to be the input file name
  // if none is given, that file name is "input"
  cout << "Opening input file ... " << endl;
  bool ic_cli=false;  // true if initial condition file is a cli argument
  string file_name="input";
  if (argc > 1)
  {
    if (string(argv[1])=="-ic")
    {
      ic_cli = true;
      param->setInputWline(argv[2]);
    }
  }
  if (argc>1)
  {
    if (string(argv[1])!="-ic")
    {
      file_name = argv[1];
      cout << "Using file name \"" << file_name << "\"." << endl;
    }
  }
  else
  {
    cout << "No input file name given. Using default \"input\"." << endl;
    file_name = "input";
  }
  
  // read and set all the parameters in the "param" object of class "Parameters"
  cout << "Reading parameters from file ... ";
  param->setMode(setup->IFind(file_name.c_str(),"mode"));
  param->setInitMethod(setup->IFind(file_name.c_str(),"initMethod"));
  param->setRunningCoupling(setup->IFind(file_name.c_str(),"runningCoupling"));
  param->setg2mu(setup->DFind(file_name.c_str(),"g2mua")); // g^2mu times a
  param->setLambdaQCD(setup->DFind(file_name.c_str(),"Lambda_QCD")); // in units of g^2mu
  param->setR(setup->DFind(file_name.c_str(),"R"));
  param->setDs(setup->DFind(file_name.c_str(),"ds"));
  param->setSize(setup->IFind(file_name.c_str(),"size"));
  param->setNc(setup->IFind(file_name.c_str(),"Nc"));
  param->setSeed(setup->IFind(file_name.c_str(),"seed"));
  param->setNy(setup->IFind(file_name.c_str(),"Ny"));
  param->setSteps(setup->IFind(file_name.c_str(),"steps"));
  param->setMeasureSteps(setup->IFind(file_name.c_str(),"measureSteps"));
  param->setMu0(setup->DFind(file_name.c_str(),"mu0"));
  param->setg(setup->DFind(file_name.c_str(),"g"));
  param->setkappa4Factor(setup->DFind(file_name.c_str(),"kappa4Factor"));
  param->setm(setup->DFind(file_name.c_str(),"m"));
  param->setL(setup->DFind(file_name.c_str(),"L"));
  if (!ic_cli)
    param->setInputWline(setup->StringFind(file_name.c_str(), "input_wline"));
  param->setOutputDir(setup->StringFind(file_name.c_str(), "output_dir"));
  
  // write the used parameters into file "usedParameters.dat" as a double check for later
  time_t rawtime;
  time ( &rawtime );
  fstream fout1("usedParameters.dat",ios::out);
  fout1 << "File created on " << ctime(&rawtime) << endl;
  fout1 << "Program run in mode " << param->getMode() << endl;
  fout1 << "Running coupling (0=no, 1=yes) " << param->getRunningCoupling() << endl;
  fout1 << "LambdaQCD " << param->getLambdaQCD() << endl;
  fout1 << "Nc " << param->getNc() << endl;
  fout1 << "size " << param->getSize() << endl;
  fout1 << "ds " << param->getDs() << endl;
  fout1 << "steps " << param->getSteps() << endl;
  fout1 << "measureSteps " << param->getMeasureSteps() << endl;
  fout1 << "initMethod " << param->getInitMethod() << endl;
  if (param->getInitMethod()==10)
    fout1 << "input_filename " << param->getInputWline() << endl;
  fout1 << "R " << param->getR() << endl;
  fout1 << "Ny " << param->getNy() << endl;
  fout1 << "g^2 mu a " << param->getg2mu() << endl;
  fout1 << "g " << param->getg() << endl;
  fout1 << "m " << param->getm() << " GeV" << endl;
  fout1 << "L " << param->getL() << " fm" << endl;
  fout1 << "mu0 " << param->getMu0() << endl;
  fout1 << "kappa4Factor " << param->getkappa4Factor() << endl;
  
  fout1.close();
  cout << "done." << endl;
  return 0;
}
