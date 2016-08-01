// Measure.cpp is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.

#include "Measure.h"

// this returns |Psi_{T,L}(z,Q^2,r)|^2
double Measure::PsiSq(double z, double QSq, double r, int selector, Parameters *param)
{
  double alpha_em = 1./137.035999679;
  double mSq[4], ef[4];
  double eps, K0, K1;
  int Nf=3; // number of flavors to include
  //mSq[0] = 0.00000576; // mass squared in GeV^2 of the up quark
  //mSq[1] = 0.00002304; // down quark
  //mSq[2] = 0.010816;   // strange quark
  mSq[0] = 0.14*0.14; //mass squared in GeV^2 of the up quark
  mSq[1] = 0.14*0.14; // down quark
  mSq[2] = 0.14*0.14; // strange quark
  mSq[3] = 1.27*1.27; // charm quark
  ef[0] = 2./3.; // charged (in units of alpha_em) of the up quark
  ef[1] = -1./3.; // down quark
  ef[2] = -1./3.; // strange quark
  ef[3] = 2./3.; // charm quark
  double f;

  f = 0; 
  
  if (selector==0) // do transverse
    {  
      for (int i=0; i<Nf; i++)
	{
	  eps = sqrt(z*(1-z)*QSq+mSq[i]);
	  K0 = gsl_sf_bessel_K0(eps*r);
	  K1 = gsl_sf_bessel_K1(eps*r);
	  f+=ef[i]*ef[i]*((z*z+(1-z)*(1-z))*eps*eps*K1*K1+mSq[i]*K0*K0);
	}
    }
  else if (selector==1)
    {
      for (int i=0; i<Nf; i++)
	{
	  eps = sqrt(z*(1-z)*QSq+mSq[i]);
	  K0 = gsl_sf_bessel_K0(eps*r);
	  f+=ef[i]*ef[i]*(4.*QSq*z*z*(1-z)*(1-z)*K0*K0);
	}      
    }
  
  f*=6.*alpha_em/(4.*param->PI*param->PI);
  
  return f;
}

// returns the integral of |Psi_{L,T}|^2 over z
double Measure::Chi(double (Measure::*f)(double, double, double, int, Parameters *param), double QSq, double r, int selector, Parameters *param)
{
  double PsiSq1, PsiSq2, PsiSq3, I;
  double deltaZ, z;
  int steps = 1000;
  deltaZ = 1./static_cast<double>(steps); // the range of integration is 0 to 1

  I = 0.;
  // using Simpson rule - can try something faster lateron ...
  PsiSq1 = (this->*f)((0)*deltaZ, QSq, r, selector, param);
  PsiSq2 = (this->*f)((0)*deltaZ+deltaZ*0.5, QSq, r, selector, param);
  PsiSq3 = (this->*f)(1*deltaZ, QSq, r, selector, param);
  I+=deltaZ/6.*(PsiSq1+4.*PsiSq2+PsiSq3);

  for (int i=2; i<=steps; i++)
    {
      PsiSq1 = PsiSq3;
      PsiSq2 = (this->*f)((i-1)*deltaZ+deltaZ*0.5, QSq, r, selector, param);
      PsiSq3 = (this->*f)(i*deltaZ, QSq, r, selector, param);
      I+=deltaZ/6.*(PsiSq1+4.*PsiSq2+PsiSq3);
  }
  return I;
}

void Measure::computeF2(double ***C, double ds, int xsteps, Parameters *param) // routine to analyze data that is read from disk  
{
  cout << "computing F_2" << endl;
  cout << "|Psi_T|^2 =" << PsiSq(0.5,1,0.2,0,param) << endl;
  cout << "|Psi_L|^2 =" << PsiSq(0.5,1,0.2,1,param) << endl;
  cout << "Chi_T=" << Chi(&Measure::PsiSq,1,0.2,0,param) << endl;
  cout << "Chi_L=" << Chi(&Measure::PsiSq,1,0.2,1,param) << endl;
  double alpha_em = 1./137.035999679;
  double hbarc = param->hbarc;
  double sigma0=2.4/hbarc/hbarc;
  double I, ChiT1, ChiT2, ChiL1, ChiL2, I2;
  int steps = param->getSize()/2;
  double deltaR = param->getg2mu()/0.56; // here put a variable that determines g^2mu in GeV (results in a[1/GeV])
  // this is delta r in GeV^-1
  double QSq;
  double F2, FL;

  QSq=45.;

  ofstream fout("F2-Q45.dat",ios::out); 
 
  for (int j=0; j<xsteps; j++)
    {
      I = 0.;
      I2 = 0.;
      ChiT1 = 0.;
      ChiT2 = Chi(&Measure::PsiSq, QSq, deltaR, 0, param)*deltaR*(*C[j][1]);
      I+=deltaR/2.*(ChiT1+ChiT2);
      ChiL1 = 0.;
      ChiL2 = Chi(&Measure::PsiSq, QSq, deltaR, 1, param)*deltaR*(*C[j][1]);
      I2+=deltaR/2.*(ChiL1+ChiL2);
      
      for (int i=2; i<steps; i++)
	{
	  ChiT1 = ChiT2;
	  ChiT2 = Chi(&Measure::PsiSq, QSq, i*deltaR, 0, param)*i*deltaR*(*C[j][i]);
	  I+=deltaR/2.*(ChiT1+ChiT2);
	  ChiL1 = ChiL2;
	  ChiL2 = Chi(&Measure::PsiSq, QSq, i*deltaR, 1, param)*i*deltaR*(*C[j][i]);
	  I2+=deltaR/2.*(ChiL1+ChiL2);
	  //cout << I << " " << I2 << endl;
	}
      
      I*=2*param->PI*sigma0;
      I2*=2*param->PI*sigma0;
      
      cout << "sigma_T=" << I << endl;
      cout << "sigma_L=" << I2 << endl;
      
      F2 = QSq/(4.*param->PI*param->PI*alpha_em)*(I+I2);
      FL = QSq/(4.*param->PI*param->PI*alpha_em)*(I2);
      
      cout << "F_2(Q^2=" << QSq << " GeV^2, x=" << exp(log(0.01)-ds*j*param->PI*param->PI) << ")="  << F2 << endl;
      cout << "F_L(Q^2=" << QSq << " GeV^2, x=" << exp(log(0.01)-ds*j*param->PI*param->PI) << ")="  << FL << endl;
      fout << QSq << " " << exp(log(0.01)-ds*j*param->PI*param->PI) << " " << F2 << endl; 
    }

  fout.close();

  QSq=6.5;

  ofstream fout2("F2-Q6.5.dat",ios::out); 
 
  for (int j=0; j<xsteps; j++)
    {
      I = 0.;
      I2 = 0.;
      ChiT1 = 0.;
      ChiT2 = Chi(&Measure::PsiSq, QSq, deltaR, 0, param)*deltaR*(*C[j][1]);
      I+=deltaR/2.*(ChiT1+ChiT2);
      ChiL1 = 0.;
      ChiL2 = Chi(&Measure::PsiSq, QSq, deltaR, 1, param)*deltaR*(*C[j][1]);
      I2+=deltaR/2.*(ChiL1+ChiL2);
      
      for (int i=2; i<steps; i++)
	{
	  ChiT1 = ChiT2;
	  ChiT2 = Chi(&Measure::PsiSq, QSq, i*deltaR, 0, param)*i*deltaR*(*C[j][i]);
	  I+=deltaR/2.*(ChiT1+ChiT2);
	  ChiL1 = ChiL2;
	  ChiL2 = Chi(&Measure::PsiSq, QSq, i*deltaR, 1, param)*i*deltaR*(*C[j][i]);
	  I2+=deltaR/2.*(ChiL1+ChiL2);
	  //cout << I << " " << I2 << endl;
	}
      
      I*=2*param->PI*sigma0;
      I2*=2*param->PI*sigma0;
      
      cout << "sigma_T=" << I << endl;
      cout << "sigma_L=" << I2 << endl;
      
      F2 = QSq/(4.*param->PI*param->PI*alpha_em)*(I+I2);
      FL = QSq/(4.*param->PI*param->PI*alpha_em)*(I2);
      
      cout << "F_2(Q^2=" << QSq << " GeV^2, x=" << exp(log(0.01)-ds*j*param->PI*param->PI) << ")="  << F2 << endl;
      cout << "F_L(Q^2=" << QSq << " GeV^2, x=" << exp(log(0.01)-ds*j*param->PI*param->PI) << ")="  << FL << endl;
      fout2 << QSq << " " << exp(log(0.01)-ds*j*param->PI*param->PI) << " " << F2 << endl; 
    }

  fout2.close();

  QSq=0.85;

  ofstream fout3("F2-Q0.85.dat",ios::out); 
 
  for (int j=0; j<xsteps; j++)
    {
      I = 0.;
      I2 = 0.;
      ChiT1 = 0.;
      ChiT2 = Chi(&Measure::PsiSq, QSq, deltaR, 0, param)*deltaR*(*C[j][1]);
      I+=deltaR/2.*(ChiT1+ChiT2);
      ChiL1 = 0.;
      ChiL2 = Chi(&Measure::PsiSq, QSq, deltaR, 1, param)*deltaR*(*C[j][1]);
      I2+=deltaR/2.*(ChiL1+ChiL2);
      
      for (int i=2; i<steps; i++)
	{
	  ChiT1 = ChiT2;
	  ChiT2 = Chi(&Measure::PsiSq, QSq, i*deltaR, 0, param)*i*deltaR*(*C[j][i]);
	  I+=deltaR/2.*(ChiT1+ChiT2);
	  ChiL1 = ChiL2;
	  ChiL2 = Chi(&Measure::PsiSq, QSq, i*deltaR, 1, param)*i*deltaR*(*C[j][i]);
	  I2+=deltaR/2.*(ChiL1+ChiL2);
	  //cout << I << " " << I2 << endl;
	}
      
      I*=2*param->PI*sigma0;
      I2*=2*param->PI*sigma0;
      
      cout << "sigma_T=" << I << endl;
      cout << "sigma_L=" << I2 << endl;
      
      F2 = QSq/(4.*param->PI*param->PI*alpha_em)*(I+I2);
      FL = QSq/(4.*param->PI*param->PI*alpha_em)*(I2);
      
      cout << "F_2(Q^2=" << QSq << " GeV^2, x=" << exp(log(0.01)-ds*j*param->PI*param->PI) << ")="  << F2 << endl;
      cout << "F_L(Q^2=" << QSq << " GeV^2, x=" << exp(log(0.01)-ds*j*param->PI*param->PI) << ")="  << FL << endl;
      fout3 << QSq << " " << exp(log(0.01)-ds*j*param->PI*param->PI) << " " << F2 << endl; 
    }

  fout3.close();


}

void Measure::analyze(Parameters *param) // routine to analyze data that is read from disk  
{
  cout << "Analysis mode" << endl;
  
  ifstream fin("corr2.dat"); 
  ifstream finS6("sixPointLine.dat"); 
  int bins, steps;
  double ds;
  char cdummy[2]; // read the ">" into this

  // read parameters from first line of the data file
  fin >> cdummy;
  fin >> bins;
  fin >> steps;
  fin >> ds;
  cout << "bins=" << bins << ", steps=" << steps << ", ds=" << ds << endl;
  // got parameters.

  double ***C;
  double dummy;
  C = new double **[steps]; 
  double rAverage[steps][bins];
  
  for (int s=0; s<steps; s++)
    {
      C[s] = new double *[bins];
      for (int i=0; i<bins; i++)
	{
	  C[s][i] = new double;
	  if (!fin.eof())
	    {  
	      fin >> rAverage[s][i];
	      fin >> *C[s][i];
	      fin >> dummy;
	    }
	  else 
	    {
	      cerr << "End of file reached prematurely. Was the job finished? Exiting." << endl;
	      exit(1);
	    }
	}
    }

  cout << endl;

  // Determine Qs
  // if correlator is 1-e^(-1/2), define that |x-y| as sqrt(2)/Q_s
  double Qs=0.;
  double R1=0.;
  double R2=0.;
  double r, Rs;
  double lambda, Qsp, Qspp;
  ofstream foutq("Qs-coord.dat",ios::out); 
    
  // compute Qs from coordinate space correlator:
  for(int s=0; s<steps; s++)
    {
      Qs=0;
      for (int j=1; j<bins; j++)
	{
	  if (Nc==2)
	    {
	      if(*C[s][j] > 1.-exp(-1./2.) && Qs==0.)
		{
		  R2 = rAverage[s][j];
		  R1 = rAverage[s][j-1];
		  r= (*C[s][j]-(1.-exp(-1./2.)))/(*C[s][j]-*C[s][j-1]);
		  Rs = R2*(1-r)+r*R1;
		  Qs = 1./Rs;
		}
	    }
	  else if (Nc==3)
	    {
	      if(*C[s][j] > 1.-exp(-1./2.) && Qs==0.)
		{
		  R2 = rAverage[s][j];
		  R1 = rAverage[s][j-1];
		  r= (*C[s][j]-(1.-exp(-1./2.)))/(*C[s][j]-*C[s][j-1]);
		  Rs = R2*(1-r)+r*R1;
		  Qs = 1./Rs;
		}
	    }
	}
      // compute the evolution speed of Qs: lambda =  2 d log(Qs)/d y (2, because the original definition is with Qs^2)
      if (s==1)
	{
	  lambda = 2.*(Qs-Qsp)/(ds*Pi*Pi)/((Qs+Qsp)/2.); // this is lambda/alpha_s
	}
      else if (s>1)
	{
	  lambda = 2.*(Qs-Qspp)/(2.*ds*Pi*Pi)/(Qsp); // this is lambda/alpha_s
	}
      else
	lambda = 0;

      cout << "alpha_s y=" << s*ds*Pi*Pi << " Qs a=" << sqrt(2.)*Qs << " R1=" << R1 << " R2=" << R2 << " Rs=" << Rs << " lambda=" << lambda 
	   << " Qs-Qsp=" << Qs-Qsp << endl;

      
      Qspp = Qsp; // set previous Qs to Qs
      Qsp = Qs; // set previous Qs to Qs

      foutq << s*ds*Pi*Pi << " " << sqrt(2.)*Qs << " " << (s-1)*ds*Pi*Pi << " " << lambda << " " << floor(R2) << endl; // sqrt(2) because Rs^2=2/Qs^2
      // this is alpha_s*y = alpha_s*1/ln(x) in column 1 and Qs/g^2mu in columns two
    }
  
  foutq.close();

  computeF2(C,ds,steps,param);
      
  // ---- finalize ----------------------------------
      
  for (int s=0; s<steps; s++)
    {
      for(int i=0; i<bins; i++)
	{
	  delete C[s][i];
	}
    }
  
  for (int s=0; s<steps; s++)
    {
      delete [] C[s];
    }
  
  delete [] C;
}

void Measure::analyzeS6(Parameters *param) // routine to analyze data that is read from disk  
{
  cout << "Analysis mode for S6" << endl;
  
  ifstream fin("sixPointLine.dat"); 
  int bins, steps;
  double ds;
  char cdummy[2]; // read the ">" into this

  double dummy;
  
  // read parameters from first line of the data file
  fin >> cdummy;
  fin >> bins;
  fin >> steps;
  fin >> ds;

  bins+=1;
  cout << "S6Line: bins=" << bins << ", steps=" << steps << ", ds=" << ds << endl;
  // got parameters.

  double ***S6;
  S6 = new double **[steps]; 
  double rAverage[steps][bins];
  double rQsAverage[steps][bins];

  for (int s=0; s<steps; s++)
    {
      S6[s] = new double *[bins];
      for (int i=0; i<bins; i++)
	{
	  S6[s][i] = new double;
	  if (!fin.eof())
	    {  
	      fin >> rAverage[s][i];
	      fin >> rQsAverage[s][i];
	      fin >> dummy;
	      fin >> *S6[s][i];
	      fin >> dummy;
	      fin >> dummy;
	      fin >> dummy;
	      fin >> dummy;
	      fin >> dummy;
	      fin >> dummy;
	      fin >> dummy;
	      fin >> dummy;
	      fin >> dummy;
	      fin >> dummy;
	      fin >> dummy;
	    }
	  else 
	    {
	      cerr << "End of file reached prematurely. Was the job finished? Exiting." << endl;
	      exit(1);
	    }
	}
    }

  ofstream fout("S6LrQs0.1.dat",ios::out); 
  double rQs=0.1;
  double frac;
  double S6v[steps];
  double speed;

  for(int s=0; s<steps; s++)
    {
      for (int j=1; j<bins; j++)
	{
	  if ((rQsAverage[s][j]>rQs && rQsAverage[s][j-1]<rQs) || (rQsAverage[s][j]>rQs && j-1==0))
	    {
	      frac=rQsAverage[s][j]-rQs;
	      S6v[s] = (1.-frac)*(*S6[s][j])+frac*(*S6[s][j-1]);
	      cout << S6v[s] << " " << *S6[s][j] << " " << *S6[s][j-1] << endl;
	    }
	}
      if(s>0)
	{
	  speed = -(S6v[s]-S6v[s-1])/(ds*Pi*Pi);
	}
      fout << ds*s*Pi*Pi << " " << speed << endl;
    }

  // Determine Qs
  // if correlator is 1-e^(-1/2), define that |x-y| as sqrt(2)/Q_s
  double Qs=0.;
  double R1=0.;
  double R2=0.;
  double r, Rs;
  double lambda, Qsp, Qspp;
  ofstream foutq("QsSL.dat",ios::out); 
    
  // compute Qs from coordinate space correlator:
  for(int s=0; s<steps; s++)
    {
      Qs=0;
      for (int j=1; j<bins; j++)
	{
	  if (Nc==2)
	    {
	      if(*S6[s][j] < exp(-1./2.) && Qs==0.)
		{
		  R2 = rAverage[s][j];
		  R1 = rAverage[s][j-1];
		  r= (*S6[s][j]-(exp(-1./2.)))/(*S6[s][j]-*S6[s][j-1]);
		  Rs = R2*(1-r)+r*R1;
		  Qs = 1./Rs;
		}
	    }
	  else if (Nc==3)
	    {
	      if(*S6[s][j] < exp(-1./2.) && Qs==0.)
		{
		  R2 = rAverage[s][j];
		  R1 = rAverage[s][j-1];
		  r= (*S6[s][j]-(exp(-1./2.)))/(*S6[s][j]-*S6[s][j-1]);
		  Rs = R2*(1-r)+r*R1;
		  Qs = 1./Rs;
		}
	    }
	}
      // compute the evolution speed of Qs: lambda =  2 d log(Qs)/d y (2, because the original definition is with Qs^2)
      if (s==1)
	{
	  lambda = 2.*(Qs-Qsp)/(ds*Pi*Pi)/((Qs+Qsp)/2.); // this is lambda
	}
      else if (s>1)
	{
	  lambda = 2.*(Qs-Qspp)/(2.*ds*Pi*Pi)/(Qsp); // this is lambda
	}
      else
	lambda = 0;

      cout << "y=" << s*ds*Pi*Pi << " Qs a=" << sqrt(2.)*Qs << " R1=" << R1 << " R2=" << R2 << " Rs=" << Rs << " lambda=" << lambda 
	   << " Qs-Qsp=" << Qs-Qsp << endl;

      
      Qspp = Qsp; // set previous Qs to Qs
      Qsp = Qs; // set previous Qs to Qs

      foutq << s*ds*Pi*Pi << " " << sqrt(2.)*Qs << " " << (s-1)*ds*Pi*Pi << " " << lambda << " " << floor(R2) << endl; // sqrt(2) because Rs^2=2/Qs^2
    }
  
  foutq.close();

      
  // ---- finalize ----------------------------------
      
  for (int s=0; s<steps; s++)
    {
      for(int i=0; i<bins; i++)
	{
	  delete S6[s][i];
	}
    }
  
  for (int s=0; s<steps; s++)
    {
      delete [] S6[s];
    }
  
  delete [] S6;

  exit(1);
}

void Measure::TwoDCorrelator(Parameters *param, Lattice *lat, int ids)
{
  int nn[2];
  nn[0]=size;
  nn[1]=size;
  Matrix UDag(Nc,0);
  int pos;

  Matrix ** U;
  U = new Matrix *[size*size];
  
  for (int i=0; i<size*size; i++)
    {
      U[i] = new Matrix(Nc,0);
      *U[i] = lat->cells[i]->getU();
    }

  
  UDag=*U[(nn[0]/2)*nn[1]+nn[1]/2]; //U at (0,0)
  UDag.conjg();
  
  double x,y;


 ofstream fout2("traceU.dat",ios::app); 
  for (int i=0; i<nn[1]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  x = lat->cells[pos]->getX();
	  y = lat->cells[pos]->getY();
	  if(Nc==2)
	    {
	      fout2 << x*param->getg2mu() << " " << y*param->getg2mu() << " " << 1.-(U[pos]->getRe(0) + U[pos]->getRe(3))/static_cast<double>(Nc) << endl;
	    }
	  else if (Nc==3)
	    {
	      fout2 <<  x*param->getg2mu() << " " << y*param->getg2mu() << " " << 1.-(U[pos]->getRe(0) + U[pos]->getRe(4) + U[pos]->getRe(8))/static_cast<double>(Nc) << endl;
	    }
	}
      fout2 << endl;
    }
  fout2.close();

  ofstream fout("2dcorr.dat",ios::app); 
  for (int i=0; i<nn[1]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  x = lat->cells[pos]->getX();
	  y = lat->cells[pos]->getY();
	  if(Nc==2)
	    {
	      *U[pos]=UDag*(*U[pos]);
	      fout << x*param->getg2mu() << " " << y*param->getg2mu() << " " << (U[pos]->getRe(0) + U[pos]->getRe(3))/static_cast<double>(Nc) << endl;
	    }
	  else if (Nc==3)
	    {
	      *U[pos]=UDag*(*U[pos]);
	      fout <<  x*param->getg2mu() << " " << y*param->getg2mu() << " " << (U[pos]->getRe(0) + U[pos]->getRe(4) + U[pos]->getRe(8))/static_cast<double>(Nc) << endl;
	    }
	}
      fout << endl;
    }
  fout.close();



  for (int i=0; i<size*size; i++)
    {
      delete U[i];
    }
  
  delete [] U; 
}


void Measure::twoPointFunctionInK(Parameters *param, Lattice *lat, int ids)
{
  int bins;
  bins = size/2;
  int nn[2];
  nn[0]=size;
  nn[1]=size;
  Matrix UDag(Nc,0);
  double dNc = static_cast<double>(Nc);

  Matrix ** U;
  U = new Matrix *[size*size];
  
  for (int i=0; i<size*size; i++)
    {
      U[i] = new Matrix(Nc,0);
      *U[i] = lat->cells[i]->getU();
    }
  
  double **C;
  C = new double *[bins]; 

  // correlator tr U U^\dag tr U U^\dag
  double **C4;
  C4 = new double *[bins]; 

  double count[bins]; // counts entries in bin i
 
  for (int i=0; i<bins; i++)
    {
      C[i] = new double;
      *C[i] = 0.;
      C4[i] = new double;
      *C4[i] = 0.;
      count[i]=0.;
    }

  double kRange = 2.*sqrt(2.);
  double step = kRange/static_cast<double>(bins);
  int position;

  fft->fftn(U,U,nn,2,1);

  int pos;
  double kt, kx, ky;
  for (int i=0; i<size*size; i++)
    {
      UDag = *U[i];
      UDag.conjg();
      *U[i] = UDag*(*U[i]); // U is now U^dag U
    }

  double sum = 0.;

  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  // continuum momentum from FFT
	  kx = 2.*Pi*(-0.5+static_cast<double>(i)/static_cast<double>(nn[0]));
	  ky = 2.*Pi*(-0.5+static_cast<double>(j)/static_cast<double>(nn[1]));
	  
	  kt = 2.*sqrt(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); // lattice momentum
	  
	  position = static_cast<int>(floor(kt/step));
	
	  if (position<bins)
	    {
	      if(Nc==2)
		{
		  *C[position]+=(U[pos]->getRe(0) + U[pos]->getRe(3))*kt*kt;
		  //*C[position]+=(U[pos]->getRe(0) + U[pos]->getRe(3))/dNc;
		  *C4[position]+=((U[pos]->getRe(0) + U[pos]->getRe(3))*(U[pos]->getRe(0) + U[pos]->getRe(3))
					-(U[pos]->getIm(0) + U[pos]->getIm(3))*(U[pos]->getIm(0) + U[pos]->getIm(3)))/dNc/dNc;
		  count[position]++;
		}
	      else if (Nc==3)
		{
		  *C[position]+=(U[pos]->getRe(0) + U[pos]->getRe(4) + U[pos]->getRe(8))*kt*kt;
		  //*C[position]+=(U[pos]->getRe(0) + U[pos]->getRe(4) + U[pos]->getRe(8))/dNc;
		  *C4[position]+=((U[pos]->getRe(0) + U[pos]->getRe(4) + U[pos]->getRe(8))*(U[pos]->getRe(0) + U[pos]->getRe(4) + U[pos]->getRe(8))
				  -(U[pos]->getIm(0) + U[pos]->getIm(4) + U[pos]->getIm(8))*(U[pos]->getIm(0) + U[pos]->getIm(4) + U[pos]->getIm(8)))/dNc/dNc;
		  count[position]++;
		}
	    }
	}
    }

  for(int i=0; i<bins; i++)
    {
      *C[i]/=static_cast<double>(count[i])*size*size; // divide by N^2, because sum_k C_k = N^2* N_c (instead of (2\pi)^2 N_c in continuum
      *C4[i]/=static_cast<double>(count[i])*size*size*size*size;
    }

  ofstream fout("k-corr.dat",ios::app); 
  for (int j=0; j<bins; j++)
    {
      fout << (j*step+step/2.)/g2mu << " "  << *C[j]/g2mu/g2mu << " " << *C4[j]/(*C[j])/(*C[j]) - 1. << endl; // scale k by g^2 mu
    }
  fout << endl;
  fout.close();

  // now Fourier transform U^dag U to coordinate space - that will give me the averaged correlator in coordinate space
  fft->fftn(U,U,nn,2,-1);

  double x,y,r;
  double rAverage[bins];
  double rstep=static_cast<double>(size/2)/static_cast<double>(bins);

  for (int i=0; i<bins; i++)
    {
      *C[i] = 0.;
      count[i] = 0;
      rAverage[i] = 0.; // average value of r in each bin
    }

  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  x = lat->cells[pos]->getX();
	  y = lat->cells[pos]->getY();
	  r = sqrt(x*x+y*y); 
	  position = static_cast<int>(floor(r/static_cast<double>(rstep)));
	  if (position<bins)
	    {
	      if(Nc==2)
		{
		  rAverage[position]+=r;
		  *C[position]+=(U[pos]->getRe(0) + U[pos]->getRe(3));
		  count[position]++;
		}
	      else if (Nc==3)
		{
		  rAverage[position]+=r;
		  *C[position]+=(U[pos]->getRe(0) + U[pos]->getRe(4) + U[pos]->getRe(8));
		  count[position]++;
		}
	    }
	}
    }

  ofstream fout2("corr2.dat",ios::app); 
  for (int j=0; j<bins; j++)
    {
      rAverage[j]/=static_cast<double>(count[j]);
      *C[j]/=size*size*count[j];
      fout2 << rAverage[j] << " " << 1.-*C[j]/static_cast<double>(Nc) << " " << 0 << endl;
    }
  fout2 << endl;
  fout2.close();

  // compute Q_s
  // if correlator is 1-e^(-1/2), define that |x-y| as sqrt(2)/Q_s
  double Qs=0.;
  double R1=0.;
  double R2=0.;
  double Rs;
 
  // compute Qs from coordinate space correlator:
  for (int j=1; j<bins; j++)
    {
      if (Nc==2)
	{
	  if(2.-*C[j] > 2.-2.*exp(-1./2.) && Qs==0.)
	    {
	      R2 = rAverage[j];
	      R1 = rAverage[j-1];
	      r=(-*C[j]+2.*exp(-1./2.))/(-*C[j]+*C[j-1]);
	      Qs = r*1./R1+(1-r)*1./R2;
	    }
	}
      else if (Nc==3)
	{
// 	  if(*C[j] < 3.*exp(-1./4.) && Qs==0.)
// 	    {
// 	      R2 = rAverage[j];
// 	      R1 = rAverage[j-1];
// 	      cout << "C(j)=" << *C[j] << " C(j-1)=" << *C[j-1] << endl;
// 	      r=(-*C[j]+3.*exp(-1./4.))/(-*C[j]+*C[j-1]);
// 	      Rs = R2*(1-r)+r*R1;
// 	      Qs = 1./Rs;
// 	    }
	  if(*C[j] < 3.*exp(-1./2.) && Qs==0.)
	    {
	      R2 = rAverage[j];
	      R1 = rAverage[j-1];
	      cout << "C(j)=" << *C[j] << " C(j-1)=" << *C[j-1] << endl;
	      r=(-*C[j]+3.*exp(-1./2.))/(-*C[j]+*C[j-1]);
	      Rs = R2*(1-r)+r*R1;
	      Qs = 1./Rs;
	    }
	}
    }

  Qs*=sqrt(2.);
  // cout << "Qs a=" << Qs << " (j=" << R2 << ") " << (R1+(1-r)*(R2-R1)) << endl;
  cout << "Qs=" << Qs/param->getg2mu() << " g^2 mu" << endl;

  //  cout << "Qs a=" << Qs << " (j=" << R2 << ") " << (R1+(1-r)*(R2-R1)) << endl;
  //cout << "Qs=" << Qs/param->getg2mu() << " g^2 mu" << endl;

  // compute the evolution speed of Qs: lambda =  2 d log(Qs)/d y (2, because the original definition is with Qs^2)
  double lambda, Qsp;
  Qsp = param->getQs();
  
  lambda = 2.*(log(Qs)-log(Qsp))/(param->getMeasureSteps()*ds*Pi*Pi); // this is lambda/alpha_s
  cout << "lambda=" << lambda << endl;
  param->setQs(Qs); //sets Qs

  double alphas;
  double Nf = 3;
  double mu0=param->getMu0();
  double Lambda2 = param->getLambdaQCD()*param->getLambdaQCD();
  double c=0.2;
  alphas = 4.*param->PI
    /((11*param->getNc()-2*Nf)/3.
      *log(pow((pow(mu0*mu0/Lambda2,1./c)+pow((Qs*Qs/(Lambda2*param->getg2mu()*param->getg2mu())),1./c)),c)));
	
  //  alphas = 4.*param->PI/((11*param->getNc()-2*Nf)/3.*log(mu0*mu0+Qs*Qs/(param->getg2mu()*param->getg2mu()*Lambda2)));
  //alphas = 4.*param->PI/((11*param->getNc()-2*Nf)/3.*log(mu0*mu0+4*Qs*Qs/(Lambda2*param->getg2mu()*param->getg2mu())));
  cout << "alphas(Qs)=" << alphas << endl;

  ofstream foutq("Qs-coord.dat",ios::app); 
  cout << ids << " " << ds << " " << Pi << endl;
  foutq << ids*ds*Pi*Pi << " " << Qs << " " << lambda << " " << R2 << endl; // sqrt(2) because Rs^2=2/Qs^2
  // this is y = 1/ln(x) in column 1 and Qs/g^2mu in columns two
  foutq.close();

  ofstream fout3("corr-rQs.dat",ios::app); 
  for (int j=0; j<bins; j++)
    {
      fout3 << rAverage[j]*Qs << " " << 1.-*C[j]/static_cast<double>(Nc) << " " << 0 << endl;
    }
  fout3 << endl;
  fout3.close();

  
  for (int i=0; i<bins; i++)
    {
      delete C[i];
      delete C4[i];
    }

  for (int i=0; i<size*size; i++)
    {
      delete U[i];
    }
  
  delete [] C;
  delete [] C4;
  delete [] U;
}

void Measure::storeU(Parameters *param, Lattice *lat, int a)
{
  int pos;
  int nn[2];
  nn[0]=size;
  nn[1]=size;

  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  lat->cells[pos]->setUy(a,lat->cells[pos]->getU());  //save current U into Uy[i] 
	}
    }
}

void Measure::fourPointFunctionInK(Parameters *param, Lattice *lat, int ids)
{
  int bins;
  bins = size/2;
  int nn[2];
  nn[0]=size;
  nn[1]=size;
  Matrix UDag(Nc,0);
  Matrix UiDag(Nc,0);
  double dNc = static_cast<double>(Nc);

  Matrix ** U;
  U = new Matrix *[size*size];
  Matrix ** Ui;
  Ui = new Matrix *[size*size];
  Matrix *** Uy;
  Uy = new Matrix **[10];


  // also add Q

  //get the U distribution from the lattice
  for (int i=0; i<size*size; i++)
    {
      U[i] = new Matrix(Nc,0);
      *U[i] = lat->cells[i]->getU();
      Ui[i] = new Matrix(Nc,0);
      *Ui[i] = lat->cells[i]->getUi();
    }

  for (int i=0; i<10; i++)
    {
      if (ids>=i*param->getMeasureSteps()*5)
	{
	  Uy[i] = new Matrix *[size*size];

	  for (int j=0; j<size*size; j++)
	    {
	      Uy[i][j] = new Matrix(Nc,0);
	      *Uy[i][j] = lat->cells[j]->getUy(i);
	    }
	}
    }
  
  double **C;
  C = new double *[bins]; 
  double **Ci;
  Ci = new double *[bins]; 
  
  // correlator tr U U^\dag tr U U^\dag
  double **C4;
  C4 = new double *[bins]; 

  double count[bins]; // counts entries in bin i
 
  double kRange = 2.*sqrt(2.);
  double step = kRange/static_cast<double>(bins);
  int position;

  // Fourier transform the U's:
  fft->fftn(U,U,nn,2,1);
  fft->fftn(Ui,Ui,nn,2,1);

  for (int i=0; i<10; i++)
    {
      if (ids>=i*param->getMeasureSteps()*5)
	{
	  fft->fftn(Uy[i],Uy[i],nn,2,1);
	}
    }

      for (int i=0; i<bins; i++)
	{
	  Ci[i] = new double;
	  *Ci[i] = 0.;
	  C[i] = new double;
	  *C[i] = 0.;
	  C4[i] = new double;
	  *C4[i] = 0.;
	  count[i]=0.;
	}

  int pos;
  double kt, kx, ky;
  for (int i=0; i<size*size; i++)
    {
      UDag = *U[i];
      UDag.conjg();
      *U[i] = UDag*(*U[i]); // U is now U^dag U
      UiDag = *Ui[i];
      UiDag.conjg();
      *Ui[i] = UiDag*(*Ui[i]); // Ui is now Ui^dag Ui
      for (int j=0; j<10; j++)
	{
	  if (ids>=j*param->getMeasureSteps()*5)
	    {
	      UiDag = *Uy[j][i];
	      UiDag.conjg();
	      *Uy[j][i] = UiDag*(*Uy[j][i]); // Uy[j] is now Uy[j]^dag Uy[j]
	    }
	}
    }

  double sum = 0.;

  for(int a=0; a<10; a++)
    {
      if(ids>=a*param->getMeasureSteps()*5)
	{
	  for (int i=0; i<bins; i++)
	    {
	      *Ci[i] = 0.;
	      *C[i] = 0.;
	      *C4[i] = 0.;
	      count[i]=0.;
	    }
	  for (int i=0; i<nn[0]; i++)
	    {
	      for (int j=0; j<nn[1]; j++)
		{
		  pos = i*nn[1]+j;
		  // continuum momentum from FFT
		  kx = 2.*Pi*(-0.5+static_cast<double>(i)/static_cast<double>(nn[0]));
		  ky = 2.*Pi*(-0.5+static_cast<double>(j)/static_cast<double>(nn[1]));
		  
		  kt = 2.*sqrt(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); // lattice momentum
		  
		  position = static_cast<int>(floor(kt/step));
		  
		  if (position<bins)
		    {
		      if(Nc==2)
			{
			  *C[position]+=(U[pos]->getRe(0) + U[pos]->getRe(3))/dNc;
			  *Ci[position]+=(Uy[a][pos]->getRe(0) + Uy[a][pos]->getRe(3))/dNc;
			  *C4[position]+=((Uy[a][pos]->getRe(0) + Uy[a][pos]->getRe(3))*(U[pos]->getRe(0) + U[pos]->getRe(3))
					  -(Uy[a][pos]->getIm(0) + Uy[a][pos]->getIm(3))*(U[pos]->getIm(0) + U[pos]->getIm(3)))/dNc/dNc;
			  count[position]++;
			}
		      else if (Nc==3)
			{
			  *C[position]+=(U[pos]->getRe(0) + U[pos]->getRe(4) + U[pos]->getRe(8))/dNc;
			  *Ci[position]+=(Uy[a][pos]->getRe(0) + Uy[a][pos]->getRe(4) + Uy[a][pos]->getRe(8))/dNc;
			  *C4[position]+=((Uy[a][pos]->getRe(0) + Uy[a][pos]->getRe(4) + Uy[a][pos]->getRe(8))
 					  *(U[pos]->getRe(0) + U[pos]->getRe(4) + U[pos]->getRe(8))
 					  -(Uy[a][pos]->getIm(0) + Uy[a][pos]->getIm(4) + Uy[a][pos]->getIm(8))
 					  *(U[pos]->getIm(0) + U[pos]->getIm(4) + U[pos]->getIm(8)))/dNc/dNc;
			  count[position]++;
			}
		    }
		}
	    }
	  
	  for(int i=0; i<bins; i++)
	    {
	      *C[i]/=static_cast<double>(count[i]*size*size); // divide by N^2, because sum_k C_k = N^2* N_c (instead of (2\pi)^2 N_c in continuum
	      *Ci[i]/=static_cast<double>(count[i]*size*size); // divide by N^2, because sum_k C_k = N^2* N_c (instead of (2\pi)^2 N_c in continuum
	      *C4[i]/=static_cast<double>(count[i]*size*size*size*size);
	    }
	  
	  double Qs = param->getQs();
	  char outname[25];
	  sprintf(outname, "k-corr-unequal-%0.5f",a*param->getMeasureSteps()*5*ds*Pi*Pi);
	  ofstream fout(outname,ios::app); 
	  for (int j=0; j<bins; j++)
	    {
	      fout << (j*step+step/2.)/Qs << " "  << *C[j] << " " << *C4[j]/(*Ci[j])/(*C[j]) - 1. << " " << *C4[j] << " " << *Ci[j] << endl; // scale k by g^2 mu
	    }
	  fout << endl;
	  fout.close();
	
	} // if a
    } // a loop
  
  for (int i=0; i<bins; i++)
    {
      delete C[i];
      delete Ci[i];
      delete C4[i];
    }

  for (int i=0; i<size*size; i++)
    {
      delete U[i];
      delete Ui[i];
    }
  
  delete [] C;
  delete [] Ci;
  delete [] C4;
  delete [] U;
  delete [] Ui;

  for (int i=0; i<10; i++)
    {
      if (ids>=i*param->getMeasureSteps()*5)
	{
	  for (int j=0; j<size*size; j++)
	    {
	      delete Uy[i][j];
	    }
	  delete [] Uy[i];
	}
    }
  delete [] Uy;
}


void Measure::dipoleOperator(Parameters *param, Lattice *lat, int ids)
{
  for (int j=0; j<size; j++)
    {
      *tempC[j] = 0.; 
    }
  Matrix Udag(Nc,0.);
  Matrix temp2(Nc,0.);
 
  int pos,pos2,pos3,counts=0;
  for (int i=0; i<size; i++)
    {  
      for (int k=0; k<size; k++)
	{
	  for (int j=0; j<size; j++)
	    {
	      pos = i*size+j;
	      pos2 = i*size+k;
	      Udag = lat->cells[pos]->getU();
	      Udag.conjg(); // now Udag is U^dagger
	      temp2 = Udag * lat->cells[pos2]->getU();
	      //cout << "U=" << lat->cells[i*size]->getU() << endl;
	      //cout << "Udag=" << Udag << endl;
	      pos3 = abs(j-k); 
	      if (pos3==0) counts++;
	      if(pos3<size/2)
		{
		  *tempC[pos3] = *tempC[pos3]+temp2;
		}
	      else 
		{
		  *tempC[size-pos3] = *tempC[size-pos3]+temp2;
		}
	      //cout << *tempC[j] << endl;

	      //	  fouta << i << " " << setw(12) << lat->cells[pos]->getU().getRe(0) << endl;
	    }
	}
    }
  
  for (int i=0; i<size; i++)
    {  
      for (int k=0; k<size; k++)
	{
	  for (int j=0; j<size; j++)
	    {
	      pos = j*size+i;
	      pos2 = k*size+i;
	      Udag = lat->cells[pos]->getU();
	      Udag.conjg(); // now Udag is U^dagger
	      temp2 = Udag * lat->cells[pos2]->getU();
	      //cout << "U=" << lat->cells[i*size]->getU() << endl;
	      //cout << "Udag=" << Udag << endl;
	      pos3 = abs(j-k); 
	      if (pos3==0) counts++;
	      if(pos3<size/2)
		{
		  *tempC[pos3] = *tempC[pos3]+temp2;
		}
	      else 
		{
		  *tempC[size-pos3] = *tempC[size-pos3]+temp2;
		}
	      //cout << *tempC[j] << endl;
	      //	  fouta << i << " " << setw(12) << lat->cells[pos]->getU().getRe(0) << endl;
	    }
	}
    }
  
  for (int j=0; j<=size/2; j++)
    {
      *tempC[j]/=static_cast<double>(counts); 
      if (j>0)
	*tempC[j]/=2.; // because for all distances except 0 I double count. 
    }
  

  ofstream foutg("history.dat",ios::app); 
  for (int j=0; j<=size/2; j++)
    {
      if (Nc==2)
	{
	  foutg << j << " "  << tempC[j]->getRe(0) << " " 
		<< (2.-(tempC[j]->getRe(0) + tempC[j]->getRe(3)))/2. << " " << (-(tempC[j]->getIm(0) + tempC[j]->getIm(3)))/2. << endl;
	}
      else if (Nc==3)
	{
	  foutg << j << " "  << tempC[j]->getRe(0) << " " 
		<< (3.-(tempC[j]->getRe(0) + tempC[j]->getRe(4) + tempC[j]->getRe(8)))/3.
		<< " " << (-(tempC[j]->getIm(0) + tempC[j]->getIm(4) + tempC[j]->getIm(8)))/3. << endl;
	}
    }
  foutg << endl;
  foutg.close();
  
}  

void Measure::fourPointFunction(Parameters *param, Lattice *lat, int ids)
{
  for (int j=0; j<size; j++)
    {
      *tempC[j] = 0.; 
    }
  Matrix Udag(Nc,0.);
  Matrix temp2(Nc,0.);
  Matrix temp3(Nc,0.);
  double A[size][size];
  double B[size][size];
  double C[size][size];
  int counts[size][size];
  
  for (int j=0; j<size; j++)
    {  
      A[j][j]=0;
      B[j][j]=0;
      C[j][j]=0;
      counts[j][j]=0;
    }

  int pos,pos2,pos3,posx,posy,posz;
  for (int j=0; j<size; j++)
    {  
      for (int l=0; l<size; l++)
	{
	  for (int m=0; m<size; m++)
	    {
	      if (m+j<size)
		{
		  posx = size*l+m+j;
		}
	      else
		{
		  posx = size*l+(m+j-size);
		}
	      
	      posy = size*l+m;
	      
	      if (l+j<size)
		{
		  posz = size*(l+j)+m;
		}
	      else
		{
		  posz = size*(l+j-size)+m;
		}
	      
	      Udag = lat->cells[posx]->getU();
	      Udag.conjg(); // now Udag is U^dagger
	      temp2 = Udag * lat->cells[posy]->getU();
	      
	      Udag = lat->cells[posy]->getU();
	      Udag.conjg(); // now Udag is U^dagger
	      temp3 = Udag * lat->cells[posz]->getU();
	      //cout << "U=" << lat->cells[i*size]->getU() << endl;
	      //cout << "Udag=" << Udag << endl;
	      pos = abs(j); 
	      pos2 = abs(j);
	      if(pos2<size/2)
		{
		  if(pos<size/2)
		    {
		      if (Nc==2)
			{
			  A[pos][pos2]+=(2.- temp2.getRe(0) - temp2.getRe(3))*(2.-temp3.getRe(0) - temp3.getRe(3));
			  B[pos][pos2]+=(2.-temp2.getRe(0) - temp2.getRe(3));
			  C[pos][pos2]+=(2.-temp3.getRe(0) - temp3.getRe(3));
			}
		      else if (Nc==3)
			{
			  A[pos][pos2]+=(3.-temp2.getRe(0) - temp2.getRe(4) - temp2.getRe(8))*(3.-temp3.getRe(0) - temp3.getRe(4) - temp3.getRe(8));
			  B[pos][pos2]+=(3.-temp2.getRe(0) - temp2.getRe(4) - temp2.getRe(8));
			  C[pos][pos2]+=(3.-temp3.getRe(0) - temp3.getRe(4) - temp3.getRe(8));
			}
		      counts[pos][pos2]++;
		    }
		  else if (pos>=size/2)
		    {
		      if (Nc==2)
			{
			  A[size-pos][pos2]+=(2.-temp2.getRe(0) - temp2.getRe(3))*(2.-temp3.getRe(0) - temp3.getRe(3));
			  B[size-pos][pos2]+=(2.-temp2.getRe(0) - temp2.getRe(3));
			  C[size-pos][pos2]+=(2.-temp3.getRe(0) - temp3.getRe(3));
			}
		      else if (Nc==3)
			{
			  A[size-pos][pos2]+=(3.-temp2.getRe(0) - temp2.getRe(4) - temp2.getRe(8))*(3.-temp3.getRe(0) - temp3.getRe(4) - temp3.getRe(8));
			  B[size-pos][pos2]+=(3.-temp2.getRe(0) - temp2.getRe(4) - temp2.getRe(8));
			  C[size-pos][pos2]+=(3.-temp3.getRe(0) - temp3.getRe(4) - temp3.getRe(8));
			}
		      counts[size-pos][pos2]++;
		    }
		}
	      else
		{
		  if(pos<size/2)
		    {
		      if (Nc==2)
			{
			  A[pos][size-pos2]+=(2.-temp2.getRe(0) - temp2.getRe(3))*(2.-temp3.getRe(0) - temp3.getRe(3));
			  B[pos][size-pos2]+=(2.-temp2.getRe(0) - temp2.getRe(3));
			  C[pos][size-pos2]+=(2.-temp3.getRe(0) - temp3.getRe(3));
			}
		      else if (Nc==3)
			{
			  A[pos][size-pos2]+=(3.-temp2.getRe(0) - temp2.getRe(4) - temp2.getRe(8))*(3.-temp3.getRe(0) - temp3.getRe(4) - temp3.getRe(8));
			  B[pos][size-pos2]+=(3.-temp2.getRe(0) - temp2.getRe(4) - temp2.getRe(8));
			  C[pos][size-pos2]+=(3.-temp3.getRe(0) - temp3.getRe(4) - temp3.getRe(8));
			}
		      counts[pos][size-pos2]++;
		    }
		  else if (pos>=size/2)
		    {
		      if (Nc==2)
			{
			  A[size-pos][size-pos2]+=(2.-temp2.getRe(0) - temp2.getRe(3))*(2.-temp3.getRe(0) - temp3.getRe(3));
			  B[size-pos][size-pos2]+=(2.-temp2.getRe(0) - temp2.getRe(3));
			  C[size-pos][size-pos2]+=(2.-temp3.getRe(0) - temp3.getRe(3));
			}
		      else if (Nc==3)
			{
			  A[size-pos][size-pos2]+=(3.-temp2.getRe(0) - temp2.getRe(4) - temp2.getRe(8))*(3.-temp3.getRe(0) - temp3.getRe(4) - temp3.getRe(8));
			  B[size-pos][size-pos2]+=(3.-temp2.getRe(0) - temp2.getRe(4) - temp2.getRe(8));
			  C[size-pos][size-pos2]+=(3.-temp3.getRe(0) - temp3.getRe(4) - temp3.getRe(8));
			}
		      counts[size-pos][size-pos2]++;
		    }
		}   
	      //cout << *tempC[j] << endl;
	      //	  fouta << i << " " << setw(12) << lat->cells[pos]->getU().getRe(0) << endl;
	    }
	}
    }
  
  
  for (int j=0; j<=size/2; j++)
    {
      A[j][j]/=static_cast<double>(counts[j][j]); 
      B[j][j]/=static_cast<double>(counts[j][j]); 
      C[j][j]/=static_cast<double>(counts[j][j]); 
    }
   
  ofstream foutg("fourPointDifference.dat",ios::app); 
  for (int j=0; j<=size/2; j++)
    {
      foutg << j << " "  << (A[j][j]-(B[j][j])*(C[j][j]))/param->getNc()/param->getNc() << " " << A[j][j] << " " << B[j][j] << " " << C[j][j] << endl;
    }
  foutg << endl;
  foutg.close();

  ofstream fout2("fourPointDifference2D.dat",ios::app); 
  for (int i=0; i<=size/2; i++)
    {
      for (int j=0; j<=size/2; j++)
	{
	  fout2 << i << " " << j << " "  << (A[i][j]-(B[i][j])*(C[i][j]))/param->getNc()/param->getNc() << endl;
	}
      fout2 << endl;
    }
  fout2.close();
} 

// ---------------------- Q and S_6 -----------------------

void Measure::fourPointFunctionLine(Parameters *param, Lattice *lat, int ids) //Q
{
  int bins;
  bins = size/2;

  for (int j=0; j<size; j++)
    {
      *tempC[j] = 0.; 
    }
  Matrix Udag(Nc,0.);
  Matrix Udagx(Nc,0.);
  Matrix Udaga(Nc,0.);
  Matrix temp2(Nc,0.);
  Matrix temp3(Nc,0.);
  Matrix tya(Nc,0.);
  Matrix txz(Nc,0.);
  Matrix taz(Nc,0.);
  Matrix tyz(Nc,0.);
  Matrix full(Nc,0.);
  double FULL[size];
  double YX[size];
  double YXYX[size];
  double c1mYX1mYX[size];
  int counts[size];

  double x,y,r;
  double rAverage[bins];
  double rstep=static_cast<double>(size/2)/static_cast<double>(bins);

  for (int j=0; j<size; j++)
    {  
      FULL[j]=0;
      YX[j]=0;
      YXYX[j]=0;
      c1mYX1mYX[j]=0;
      counts[j]=0;
      rAverage[j] = 0.; // average value of r in each bin
    }

  int pos,pos2,pos3,posx,posy,posz,posa;
  // measure on a line:
  // y---------x

  for (int j=0; j<size; j++)
    {  
      for (int l=0; l<size; l++)
	{
	  for (int m=0; m<size; m++)
	    {
	      if (m+j<size)
		{
		  posx = size*l+m+j;
		}
	      else
		{
		  posx = size*l+(m+j-size);
		}
	      
	      posy = size*l+m;
	      
	      Udagx = lat->cells[posx]->getU();
	      Udagx.conjg(); // now Udag is U^dagger
	      temp2 = lat->cells[posy]->getU()*Udagx; // this is U_y U_x^dag
	 
	      full = temp2*temp2; // this is U_y U_x^dag U_y U_x^dag
	  
	      pos = abs(j); 
	      
	      //	      r = pos; 
	      //position = static_cast<int>(floor(r/static_cast<double>(rstep)));
	
	      if(pos<size/2)
		{
		  if (Nc==2)
		    {
		      FULL[pos]+=(full.getRe(0) + full.getRe(3));
		      YX[pos]+=(temp2.getRe(0) + temp2.getRe(3));
		      YXYX[pos]+=(temp2.getRe(0) + temp2.getRe(3))*(temp2.getRe(0) + temp2.getRe(3));
		      c1mYX1mYX[pos]+=(2.-(temp2.getRe(0) + temp2.getRe(3)))*(2.-(temp2.getRe(0) + temp2.getRe(3)));
		      //      rAverage[pos]+=r;
		    }
		  else if (Nc==3)
		    {
		      FULL[pos]+=(full.getRe(0) + full.getRe(4) + full.getRe(8));	
		      YX[pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8));
		      YXYX[pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8))*(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8));
		      c1mYX1mYX[pos]+=(3.-(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8)))*(3.-(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8)));
		      //rAverage[pos]+=r;
		    }
		  counts[pos]++;
		}
	      else if (pos>=size/2)
		{
		  if (Nc==2)
		    {
		      FULL[size-pos]+=(full.getRe(0) + full.getRe(3));
		      YX[size-pos]+=(temp2.getRe(0) + temp2.getRe(3));
		      YXYX[size-pos]+=(temp2.getRe(0) + temp2.getRe(3))*(temp2.getRe(0) + temp2.getRe(3));
		      c1mYX1mYX[size-pos]+=(2.-(temp2.getRe(0) + temp2.getRe(3)))*(2.-(temp2.getRe(0) + temp2.getRe(3)));
		      //rAverage[size-pos]+=r;
		    }
		  else if (Nc==3)
		    {
		      FULL[size-pos]+=(full.getRe(0) + full.getRe(4) + full.getRe(8));
		      YX[size-pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8));
		      YXYX[size-pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8))*(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8));
		      c1mYX1mYX[size-pos]+=(3.-(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8)))*(3.-(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8)));
		      //rAverage[size-pos]+=r;
		    }
		  counts[size-pos]++;
		}
	      //cout << *tempC[j] << endl;
	      //	  fouta << i << " " << setw(12) << lat->cells[pos]->getU().getRe(0) << endl;
	    }
	}
    }
  
  
  for (int j=0; j<=size/2; j++)
    {
      FULL[j]/=static_cast<double>(counts[j]);
      YX[j]/=static_cast<double>(counts[j]); 
      YXYX[j]/=static_cast<double>(counts[j]); 
      c1mYX1mYX[j]/=static_cast<double>(counts[j]); 
      //  rAverage[j]/=static_cast<double>(counts[j]);   
    }
   
  double dNc = static_cast<double>(Nc);
  ofstream foutg("fourPointLine.dat",ios::app); 
  for (int j=0; j<=size/2; j++)
    {
      //      cout << j*param->getQs() << " " << static_cast<double>(size)/Pi
      //*sqrt(sin(Pi*static_cast<double>(j)/static_cast<double>(size))*sin(Pi*static_cast<double>(j)/static_cast<double>(size)))*param->getQs() << endl; 
      foutg << j << " " << j*param->getQs() << " " << FULL[j]/static_cast<double>(Nc) << " " 
	    << YX[j]*YX[j]/static_cast<double>(Nc)/static_cast<double>(Nc) << " " 
	    << (dNc+1.)/2.*pow(YX[j]/dNc,2.*(dNc+2.)/(dNc+1.)) - (dNc-1.)/2.*pow(YX[j]/dNc,2.*(dNc-2.)/(dNc-1.)) 
	    << " " << YX[j]*YX[j]/dNc/dNc*(1.+2.*log(YX[j]/dNc)) << " " << YX[j]/dNc << " " << YXYX[j]/dNc/dNc << " " << c1mYX1mYX[j]/dNc/dNc << endl;
    }
  foutg << endl;
  foutg.close();

  // compute Q_sQL
  // if correlator is 1-e^(-1/2), define that |x-y| as sqrt(2)/Q_s

  double Qs=0.;
  double R1=0.;
  double R2=0.;
  double Rs;
 
  // compute Qs from coordinate space correlator:
  for (int j=1; j<bins; j++)
    {
      if (Nc==2)
	{
	  if(2.-FULL[j] > 2.-2.*exp(-1./2.) && Qs==0.)
	    {
	      R2 = j;
	      R1 = j-1;
	      r=(-FULL[j]+2.*exp(-1./2.))/(-FULL[j]+FULL[j-1]);
	      Qs = r*1./R1+(1-r)*1./R2;
	    }
	}
      else if (Nc==3)
	{
	  if(FULL[j] < 3.*exp(-1./2.) && Qs==0.)
	    {
	      R2 = j;
	      R1 = j-1;
	      r=(-FULL[j]+3.*exp(-1./2.))/(-FULL[j]+FULL[j-1]);
	      Rs = R2*(1-r)+r*R1;
	      Qs = 1./Rs;
	    }
	}
    }

  Qs*=sqrt(2.);
  // cout << "Qs a=" << Qs << " (j=" << R2 << ") " << (R1+(1-r)*(R2-R1)) << endl;
  cout << "QsQL=" << Qs/param->getg2mu() << " g^2 mu" << endl;
  cout << "QsQL a=" << Qs << endl;
  //  cout << "Qs a=" << Qs << " (j=" << R2 << ") " << (R1+(1-r)*(R2-R1)) << endl;
  //cout << "Qs=" << Qs/param->getg2mu() << " g^2 mu" << endl;

  ofstream foutq("QsQL.dat",ios::app); 
  foutq << ids*ds*Pi*Pi << " " << Qs << endl;
  foutq.close();
} 

void Measure::fourPointFunctionSquare(Parameters *param, Lattice *lat, int ids)
{
  int bins;
  bins = size/2.;
  
  for (int j=0; j<size; j++)
    {
      *tempC[j] = 0.; 
    }
  Matrix Udag(Nc,0.);
  Matrix Udagx(Nc,0.);
  Matrix Udaga(Nc,0.);
  Matrix temp2(Nc,0.);
  Matrix temp3(Nc,0.);
  Matrix tya(Nc,0.);
  Matrix txz(Nc,0.);
  Matrix taz(Nc,0.);
  Matrix tyz(Nc,0.);
  Matrix full(Nc,0.);
  double FULL[size];
  double YX[size];
  double ZA[size];
  double AZ[size];
  double YA[size];
  double YZ[size];
  double XZ[size];
  int counts[size];
  
  for (int j=0; j<size; j++)
    {  
      FULL[j]=0;
      YX[j]=0;
      ZA[j]=0;
      AZ[j]=0;
      YA[j]=0;
      YZ[j]=0;
      XZ[j]=0;
      counts[j]=0;
    }

  int pos,pos2,pos3,posx,posy,posz,posa;
  // measure on a square: (z=u, v=a)
  // z---------a
  // |         |
  // |         |
  // |         |
  // y---------x

  for (int j=0; j<size; j++)
    {  
      for (int l=0; l<size; l++)
	{
	  for (int m=0; m<size; m++)
	    {
	      if (m+j<size)
		{
		  posx = size*l+m+j;
		}
	      else
		{
		  posx = size*l+(m+j-size);
		}
	      
	      posy = size*l+m;
	      
	      if (l+j<size)
		{
		  posz = size*(l+j)+m;
		}
	      else
		{
		  posz = size*(l+j-size)+m;
		}

	      if (m+j<size)
		{
		  if (l+j<size)
		    {
		      posa = size*(l+j)+m+j;
		    }
		  else
		    {
		      posa = size*(l+j-size)+m+j;
		    }
		}
	      else
		{
		  if (l+j<size)
		    {
		      posa = size*(l+j)+(m+j-size);
		    }
		  else
		    {
		      posa = size*(l+j-size)+(m+j-size);
		    }
		}
	      
	      Udagx = lat->cells[posx]->getU();
	      Udagx.conjg(); // now Udag is U^dagger
	      temp2 = lat->cells[posy]->getU()*Udagx; // this is U_y U_x^dag
	
	      Udag = lat->cells[posz]->getU();
	      Udag.conjg(); // now Udag is U^dagger
	      taz = lat->cells[posa]->getU()*Udag; // this is U_a U_z^dag

	      full = temp2*taz; // this is U_y U_x^dag U_a U_z^dag

	      pos = abs(j); 
	      
	      if(pos<size/2)
		{
		  if (Nc==2)
		    {
		      AZ[pos]+=(taz.getRe(0) + taz.getRe(3));
		      FULL[pos]+=(full.getRe(0) + full.getRe(3));
		      YX[pos]+=(temp2.getRe(0) + temp2.getRe(3));
		    }
		  else if (Nc==3)
		    {
		      AZ[pos]+=(taz.getRe(0) + taz.getRe(4) + taz.getRe(8));
		      FULL[pos]+=(full.getRe(0) + full.getRe(4) + full.getRe(8));	
		      YX[pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8));
		    }
		  counts[pos]++;
		}
	      else if (pos>=size/2)
		{
		  if (Nc==2)
		    {
		      AZ[size-pos]+=(taz.getRe(0) + taz.getRe(3));
		      FULL[size-pos]+=(full.getRe(0) + full.getRe(3)); // full correlator
		      YX[size-pos]+=(temp2.getRe(0) + temp2.getRe(3));
		    }
		  else if (Nc==3)
		    {
		      AZ[size-pos]+=(taz.getRe(0) + taz.getRe(4) + taz.getRe(8));
		      FULL[size-pos]+=(full.getRe(0) + full.getRe(4) + full.getRe(8)); // full correlator
		      YX[size-pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8));
		    }
		  counts[size-pos]++;
		}
	      //cout << *tempC[j] << endl;
	      //	  fouta << i << " " << setw(12) << lat->cells[pos]->getU().getRe(0) << endl;
	    }
	}
    }
  
  
  for (int j=0; j<=size/2; j++)
    {
      FULL[j]/=static_cast<double>(counts[j]);
      YX[j]/=static_cast<double>(counts[j]); 
      AZ[j]/=static_cast<double>(counts[j]); 
    }
  
  double dNc = static_cast<double>(Nc);
  double sqrtpos,S2sqrt, diff;
  int p1, p2;
  ofstream foutg("fourPointSquare.dat",ios::app); 
  for (int j=0; j<=size/2; j++)
    {
      sqrtpos = sqrt(2.)*j;
      if(sqrtpos<size/2)
	{
	  p1 = floor(sqrtpos);
	  p2 = p1+1;
	  diff = sqrtpos-static_cast<double>(p1);
	  S2sqrt = (YX[p1]*(1-diff)+YX[p2]*diff+AZ[p1]*(1-diff)+AZ[p2]*diff)/2.;
	  //	  if (S2sqrt<0.)
	  // cout << YX[p1] << ", " << YX[p2] << ", " << AZ[p1] << ", " << AZ[p2] << ", diff=" << diff << endl; 
	  
	}
      else
	continue;
      // cell, Qs*r, full, naive, Gaussian, Gaussian large Nc, S_2(r), S_2(\sqrt{2} r)
      foutg << j << " " 
	    << param->getQs()*j << " "  
	    << FULL[j]/static_cast<double>(Nc) << " " 
	    << YX[j]*AZ[j]/dNc/dNc << " " 
	    << pow((AZ[j]+YX[j])/2./dNc,2.) * ((dNc+1.)/2.*pow((AZ[j]+YX[j])/2./S2sqrt,2./(dNc+1))-(dNc-1.)/2.*pow(S2sqrt/((AZ[j]+YX[j])/2.),2./(dNc-1))) 
	    << " " 
	    << pow((AZ[j]+YX[j])/2./dNc,2.) * ( 1.+ 2.*log((AZ[j]+YX[j])/2./(S2sqrt)) ) << " " 
	    << (AZ[j]+YX[j])/2./dNc << " " 
	    << S2sqrt/dNc << endl;
      //cout << pow((AZ[j]+YX[j])/2./dNc,2.) << ", " << S2sqrt << endl;
    }
  foutg << endl;
  foutg.close();


  // compute Q_sQS

  double Qs=0.;
  double R1=0.;
  double R2=0.;
  double Rs;
  double r;

  // compute Qs from coordinate space correlator:
  for (int j=1; j<bins; j++)
    {
      if (Nc==2)
	{
	  if(2.-FULL[j] > 2.-2.*exp(-1./2.) && Qs==0.)
	    {
	      R2 = j;
	      R1 = j-1;
	      r=(-FULL[j]+2.*exp(-1./2.))/(-FULL[j]+FULL[j-1]);
	      Qs = r*1./R1+(1-r)*1./R2;
	    }
	}
      else if (Nc==3)
	{
	  if(FULL[j] < 3.*exp(-1./2.) && Qs==0.)
	    {
	      R2 = j;
	      R1 = j-1;
	      r=(-FULL[j]+3.*exp(-1./2.))/(-FULL[j]+FULL[j-1]);
	      Rs = R2*(1-r)+r*R1;
	      Qs = 1./Rs;
	    }
	}
    }

  Qs*=sqrt(2.);
  // cout << "Qs a=" << Qs << " (j=" << R2 << ") " << (R1+(1-r)*(R2-R1)) << endl;
  cout << "QsQS=" << Qs/param->getg2mu() << " g^2 mu" << endl;
  cout << "QsQS a=" << Qs << endl;
  //  cout << "Qs a=" << Qs << " (j=" << R2 << ") " << (R1+(1-r)*(R2-R1)) << endl;
  //cout << "Qs=" << Qs/param->getg2mu() << " g^2 mu" << endl;

  ofstream foutq("QsQS.dat",ios::app); 
  foutq << ids*ds*Pi*Pi << " " << Qs << endl;
  foutq.close();
} 



void Measure::sixPointFunctionLine(Parameters *param, Lattice *lat, int ids)
{
  for (int j=0; j<size; j++)
    {
      *tempC[j] = 0.; 
    }
  Matrix Udag(Nc,0.);
  Matrix Udagx(Nc,0.);
  Matrix Udaga(Nc,0.);
  Matrix temp2(Nc,0.);
  Matrix temp3(Nc,0.);
  Matrix txy(Nc,0.);
  Matrix full(Nc,0.);
  double FULL[size];
  double YX[size];
  double XY[size];
  int counts[size];
  
  for (int j=0; j<size; j++)
    {  
      FULL[j]=0;
      YX[j]=0;
      XY[j]=0;
      counts[j]=0;
    }

  int pos,pos2,pos3,posx,posy,posz,posa;
  // measure on a line:
  // y---------x

  for (int j=0; j<size; j++)
    {  
      for (int l=0; l<size; l++)
	{
	  for (int m=0; m<size; m++)
	    {
	      if (m+j<size)
		{
		  posx = size*(l)+m+j;
		}
	      else
		{
		  posx = size*(l)+m+j-size;
		}
	      
	      posy = size*l+m;
	      
	      Udagx = lat->cells[posx]->getU();
	      Udagx.conjg(); // now Udag is U^dagger
	      temp2 = lat->cells[posy]->getU()*Udagx; // this is U_y U_x^dag
	    
	      Udag = lat->cells[posy]->getU();
	      Udag.conjg(); // now Udag is U^dagger
	      txy = lat->cells[posx]->getU()*Udag; // this is U_x U_y^dag
	    
	      full = temp2*temp2; // this is U_y U_x^dag U_y U_x^dag

	      pos = abs(j); 
	      
	      if(pos<size/2)
		{
		  if (Nc==2)
		    {
		      FULL[pos]+=(full.getRe(0) + full.getRe(3))*(txy.getRe(0) + txy.getRe(3))
			-1./static_cast<double>(Nc)*(temp2.getRe(0) + temp2.getRe(3)); // full correlator
		      YX[pos]+=(temp2.getRe(0) + temp2.getRe(3));
		    }
		  else if (Nc==3)
		    {
		      FULL[pos]+=
			//	(full.getRe(0) + full.getRe(4) + full.getRe(8))*(txy.getRe(0) + txy.getRe(4) + txy.getRe(8));
			//(txy.getRe(0) + txy.getRe(4) + txy.getRe(8));
			(full.getRe(0) + full.getRe(4) + full.getRe(8))*(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))
			-(full.getIm(0) + full.getIm(4) + full.getIm(8))*(txy.getIm(0) + txy.getIm(4) + txy.getIm(8))
			-1./static_cast<double>(Nc)*(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8)); // full correlator
		      YX[pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8));
		    }
		  counts[pos]++;
		}
	      else if (pos>=size/2)
		{
		  if (Nc==2)
		    {
		      FULL[size-pos]+=(full.getRe(0) + full.getRe(3))*(txy.getRe(0) + txy.getRe(3))
			-1./static_cast<double>(Nc)*(temp2.getRe(0) + temp2.getRe(3)); // full correlator
		      YX[size-pos]+=(temp2.getRe(0) + temp2.getRe(3));
		    }
		  else if (Nc==3)
		    {
		      FULL[size-pos]+=
			//	(full.getRe(0) + full.getRe(4) + full.getRe(8))*(txy.getRe(0) + txy.getRe(4) + txy.getRe(8));
			//(txy.getRe(0) + txy.getRe(4) + txy.getRe(8));
			(full.getRe(0) + full.getRe(4) + full.getRe(8))*(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))
			-(full.getIm(0) + full.getIm(4) + full.getIm(8))*(txy.getIm(0) + txy.getIm(4) + txy.getIm(8))
			-1./static_cast<double>(Nc)*(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8)); // full correlator
		      YX[size-pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8));
		    }
		  counts[size-pos]++;
		}
	      //cout << *tempC[j] << endl;
	      //	  fouta << i << " " << setw(12) << lat->cells[pos]->getU().getRe(0) << endl;
	    }
	}
    }
  
  
  for (int j=0; j<=size/2; j++)
    {
      FULL[j]/=static_cast<double>(counts[j]);
      YX[j]/=static_cast<double>(counts[j]); 
    }
   
  double dNc=static_cast<double>(Nc);
  ofstream foutg("sixPointLine.dat",ios::app); 
  // cell, Qs*r, full, naive, Gaussian.

  for (int j=0; j<=size/2; j++)
    {
      foutg << j << " "  << param->getQs()*j << " "  << FULL[j]/(dNc*dNc-1.) << " " 
	    << YX[j]*YX[j]*YX[j]/(dNc*dNc*dNc*2.)+YX[j]*YX[j]*YX[j]/(dNc*dNc*dNc*2.) << " "
	    << dNc*dNc/(dNc*dNc-1)*((dNc+1.)/2.*pow(YX[j]/dNc,2.*(dNc+2.)/(dNc+1.)) - (dNc-1.)/2.*pow(YX[j]/dNc,2.*(dNc-2.)/(dNc-1.)))*YX[j]/dNc-1/(dNc*dNc-1)*YX[j]/dNc << endl;
    }
  foutg << endl;
  foutg.close();
} 

void Measure::sixPointFunctionSquare(Parameters *param, Lattice *lat, int ids)
{
  for (int j=0; j<size; j++)
    {
      *tempC[j] = 0.; 
    }
  Matrix Udag(Nc,0.);
  Matrix Udagx(Nc,0.);
  Matrix Udaga(Nc,0.);
  Matrix temp2(Nc,0.);
  Matrix temp3(Nc,0.);
  Matrix tya(Nc,0.);
  Matrix txz(Nc,0.);
  Matrix tzx(Nc,0.);
  Matrix taz(Nc,0.);
  Matrix tyz(Nc,0.);
  Matrix full(Nc,0.);
  double FULL[size];
  double YX[size];
  double ZA[size];
  double AZ[size];
  double YA[size];
  double YZ[size];
  double XZ[size];
  double ZX[size];
  int counts[size];
  
  for (int j=0; j<size; j++)
    {  
      FULL[j]=0;
      YX[j]=0;
      ZA[j]=0;
      AZ[j]=0;
      YA[j]=0;
      YZ[j]=0;
      XZ[j]=0;
      ZX[j]=0;
      counts[j]=0;
    }

  int pos,pos2,pos3,posx,posy,posz,posa;
  // measure on a square:
  // z---------a
  // |         |
  // |         |
  // |         |
  // y---------x

  for (int j=0; j<size; j++)
    {  
      for (int l=0; l<size; l++)
	{
	  for (int m=0; m<size; m++)
	    {
	      if (m+j<size)
		{
		  posx = size*l+m+j;
		}
	      else
		{
		  posx = size*l+(m+j-size);
		}
	      
	      posy = size*l+m;
	      
	      if (l+j<size)
		{
		  posz = size*(l+j)+m;
		}
	      else
		{
		  posz = size*(l+j-size)+m;
		}

	      if (m+j<size)
		{
		  if (l+j<size)
		    {
		      posa = size*(l+j)+m+j;
		    }
		  else
		    {
		      posa = size*(l+j-size)+m+j;
		    }
		}
	      else
		{
		  if (l+j<size)
		    {
		      posa = size*(l+j)+(m+j-size);
		    }
		  else
		    {
		      posa = size*(l+j-size)+(m+j-size);
		    }
		}
	      
	      Udagx = lat->cells[posx]->getU();
	      Udagx.conjg(); // now Udag is U^dagger
	      temp2 = lat->cells[posy]->getU()*Udagx; // this is U_y U_x^dag
	      txz = lat->cells[posz]->getU() *Udagx; // this is U_z U_x^dag
    
	      Udag = lat->cells[posz]->getU();
	      Udag.conjg(); // now Udag is U^dagger
	      tzx = lat->cells[posx]->getU() * Udag; // this is U_x U_z^dag
    
	      Udaga = lat->cells[posa]->getU();
	      Udaga.conjg(); // now Udag is U^dagger
	      temp3 = lat->cells[posz]->getU() * Udaga; // this is U_z U_a^dag
	      tya = lat->cells[posy]->getU() * Udaga; // this is U_y U_a^dag

	      Udag = lat->cells[posz]->getU();
	      Udag.conjg(); // now Udag is U^dagger
	      taz = lat->cells[posa]->getU() * Udag; // this is U_a U_z^dag

	  //     Udag = lat->cells[posy]->getU();
// 	      Udag.conjg(); // now Udag is U^dagger
// 	      tyz = lat->cells[posz]->getU() * Udag; // this is U_z U_y^dag
	      
	      full = temp2*taz; // this is U_y U_x^dag U_z U_a^dag

	      pos = abs(j); 
	      
	      if(pos<size/2)
		{
		  if (Nc==2)
		    {
		      AZ[pos]+=(taz.getRe(0) + taz.getRe(3));
		      FULL[pos]+=(full.getRe(0) + full.getRe(3))*(temp3.getRe(0) + temp3.getRe(3))
			-1./static_cast<double>(Nc)*(temp2.getRe(0) + temp2.getRe(3)); // full correlator
		      YX[pos]+=(temp2.getRe(0) + temp2.getRe(3));
		      ZA[pos]+=(temp3.getRe(0) + temp3.getRe(3));
		      YA[pos]+=(tya.getRe(0) + tya.getRe(3));
		      XZ[pos]+=(txz.getRe(0) + txz.getRe(3));
		      ZX[pos]+=(tzx.getRe(0) + tzx.getRe(3));
		      //      YZ[pos]+=(tyz.getRe(0) + tyz.getRe(3));
		    }
		  else if (Nc==3)
		    {
		      AZ[pos]+=(taz.getRe(0) + taz.getRe(4) + taz.getRe(8));
		      FULL[pos]+=(full.getRe(0) + full.getRe(4) + full.getRe(8))*(temp3.getRe(0) + temp3.getRe(4) + temp3.getRe(8))	
			-1./static_cast<double>(Nc)*(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8)); // full correlator
		      YX[pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8));
		      ZA[pos]+=(temp3.getRe(0) + temp3.getRe(4) + temp3.getRe(8));
		      YA[pos]+=(tya.getRe(0) + tya.getRe(4) + tya.getRe(8));
		      XZ[pos]+=(txz.getRe(0) + txz.getRe(4) + txz.getRe(8));
		      ZX[pos]+=(tzx.getRe(0) + tzx.getRe(4) + tzx.getRe(8));
		      //YZ[pos]+=(tyz.getRe(0) + tyz.getRe(4) + tyz.getRe(8));
		    }
		  counts[pos]++;
		}
	      else if (pos>=size/2)
		{
		  if (Nc==2)
		    {
		      AZ[size-pos]+=(taz.getRe(0) + taz.getRe(3));
		      FULL[size-pos]+=(full.getRe(0) + full.getRe(3))*(temp3.getRe(0) + temp3.getRe(3))
			-1./static_cast<double>(Nc)*(temp2.getRe(0) + temp2.getRe(3)); // full correlator
		      YX[size-pos]+=(temp2.getRe(0) + temp2.getRe(3));
		      ZA[size-pos]+=(temp3.getRe(0) + temp3.getRe(3));
		      YA[size-pos]+=(tya.getRe(0) + tya.getRe(3));
		      XZ[size-pos]+=(txz.getRe(0) + txz.getRe(3));
		      ZX[size-pos]+=(tzx.getRe(0) + tzx.getRe(3));
		      //YZ[size-pos]+=(tyz.getRe(0) + tyz.getRe(3));
		    }
		  else if (Nc==3)
		    {
		      AZ[size-pos]+=(taz.getRe(0) + taz.getRe(4) + taz.getRe(8));
		      FULL[size-pos]+=(full.getRe(0) + full.getRe(4) + full.getRe(8))*(temp3.getRe(0) + temp3.getRe(4) + temp3.getRe(8))	
			-1./static_cast<double>(Nc)*(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8)); // full correlator
		      YX[size-pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8));
		      ZA[size-pos]+=(temp3.getRe(0) + temp3.getRe(4) + temp3.getRe(8));
		      YA[size-pos]+=(tya.getRe(0) + tya.getRe(4) + tya.getRe(8));
		      XZ[size-pos]+=(txz.getRe(0) + txz.getRe(4) + txz.getRe(8));
		      ZX[size-pos]+=(tzx.getRe(0) + tzx.getRe(4) + tzx.getRe(8));
		      //YZ[size-pos]+=(tyz.getRe(0) + tyz.getRe(4) + tyz.getRe(8));
		    }
		  counts[size-pos]++;
		}
	      //cout << *tempC[j] << endl;
	      //	  fouta << i << " " << setw(12) << lat->cells[pos]->getU().getRe(0) << endl;
	    }
	}
    }
  
  
  for (int j=0; j<=size/2; j++)
    {
      FULL[j]/=static_cast<double>(counts[j]);
      YX[j]/=static_cast<double>(counts[j]); 
      ZA[j]/=static_cast<double>(counts[j]); 
      AZ[j]/=static_cast<double>(counts[j]); 
      YA[j]/=static_cast<double>(counts[j]); 
      XZ[j]/=static_cast<double>(counts[j]); 
      ZX[j]/=static_cast<double>(counts[j]); 
      //YZ[j]/=static_cast<double>(counts[j]); 
    }
   
  double dNc=static_cast<double>(Nc);
  double sqrtpos,S2sqrt, diff;
  int p1, p2;
  ofstream foutg("sixPointSquare.dat",ios::app); 
  for (int j=0; j<=size/2; j++)
    {
      sqrtpos = sqrt(2.)*j;
      if(sqrtpos<size/2)
	{
	  p1 = floor(sqrtpos);
	  p2 = p1+1;
	  diff = sqrtpos-static_cast<double>(p1);
	  S2sqrt = (YX[p1]*(1-diff)+YX[p2]*diff+AZ[p1]*(1-diff)+AZ[p2]*diff)/2.;
	}
      else
	continue;
      // cell, Qs*r, full, naive, Gaussian.
          foutg << j << " "  << param->getQs()*j << " "  << FULL[j]/(dNc*dNc-1.) << " " 
	    << YA[j]*XZ[j]*AZ[j]/(dNc*dNc*dNc*2.)+YX[j]*ZA[j]*AZ[j]/(dNc*dNc*dNc*2.) << " "
	    << dNc*dNc/(dNc*dNc-1)*(pow((AZ[j]+YX[j])/2./dNc,2.) * ((dNc+1.)/2.*pow((AZ[j]+YX[j])/2./S2sqrt,2./(dNc+1))-(dNc-1.)/2.*pow(S2sqrt/((AZ[j]+YX[j])/2.),2./(dNc-1))))*(AZ[j]+YX[j])/2./dNc-1./(dNc*dNc-1)*(AZ[j]+YX[j])/2./dNc << endl;
    }
  foutg << endl;
  foutg.close();
} 


void Measure::twoPointFunctions(Parameters *param, Lattice *lat, int ids)
{
  for (int j=0; j<size; j++)
    {
      *tempC[j] = 0.; 
    }
  Matrix Udag(Nc,0.);
  Matrix Udagx(Nc,0.);
  Matrix Udaga(Nc,0.);
  Matrix temp2(Nc,0.);
  Matrix temp3(Nc,0.);
  Matrix tya(Nc,0.);
  Matrix tyx(Nc,0.);
  Matrix txy(Nc,0.);
  Matrix txyxy(Nc,0.);
  Matrix txyxyxy(Nc,0.);
  Matrix txz(Nc,0.);
  Matrix tzx(Nc,0.);
  Matrix taz(Nc,0.);
  Matrix tyz(Nc,0.);
  Matrix full(Nc,0.);
  double FULL[size];
  double S12[size];
  double S12S21[size];
  double S12S23[size];
  double Q1212S12[size];
  double X121212[size];
  double S12S12S12[size];
  double S21[size];
  double S23[size];
  double S13[size];
 
  double YX[size];
  double ZA[size];
  double AZ[size];
  double YA[size];
  double YZ[size];
  double XZ[size];
  double ZX[size];
  int counts[size];
  double dNc=static_cast<double>(Nc);
  
  for (int j=0; j<size; j++)
    {  
      FULL[j]=0;
      YX[j]=0;
      ZA[j]=0;
      AZ[j]=0;
      YA[j]=0;
      YZ[j]=0;
      XZ[j]=0;
      ZX[j]=0;
      counts[j]=0;
      S12[j]=0;
      S12S21[j]=0;
      S12S23[j]=0;
      S12S12S12[j]=0;
      Q1212S12[j]=0;
      X121212[j]=0;
      S13[j]=0;
      S21[j]=0;
      S23[j]=0;
    }

  int pos,pos2,pos3,posx,posy,posz,posa;

  for (int j=0; j<size; j++)
    {  
      for (int l=0; l<size; l++)
	{
	  for (int m=0; m<size; m++)
	    {
	      if (m+j<size)
		{
		  posx = size*l+m+j;
		}
	      else
		{
		  posx = size*l+(m+j-size);
		}
	      
	      posy = size*l+m;
	      
	      if (l+j<size)
		{
		  posz = size*(l+j)+m;
		}
	      else
		{
		  posz = size*(l+j-size)+m;
		}

	      if (m+j<size)
		{
		  if (l+j<size)
		    {
		      posa = size*(l+j)+m+j;
		    }
		  else
		    {
		      posa = size*(l+j-size)+m+j;
		    }
		}
	      else
		{
		  if (l+j<size)
		    {
		      posa = size*(l+j)+(m+j-size);
		    }
		  else
		    {
		      posa = size*(l+j-size)+(m+j-size);
		    }
		}
	      
	      // for S12
	      Udag = lat->cells[posx]->getU();
	      Udag.conjg(); // now Udag is U^dagger
	      txy = Udag*lat->cells[posy]->getU(); // this is U_x^dag U_y
	    
	      // for Q1212
	      txyxy = txy*txy;
	      // for X121212
	      txyxyxy = txyxy*txy;
  
	      // for S23
	      Udag = lat->cells[posy]->getU();
	      Udag.conjg(); // now Udag is U^dagger
	      tyz = Udag*lat->cells[posz]->getU(); // this is U_y^dag U_z
    
	      // for S13
	      Udag = lat->cells[posx]->getU();
	      Udag.conjg(); // now Udag is U^dagger
	      txz = lat->cells[posz]->getU() * Udag; // this is U_x^dag U_z
    
	      // for S21
	      Udag = lat->cells[posy]->getU();
	      Udag.conjg(); // now Udag is U^dagger
	      tyx = Udag*lat->cells[posx]->getU(); // this is U_y^dag U_x

	      pos = abs(j); 
	      
	      if(pos<size/2)
		{
		  if (Nc==2)
		    {
		      S12[pos]+=(txy.getRe(0) + txy.getRe(3))/dNc;
		      S12S21[pos]+=(txy.getRe(0) + txy.getRe(3))/dNc*(tyx.getRe(0) + tyx.getRe(3))/dNc
			-(txy.getIm(0) + txy.getIm(3))/dNc*(tyx.getIm(0) + tyx.getIm(3))/dNc;
		      S12S23[pos]+=(txy.getRe(0) + txy.getRe(3))/dNc*(tyz.getRe(0) + tyz.getRe(3))/dNc
			-(txy.getIm(0) + txy.getIm(3))/dNc*(tyz.getIm(0) + tyz.getIm(3))/dNc;
		      S21[pos]+=(tyx.getRe(0) + tyx.getRe(3))/dNc;
		      S23[pos]+=(tyz.getRe(0) + tyz.getRe(3))/dNc;
		      S13[pos]+=(txz.getRe(0) + txz.getRe(3))/dNc;
		      Q1212S12[pos]+=(txyxy.getRe(0) + txyxy.getRe(3))*(txy.getRe(0) + txy.getRe(3))/dNc/dNc
			-(txyxy.getIm(0) + txyxy.getIm(3))*(txy.getIm(0) + txy.getIm(3))/dNc/dNc;
		      X121212[pos]+=(txyxyxy.getRe(0) + txyxyxy.getRe(3))/dNc;
		      S12S12S12[pos]+=(txy.getRe(0) + txy.getRe(3))*(txy.getRe(0) + txy.getRe(3))*(txy.getRe(0) + txy.getRe(3))/dNc/dNc/dNc
			-3.*(txy.getIm(0) + txy.getIm(3))*(txy.getIm(0) + txy.getIm(3))*(txy.getRe(0) + txy.getRe(3))/dNc/dNc/dNc;
		    }
		  else if (Nc==3)
		    {
		      S12[pos]+=(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))/dNc;
		      S12S21[pos]+=(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))/dNc*(tyx.getRe(0) + tyx.getRe(4) + tyx.getRe(8))/dNc
			-(txy.getIm(0) + txy.getIm(4) + txy.getIm(8))/dNc*(tyx.getIm(0) + tyx.getIm(4) + tyx.getIm(8))/dNc;
		      S12S23[pos]+=(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))/dNc*(tyz.getRe(0) + tyz.getRe(4) + tyz.getRe(8))/dNc
			-(txy.getIm(0) + txy.getIm(4) + txy.getIm(8))/dNc*(tyz.getIm(0) + tyz.getIm(4) + tyz.getIm(8))/dNc;
		      S21[pos]+=(tyx.getRe(0) + tyx.getRe(4) + tyx.getRe(8))/dNc;
		      S23[pos]+=(tyz.getRe(0) + tyz.getRe(4) + tyz.getRe(8))/dNc;
		      S13[pos]+=(txz.getRe(0) + txz.getRe(4) + txz.getRe(8))/dNc;
		      Q1212S12[pos]+=(txyxy.getRe(0) + txyxy.getRe(4) + txyxy.getRe(8))*(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))/dNc/dNc
			-(txyxy.getIm(0) + txyxy.getIm(4) + txyxy.getIm(8))*(txy.getIm(0) + txy.getIm(4) + txy.getIm(8))/dNc/dNc;
		      X121212[pos]+=(txyxyxy.getRe(0) + txyxyxy.getRe(4) + txyxyxy.getRe(8))/dNc;
		      S12S12S12[pos]+=(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))*(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))
			*(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))/dNc/dNc/dNc
			-3.*(txy.getIm(0) + txy.getIm(4) + txy.getIm(8))*(txy.getIm(0) + txy.getIm(4) + txy.getIm(8))
			*(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))/dNc/dNc/dNc;
		    }
		  counts[pos]++;
		}
	      else if (pos>=size/2)
		{
		  if (Nc==2)
		    {
		      S12[size-pos]+=(txy.getRe(0) + txy.getRe(3))/dNc;
		      S12S21[size-pos]+=(txy.getRe(0) + txy.getRe(3))/dNc*(tyx.getRe(0) + tyx.getRe(3))/dNc
		      -(txy.getIm(0) + txy.getIm(3))/dNc*(tyx.getIm(0) + tyx.getIm(3))/dNc;
		      S12S23[size-pos]+=(txy.getRe(0) + txy.getRe(3))/dNc*(tyz.getRe(0) + tyz.getRe(3))/dNc
			-(txy.getIm(0) + txy.getIm(3))/dNc*(tyz.getIm(0) + tyz.getIm(3))/dNc;
		      S21[size-pos]+=(tyx.getRe(0) + tyx.getRe(3))/dNc;
		      S23[size-pos]+=(tyz.getRe(0) + tyz.getRe(3))/dNc;
		      S13[size-pos]+=(txz.getRe(0) + txz.getRe(3))/dNc;
		      Q1212S12[size-pos]+=(txyxy.getRe(0) + txyxy.getRe(3))*(txy.getRe(0) + txy.getRe(3))/dNc/dNc
			-(txyxy.getIm(0) + txyxy.getIm(3))*(txy.getIm(0) + txy.getIm(3))/dNc/dNc;
		      X121212[size-pos]+=(txyxyxy.getRe(0) + txyxyxy.getRe(3))/dNc;
		      S12S12S12[size-pos]+=(txy.getRe(0) + txy.getRe(3))*(txy.getRe(0) + txy.getRe(3))*(txy.getRe(0) + txy.getRe(3))/dNc/dNc/dNc
			-3.*(txy.getIm(0) + txy.getIm(3))*(txy.getIm(0) + txy.getIm(3))*(txy.getRe(0) + txy.getRe(3))/dNc/dNc/dNc;
		    }
		  else if (Nc==3)
		    {
		      S12[size-pos]+=(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))/dNc;
		      S12S21[size-pos]+=(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))/dNc*(tyx.getRe(0) + tyx.getRe(4) + tyx.getRe(8))/dNc
			-(txy.getIm(0) + txy.getIm(4) + txy.getIm(8))/dNc*(tyx.getIm(0) + tyx.getIm(4) + tyx.getIm(8))/dNc;
		      S12S23[size-pos]+=(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))/dNc*(tyz.getRe(0) + tyz.getRe(4) + tyz.getRe(8))/dNc
			-(txy.getIm(0) + txy.getIm(4) + txy.getIm(8))/dNc*(tyz.getIm(0) + tyz.getIm(4) + tyz.getIm(8))/dNc;
		      S21[size-pos]+=(tyx.getRe(0) + tyx.getRe(4) + tyx.getRe(8))/dNc;
		      S23[size-pos]+=(tyz.getRe(0) + tyz.getRe(4) + tyz.getRe(8))/dNc;
		      S13[size-pos]+=(txz.getRe(0) + txz.getRe(4) + txz.getRe(8))/dNc;
		      Q1212S12[size-pos]+=(txyxy.getRe(0) + txyxy.getRe(4) + txyxy.getRe(8))*(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))/dNc/dNc
			-(txyxy.getIm(0) + txyxy.getIm(4) + txyxy.getIm(8))*(txy.getIm(0) + txy.getIm(4) + txy.getIm(8))/dNc/dNc;
		      X121212[size-pos]+=(txyxyxy.getRe(0) + txyxyxy.getRe(4) + txyxyxy.getRe(8))/dNc;
		      S12S12S12[size-pos]+=(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))*(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))
			*(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))/dNc/dNc/dNc
			-3.*(txy.getIm(0) + txy.getIm(4) + txy.getIm(8))*(txy.getIm(0) + txy.getIm(4) + txy.getIm(8))
			*(txy.getRe(0) + txy.getRe(4) + txy.getRe(8))/dNc/dNc/dNc;
		    }
		  counts[size-pos]++;
		}
	      //cout << *tempC[j] << endl;
	      //	  fouta << i << " " << setw(12) << lat->cells[pos]->getU().getRe(0) << endl;
	    }
	}
    }
    
  for (int j=0; j<=size/2; j++)
    {
      S12[j]/=static_cast<double>(counts[j]); 
      S12S12S12[j]/=static_cast<double>(counts[j]); 
      S12S21[j]/=static_cast<double>(counts[j]); 
      S12S23[j]/=static_cast<double>(counts[j]); 
      S21[j]/=static_cast<double>(counts[j]); 
      S23[j]/=static_cast<double>(counts[j]); 
      S13[j]/=static_cast<double>(counts[j]); 
      Q1212S12[j]/=static_cast<double>(counts[j]); 
      X121212[j]/=static_cast<double>(counts[j]); 
    }
   
  double sqrtpos,S2sqrt, diff;
  int p1, p2;
  ofstream foutg("dionysis.dat",ios::app); 
  for (int j=0; j<=size/2; j++)
    {
      foutg << j << " "  << param->getQs()*j << " "       //1 2 3
	    << S12S23[j]/S12[j]/S23[j] << " "     //3     4 5
	    << dNc*dNc*S12S23[j]/S13[j] << " "    //4     6 7
	    << S12S21[j] << " "                   //5     8 9
	    << X121212[j] << " "                  //6     10 11
	    << Q1212S12[j] << " "                 //7     12 13
	    << S12S12S12[j] << " "                //8     14 15
	    << S12[j] << " "                      //9     16 17
	    << S23[j] << " "                      //10    18 19
	    << S13[j] << " "                      //11    20 21
	    << S12S23[j] << endl;                 //12    22 23
    }
  foutg << endl;
  foutg.close();
} 


// -- unequal rapidity correlations
void Measure::fourPointFunctionLineUnequalY(Parameters *param, Lattice *lat, int ids) //Q
{
  int bins;
  bins = size/2;

  for (int j=0; j<size; j++)
    {
      *tempC[j] = 0.; 
    }
  Matrix Udag(Nc,0.);
  Matrix Uidag(Nc,0.);
  Matrix Udagx(Nc,0.);
  Matrix Uidagx(Nc,0.);
  Matrix Udaga(Nc,0.);
  Matrix initial(Nc,0.);
  Matrix temp2(Nc,0.);
  Matrix temp3(Nc,0.);
  Matrix tya(Nc,0.);
  Matrix txz(Nc,0.);
  Matrix taz(Nc,0.);
  Matrix tyz(Nc,0.);
  Matrix full(Nc,0.);
  double FULL[size];
  double DtrYXtrYXi[size];
  double trYXtrYXi[size];
  double YX[size];
  double YXi[size];
  double YXni[size];
  int counts[size];
  double dNc = static_cast<double>(Nc);

  double x,y,r;
  double rAverage[bins];
  double rstep=static_cast<double>(size/2)/static_cast<double>(bins);

  for (int j=0; j<size; j++)
    {  
      FULL[j]=0;
      YX[j]=0;
      YXi[j]=0;
      YXni[j]=0;
      DtrYXtrYXi[j]=0;
      trYXtrYXi[j]=0;
      counts[j]=0;
      rAverage[j] = 0.; // average value of r in each bin
    }

  int pos,pos2,pos3,posx,posy,posz,posa;
  // measure on a line:
  // y---------x

  for (int j=0; j<size; j++)
    {  
      for (int l=0; l<size; l++)
	{
	  for (int m=0; m<size; m++)
	    {
	      if (m+j<size)
		{
		  posx = size*l+m+j;
		}
	      else
		{
		  posx = size*l+(m+j-size);
		}
	      
	      posy = size*l+m;
	      
	      Udagx = lat->cells[posx]->getU();
	      Udagx.conjg(); // now Udag is U^dagger
	      temp2 = lat->cells[posy]->getU()*Udagx; // this is U_y U_x^dag
	      
	      Uidagx = lat->cells[posx]->getUi();
	      Uidagx.conjg(); // now Uidag is Ui^dagger
	      initial = lat->cells[posy]->getUi()*Uidagx; // this is U_y^initial U_x^dag,initial
	 
	      Uidagx = lat->cells[posx]->getUi();
	      Uidagx.conjg(); // now Uidag is Ui^dagger
	      temp3 = lat->cells[posy]->getU()*Uidagx; // this is U_y U_x^dag,initial
	   

	      full = initial*temp2; // this is U_y,init U_x^dag,init U_y U_x^dag
	  
	      pos = abs(j); 
	      
	      if(pos<size/2)
		{
		  if (Nc==2)
		    {
		      FULL[pos]+=(full.getRe(0) + full.getRe(3));
		      DtrYXtrYXi[pos]+=(1.-(temp2.getRe(0) + temp2.getRe(3))/dNc)*(1.-(initial.getRe(0) + initial.getRe(3))/dNc)
			-(temp2.getIm(0) + temp2.getIm(3))/dNc*(initial.getIm(0) + initial.getIm(3))/dNc;
		      trYXtrYXi[pos]+=(temp2.getRe(0) + temp2.getRe(3))/dNc*(initial.getRe(0) + initial.getRe(3))/dNc
			-(temp2.getIm(0) + temp2.getIm(3))/dNc*(initial.getIm(0) + initial.getIm(3))/dNc;
		      YX[pos]+=(temp2.getRe(0) + temp2.getRe(3));
		      YXni[pos]+=(temp3.getRe(0) + temp3.getRe(3));
		      YXi[pos]+=(initial.getRe(0) + initial.getRe(3));
		    }
		  else if (Nc==3)
		    {
		      FULL[pos]+=(full.getRe(0) + full.getRe(4) + full.getRe(8));	
		      DtrYXtrYXi[pos]+=(1.-(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8))/dNc)
			*(1.-(initial.getRe(0) + initial.getRe(4) + initial.getRe(8))/dNc)
			-(temp2.getIm(0) + temp2.getIm(4) + temp2.getIm(8))/dNc
			*(initial.getIm(0) + initial.getIm(4) + initial.getIm(8))/dNc;
		      trYXtrYXi[pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8))/dNc
			*(initial.getRe(0) + initial.getRe(4) + initial.getRe(8))/dNc
			-(temp2.getIm(0) + temp2.getIm(4) + temp2.getIm(8))/dNc
			*(initial.getIm(0) + initial.getIm(4) + initial.getIm(8))/dNc;
		      YX[pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8));
		      YXni[pos]+=(temp3.getRe(0) + temp3.getRe(4) + temp3.getRe(8));
		      YXi[pos]+=(initial.getRe(0) + initial.getRe(4) + initial.getRe(8));
		    }
		  counts[pos]++;
		}
	      else if (pos>=size/2)
		{
		  if (Nc==2)
		    {
		      FULL[size-pos]+=(full.getRe(0) + full.getRe(3));
		      DtrYXtrYXi[size-pos]+=(1.-(temp2.getRe(0) + temp2.getRe(3))/dNc)*(1.-(initial.getRe(0) + initial.getRe(3))/dNc)
			-(temp2.getIm(0) + temp2.getIm(3))/dNc*(initial.getIm(0) + initial.getIm(3))/dNc;
		      trYXtrYXi[size-pos]+=(temp2.getRe(0) + temp2.getRe(3))/dNc*(initial.getRe(0) + initial.getRe(3))/dNc
			-(temp2.getIm(0) + temp2.getIm(3))/dNc*(initial.getIm(0) + initial.getIm(3))/dNc;
		      YX[size-pos]+=(temp2.getRe(0) + temp2.getRe(3));
		      YXni[size-pos]+=(temp3.getRe(0) + temp3.getRe(3));
		      YXi[size-pos]+=(initial.getRe(0) + initial.getRe(3));
		    }
		  else if (Nc==3)
		    {
		      FULL[size-pos]+=(full.getRe(0) + full.getRe(4) + full.getRe(8));
		      DtrYXtrYXi[size-pos]+=(1.-(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8))/dNc)
			*(1.-(initial.getRe(0) + initial.getRe(4) + initial.getRe(8))/dNc)
			-(temp2.getIm(0) + temp2.getIm(4) + temp2.getIm(8))/dNc
			*(initial.getIm(0) + initial.getIm(4) + initial.getIm(8))/dNc;
		      trYXtrYXi[size-pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8))/dNc
			*(initial.getRe(0) + initial.getRe(4) + initial.getRe(8))/dNc
			-(temp2.getIm(0) + temp2.getIm(4) + temp2.getIm(8))/dNc
			*(initial.getIm(0) + initial.getIm(4) + initial.getIm(8))/dNc;
		      YX[size-pos]+=(temp2.getRe(0) + temp2.getRe(4) + temp2.getRe(8));
		      YXni[size-pos]+=(temp3.getRe(0) + temp3.getRe(4) + temp3.getRe(8));
		      YXi[size-pos]+=(initial.getRe(0) + initial.getRe(4) + initial.getRe(8));
		    }
		  counts[size-pos]++;
		}
	      //cout << *tempC[j] << endl;
	      //	  fouta << i << " " << setw(12) << lat->cells[pos]->getU().getRe(0) << endl;
	    }
	}
    }
  
  
  for (int j=0; j<=size/2; j++)
    {
      //      cout << "counts[" << j << "]=" << counts[j] << endl;
      FULL[j]/=static_cast<double>(counts[j]);
      DtrYXtrYXi[j]/=static_cast<double>(counts[j]);
      trYXtrYXi[j]/=static_cast<double>(counts[j]);
      YX[j]/=static_cast<double>(counts[j]); 
      YXi[j]/=static_cast<double>(counts[j]); 
      YXni[j]/=static_cast<double>(counts[j]); 
    }
   
  ofstream foutg("fourPointLineUnequalY.dat",ios::app); 
  for (int j=0; j<=size/2; j++)
    {
      foutg << j << " " << j*param->getQs() << " " << FULL[j]/dNc << " " << trYXtrYXi[j] << " " 
	    << YX[j]/dNc << " " << YXi[j]/dNc << " " << YXni[j]/dNc << " " << DtrYXtrYXi[j] << endl;
    }
  foutg << endl;
  foutg.close();
} 
