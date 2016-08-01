// Infrared.cpp is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.
#include "Infrared.h"

//**************************************************************************
// Infrared regulator class.

void Infrared::regulate(Lattice *lat, Group *group, Parameters *param, Random *random, int ids)
{
  // initialize rho as Gaussian noise
  int N = param->getSize();
  int Ny=param->getNy();
  int Nc2m1;
  int Nc = param->getNc();
  Nc2m1 = Nc*Nc-1;
  double m = param->getm();
  double L = param->getL();
  double a = L/static_cast<double>(N); // lattice spacing in fm

  m=m*a/param->hbarc;
  //  cout << "m_lat =" << m << endl; 

  int nn[2];
  nn[0]=N;
  nn[1]=N;
    
  int pos;
  double kt2, kx, ky;

  Matrix temp(Nc,1.);
  Matrix U(Nc,1.);
  Matrix UD(Nc,1.);
  Matrix temp2(Nc,0.);
  double temp3;

  Matrix **A;

  A = new Matrix*[nn[0]*nn[1]];
  for(int i=0; i<nn[0]*nn[1]; i++)
    {
      A[i] = new Matrix(Nc,0.);
    }
  
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  U = lat->cells[pos]->getU();
	  U.logm();
	  *A[pos] = U;
	}
    }

  // Fourier transform A
  fft->fftn(A,A,nn,2,1);


  // then filter out infrared tail - current version: compute \rho (no filter)
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  kx = 2.*param->PI*(-0.5+static_cast<double>(i)/static_cast<double>(nn[0]));
	  ky = 2.*param->PI*(-0.5+static_cast<double>(j)/static_cast<double>(nn[1]));
	  //kt2 = kx*kx+ky*ky;
	  kt2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); //lattice momentum
	  //divide by i:
	  for(int nc=0; nc<Nc*Nc; nc++)
	    {
	      temp3 = A[pos]->getRe(nc); 
	      A[pos]->setRe(nc,A[pos]->getIm(nc));
	      A[pos]->setIm(nc,-temp3);
	    }
	  
	  *A[pos] = *A[pos]*kt2;
	}
    }
      
  // Fourier transform back rho (A contains rho(x_\perp) here)
  fft->fftn(A,A,nn,2,-1);

  complex<double> rho[8];
  Matrix rho0t0(Nc,0.);
  Matrix rho1t1(Nc,0.);
  Matrix rho2t2(Nc,0.);
  Matrix rho3t3(Nc,0.);
  Matrix rho4t4(Nc,0.);
  complex<double> corr2[Nc2m1][Nc2m1];
  complex<double> corr4[Nc2m1][Nc2m1][Nc2m1][Nc2m1];
  
  for(int i=0; i<Nc2m1; i++)
    for(int j=0; j<Nc2m1; j++)
      {
	corr2[i][j] = 0;
	for(int k=0; k<Nc2m1; k++)
	  for(int l=0; l<Nc2m1; l++)
	    {
	      corr4[i][j][k][l] = 0;
	    }
      }


  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  
	  for (int k=0; k<Nc2m1; k++)
	    {
	      rho[k] = 2.*(group->getT(k)*(*A[pos])).trace(); // 2 is here because tr(t^a t^b) = 1/2.*\delta^{ab}
	    }

	  for (int k=0; k<Nc2m1; k++)
	    {
	  //     if (rho[k].real()>4.)
// 		{
// 		  cout << "rho(" << k << ")=" << rho[k] << endl;
// 		  cout << "at " << i << " " << j << endl;
// 		}
	      for (int l=0; l<Nc2m1; l++)
		{
		  corr2[k][l] += rho[k]*rho[l];
		  if(k!=0 && l!=1)
		    corr4[k][k][l][l] += rho[k]*rho[k]*rho[l]*rho[l];
		  if(k!=l)
		    corr4[0][0][1][1] += rho[k]*rho[k]*rho[l]*rho[l];
		}
	    }
	      
	  corr4[0][1][2][3] += rho[0]*rho[1]*rho[2]*rho[3];
	}
    }
 
  for (int k=0; k<Nc2m1; k++)
    {
      for (int l=0; l<Nc2m1; l++)
	{
	  corr2[k][l]/=static_cast<double>(nn[0]*nn[1]);
	}
    }
  
  corr4[0][1][2][3]/=static_cast<double>(nn[0]*nn[1]);
  corr4[0][0][1][1]/=static_cast<double>(nn[0]*nn[1])*(Nc2m1*(Nc2m1-1));

  cout << "<rho0 rho0>=" << corr2[0][0] << endl;
  cout << "<rho0 rho1>=" << corr2[0][1] << endl;
  cout << "<rho1 rho1>=" << corr2[1][1] << endl;

  cout << "<rho0 rho0 rho1 rho1>=" << corr4[0][0][1][1] << endl;
  cout << "<rho0 rho0 rho1 rho1> - Gaussians =" 
       << (corr4[0][0][1][1].real()-2.*corr2[0][1].real()*corr2[0][1].real()-corr2[0][0].real()*corr2[1][1].real())/corr4[0][0][1][1].real() << endl;
  

  double ds = param->getDs();
  double Pi = param->PI;
     
      
  ofstream fout("rhorho.dat",ios::app); 
  fout << ids*ds*Pi*Pi << " " <<  (corr4[0][0][1][1].real()-2.*corr2[0][1].real()*corr2[0][1].real()-corr2[0][0].real()*corr2[1][1].real())/corr4[0][0][1][1].real() << endl;
  fout.close();

// no filtering -> no rewriting U: (uncommment if filter is desired)
  
// // recompute U
//   for (int i=0; i<nn[0]; i++)
//     {
//       for (int j=0; j<nn[1]; j++)
// 	{
// 	  pos = i*nn[1]+j;
// 	  temp2 = *A[pos];
// 	  temp2.expm();
// 	  // set U
// 	  lat->cells[pos]->setU(temp2);
// 	}
//     }
        
  for(int i=0; i<nn[0]*nn[1]; i++) 
    {
      delete A[i];
    }
  
  delete[] A;
      
  // done. 
}
