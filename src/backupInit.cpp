// Init.cpp is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.
#include "Init.h"

//**************************************************************************
// Init class.


//**************************************************************************

//computes different values for A in all n steps
void Init::computeAx(Lattice *lat, Group *group, Parameters *param, Random *random)
{
  //compute C_k^A
  int size = lat->getSize();
  int Nc = param->getNc();
  int Nc2;
  Nc2 = param->getNc()*param->getNc()-1;
  int nn[2];
  int pos;
  nn[0]=param->getSize();
  nn[1]=param->getSize();

  Matrix** Cx;
  Cx = new Matrix*[param->getSize()*param->getSize()];

   for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      Cx[i] = new Matrix(1);
    }
 
  double x, y, R;

  // initialize C_x
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  x = lat->cells[pos]->getX();
	  y = lat->cells[pos]->getY();
	  //	  Cx[pos]->setRe(0,((x*x+y*y)/(4.*param->getR()*param->getR()))); // this is -log(C_x), x*x+y*y is \vec{x}^2, C_x=\vec{x}^2/(4R^2)
	  Cx[pos]->setRe(0,y*y); // this is -log(C_x), x*x+y*y is \vec{x}^2, C_x=\vec{x}^2/(4R^2)
	  Cx[pos]->setIm(0,0.);  
	}
    }
	  //	  fouta << j << " " << Aa[0][pos]->getRe(0)*(Aa[0][pos]->getRe(0)) << endl;

//  for (int i=0; i<nn[0]; i++)
//     {
//       for (int j=0; j<nn[1]; j++)
// 	{
// 	  pos = i*nn[1]+j;
// 	  x = lat->cells[pos]->getX();
// 	  y = lat->cells[pos]->getY();
// 	  cout << x << " " << y << " " << *Cx[pos] << endl;  
// 	}
//     }
//   cout << endl;

//   fstream fouta("Cx",ios::out); 
//   for (int i=0; i<nn[0]; i++)
//     {
//       for (int j=0; j<nn[1]; j++)
// 	{
// 	  pos = i*nn[1]+j;
// 	  x = lat->cells[pos]->getX();
// 	  y = lat->cells[pos]->getY();
// 	  fouta << x << " " << y << " " << Cx[pos]->getRe(0) << " " << Cx[pos]->getIm(0) << endl;  
// 	}
//       fouta << endl;
//     }

//   cout << endl;
//   fouta.close();

  
  fft->fftn(Cx,nn,2,1);

//   fstream foutb("Ck",ios::out); 
//   for (int i=0; i<nn[0]; i++)
//     {
//       for (int j=0; j<nn[1]; j++)
// 	{
// 	  pos = i*nn[1]+j;
// 	  x = lat->cells[pos]->getX();
// 	  y = lat->cells[pos]->getY();
// 	  if(x==0 && y==0)
// 	    {
// 	      Cx[pos]->setRe(0,0.);
// 	      Cx[pos]->setIm(0,0.);
// 	    }	 
// 	  foutb << 2.*param->PI*(-0.5+static_cast<double>(i)/static_cast<double>(nn[0])) << " " 
// 		<< 2.*param->PI*(-0.5+static_cast<double>(j)/static_cast<double>(nn[1])) << " " << Cx[pos]->getRe(0) << " " << Cx[pos]->getIm(0) << endl;  
// 	  cout << 2.*param->PI*(-0.5+static_cast<double>(i)/static_cast<double>(nn[0])) << " " 
// 	       << 2.*param->PI*(-0.5+static_cast<double>(j)/static_cast<double>(nn[1])) << " " << Cx[pos]->getRe(0) << " " << Cx[pos]->getIm(0) << endl;  
// 	}
//       foutb << endl;
//     }

//   cout << endl;
//   foutb.close();
  

  // now Cx contains the C_k^A, they're real, so take the square root trivially:
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  if(Cx[pos]->getIm(0)>pow(10.,-12))
	    cout << "Im(C_k^A)=" << Cx[pos]->getIm(0) << endl;
	  //cout << "C_k=" << *Cx[pos] << endl;
	  x = Cx[pos]->getRe(0);
	  y = Cx[pos]->getIm(0);
	  
	  Cx[pos]->setRe(0,1/sqrt(2.)*sqrt(sqrt(x*x+y*y)+x));
	    
	  if(y<0)
	    Cx[pos]->setIm(0,-1/sqrt(2.)*sqrt(sqrt(x*x+y*y)-x));
	  else
	    Cx[pos]->setIm(0,1/sqrt(2.)*sqrt(sqrt(x*x+y*y)-x));
	  
	  //cout << "sqrt(C_k)=" << *Cx[pos] << endl;
	}
    }
  
  // -----------------------------------------------------------
  // compute A_k
  // first, define Nc^2-1 coefficients A_a:
  
  Matrix*** Aa;
  Aa = new Matrix**[Nc2];
  
  for(int n=0; n<Nc2; n++)
    {
      Aa[n] = new Matrix*[param->getSize()*param->getSize()]; // the coefficients in Aa*t^a
    }
  
  for(int n=0; n<Nc2; n++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  Aa[n][i] = new Matrix(1); // the coefficients in Aa*t^a
	}
    }
  
  double xir[param->getSize()*param->getSize()][Nc2];
  double xii[param->getSize()*param->getSize()][Nc2];
  
  int pos2, j2, i2;

  for(int n=0; n<Nc2; n++)
    {
      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<nn[1]; j++)
	    {
	      pos = i*nn[1]+j;
	      xir[pos][n]=random->Gauss(); //Re(xi)
	      xii[pos][n]=random->Gauss(); //Im(xi)
	      Aa[n][pos]->setRe(0,Cx[pos]->getRe(0)*xir[pos][n]-Cx[pos]->getIm(0)*xii[pos][n]); // Aa=sqrt(C_k)*xi_k
	      Aa[n][pos]->setIm(0,Cx[pos]->getIm(0)*xir[pos][n]+Cx[pos]->getRe(0)*xii[pos][n]); 
	    }
	}
      // here I put in the correct mirroring so that A_x will be real.
//    for (int i=1; i<nn[0]/2; i++)
// 	{
// 	  for (int j=1; j<nn[1]; j++)
// 	    {
// 	      pos = i*nn[1]+j;
// 	      i2 = nn[0]-i;
// 	      j2 = nn[1]-j;
// 	      pos2 = i2*nn[1]+j2;
// 	      Aa[n][pos2]->setRe(0,Aa[n][pos]->getRe(0)); // Aa=sqrt(C_k)*xi_k
// 	      Aa[n][pos2]->setIm(0,-Aa[n][pos]->getIm(0)); // Aa=sqrt(C_k)*xi_k
// 	    }
// 	}
//       for (int i=nn[0]/2+1; i<nn[0]; i++)
// 	{
// 	  for (int j=0; j<1; j++)
// 	    {
// 	      pos = i*nn[1]+j;
// 	      xir[pos][n]=random->Gauss(); //Re(xi)
// 	      xii[pos][n]=random->Gauss(); //Im(xi)
// 	      Aa[n][pos]->setRe(0,Cx[pos]->getRe(0)*xir[pos][n]-Cx[pos]->getIm(0)*xii[pos][n]); // Aa=sqrt(C_k)*xi_k
// 	      Aa[n][pos]->setIm(0,Cx[pos]->getIm(0)*xir[pos][n]+Cx[pos]->getRe(0)*xii[pos][n]); 
// 	    }
// 	}
//       for (int i=nn[0]/2; i<nn[0]/2+1; i++)
// 	{
// 	  for (int j=nn[1]/2+1; j<nn[1]; j++)
// 	    {
// 	      pos = i*nn[1]+j;
// 	      j2 = nn[1]-j;
// 	      pos2 = i*nn[1]+j2;
// 	      Aa[n][pos2]->setRe(0,Aa[n][pos]->getRe(0)); // Aa=sqrt(C_k)*xi_k
// 	      Aa[n][pos2]->setIm(0,-Aa[n][pos]->getIm(0)); // Aa=sqrt(C_k)*xi_k
// 	    }
// 	}
    }
  // -----------------------------------------------------------

//   cout.precision(4);
//   for (int i=0; i<nn[0]; i++)
//     {
//       for (int j=0; j<nn[1]; j++)
// 	{
// 	  pos = i*nn[1]+j;
// 	  cout << setw(12) << *Aa[0][pos];
// 	}
//       cout << endl;
//       cout << endl;
//     }
//   cout << endl;
//   cout << endl;


  // -----------------------------------------------------------
  // get A_x^a by Fourier transform
  for(int n=0; n<Nc2; n++)
    {
      fft->fftn(Aa[n],nn,2,-1);
    }
  // -----------------------------------------------------------

//   for(int n=0; n<Nc2; n++)
//     {
//       for (int i=0; i<nn[0]; i++)
// 	{
// 	  for (int j=0; j<nn[1]; j++)
// 	    {
// 	      pos = i*nn[1]+j;
// 	      if (Aa[n][pos]->getIm(0)<pow(10.,-15))
// 		Aa[n][pos]->setIm(0,0.);
// 	      //cout << i << " " << j << " " << *Aa[n][pos] << endl;
// 	    }
// 	}
//     }

//    cout.precision(4);
//    for (int i=0; i<nn[0]; i++)
//      {
//        for (int j=0; j<nn[1]; j++)
//  	{
//  	  pos = i*nn[1]+j;
//  	  cout << setw(12) << *Aa[0][pos];
//  	}
//        cout << endl;
//        cout << endl;
//      }
//    cout << endl;
//    cout << endl;


  
  // -----------------------------------------------------------
  // compute A_x=A_x^a*t^a
  
  for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      A[i] = new Matrix(param->getNc()); // sets all entries 0.
    }
  
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  for(int n=0; n<Nc2; n++)
	    {
	      *A[pos]+=(Aa[n][pos]->getRe(0))*group->getT(n);
	    }
	}
    }

  for (int i=nn[0]/2; i<nn[0]/2+1; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  As[j] += Aa[0][pos]->getRe(0)*Aa[0][nn[1]/4]->getRe(0);
	}
    }



  // -----------------------------------------------------------
  // finalize --------------------------------------------------
  // -----------------------------------------------------------
  for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      delete Cx[i];
    }
  
  for(int n=0; n<Nc2; n++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  delete Aa[n][i];
	}
    }
  for(int n=0; n<Nc2; n++)
    {
      delete [] Aa[n];
    }	

  delete[] Cx;
  delete[] Aa;
}

//computes different values for A in all n steps
void Init::computeAx2(Lattice *lat, Group *group, Parameters *param, Random *random)
{
  //compute C_k^A
  int size = lat->getSize();
  int Nc = param->getNc();
  int Nc2;
  Nc2 = param->getNc()*param->getNc()-1;
  int nn[2];
  int pos;
  nn[0]=param->getSize();
  nn[1]=param->getSize();

  Matrix** Cx;
  Cx = new Matrix*[param->getSize()*param->getSize()];

   for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      Cx[i] = new Matrix(1);
    }

  double x, y, R;

  // initialize C_x
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  x = lat->cells[pos]->getX();
	  y = lat->cells[pos]->getY();
	  Cx[pos]->setRe(0,((x*x+y*y)/(4.*param->getR()*param->getR()))); // this is -log(C_x), x*x+y*y is \vec{x}^2, C_x=\vec{x}^2/(4R^2)
	  Cx[pos]->setIm(0,0.);  
	}
    }
  
  fft->fftn(Cx,nn,2,1);

  // generate uncorrelated noise
  Matrix*** xi;
  xi = new Matrix**[Nc2];
  Matrix*** xi2;
  xi2 = new Matrix**[Nc2];
  
  for(int n=0; n<Nc2; n++)
    {
      xi[n] = new Matrix*[param->getSize()*param->getSize()];
      xi2[n] = new Matrix*[param->getSize()*param->getSize()];
    }
  
  for(int n=0; n<Nc2; n++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  xi[n][i] = new Matrix(1);
	  xi2[n][i] = new Matrix(1);
	}
    }

  int pos2, j2, i2;

  for(int n=0; n<Nc2; n++)
    {
      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<nn[1]; j++)
	    {
	      pos = i*nn[1]+j;
	      xi[n][pos]->setRe(0,random->Gauss()); //Re(xi)
	      xi[n][pos]->setIm(0,random->Gauss()); //Im(xi)
	      xi2[n][pos]->setRe(0,xi[n][pos]->getRe(0)); //Re(xi)
	      xi2[n][pos]->setIm(0,xi[n][pos]->getIm(0)); //Im(xi)
	    }
	}
    }

  for(int n=0; n<Nc2; n++)
    {
      fft->fftn(xi[n],nn,2,1);
    }
  
  //square root of C_k
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  x = Cx[pos]->getRe(0);
	  y = Cx[pos]->getIm(0);
	  Cx[pos]->setRe(0,1/sqrt(2.)*sqrt(sqrt(x*x+y*y)+x));
	    
	  if(y<0)
	    Cx[pos]->setIm(0,-1/sqrt(2.)*sqrt(sqrt(x*x+y*y)-x));
	  else
	    Cx[pos]->setIm(0,1/sqrt(2.)*sqrt(sqrt(x*x+y*y)-x));
	}
    }
  
  Matrix*** B;
  B = new Matrix**[Nc2];
  
  for(int n=0; n<Nc2; n++)
    {
      B[n] = new Matrix*[param->getSize()*param->getSize()];
    }
  
  for(int n=0; n<Nc2; n++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  B[n][i] = new Matrix(1);
	}
    }
  
  double absXi;
  for(int n=0; n<Nc2; n++)
    {
      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<nn[1]; j++)
	    {
	      pos = i*nn[1]+j;
	      absXi = sqrt(xi[n][pos]->getRe(0)*xi[n][pos]->getRe(0)+xi[n][pos]->getIm(0)*xi[n][pos]->getIm(0));
	      B[n][pos]->setRe(0,Cx[pos]->getRe(0)/absXi); 
	      B[n][pos]->setIm(0,Cx[pos]->getIm(0)/absXi); 
	    }
	}
    }

  for(int n=0; n<Nc2; n++)
    {
      fft->fftn(B[n],nn,2,-1);
    }
 
  // -----------------------------------------------------------
  // compute A_k
  // first, define Nc^2-1 coefficients A_a:
  
  Matrix*** Aa;
  Aa = new Matrix**[Nc2];
  
  for(int n=0; n<Nc2; n++)
    {
      Aa[n] = new Matrix*[param->getSize()*param->getSize()]; // the coefficients in Aa*t^a
    }
  
  for(int n=0; n<Nc2; n++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  Aa[n][i] = new Matrix(1); // the coefficients in Aa*t^a
	}
    }

  double t1, t2;
  int pos3;
  for(int n=0; n<Nc2; n++)
    {
      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<nn[1]; j++)
	    {
	      pos = i*nn[1]+j;
	      t1=0.;
	      t2=0.;
 	      for (int l=0; l<nn[0]; l++)
		{
		  for (int m=0; m<nn[1]; m++)
		    {
		      if(i-l<0)
			{
			  if(j-m<0)
			    pos2 = (i-l+nn[0])*nn[1]+(j-m+nn[1]);
			  else
			    pos2 = (i-l+nn[0])*nn[1]+(j-m);
			}
		      else
			{
			  if(j-m<0)
			    pos2 = (i-l)*nn[1]+(j-m+nn[1]);
			  else
			    pos2 = (i-l)*nn[1]+(j-m);
			}
		      pos3 = l*nn[1]+m;
		      
		      t1+=B[n][pos2]->getRe(0)*xi2[n][pos3]->getRe(0)-B[n][pos2]->getIm(0)*xi2[n][pos3]->getIm(0); 
		      t2+=B[n][pos2]->getIm(0)*xi2[n][pos3]->getRe(0)+B[n][pos]->getRe(0)*xi2[n][pos]->getIm(0); 
		      //cout << t1 << endl;
		    }
		}
	      Aa[n][pos]->setRe(0,t1); 
	      Aa[n][pos]->setIm(0,t2); 
	    }
	}
    }

  // -----------------------------------------------------------
  // compute A_x=A_x^a*t^a
  
  for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      A[i] = new Matrix(param->getNc()); // sets all entries 0.
    }
  
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  for(int n=0; n<Nc2; n++)
	    {
	      *A[pos]+=(Aa[n][pos]->getRe(0))*group->getT(n);
	    }
	}
    }

  for (int i=nn[0]/2; i<nn[0]/2+1; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  As[j] += Aa[0][pos]->getRe(0)*Aa[0][0]->getRe(0);
	}
    }



  // -----------------------------------------------------------
  // finalize --------------------------------------------------
  // -----------------------------------------------------------
  for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      delete Cx[i];
    }
  
  for(int n=0; n<Nc2; n++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  delete Aa[n][i];
	}
    }
  for(int n=0; n<Nc2; n++)
    {
      delete [] Aa[n];
    }	

  delete[] Cx;
  delete[] Aa;
}


void Init::initU(Lattice *lat, Group *group, Parameters *param, Random *random)
{
  int size = lat->getSize();
  int Nc = param->getNc();
  int Nc2;
  Nc2 = param->getNc()*param->getNc()-1;
  int pos;
  int nn[2];
  nn[0]=param->getSize();
  nn[1]=param->getSize();
  cout << "Initializing SU(" << Nc << ") U-fields in " << size 
       << " cells, using initialization method " << param->getInitMethod() << " ..." << endl;
  
  Matrix temp(Nc,1.);
  Matrix temp2(Nc,0);
  Matrix **tempC;
  
  //  cout << "temp=\n" << temp << endl;
  double xil[Nc];
    
  A = new Matrix*[param->getSize()*param->getSize()];
  tempC = new Matrix*[param->getSize()];
  
  for(int i=0; i<param->getSize(); i++)
    {
      tempC[i] = new Matrix(Nc,0.);
    }
 

  // test random distribution
  //   const int bins=200;
  //   const int samples=1000;
  //   int bin[bins];
  //   double length = 10.;
  //   double step = length/bins;
  //   fstream fout("graph",ios::out); 
  //   for(int i=0; i<bins; i++)
  //     {
  //       bin[i]=0; 
  //     }
  
  //   double mean = 1.;
  //   double width = 2.;
  
  //   for(int i=0; i<samples; i++)
  //     {
  //       for(int ci=0; ci<Nc; ci++)
  // 	{
  // 	  xil[ci] = random->Gauss(1,2); 
  // 	}
  
  //       pos = floor((xil[1]+length/2.)/step);
  //       if(pos<bins && pos>0) bin[pos]+=1;
  //     }
  
  //   for(int i=0; i<bins; i++)
  //     {
  //       fout << i*step-length/2.+step/2. << " " << static_cast<double>(bin[i])/samples/step << endl; 
  //     }
  //   fout.close();
  // ------------------------
  
  
  
  
  // -----------------------------------------------------------
  // compute U
  double temp3;
  int runs = 1;  


  for (int i=0; i<1; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  As[pos] = 0.;
	}
    }
  
  for (int r=0; r<runs; r++)
    {
      
      int n = 500; // steps in U_x <- e^{iA_x/\sqrt{n}} U_x
      
      for(int in=0; in<n; in++)
	{
	  // get a new A_x
	  computeAx(lat, group, param, random);
	  for (int i=0; i<nn[0]; i++)
	    {
	      for (int j=0; j<nn[1]; j++)
		{
		  pos = i*nn[1]+j;
		  //  cout << "A=" << *A[pos] << endl << endl;
		  // multiply by i/sqrt{n}:
		  for(int nc=0; nc<Nc*Nc; nc++)
		    {
		      temp3 = A[pos]->getRe(nc);
		      A[pos]->setRe(nc,-A[pos]->getIm(nc)/sqrt(n));
		      A[pos]->setIm(nc,temp3/sqrt(n));
		    }
		  // -------------
		  // get U and multiply it by exp(iA/sqrt(n))
		  //cout << "iA/sqrt(n)=" << *A[pos] << endl << endl;
		  temp2 = A[pos]->expm();
		  //cout << "exp(iA/sqrt(n))=" << temp2 << endl << endl;
		  temp = temp2 * lat->cells[pos]->getU();
		  
		  // set U
		  lat->cells[pos]->setU(temp);

		  if(i==0 && j ==0)
		    {
		      // cout << "U=" << lat->cells[pos]->getU() << endl << endl;
		      //cout << "U^dag U=\n" << lat->cells[pos]->getU().conjg()*lat->cells[pos]->getU() << endl << endl;
		    }
		}
	    }
	}   
      
      
      // -----------------------------------------------------------
      
      fstream fouta("A",ios::out); 
      for (int i=nn[0]/2; i<nn[0]/2+1; i++)
	{
	  for (int j=0; j<nn[1]; j++)
	    {
	      pos = j;
	      As[pos] /= n;
	      fouta << j << " " << As[pos] << endl;
	    }
	}
      fouta.close();

      
      // check correlation between U's
      
      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<1; j++)
	    {
	      pos = i*nn[1]+j;
	      temp2 = lat->cells[pos]->getU().conjg() * lat->cells[0]->getU();
	      *tempC[i] += temp2; 
	      //	  fouta << i << " " << setw(12) << lat->cells[pos]->getU().getRe(0) << endl;
	    }
	}
      
    }
  
      fstream fouta("corr.dat",ios::out); 
      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<1; j++)
	    {
	      *tempC[i] /= static_cast<double>(runs);
	      fouta << i << " " << setw(12) << tempC[i]->getRe(0) + tempC[i]->getRe(2) << " " 
		    << tempC[i]->getIm(0)+tempC[i]->getIm(2) << endl;
	      //	  fouta << i << " " << setw(12) << lat->cells[pos]->getU().getRe(0) << endl;
	    }
	}
      
      cout << endl << endl;
      fouta.close();
      
      
  
  for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      delete A[i];
    }
  
  delete[] A;

//   for(int i=0; i<param->getSize(); i++)
//     {
//       delete tempC[i];
//     }
  
//   delete[] A;
  delete[] tempC;
  cout << "done initU" << endl;
 }



