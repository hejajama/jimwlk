// Init.cpp is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.
#include "Init.h"
#include <complex>

using namespace std;
//**************************************************************************
// Init class.


//**************************************************************************
//computes different values for A in all n steps
void Init::computeAx(Lattice *lat, Group *group, Parameters *param, Random *random)
{
  //compute C_k^A
  int size = lat->getSize();
  int Nc = param->getNc();
  int Nc2m1;
  Nc2m1 = param->getNc()*param->getNc()-1;
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
	  // this is -log(C_x), x*x+y*y is \vec{x}^2, C_x=\vec{x}^2/(4R^2):
	  Cx[pos]->setRe(0,-(x*x+y*y)/(4.*param->getR()*param->getR())); 
	  //Cx[pos]->setIm(0,(x*x+y*y)/(4.*param->getR()*param->getR())); 
	  //Cx[pos]->setRe(0,(exp(-(x*x+y*y)/(4.*param->getR()*param->getR())))); 
	}
    }
  

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

//    fstream fouta("Cx",ios::out); 
// //    for (int i=0; i<nn[0]; i++)
// //      {
// //         for (int j=0; j<nn[1]; j++)
// //   	{
// //   	  pos = i*nn[1]+j;
// //   	  x = lat->cells[pos]->getX();
// //   	  y = lat->cells[pos]->getY();
// //   	  fouta << x << " " << y << " " << Cx[pos]->getRe(0) << " " << Cx[pos]->getIm(0) << endl;  
// //   	}
// //         fouta << endl;
// //       }
// //     fouta.close();

//    for (int j=0; j<nn[1]; j++)
//      {
//        pos = nn[0]/2*nn[1]+j;
//        x = lat->cells[pos]->getX();
//        y = lat->cells[pos]->getY();
//        fouta <<  y << " " << Cx[pos]->getRe(0) << " " << Cx[pos]->getIm(0) << endl;  
//      }
//    fouta.close();
   
  fft->fftn(Cx,Cx,nn,2,1);


//    cout.precision(4);
//    for (int i=0; i<nn[0]; i++)
//      {
//        for (int j=0; j<nn[1]; j++)
// 	 {
// 	   pos = i*nn[1]+j;
// 	   cout << setw(12) << *Cx[pos];
// 	 }
//        cout << endl;
//        cout << endl;
//      }
//    cout << endl;
//    cout << " ----------------------------------------------------------------------- " << endl;
   

//     Cx[nn[0]/2*nn[1]+nn[1]/2]->setRe(0,1000+Cx[nn[0]/2*nn[1]+nn[1]/2]->getRe(0));
//     //Cx[nn[0]/2*nn[1]+nn[1]/2]->setIm(0,0.);
//     fstream foutb("Ck",ios::out); 
//     for (int i=0; i<nn[0]; i++)
//       {
// 	for (int j=0; j<nn[1]; j++)
// 	{
// 	  pos = i*nn[1]+j;
// 	  x = lat->cells[pos]->getX();
// 	  y = lat->cells[pos]->getY();
// 	  foutb << 2.*param->PI*(-0.5+static_cast<double>(i)/static_cast<double>(nn[0])) << " " 
// 		<< 2.*param->PI*(-0.5+static_cast<double>(j)/static_cast<double>(nn[1])) << " " << Cx[pos]->getRe(0) << " " << Cx[pos]->getIm(0) << endl;  
// 	  //	  cout << 2.*param->PI*(-0.5+static_cast<double>(i)/static_cast<double>(nn[0])) << " " 
// 	  //       << 2.*param->PI*(-0.5+static_cast<double>(j)/static_cast<double>(nn[1])) << " " << Cx[pos]->getRe(0) << " " << Cx[pos]->getIm(0) << endl;  
// 	}
// 	foutb << endl;
//       }
    
//     //  cout << endl;
//     foutb.close();

//     fft->fftn(Cx,nn,2,-1);
    
//     fstream fouts("Cx2",ios::out); 
//     for (int j=0; j<nn[1]; j++)
//       {
// 	pos = nn[0]/2*nn[1]+j;
// 	x = lat->cells[pos]->getX();
// 	y = lat->cells[pos]->getY();
// 	fouts <<  y << " " << Cx[pos]->getRe(0) << " " << Cx[pos]->getIm(0) << endl;  
//       }
//     fouts.close();


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  // now Cx contains the C_k^A, take the square root:
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  //if(Cx[pos]->getIm(0)>pow(10.,-12))
	  //  cout << "Im(C_k^A)=" << Cx[pos]->getIm(0) << endl;
	  //cout << "C_k=" << *Cx[pos] << endl;
	  x = Cx[pos]->getRe(0);
	  y = Cx[pos]->getIm(0);
	  //cout << "y=" << y << endl;
	  y = 0.;
	  
	  Cx[pos]->setRe(0,1/sqrt(2.)*sqrt(sqrt(x*x+y*y)+x));
	    
 	  if(y<0)
  	    Cx[pos]->setIm(0,-1/sqrt(2.)*sqrt(sqrt(x*x+y*y)-x));
  	  else
  	    Cx[pos]->setIm(0,1/sqrt(2.)*sqrt(sqrt(x*x+y*y)-x));
  

	}
    }
  
  // -----------------------------------------------------------
  // compute A_k
  // first, define Nc^2-1 coefficients A_a:
  
  Matrix*** Aa;
  Matrix*** xi;
  Aa = new Matrix**[Nc2m1];
  xi = new Matrix**[Nc2m1];
  
  for(int n=0; n<Nc2m1; n++)
    {
      Aa[n] = new Matrix*[param->getSize()*param->getSize()]; // the coefficients in Aa*t^a
    }
  
  for(int n=0; n<Nc2m1; n++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  Aa[n][i] = new Matrix(1); // the coefficients in Aa*t^a
	}
    }

  for(int n=0; n<Nc2m1; n++)
    {
      xi[n] = new Matrix*[param->getSize()*param->getSize()];
    }

  for(int n=0; n<Nc2m1; n++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  xi[n][i] = new Matrix(1);
	}
    }

  double xir[param->getSize()*param->getSize()][Nc2m1];
  double xii[param->getSize()*param->getSize()][Nc2m1];
  
  int pos2, j2, i2;

//   for(int n=0; n<Nc2m1; n++)
//     {
//       for (int i=0; i<nn[0]; i++)
// 	{
// 	  for (int j=0; j<nn[1]; j++)
// 	    {
// 	      pos = i*nn[1]+j;
// 	      xir[pos][n]=random->Gauss(); //Re(xi)
// 	      xii[pos][n]=0.; //Im(xi)
// 	      xi[n][pos]->setRe(0,xir[pos][n]);
// 	      xi[n][pos]->setIm(0,xii[pos][n]);
// 	    }
// 	}
//     }

//   for(int n=0; n<Nc2m1; n++)
//     {
//       fft->fftn(xi[n],nn,2,1);
//     }


  for(int n=0; n<Nc2m1; n++)
    {
      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<nn[1]; j++)
	    {
	      pos = i*nn[1]+j;

	      xir[pos][n]=random->Gauss(); //Re(xi)
	      xii[pos][n]=random->Gauss(); //Im(xi)
	      xi[n][pos]->setRe(0,xir[pos][n]);
	      xi[n][pos]->setIm(0,xii[pos][n]);
	      
	      Aa[n][pos]->setRe(0,(param->getSize()*(Cx[pos]->getRe(0)*xi[n][pos]->getRe(0)-Cx[pos]->getIm(0)*xi[n][pos]->getIm(0)))); //Aa=sqrt(C_k)*xi_k
	      Aa[n][pos]->setIm(0,(param->getSize()*(Cx[pos]->getIm(0)*xi[n][pos]->getRe(0)+Cx[pos]->getRe(0)*xi[n][pos]->getIm(0)))); 
	      //param->getSize() is the normalization factor \alpha on page 222 of Mecke (ed.) 
	      //Hilfer: Statistical physics and spatial statistics 
		

	      if(i!=nn[0]/2 && j!=nn[1]/2)
		{
		  Aa[n][pos]->setRe(0,0.); 
		  Aa[n][pos]->setIm(0,0.); 
		}
	      if(i==0 || j==0)
		{
		  Aa[n][pos]->setRe(0,0.); 
		  Aa[n][pos]->setIm(0,0.); 
		}

	    }
	}

      //here I put in the correct mirroring so that A_x will be real.
      for (int i=1; i<=nn[0]/2; i++)
	{
	  for (int j=1; j<nn[1]; j++)
	    {
	      pos = i*nn[1]+j;
	      i2 = nn[0]-i;
	      j2 = nn[1]-j;
	      pos2 = i2*nn[1]+j2;
	      //cout << "pos=" << pos << " pos2=" << pos2 << endl;
	      //cout << "i=" << i << " i2=" << i2 << endl;
	      //cout << "j=" << j << " j2=" << j2 << endl;
	      Aa[n][pos2]->setRe(0,Aa[n][pos]->getRe(0));
	      Aa[n][pos2]->setIm(0,-Aa[n][pos]->getIm(0));

	      if(i==nn[0]/2 && j==nn[1]/2)
		{
		  //Aa[n][pos]->setRe(0,0.); 
		  Aa[n][pos]->setIm(0,0.); 
		}
	      // actually just the imaginary part should be set to zero here
	    }
	}

// -----


//       for (int i=0; i<nn[0]; i++)
// 	{
// 	  for (int j=0; j<nn[1]; j++)
// 	    {
// 	      pos = i*nn[1]+j;
// 	      if(i==nn[0]/2-4)
// 		{
// 		  Aa[n][pos]->setRe(0,Aa[n][nn[0]/2*nn[1]+j]->getRe(0)); 
// 		  Aa[n][pos]->setIm(0,Aa[n][nn[0]/2*nn[1]+j]->getIm(0)); 
// 		}
// 	      if(j==nn[1]/2-4)
// 		{
// 		  Aa[n][pos]->setRe(0,Aa[n][i*nn[1]+nn[1]/2]->getRe(0)); 
// 		  Aa[n][pos]->setIm(0,Aa[n][i*nn[1]+nn[1]/2]->getIm(0)); 
// 		}
	      
// 	      if(i!=nn[0]/2-4 && j!=nn[1]/2-4)
// 		{
// 		  Aa[n][pos]->setRe(0,0.); 
// 		  Aa[n][pos]->setIm(0,0.); 
// 		}
// 	    }
// 	}




//     for (int i=nn[0]/2+1; i<nn[0]; i++)
// 	{
// 	  for (int j=0; j<1; j++)
// 	    {
// 	      pos = i*nn[1]+j;
// 	      xir[pos][n]=1.;//random->Gauss(); //Re(xi)
// 	      xii[pos][n]=0.;//random->Gauss(); //Im(xi)
// 	      Aa[n][pos]->setRe(0,param->getSize()*(Cx[pos]->getRe(0)*xir[pos][n]-Cx[pos]->getIm(0)*xii[pos][n])); // Aa=sqrt(C_k)*xi_k
// 	      Aa[n][pos]->setIm(0,param->getSize()*(Cx[pos]->getIm(0)*xir[pos][n]+Cx[pos]->getRe(0)*xii[pos][n])); 
// 	    }
// 	}
 //      for (int i=nn[0]/2; i<nn[0]/2+1; i++)
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
    //   for (int i=0; i<1; i++)
// 	{
// 	  for (int j=0; j<nn[1]; j++)
// 	    {
// 	      pos = i*nn[1]+j;
// 	      pos2 = (i+1)*nn[1]+j;
// 	      Aa[n][pos]->setRe(0,Aa[n][pos2]->getRe(0)); // Aa=sqrt(C_k)*xi_k
// 	      Aa[n][pos]->setIm(0,-Aa[n][pos2]->getIm(0)); // Aa=sqrt(C_k)*xi_k
// 	    }
// 	}
//       for (int i=0; i<nn[0]; i++)
// 	{
// 	  for (int j=0; j<1; j++)
// 	    {
// 	      pos = i*nn[1]+j;
// 	      pos2= i*nn[1]+(j+1);
// 	      Aa[n][pos]->setRe(0,Aa[n][pos2]->getRe(0)); // Aa=sqrt(C_k)*xi_k
// 	      Aa[n][pos]->setIm(0,-Aa[n][pos2]->getIm(0)); // Aa=sqrt(C_k)*xi_k
// 	    }
// 	}
 
      // Aa[n][nn[0]/2*nn[1]+nn[1]/2]->setRe(0,0.);
      // Aa[n][nn[0]/2*nn[1]+nn[1]/2]->setIm(0,0.);
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
  for(int n=0; n<Nc2m1; n++)
    {
      fft->fftn(Aa[n],Aa[n],nn,2,-1);
    }
  // -----------------------------------------------------------

  //  cout << "after AA FFT" << endl;
//   for(int n=0; n<Nc2m1; n++)
//     {
//       for (int i=0; i<nn[0]; i++)
// 	{
// 	  for (int j=0; j<nn[1]; j++)
// 	    {
// 	      pos = i*nn[1]+j;
// 	      //   if (Aa[n][pos]->getIm(0)<pow(10.,-15))
// 	      //cout << i << " " << j << " " << *Aa[n][pos] << endl;
// 	      Aa[n][pos]->setIm(0,0.);
// 	    }
// 	}
//     }


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
      for (int n=0; n<Nc*Nc; n++)
	{
	  A[i]->setRe(n,0.); // sets all entries 0.
	  A[i]->setIm(n,0.); // sets all entries 0.
	}
    }
  
  int postemp;
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
 	  postemp = i-nn[0]/2;
 	  if (postemp<0)
 	    {
	      postemp+=nn[0];
 	    }
 	  pos2 = (postemp)*nn[1]+j-nn[1]/2;
 	  if (pos2<postemp*nn[1])
 	    pos2 += nn[1];
	  for(int n=0; n<Nc2m1; n++)
	    {
	      *A[pos]+=(Aa[n][pos]->getRe(0))*group->getT(n);
	    }
	}
    }

//  for (int i=nn[0]/2; i<nn[0]/2+1; i++)
//     {
//       for (int j=0; j<nn[1]; j++)
// 	{
// 	  pos = i*nn[1]+j;
// 	  postemp = i-nn[0]/2;
// 	  if (postemp<0)
// 	    {
// 	      postemp+=nn[0];
// 	    }
// 	  pos2 = (postemp)*nn[1]+j-nn[1]/2;
// 	  if (pos2<postemp*nn[1])
// 	    pos2 += nn[1];
// 	  As[j] += Aa[0][pos2]->getRe(0)*Aa[0][nn[0]/2*nn[1]]->getRe(0);
// 	}
//    }

  for (int i=1; i<2; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  postemp = i-nn[0]/2;
	  if (postemp<0)
	    {
	      postemp+=nn[0];
	    }
	  pos2 = (postemp)*nn[1]+j-nn[1]/2;
	  if (pos2<postemp*nn[1])
	    pos2 += nn[1];
	  //	  As[i] += Aa[0][pos2]->getRe(0)*Aa[0][nn[1]/2]->getRe(0);
	  As[j] += (A[pos]->getRe(0)*A[0]->getRe(0));
	}
    }



  // -----------------------------------------------------------
  // finalize --------------------------------------------------
  // -----------------------------------------------------------
  for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      delete Cx[i];
    }
  
  for(int n=0; n<Nc2m1; n++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  delete Aa[n][i];
	}
    }
  for(int n=0; n<Nc2m1; n++)
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
  int Nc2m1;
  Nc2m1 = param->getNc()*param->getNc()-1;
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
	  Cx[pos]->setRe(0,(x*x+y*y)/(4.*param->getR()*param->getR())); // this is -log(C_x), x*x+y*y is \vec{x}^2, C_x=\vec{x}^2/(4R^2)
	  //Cx[pos]->setRe(0,exp(-(x*x+y*y)/(4.*param->getR()*param->getR()))); // this is -log(C_x), x*x+y*y is \vec{x}^2, C_x=\vec{x}^2/(4R^2)
	  //Cx[pos]->setRe(0,(4.*exp(-(x*x+y*y)/(4.*param->getR()*param->getR())))); 
	  Cx[pos]->setIm(0,0.); // this is -log(C_x), x*x+y*y is \vec{x}^2, C_x=\vec{x}^2/(4R^2)
	}
    }
  
  fft->fftn(Cx,Cx,nn,2,1);

  // generate uncorrelated noise
  Matrix*** xi;
  xi = new Matrix**[Nc2m1];
  Matrix*** xi2;
  xi2 = new Matrix**[Nc2m1];
  
  for(int n=0; n<Nc2m1; n++)
    {
      xi[n] = new Matrix*[param->getSize()*param->getSize()];
      xi2[n] = new Matrix*[param->getSize()*param->getSize()];
    }
  
  for(int n=0; n<Nc2m1; n++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  xi[n][i] = new Matrix(1);
	  xi2[n][i] = new Matrix(1);
	}
    }

  int pos2, j2, i2;

  for(int n=0; n<Nc2m1; n++)
    {
      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<nn[1]; j++)
	    {
	      pos = i*nn[1]+j;
	      //   xi[n][pos]->setRe(0,random->genrand64_real1()); //Re(xi)
	      //   xi[n][pos]->setIm(0,random->genrand64_real1()); //Im(xi)
	      xi[n][pos]->setRe(0,random->Gauss()); //Re(xi)
	      //xi[n][pos]->setIm(0,random->Gauss()); //Im(xi)
	      xi[n][pos]->setIm(0,0.); //Im(xi)
	      xi2[n][pos]->setRe(0,xi[n][pos]->getRe(0)); //Re(xi)
	      xi2[n][pos]->setIm(0,xi[n][pos]->getIm(0)); //Im(xi)
	    }
	}
    }

  for(int n=0; n<Nc2m1; n++)
    {
      fft->fftn(xi[n],xi[n],nn,2,1);
    }
  
  //square root of C_k
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  x = Cx[pos]->getRe(0);
	  y = 0.;//Cx[pos]->getIm(0);
	  Cx[pos]->setRe(0,1/sqrt(2.)*sqrt(sqrt(x*x+y*y)+x));
	    
	  if(y<0)
	    Cx[pos]->setIm(0,-1/sqrt(2.)*sqrt(sqrt(x*x+y*y)-x));
	  else
	    Cx[pos]->setIm(0,1/sqrt(2.)*sqrt(sqrt(x*x+y*y)-x));
	}
    }
  
  Matrix*** B;
  B = new Matrix**[Nc2m1];
  
  for(int n=0; n<Nc2m1; n++)
    {
      B[n] = new Matrix*[param->getSize()*param->getSize()];
    }
  
  for(int n=0; n<Nc2m1; n++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  B[n][i] = new Matrix(1);
	}
    }
  
  // B=sqrt(C_k) normalized
  double absXi;
  for(int n=0; n<Nc2m1; n++)
    {
      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<nn[1]; j++)
	    {
	      pos = i*nn[1]+j;
	      absXi =sqrt(xi[n][pos]->getRe(0)*xi[n][pos]->getRe(0)+xi[n][pos]->getIm(0)*xi[n][pos]->getIm(0));
	      B[n][pos]->setRe(0,param->getSize()*Cx[pos]->getRe(0)/absXi); 
	      B[n][pos]->setIm(0,param->getSize()*Cx[pos]->getIm(0)/absXi); 
	    }
	}
    }

  for(int n=0; n<Nc2m1; n++)
    {
      fft->fftn(B[n],B[n],nn,2,-1);
    }
 
  // -----------------------------------------------------------
  // compute A_k
  // first, define Nc^2-1 coefficients A_a:
  
  Matrix*** Aa;
  Aa = new Matrix**[Nc2m1];
  
  for(int n=0; n<Nc2m1; n++)
    {
      Aa[n] = new Matrix*[param->getSize()*param->getSize()]; // the coefficients in Aa*t^a
    }
  
  for(int n=0; n<Nc2m1; n++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  Aa[n][i] = new Matrix(1); // the coefficients in Aa*t^a
	}
    }

  double t1, t2;
  int pos3;
  for(int n=0; n<Nc2m1; n++)
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
		      t2+=B[n][pos2]->getIm(0)*xi2[n][pos3]->getRe(0)+B[n][pos2]->getRe(0)*xi2[n][pos3]->getIm(0); 
		      //cout << t1 << endl;
		    }
		}
	      Aa[n][pos]->setRe(0,t1); 
	      Aa[n][pos]->setIm(0,t2); 
	      // cout << Aa[n][pos]->getIm(0) << endl;
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
	  for(int n=0; n<Nc2m1; n++)
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
  
  for(int n=0; n<Nc2m1; n++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  delete Aa[n][i];
	}
    }
  for(int n=0; n<Nc2m1; n++)
    {
      delete [] Aa[n];
    }	

  delete[] Cx;
  delete[] Aa;
}

void Init::initU1(Lattice *lat, Group *group, Parameters *param, Random *random)
{
  int size = lat->getSize();
  int Nc = param->getNc();
  int Nc2m1;
  Nc2m1 = param->getNc()*param->getNc()-1;
  int pos;
  int nn[2];
  nn[0]=param->getSize();
  nn[1]=param->getSize();
  Matrix temp(Nc,1.);
  Matrix temp2(Nc,0.);
  Matrix Udag(Nc,0.);
  Matrix **tempC;
  Matrix **temptot;
  double xil[Nc];
    
  A = new Matrix*[param->getSize()*param->getSize()];
  for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      A[i] = new Matrix(Nc,0.);
    }

  tempC = new Matrix*[param->getSize()];
  
  for(int i=0; i<param->getSize(); i++)
    {
      tempC[i] = new Matrix(Nc,0.);
    }

  temptot = new Matrix*[param->getSize()*param->getSize()];
  
  for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      temptot[i] = new Matrix(Nc,0.);
    }
  
  // -----------------------------------------------------------
  // compute U
  double temp3;
  int runs = 1;
  Matrix unit(Nc,1.);

  for (int j=0; j<nn[1]; j++)
    {
      As[j] = 0.;
    }
  
  for (int r=0; r<runs; r++)
    {
      int n = 50; // steps in U_x <- e^{iA_x/\sqrt{n}} U_x
      
       if(Nc==2)
 	{
 	  temp.set(0,1.);
 	  temp.set(1,0.);
 	  temp.set(2,0.);
 	  temp.set(3,1.);
 	}
       else if (Nc==3)
 	{
 	  temp.set(0,1.);
 	  temp.set(1,0.);
 	  temp.set(2,0.);
 	  temp.set(3,0.);
 	  temp.set(4,1.);
 	  temp.set(5,0.);
 	  temp.set(6,0.);
 	  temp.set(7,0.);
 	  temp.set(8,1.);
 	}

      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<nn[1]; j++)
	    {
	      pos = i*nn[1]+j;
	      lat->cells[pos]->setU(temp);
	    }
	}
      for(int in=0; in<n; in++)
	{
	  // get a new A_x
	  computeAx(lat, group, param, random);
	  for (int i=0; i<nn[0]; i++)
	    {
	      for (int j=0; j<nn[1]; j++)
		{
		  pos = i*nn[1]+j;
		  //cout << "1 A=" << *A[pos] << endl << endl;
		  //multiply by i/sqrt{n}:
		  for(int nc=0; nc<Nc*Nc; nc++)
		    {
		      temp3 = A[pos]->getRe(nc);
		      A[pos]->setRe(nc,-A[pos]->getIm(nc)/sqrt(static_cast<double>(n)));
		      A[pos]->setIm(nc,temp3/sqrt(static_cast<double>(n)));
		    }
		  temp2 = *A[pos];
		  temp2.expm();
		  temp = temp2 * lat->cells[pos]->getU();
		  //temp = (unit-*A[pos]) * lat->cells[pos]->getU();
		  // set U
		  lat->cells[pos]->setU(temp);
		  lat->cells[pos]->setUi(temp);
		}
	    }
	}
      
      // -----------------------------------------------------------
      
      fstream fouta("A",ios::out); 
      for (int j=0; j<nn[1]; j++)
	{
	  pos = j;
	  As[pos] /= static_cast<double>(n);
	  fouta << j << " " << As[pos] << endl;
	}
      fouta.close();

         
      // check correlation between U's
      
      int pos2,pos3,counts=0;
      // compute correlation
      for (int i=0; i<nn[0]; i++)
	{  
	  for (int k=0; k<nn[1]; k++)
	    {
	      for (int j=0; j<nn[1]; j++)
		{
		  pos = i*nn[1]+j;
		  pos2 = i*nn[1]+k;
		  Udag = lat->cells[pos]->getU();
		  Udag.conjg(); // now Udag is U^dagger
		  temp2 = Udag * lat->cells[pos2]->getU();
		  //cout << "U=" << lat->cells[i*nn[1]]->getU() << endl;
		  //cout << "Udag=" << Udag << endl;
		  pos3 = abs(j-k); 
		  if (pos3==0) counts++;
		  if(pos3<=nn[1]/2)
		    {
		      *tempC[pos3] = *tempC[pos3]+temp2;
		    }
		  else 
		    {
		      *tempC[nn[1]-pos3] = *tempC[nn[1]-pos3]+temp2;
		    }
		  //cout << *tempC[j] << endl;
		  //	  fouta << i << " " << setw(12) << lat->cells[pos]->getU().getRe(0) << endl;
		}
	    }
	}
      
      for (int i=0; i<nn[1]; i++)
	{  
	  for (int k=0; k<nn[0]; k++)
	    {
	      for (int j=0; j<nn[0]; j++)
		{
		  pos = j*nn[1]+i;
		  pos2 = k*nn[1]+i;
		  Udag = lat->cells[pos]->getU();
		  Udag.conjg(); // now Udag is U^dagger
		  temp2 = Udag * lat->cells[pos2]->getU();
		  //cout << "U=" << lat->cells[i*nn[1]]->getU() << endl;
		  //cout << "Udag=" << Udag << endl;
		  pos3 = abs(j-k); 
		  if (pos3==0) counts++;
		  if(pos3<=nn[1]/2)
		    {
		      *tempC[pos3] = *tempC[pos3]+temp2;
		    }
		  else 
		    {
		      *tempC[nn[1]-pos3] = *tempC[nn[1]-pos3]+temp2;
		    }
		  //cout << *tempC[j] << endl;
		  //	  fouta << i << " " << setw(12) << lat->cells[pos]->getU().getRe(0) << endl;
		}
	    }
	}
    
      fstream foutc("corr.dat",ios::out); 
      for (int j=0; j<nn[1]; j++)
	{
	  *tempC[j] /= (static_cast<double>(r)+1);
	  foutc << j << " "  << tempC[j]->getRe(0) << " " 
		<< (2.-(tempC[j]->getRe(0) + tempC[j]->getRe(3)))/2. << " " << (-(tempC[j]->getIm(0) + tempC[j]->getIm(3)))/2. << endl;
	  *tempC[j] *= (static_cast<double>(r)+1);
	  
	}
      foutc.close();
 
      for (int j=0; j<=nn[1]/2; j++)
	{
	  *tempC[j]/=static_cast<double>(counts); 
	  if (j>0)
	    *tempC[j]/=2.; // because for all distances except 0 I double count. 
	}      
    
    }
  
  fstream foutc("corr.dat",ios::out); 
  for (int j=0; j<nn[1]; j++)
    {
      *tempC[j] /= static_cast<double>(runs);
      foutc << j << " "  << tempC[j]->getRe(0) << " " 
	    << (2.-(tempC[j]->getRe(0) + tempC[j]->getRe(3)))/2. << " " << (-(tempC[j]->getIm(0) + tempC[j]->getIm(3)))/2. << endl;
    }
  foutc.close();
    
  for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      delete A[i];
    }
  
  delete[] A;

  for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      delete temptot[i];
    }
  
  delete[] temptot;

   for(int i=0; i<param->getSize(); i++)
     {
       delete tempC[i];
     }
  
  delete[] tempC;
  cout << "done initU" << endl;
  //exit(1);
}

// ------------------------------------------------------------------------------


//SPATIAL CHARGE DISTRIBUTION
double Init::SpatialDistribution(int x,int y,double R, int N_T)
{
  //  return 10.*(exp(-((x-N_T/2)*(x-N_T/2)+(y-N_T/2)*(y-N_T/2))/(R*R)) + exp(-(((x+20)-N_T/2)*((x+20)-N_T/2)+((y+20)-N_T/2)*((y+20)-N_T/2))/(R*R))
  //	       + exp(-(((x-20)-N_T/2)*((x-20)-N_T/2)+((y+20)-N_T/2)*((y+20)-N_T/2))/(R*R)));        
  return exp(-((x-N_T/2)*(x-N_T/2)+(y-N_T/2)*(y-N_T/2))/(R*R));        
}


// MV initial conditions
void Init::initU2(Lattice *lat, Group *group, Parameters *param, Random *random)
{
  // initialize rho as Gaussian noise
  int Ny=param->getNy();
  int Nc2m1;
  Nc2m1 = param->getNc()*param->getNc()-1;
  double R = param->getR();
  double g2mu = param->getg2mu(); //g^2 mu
  double norm = g2mu*g2mu/static_cast<double>(Ny);
  int Nc = param->getNc();

  int nn[2];
  nn[0]=param->getSize();
  nn[1]=param->getSize();

  Matrix** rho;
    
  //  rho = new Matrix**[Ny]; // a rho for every longitudinal 'cell'
  // keep in mind that the normalization is such that \rho is actually g\rho (\rho is prop to g^2\mu here)
  int pos;
  int L=param->getSize();
  double fac = sqrt(norm);
      
  double m = param->getm();
  double length = param->getL();
  double a = length/static_cast<double>(nn[0]); // lattice spacing in fm

  m=m*a/param->hbarc;

  int pos2, pos3;
  int counts;
  double corr[L];
  
  double kt2, kx, ky;

  Matrix temp(Nc,1.);
  Matrix UD(Nc,1.);
  Matrix temp2(Nc,0.);
  double temp3;
   
  double correlator[nn[0]];
  double correlator4[nn[0]];
  double val[nn[0]][nn[0]];

//       // test correlator

//       for(int j=0; j<nn[0]; j++)
// 	{
// 	  for(int i=0; i<nn[0]; i++)
// 	    {
// 	      val[j][i]=(random->Gauss(0,2));
// 	    }
// 	}

//       for(int j=0; j<nn[0]; j++)
// 	{
// 	  for(int i=0; i<nn[0]; i++)
// 	    {
// 	      correlator[i]+=val[j][0]*val[j][i];
// 	      correlator4[i]+=val[j][0]*val[j][0]*val[j][i]*val[j][i];
// 	    }
// 	}


//       for(int j=0; j<5; j++)
// 	{
// 	  correlator[j]/=nn[0];
// 	  correlator4[j]/=nn[0];
// 	  cout << j << " " << correlator[j] << " " << correlator4[j] << endl;
// 	}

//       // end test


  for(int k=0; k<Ny; k++)
    {
      rho = new Matrix*[nn[0]*nn[1]];
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  rho[i] = new Matrix(param->getNc(),0.);
	}

      
      // compute \rho
      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<nn[1]; j++)
	    {
	      pos = i*nn[1]+j;
	      for(int n=0; n<Nc2m1; n++)
		{
		  *rho[pos]+=sqrt(SpatialDistribution(i,j,R,nn[0]))*fac*random->Gauss()*group->getT(n);
		  //		  *rho[k][pos]+=rho_a[k][n][pos]*group->getT(n);
		}
	    }
	}
      
      // Fourier transform rho
      fft->fftn(rho,rho,nn,2,1);

      // then compute A^+
      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<nn[1]; j++)
	    {
	      pos = i*nn[1]+j;
	      kx = 2.*param->PI*(-0.5+static_cast<double>(i)/static_cast<double>(nn[0]));
	      ky = 2.*param->PI*(-0.5+static_cast<double>(j)/static_cast<double>(nn[1]));
	      //kt2 = kx*kx+ky*ky;
	      //kt2 = 2.*sqrt(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.))*2.*sqrt(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); //lattice momentum
	      kt2 = 4.*(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); //lattice momentum
	      if (kt2!=0) // if m is set to zero this is only true if kt2!=0, otherwise it is always true
		*rho[pos] = *rho[pos]*(1./(kt2+m*m)); // rho contains A to save memory
	      else 
		*rho[pos] = *rho[pos]*(0.); // rho contains A to save memory
	    }
	}
      
      // Fourier transform back A^+
      fft->fftn(rho,rho,nn,2,-1);
	
      // A is saved in rho now.
    
      // compute U
      
      for (int i=0; i<nn[0]; i++)
	{
	  for (int j=0; j<nn[1]; j++)
	    {
	      pos = i*nn[1]+j;
	      //multiply by -i:
	      for(int nc=0; nc<Nc*Nc; nc++)
		{
		  temp3 = rho[pos]->getRe(nc); // rho contains A!
		  rho[pos]->setRe(nc,rho[pos]->getIm(nc));
		  rho[pos]->setIm(nc,-temp3);
		}
	      temp2 = *rho[pos];
	      //cout << temp2 << endl;
	      temp2.expm();
	      temp = temp2 * lat->cells[pos]->getU();
	      //temp = (unit-*A[pos]) * lat->cells[pos]->getU();
	      // set U
	      lat->cells[pos]->setU(temp);
	      lat->cells[pos]->setUi(temp); //this to keep the initial U's at every rapidity step (for unequal rapidity correlations)
	    }
	}
        
      //      UD=lat->cells[pos]->getU();
      //UD.conjg();

      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  delete rho[i];
	}
      
      delete[] rho;
    }


   // done. 

  // -----------------------------------------------------------------------------
  // finish
  // -----------------------------------------------------------------------------

}


// MV initial conditions with quartic term in the action
void Init::initU3(Lattice *lat, Group *group, Parameters *param, Random *random)
{
  // initialize rho as Gaussian noise
  int Ny=param->getNy();
  int Nc2m1;
  Nc2m1 = param->getNc()*param->getNc()-1;
  
  double g2mu = param->getg2mu(); //g^2 mu a
  double norm = g2mu*g2mu/static_cast<double>(Ny);
  int Nc = param->getNc();
  
  int nn[2];
  nn[0]=param->getSize();
  nn[1]=param->getSize();
  
  Matrix** rho;
  double*** rhoa;
  
  //  rho = new Matrix**[Ny]; // a rho for every longitudinal 'cell'
  
  int pos;
  int L=param->getSize();
  double fac = sqrt(norm);
  
  int pos2, pos3;
  int counts;
  double corr[L];
  
  double kt2, kx, ky;
  
  Matrix temp(Nc,1.);
  Matrix UD(Nc,1.);
  Matrix temp2(Nc,0.);
  double temp3;
  
  double correlator[nn[0]];
  double correlator4[nn[0]];
  double val[nn[0]][nn[0]];
  
  // keep in mind that the normalization was such that \rho is actually g\rho (\rho is prop to g^2\mu here)
  // in this scenario we have to pick g and \mu separately.
  // now define \rho without the g:
  double g=param->getg();
  fac/=g;
  
  // compute \rho
  // here check the value of S_4=1/kappa_4 3 \int d^2x \rho^a(x)\rho^a(x)\rho^b(x)\rho^b(x) (sum over a and b)
  double SG=0.;
  double S4=0.;
  double mu=g2mu/g/g; // this is mu times a
  double kappa4=param->getkappa4Factor()*pow(mu,6.)*pow(g,4.)/pow(g,2.);
  // kappa_4 ~ a^6 (good, since rho~a and the d^2x gives another a^2, such that all a's cancel in s_4)
  // the g^4 comes in because my rho's are proportional to g (and I multiply four of them in the 1/kappa_4 term)
  // I divide by another g^2 because of the relation kappa_4 ~ 144 mu^6/g^2, coming from
  // mu^2 = g^2 A/ (2\piR^2) and kappa_4 = 18 g^4A^3/(\pi R^2)^3 (1105.4155, (15) and (17))
  
  cout << "kappa4=" << kappa4 << endl;
  kappa4/=static_cast<double>(Ny); // divide kappa_4 by N_y to compensate for the additional Ny that comes from \rho\rho\rho\rho.
  cout << "kappa4/Ny=" << kappa4 << endl;
  double y;
  
  cout << "g " << g << " mu " << mu << " g2mu " << g2mu << " rho_normalization " << mu*g/sqrt(static_cast<double>(Ny)) << endl;
  exit(1);
  
  int accepted;
  double S4avg;
  int rejections=0;
  
  for(int k=0; k<Ny; k++)
  {
    rejections=0;
    rho = new Matrix*[nn[0]*nn[1]];
    rhoa = new double**[nn[0]*nn[1]];
    
    for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      rho[i] = new Matrix(param->getNc(),0.);
      rhoa[i] = new double*[Nc2m1];
    }
    
    for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      for(int n=0; n<Nc2m1; n++)
      {
        rhoa[i][n] = new double;
      }
    }
    
    S4avg=0.;
    
    for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
      {
        pos = i*nn[1]+j;
        accepted=0;
        SG=0.;
        S4=0.;
        while(accepted==0)
        {
          SG=0.;
          S4=0.;
          for(int n=0; n<Nc2m1; n++)
          {
            *rhoa[pos][n] = mu*g*random->Gauss()/sqrt(static_cast<double>(Ny)); // generate Gaussian random number
          }
          // compute Gaussian action and the quartic part of the action
          for(int a=0; a<Nc2m1; a++)
          {
            SG += (*rhoa[pos][a])*(*rhoa[pos][a])/2./g/g/mu/mu;
            for(int b=0; b<Nc2m1; b++)
            {
              S4 += 3./kappa4*(*rhoa[pos][a])*(*rhoa[pos][a])*(*rhoa[pos][b])*(*rhoa[pos][b]);
            }
          }
          // 		  cout << "cell " << pos << endl;
          //		  cout << "SG=" << SG << endl;
          // 		  cout << "S4=" << S4 << endl;
          // 		  cout << "exp(-S4) = " << exp(-S4) << endl;
          // 		  cout << "kappa4=" << kappa4 << endl;
          // 		  cout << "kappa4/(mu^2)^3=" << kappa4/pow(mu,6)/pow(g,6.) << endl;
          
          //accept rho^a's with probability exp(-S4), otherwise reject
          y=random->genrand64_real1();
          if (y<exp(-S4))
          {
            //cout << "*kappa4/(kappa4-3.*pow(mu,4.))=" << kappa4/(kappa4-3.*pow(mu,4.)) << endl;
            S4avg+=S4;
            accepted=1;
            //cout << " ***** Values at (" << i << "," << j << ") accepted with y=" << y << ", exp(-S4)=" << exp(-S4)
            //<< ", (1-S4)=" << (1-S4) << endl;
            
            //cout << "SG=" << SG << endl;
            //cout << "S4=" << S4 << endl;
            for(int n=0; n<Nc2m1; n++)
            {
              *rho[pos]+=(*rhoa[pos][n])*group->getT(n);
            }
          }
          else
          {
            //		cout << "Values at (" << i << "," << j << ") rejected with y=" << y << ", exp(-S4)=" << exp(-S4) << endl;
            rejections+=1;
          }
        }
      }
    }
    
    cout << rejections << " rejections." << endl;
    cout << "average S4=" << S4avg/nn[0]/nn[1] << endl;
    
    double C2[Nc2m1][Nc2m1];
    double C4[Nc2m1][Nc2m1];
    double C4m[Nc2m1][Nc2m1];
    
    for(int a=0; a<Nc2m1; a++)
    {
      for(int b=0; b<Nc2m1; b++)
      {
        C2[a][b]=0.;
        C4[a][b]=0.;
        C4m[a][b]=0.;
      }
    }
    
    int pos2;
    int C4mcount=0;
    
    for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
      {
        pos = i*nn[1]+j;
        pos2 = (nn[0]-1-i)*nn[1]+(nn[1]-1-j);
        for(int a=0; a<Nc2m1; a++)
        {
          for(int b=0; b<Nc2m1; b++)
          {
            //cout << i << " " << j << " " << a << " " << b << " " << (*rhoa[pos][a]) << " " << (*rhoa[pos][b]) << endl;
            C2[a][b] += (*rhoa[pos][a])*(*rhoa[pos][b]);
            C4[a][b] += (*rhoa[pos][a])*(*rhoa[pos][a])*(*rhoa[pos][b])*(*rhoa[pos][b]);
            if (pos!=pos2)
            {
              C4m[a][b] += (*rhoa[pos][a])*(*rhoa[pos][a])*(*rhoa[pos2][b])*(*rhoa[pos2][b]);
            }
          }
        }
        C4mcount+=1;
      }
    }
    
    for(int a=0; a<Nc2m1; a++)
    {
      for(int b=0; b<Nc2m1; b++)
      {
        C2[a][b]/=(nn[0]*nn[1]);
        C4[a][b]/=(nn[0]*nn[1]);
        C4m[a][b]/=C4mcount;
      }
    }
    cout << "mug=" << mu*g << endl;
    cout << "fac=" << fac << endl;
    cout << "Ny=" << Ny << endl;
    cout << "g^2mu^2=" << g*g*mu*mu << endl;
    cout << "kappa4=" << kappa4 << endl;
    cout << " ------------- " << endl;
    cout << "new g^2mu^2=" << g*g*mu*mu*(1-4.*mu*mu*mu*mu*g*g*g*g/(kappa4/3.)*(Nc*Nc+1.)) << endl;
    cout << "avg C2[a][a]=" << 1./static_cast<double>(Nc2m1)*(C2[0][0]+C2[1][1]+C2[2][2]+C2[3][3]+C2[4][4]+C2[5][5]+C2[6][6]+C2[7][7]) << endl;
    cout << "C2[1][1]=" << C2[1][1] << endl;
    cout << "C2[0][1]=" << C2[0][1] << endl;
    cout << " ------------- " << endl;
    
    
    cout << "perturbative C4(0000)="
	   << pow(g*g*mu*mu*(1-4.*mu*mu*mu*mu*g*g*g*g/(kappa4/3.)*(Nc*Nc+1)),2.)*(3*(1-8*pow(g*g*mu*mu*(1-4.*mu*mu*mu*mu*g*g*g*g/(kappa4/3.)*(Nc*Nc+1)),2.)/(kappa4/3.))) << endl;
    cout << "C4[0][0]=" << C4[0][0] << endl;
    cout << "perturbative C4(0011)="
	   << pow(g*g*mu*mu*(1-4.*mu*mu*mu*mu*g*g*g*g/(kappa4/3.)*(Nc*Nc+1)),2.)*((1-8*pow(g*g*mu*mu*(1-4.*mu*mu*mu*mu*g*g*g*g/(kappa4/3.)*(Nc*Nc+1)),2.)/(kappa4/3.))) << endl;
    cout << "C4[0][1]=" << C4[0][1] << endl;
    cout << "C4[2][4]=" << C4[2][4] << endl;
    cout << "perturbative C4(0x0x0u0u)="
	   << pow(g*g*mu*mu*(1-4.*mu*mu*mu*mu*g*g*g*g/(kappa4/3.)*(Nc*Nc+1)),2.) << endl;
    cout << "avg C4m[0][0]=" << 1./static_cast<double>(Nc2m1)*(C4m[0][0]+C4m[1][1]+C4m[2][2]+C4m[3][3]+C4m[4][4]+C4m[5][5]+C4m[6][6]+C4m[7][7]) << endl;
    
    // Fourier transform rho
    fft->fftn(rho,rho,nn,2,1);
    
    // then compute A^+
    for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
      {
        pos = i*nn[1]+j;
        kx = 2.*param->PI*(-0.5+static_cast<double>(i)/static_cast<double>(nn[0]));
        ky = 2.*param->PI*(-0.5+static_cast<double>(j)/static_cast<double>(nn[1]));
        //kt2 = kx*kx+ky*ky;
        kt2 = 2.*sqrt(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.))*2.*sqrt(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); //lattice momentum
        if (kt2!=0)
        *rho[pos] = *rho[pos]*(1./kt2); // rho contains A to save memory
        else
        *rho[pos] = *rho[pos]*(0.); // rho contains A to save memory
      }
    }
    
    // Fourier transform back A^+
    fft->fftn(rho,rho,nn,2,-1);
    
    // A is saved in rho now.
    
    // compute U
    
    for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
      {
        pos = i*nn[1]+j;
        //multiply by -i:
        for(int nc=0; nc<Nc*Nc; nc++)
        {
          temp3 = rho[pos]->getRe(nc); // rho contains A!
          rho[pos]->setRe(nc,rho[pos]->getIm(nc));
          rho[pos]->setIm(nc,-temp3);
        }
        temp2 = *rho[pos]*g; //   here rho is prop to g\mu, so we need another g here (unlike in initU2)
        //cout << temp2 << endl;
        temp2.expm();
        temp = temp2 * lat->cells[pos]->getU();
        //temp = (unit-*A[pos]) * lat->cells[pos]->getU();
        // set U
        lat->cells[pos]->setU(temp);
        lat->cells[pos]->setUi(temp);
      }
    }
    
    //      UD=lat->cells[pos]->getU();
    //UD.conjg();
    
    for(int i=0; i<param->getSize()*param->getSize(); i++)
    {
      delete rho[i];
    }
    
    delete[] rho;
  }
  
  
  // done. 
  
  // -----------------------------------------------------------------------------
  // finish
  // -----------------------------------------------------------------------------
  
}

// MV initial conditions with quartic term in the action and correlations in x^-
void Init::initU4(Lattice *lat, Group *group, Parameters *param, Random *random)
{
  // initialize rho as Gaussian noise
  int Ny=param->getNy();
  int Nc2m1;
  Nc2m1 = param->getNc()*param->getNc()-1;
  
  double g2mu = param->getg2mu(); //g^2 mu a
  double norm = g2mu*g2mu/static_cast<double>(Ny);
  int Nc = param->getNc();

  int nn[2];
  nn[0]=param->getSize();
  nn[1]=param->getSize();

  Matrix*** rho;
  double**** rhoa;
    
  //  rho = new Matrix**[Ny]; // a rho for every longitudinal 'cell'

  int pos;
  int L=param->getSize();
  double fac = sqrt(norm);
      
  int pos2, pos3;
  int counts;
  double corr[L];
  
  double kt2, kx, ky;

  Matrix temp(Nc,1.);
  Matrix UD(Nc,1.);
  Matrix temp2(Nc,0.);
  double temp3;
   
  double correlator[nn[0]];
  double correlator4[nn[0]];
  double val[nn[0]][nn[0]];

  // keep in mind that the normalization was such that \rho is actually g\rho (\rho is prop to g^2\mu here)
  // in this scenario we have to pick g and \mu separately.
  // now define \rho without the g:
  double g=param->getg();
  fac/=g;
  
  // compute \rho
  // here check the value of S_4=1/kappa_4 3 \int d^2x \rho^a(x)\rho^a(x)\rho^b(x)\rho^b(x) (sum over a and b)
  double SG=0.;
  double SGNew=0.;
  double S4=0.;
  double mu=g2mu/g/g; // this is mu times a
  double muNew; // shifted mu to improve rejection method
  double kappa4=param->getkappa4Factor();
  //param->getkappa4Factor()*pow(mu,6.)*pow(g,4.)/pow(g,2.); 
  // kappa_4 ~ a^6 (good, since rho~a and the d^2x gives another a^2, such that all a's cancel in s_4)
  // the g^4 comes in because my rho's are proportional to g (and I multiply four of them in the 1/kappa_4 term)
  // I divide by another g^2 because of the relation kappa_4 ~ 144 mu^6/g^2, coming from 
  // mu^2 = g^2 A/ (2\piR^2) and kappa_4 = 18 g^4A^3/(\pi R^2)^3 (1105.4155, (15) and (17))
  
  //  cout << "kappa4=" << kappa4 << endl;
  //kappa4/=static_cast<double>(Ny); // divide kappa_4 by N_y to compensate for the additional Ny that comes from \rho\rho\rho\rho.
  //cout << "kappa4/Ny=" << kappa4 << endl;
  
  double y;
  
  int accepted;
  double S4avg;
  int rejections=0;

  rejections=0;
  rho = new Matrix**[Ny];
  rhoa = new double***[Ny];

  for(int k=0; k<Ny; k++)
    {
      rho[k] = new Matrix*[nn[0]*nn[1]];
      rhoa[k] = new double**[nn[0]*nn[1]];
    }
  
  for(int k=0; k<Ny; k++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  rho[k][i] = new Matrix(param->getNc(),0.);
	  rhoa[k][i] = new double*[Nc2m1];
	}
    }
  
  for(int k=0; k<Ny; k++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  for(int n=0; n<Nc2m1; n++)
	    {
	      rhoa[k][i][n] = new double;
	    }
	}
    }    
  
  S4avg=0.;

 //  //shift mu to be closer to full action:
  
//   muNew = pow(1./pow(mu,2.)+3.*2.*pow(mu,2.)/Ny/Ny/kappa4*(Nc2m1),-0.5);
  
//   cout << "mu=" << mu << endl;
//   cout << "muNew=" << muNew << endl;

  
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  accepted=0;
	  SG=0.;
	  SGNew=0.;
	  S4=0.;
	  while(accepted==0)
	    {
	      SG=0.;
	      SGNew=0.;
	      S4=0.;
	      // save all the rhos at all the x^- positions in array
	      for(int k=0; k<Ny; k++)
		{
		  for(int n=0; n<Nc2m1; n++)
		    {
		      *rhoa[k][pos][n] = mu*g*random->Gauss()/sqrt(static_cast<double>(Ny)); // generate Gaussian random number 
		    }
		}
	      // compute Gaussian action and the quartic part of the action
	      for(int a=0; a<Nc2m1; a++)
		{
		  for(int k=0; k<Ny; k++)
		    {
		      SG += (*rhoa[k][pos][a])*(*rhoa[k][pos][a])/2./g/g/mu/mu;
		      //		      SGNew += (*rhoa[k][pos][a])*(*rhoa[k][pos][a])/2./g/g/muNew/muNew;
		    }
		  for(int b=0; b<Nc2m1; b++)
		    {
		      for(int k=0; k<Ny; k++) // two different x^- positions to sum over (the 1/Ny^2 from the 4 rho's takes care of \Delta x^-\Delta x^-)
			{
			  for(int l=0; l<Ny; l++)
			    {
			      S4 += 3./kappa4*(*rhoa[k][pos][a])*(*rhoa[k][pos][a])*(*rhoa[l][pos][b])*(*rhoa[l][pos][b]);
			    }
			}
		    }
		}

	      // 		  cout << "cell " << pos << endl;
	      //		  cout << "SG=" << SG << endl;
	      // 		  cout << "S4=" << S4 << endl;
	      // 		  cout << "exp(-S4) = " << exp(-S4) << endl;
	      // 		  cout << "kappa4=" << kappa4 << endl;
	      // 		  cout << "kappa4/(mu^2)^3=" << kappa4/pow(mu,6)/pow(g,6.) << endl;
	      
	      //accept rho^a's with probability exp(-Sfull)/exp(-SGnew), otherwise reject
	      y=random->genrand64_real1();
	      if (y<exp(-S4))
		{
		  cout << "*kappa4/(kappa4-3.*pow(mu,4.))=" << kappa4/(kappa4-3.*pow(mu,4.)) << endl;
		  S4avg+=S4;
		  accepted=1;
		  cout << " ***** Values at (" << i << "," << j << ") accepted with y=" << y << ", exp(-S4)=" << exp(-S4)
		       << ", exp(-SG-S4)/exp(-SGNew)=" << exp(-SG-S4)/exp(-SGNew) << endl;
		  
		  cout << "SG=" << SG << endl;
		  cout << "SGNew=" << SGNew << endl;
		  cout << "S4=" << S4 << endl;
		  for(int n=0; n<Nc2m1; n++)
		    {
		      for(int k=0; k<Ny; k++)
			{
			  *rho[k][pos]+=(*rhoa[k][pos][n])*group->getT(n);
			}
		    }
		}
	      else
		{
		  //		cout << "Values at (" << i << "," << j << ") rejected with y=" << y << ", exp(-S4)=" << exp(-S4) << endl;
		  rejections+=1;
		}
	    }
	}
    }
	  
  cout << rejections << " rejections." << endl;
  cout << "average S4=" << S4avg/nn[0]/nn[1] << endl;
  
  double C2[Nc2m1][Nc2m1];
  double C4[Nc2m1][Nc2m1];
  double C4m[Nc2m1][Nc2m1];
  
  for(int a=0; a<Nc2m1; a++)
    {
      for(int b=0; b<Nc2m1; b++)
	{
	  C2[a][b]=0.;
	  C4[a][b]=0.;
	  C4m[a][b]=0.;
	}
    }

  int C4mcount=0;

  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  pos2 = (nn[0]-1-i)*nn[1]+(nn[1]-1-j);
	  // for(int k=0; k<Ny; k++)
	  //{
	      for(int a=0; a<Nc2m1; a++)
		{
		  for(int b=0; b<Nc2m1; b++)
		    {
		      //cout << i << " " << j << " " << a << " " << b << " " << (*rhoa[pos][a]) << " " << (*rhoa[pos][b]) << endl;
		      C2[a][b] += (*rhoa[0][pos][a])*(*rhoa[0][pos][b]);
		      
		      
		      //      for(int l=0; l<Ny; l++)
		      //{
			  C4[a][b] += (*rhoa[0][pos][a])*(*rhoa[0][pos][a])*(*rhoa[0][pos][b])*(*rhoa[0][pos][b]);
			  if (pos!=pos2)
			    {
			      //	      C4m[a][b] += (*rhoa[0][pos][a])*(*rhoa[0][pos][a])*(*rhoa[1][pos2][b])*(*rhoa[1][pos2][b]);
			      C4mcount+=1;
			    }
			  //}
		    }
		}
	      //}
	}
    }
  
  for(int a=0; a<Nc2m1; a++)
    {
      for(int b=0; b<Nc2m1; b++)
	{
	  //	  C2[a][b]/=(nn[0]*nn[1])*Ny;
	  //C4[a][b]/=(nn[0]*nn[1])*Ny;
	  C2[a][b]/=(nn[0]*nn[1]);
	  C4[a][b]/=(nn[0]*nn[1]);
	  C4m[a][b]/=C4mcount;
	}
    }
  cout << "mug=" << mu*g << endl;
  cout << "fac=" << fac << endl;
  cout << "Ny=" << Ny << endl;
  cout << "g^2mu^2=" << g*g*mu*mu << endl;
  cout << "kappa4=" << kappa4 << endl;
  cout << " ------------- " << endl;
  cout << "new g^2mu^2=" << g*g*mu*mu*(1-4.*mu*mu*mu*mu*g*g*g*g/(kappa4/3.)*(Nc*Nc+1.)) << endl;
  cout << "avg C2[a][a]=" << 1./static_cast<double>(Nc2m1)*(C2[0][0]+C2[1][1]+C2[2][2]+C2[3][3]+C2[4][4]+C2[5][5]+C2[6][6]+C2[7][7]) << endl;
  cout << "C2[1][1]=" << C2[1][1] << endl;
  cout << "C2[0][1]=" << C2[0][1] << endl;
  cout << " ------------- " << endl;
  
  
  cout << "perturbative C4(0000)=" 
       << pow(g*g*mu*mu*(1-4.*mu*mu*mu*mu*g*g*g*g/(kappa4/3.)*(Nc*Nc+1)),2.)*(3*(1-8*pow(g*g*mu*mu*(1-4.*mu*mu*mu*mu*g*g*g*g/(kappa4/3.)*(Nc*Nc+1)),2.)/(kappa4/3.))) << endl;
  cout << "C4[0][0]=" << C4[0][0] << endl;
  cout << "perturbative C4(0011)=" 
       << pow(g*g*mu*mu*(1-4.*mu*mu*mu*mu*g*g*g*g/(kappa4/3.)*(Nc*Nc+1)),2.)*((1-8*pow(g*g*mu*mu*(1-4.*mu*mu*mu*mu*g*g*g*g/(kappa4/3.)*(Nc*Nc+1)),2.)/(kappa4/3.))) << endl;
  cout << "C4[0][1]=" << C4[0][1] << endl;
  cout << "C4[2][4]=" << C4[2][4] << endl;
  cout << "perturbative C4(0x0x0u0u)=" 
       << pow(g*g*mu*mu*(1-4.*mu*mu*mu*mu*g*g*g*g/(kappa4/3.)*(Nc*Nc+1)),2.) << endl;
  cout << "avg C4m[0][0]=" << 1./static_cast<double>(Nc2m1)*(C4m[0][0]+C4m[1][1]+C4m[2][2]+C4m[3][3]+C4m[4][4]+C4m[5][5]+C4m[6][6]+C4m[7][7]) << endl;
  
  // Fourier transform rho
  for(int k=0; k<Ny; k++)
    {
      fft->fftn(rho[k],rho[k],nn,2,1);
    }
  
  // then compute A^+
  for (int i=0; i<nn[0]; i++)
    {
      for (int j=0; j<nn[1]; j++)
	{
	  pos = i*nn[1]+j;
	  kx = 2.*param->PI*(-0.5+static_cast<double>(i)/static_cast<double>(nn[0]));
	  ky = 2.*param->PI*(-0.5+static_cast<double>(j)/static_cast<double>(nn[1]));
	  //kt2 = kx*kx+ky*ky;
	  kt2 = 2.*sqrt(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.))*2.*sqrt(sin(kx/2.)*sin(kx/2.)+sin(ky/2.)*sin(ky/2.)); //lattice momentum
	  for(int k=0; k<Ny; k++)
	    {
	      if (kt2!=0)
		*rho[k][pos] = *rho[k][pos]*(1./kt2); // rho contains A to save memory
	      else 
		*rho[k][pos] = *rho[k][pos]*(0.); // rho contains A to save memory
	    }
	}
    }
  // Fourier transform back A^+
   for(int k=0; k<Ny; k++)
    {
      fft->fftn(rho[k],rho[k],nn,2,-1);
    } 
  
   // A is saved in rho now.
   
   // compute U
   
   for (int i=0; i<nn[0]; i++)
     {
       for (int j=0; j<nn[1]; j++)
	 {
	   pos = i*nn[1]+j;
	   //multiply by -i:
	   for(int k=0; k<Ny; k++)
	     {
	       for(int nc=0; nc<Nc*Nc; nc++)
		 {
		   temp3 = rho[k][pos]->getRe(nc); // rho contains A!
		   rho[k][pos]->setRe(nc,rho[k][pos]->getIm(nc));
		   rho[k][pos]->setIm(nc,-temp3);
		 }
	       temp2 = *rho[k][pos]*g; //   here rho is prop to g\mu, so we need another g here (unlike in initU2)
	       //cout << temp2 << endl;
	       temp2.expm();
	       temp = temp2 * lat->cells[pos]->getU();
	       //temp = (unit-*A[pos]) * lat->cells[pos]->getU();
	       // set U
	       lat->cells[pos]->setU(temp);
	       lat->cells[pos]->setUi(temp);
	     }
	 }
     }
   //      UD=lat->cells[pos]->getU();
   //UD.conjg();

  for(int k=0; k<Ny; k++)
    {
      for(int i=0; i<param->getSize()*param->getSize(); i++)
	{
	  for(int n=0; n<Nc2m1; n++)
	    {
	      delete rhoa[k][i][n];
	    }
	}
    }    

   for(int k=0; k<Ny; k++)
     {
       for(int i=0; i<param->getSize()*param->getSize(); i++)
	 {
	   delete rho[k][i];
	   delete rhoa[k][i];
	 }
     }
  
   for(int k=0; k<Ny; k++)
     {
       delete rho[k];
       delete rhoa[k];
     }
   
   delete[] rho;
   delete[] rhoa;
   
   // done. 

  // -----------------------------------------------------------------------------
  // finish
  // -----------------------------------------------------------------------------

}


/// Init from file generated by IP-Glasma code
// H.M. 201608
void Init::initFromData(Lattice *lat, Group *group, Parameters *param, Random *random)
{
  cout << "Reading initial condition from file " << param->getInputWline() << endl;
  
    // TODO: should check that parameters are consistent with the datafile, now assuming it
  ifstream input(param->getInputWline().c_str());
  if (! input.good())
  {
	cerr << "File does not exist!" << endl;
	exit(1) ; 
}
  int size=param->getSize();
  //cout << "Reading " << size << " lines..." << endl;
  int i=0;
  while(!input.eof() )
  {
    Matrix m(3);
    string line;
    getline(input, line);
    if (line[0]=='#' or line.length() < 10)   // Comment line, or empty
    {
      continue;
    }
    
    stringstream l(line);
    // first two numbers are y and x coordinates,
    // followed by matrix elements (real,imag)
    int tmp; l>>tmp; l >> tmp;
    for (int j=0; j<9; j++)
    {
      double re,im;
      // real and imag parts
      l >> re;
      l >> im;
      
      complex<double> element(re,im);
      m.set(j, element);
      
    }
    lat->cells[i]->setU(m);
    lat->cells[i]->setUi(m);
    
    i++;
    
  }
  
  if (i != size*size)
  {
    cerr << "Expected " << size*size << " Wilson lines as an input, read " << i << ", check input file!" << endl;
    exit(1);
  }
    
}

void Init::initFromBinaryData(Lattice *lat, Group *group, Parameters *param, Random *random)
{
  cout << "Reading initial condition from file " << param->getInputWline() << endl;

  std::ifstream InStream;
  InStream.precision(10);
  InStream.open(param->getInputWline().c_str(), std::ios::in | std::ios::binary);
  int N;
  int Nc;
  double L,a;

  if(InStream.is_open())
  {
            // READING IN PARAMETERS -----------------------------------------------------------------//
      double temp;
      InStream.read(reinterpret_cast<char*>(&N), sizeof(int));
      InStream.read(reinterpret_cast<char*>(&Nc), sizeof(int));
      InStream.read(reinterpret_cast<char*>(&L), sizeof(double));
      InStream.read(reinterpret_cast<char*>(&a), sizeof(double));
      InStream.read(reinterpret_cast<char*>(&temp), sizeof(double));
                
      std::cout << "# BINARY Size is " << N << ", Nc " << Nc << ", length is [fm] " << L << ", a is [fm]" << a << ", yeff is " << temp  << std::endl;
                
                
      if(N != param->getSize())
      {
        std::cerr << "#---------------------------\n# ERROR wrong lattice size, data is " << N
         << " but you have specified " << param->getSize() << "\n#--------------------------" << std::endl;
         exit(0);
      }
                
                
      // READING ACTUAL DATA --------------------------------------------------------------------//
      double ValueBuffer;
      int INPUT_CTR=0;
      double re,im;
                
      while( InStream.read(reinterpret_cast<char*>(&ValueBuffer), sizeof(double)))
      {
          if(INPUT_CTR%2==0)              //this is the real part
          {
               re=ValueBuffer;
          }
          else                            // this is the imaginary part, write then to variable //
          {
              im=ValueBuffer;
                        
              int TEMPINDX=((INPUT_CTR-1)/2);
              int PositionIndx = TEMPINDX / 9;
                        
              int ix = PositionIndx / N;
              int iy = PositionIndx - N*ix;
                      
                        
              int MatrixIndx=TEMPINDX - PositionIndx*9;
              int j=MatrixIndx/3;
              int k=MatrixIndx-j*3;
                   
              int indx = N*iy + ix;
              lat->cells[indx]->getU().set(j,k, complex<double> (re,im)); 
              lat->cells[indx]->getUi().set(j,k, complex<double> (re,im)); 

          }
          INPUT_CTR++;
    }
  }

  else
  {
      std::cerr << "ERROR COULD NOT OPEN FILE";
      exit(0);
   }
   InStream.close();
    
        
   param->setSize(N);
   param->setNc(Nc);
   param->setL(L); 
}

void Init::initU(Lattice *lat, Group *group, Parameters *param, Random *random)
{
  cout << "Initializing SU(" << param->getNc() << ") U-fields in " << param->getSize()*param->getSize() 
       << " cells, using initialization method " << param->getInitMethod() << " ..." << endl;
 
  if(param->getInitMethod()==1)
    {
      initU1(lat, group, param, random);
    }
  else if(param->getInitMethod()==2)
    {
      initU2(lat, group, param, random);
    }
  else if(param->getInitMethod()==3)
    {
      initU3(lat, group, param, random);
    }
  else if(param->getInitMethod()==4)
    {
      initU4(lat, group, param, random);
    }
  else if (param->getInitMethod()==10)
    initFromData(lat, group, param, random);
  else if (param->getInitMethod()==11)
    initFromBinaryData(lat, group, param, random);

  else
    {
      cout << "[Init::initU]: No initialization method " << param->getInitMethod() 
	   << " found. Only 1(Gauss), 2(MV), and 3(MV+quartic) are available. Exiting." << endl;
      exit(1);
    }
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
  

