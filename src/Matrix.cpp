#include "Matrix.h"
#include <sstream>

//constructor if just dimension is given
Matrix::Matrix(int n)
{
  ndim = n;
  nn = ndim*ndim;
  //e = new complex<double> [nn];
  //for(int i=0; i<nn; i++) e[i] = complex<double>(0.0,0.0);
    for(int i=0; i<nn; i++) e.push_back(complex<double>(0.0,0.0));
}

//constructor if value for a and dimensions are given (a is the real value on the diagonal)
Matrix::Matrix(int n, double a)
{
  ndim = n;
  nn = ndim*ndim;
  //e = new complex<double> [nn];
    for(int i=0; i<nn; i++) e.push_back(complex<double>(0.0,0.0));
    /*
  if(e==0) 
    {
      cout << "Matrix: cannot allocate memory (matrix:) e= " << e << endl;
      abort();
    }*/
  
  //for(int i=0; i<nn; i++) e[i] = complex<double>(0.0,0.0);
    
   for(int i=0; i<ndim; i++) e[i*ndim+i] = complex<double>(a,0.0);
}


//operators:

Matrix operator * (const Matrix& a, const Matrix& b)
{
  int n = a.getNDim();
  Matrix c(n);
  if (n==2)
    {
      c.set(0,0,a(0,0)*b(0,0)+a(0,1)*b(1,0));
      c.set(0,1,a(0,0)*b(0,1)+a(0,1)*b(1,1));
      c.set(1,0,a(1,0)*b(0,0)+a(1,1)*b(1,0));
      c.set(1,1,a(1,0)*b(0,1)+a(1,1)*b(1,1));
      return c;
    }
  else if (n==3)
    {
      c.set(0,0,a(0,0)*b(0,0)+a(0,1)*b(1,0)+a(0,2)*b(2,0));
      c.set(0,1,a(0,0)*b(0,1)+a(0,1)*b(1,1)+a(0,2)*b(2,1));
      c.set(0,2,a(0,0)*b(0,2)+a(0,1)*b(1,2)+a(0,2)*b(2,2));
      c.set(1,0,a(1,0)*b(0,0)+a(1,1)*b(1,0)+a(1,2)*b(2,0));
      c.set(1,1,a(1,0)*b(0,1)+a(1,1)*b(1,1)+a(1,2)*b(2,1));
      c.set(1,2,a(1,0)*b(0,2)+a(1,1)*b(1,2)+a(1,2)*b(2,2));
      c.set(2,0,a(2,0)*b(0,0)+a(2,1)*b(1,0)+a(2,2)*b(2,0));
      c.set(2,1,a(2,0)*b(0,1)+a(2,1)*b(1,1)+a(2,2)*b(2,1));
      c.set(2,2,a(2,0)*b(0,2)+a(2,1)*b(1,2)+a(2,2)*b(2,2));
      return c;
    }
  else if (n==8)
    {
      for(int i=0; i<n; i++)
	{
	  for(int j=0; j<n; j++)
	    {
	      c.set(i,j,a(i,0)*b(0,j)+a(i,1)*b(1,j)+a(i,2)*b(2,j)+a(i,3)*b(3,j)+a(i,4)*b(4,j)+a(i,5)*b(5,j)+a(i,6)*b(6,j)+a(i,7)*b(7,j));
	    }
	}
      return c;
    }
  else
    {
      cout << "[Matrix::operator *]: Matrix product is only defined for 2x2, 3x3, and 8x8 matrixes. You gave me " << n << "x" << n << ". Exiting." << endl;
    }
}

//-
Matrix operator - (const Matrix& a, const Matrix& b)
{
  Matrix aa(a.getNDim());
  for(int i=0; i<a.getNN(); i++) aa.set(i,a(i) -b(i));
  return aa;
}

//+
Matrix operator + (const Matrix& a, const Matrix& b)
{
  Matrix aa(a.getNDim());
  for(int i=0; i<a.getNN(); i++) aa.set(i,a(i) + b(i));
  return aa;
}

//* multiply by a real scalar
Matrix operator * (const Matrix& a, const double s)
{
  Matrix aa(a.getNDim());
    for(int i=0; i<a.getNN(); i++){
      aa.set(i,a(i) * s);
    }
    return aa;
}
Matrix operator * (const double s, const Matrix& a)
{
  Matrix aa(a.getNDim());
    for(int i=0; i<a.getNN(); i++){
      aa.set(i,a(i) * s);
    }
    return aa;
}

//* multiply by a complex number
Matrix operator * (const complex<double> s,const Matrix& a)
{
  Matrix aa(a.getNDim());
    for(int i=0; i<a.getNN(); i++){
      aa.set(i,a(i) * s);
    }
    return aa;
}

// / division by scalar
Matrix operator / (const Matrix& a, const double s)
{
  Matrix aa(a.getNDim());
    for(int i=0; i<a.getNN(); i++) aa.set(i,a(i)/s);
    return aa;
}

Matrix& Matrix::conjg()
{
  // complex<double> temp[nn];
  vector < complex<double> > temp;

  for (int i=0; i<ndim; i++)
    for (int j=0; j<ndim; j++)
      {
	temp.push_back(0.);
      }  
  //save half the matrix
  for (int i=0; i<ndim; i++)
    for (int j=i; j<ndim; j++)
      {
	// 	temp[i*ndim+j]=conj(e[i*ndim+j]);
	temp.at(i*ndim+j) = (conj(e[i*ndim+j]));
      }
  //transpose half the matrix
  for (int i=0; i<ndim; i++)
    for (int j=0; j<i; j++)
      {
	e[j*ndim+i]=conj(e[i*ndim+j]);
      }
  //transpose the other half using the saved values
  for (int i=0; i<ndim; i++)
    for (int j=i; j<ndim; j++)
      {
	e[j*ndim+i]=temp[i*ndim+j];
      }
  return *this;
}

Matrix& Matrix::imag()
{
	if(ndim == 2) {
	    complex<double> e0 = e[0];
	    complex<double> e1 = e[1];
	    complex<double> e2 = e[2];
	    complex<double> e3 = e[3];
	    e[0] -= conj(e0);
	    e[1] -= conj(e2);
	    e[2] -= conj(e1);
	    e[3] -= conj(e3);
	} else {
	    cerr << " (Matrix::) invaid dimension nn= " << nn << endl;
	    exit(1);
	}
	return *this;
}


//matrix exponential using Pade approximant
//t is a scalar that multiplies the matrix (default: t=1) and p is the order in the Pade approximant (default: p=6)
Matrix& Matrix::expm(double t, const int p)
{
   const int n = this->getNDim();
   const Matrix I(n,1.);
   Matrix U(n),H2(n),P(n),Q(n);
   double norm = 0.0;
   // Calculate Pade coefficients
   if(p<6)
     {
       cout << "Matrix::expm: p should be at least 6. Exiting." << endl;
       exit(0);
     }
   // hard coded values for speed
   double c[p+1];
   c[0] = 1.;
   c[1] = 0.5;
   c[2] = 0.1136363636;
   c[3] = 0.01515151515;
   c[4] = 0.001262626263;
   c[5] = 6.313131313e-05;
   c[6] = 1.503126503e-06;
   if(p>6)
     {
       for(int i = 6; i < p; ++i) 
 	{	    
 	  c[i+1] = c[i] * ((p - i)/((i + 1.0) * (2.0 * p - i)));
 	  //std::cout << "c(" << i+1 << ")=" << c(i+1) << endl;
 	}
     }
   // Calculate the infinty norm of e, which is defined as the largest row sum of a matrix
   for(int i=0; i<n; ++i) 
     {
       double temp = 0.0;
       for(int j = 0; j < n; j++)
 	temp += abs((*this)(i,j)); 
       norm = t * max<double>(norm, temp);
     }
   // If norm = 0, and all H elements are not nan or infinity but zero, 
   // then U should be identity.
   if (norm == 0.0) 
     {
       bool all_H_are_zero = 1;
       for(int i = 0; i < n; i++)
	 {
	   for(int j = 0; j < n; j++)
	     {
	       if( (*this)(i,j) != 0.0 ) 
		 {
		   all_H_are_zero = 0;
		 }
	     }
	 }
       if( all_H_are_zero )
	 {
	   *this = I;
	   return *this;
	 }
       else
	 {
	   //	    Some error happens, H has elements which are NaN or infinity. 
	   cerr<<"Null input error in the template expm_pad.\n";
	   cout << "Null INPUT : " << *this << "\n";
	   exit(0);
	 }
     }

  // Scaling, seek s such that || e*2^(-s) || < 1/2, and set scale = 2^(-s)
  int s = 0;
  double scale = 1.0;
  if(norm > 0.5) 
    {
      s = max<int>(0, static_cast<int>((log(norm) / log(2.0) + 2.0)));
      scale /= double(pow(2.0, s)); // U <- this/2^s so ||U||_\infinity \approx 1
      U = (scale * t) * (*this); // Here U is used as temp value due to that H is const
    }
  else
    U = *this;
  
  // Horner evaluation of the irreducible fraction.
  // Initialize P (numerator) and Q (denominator) 
  H2 = U*U;
  Q = c[p]*I;
  P = c[p-1]*I;
  int odd = 1;
  
  for( int k = p - 1; k > 0; --k) 
    {
      if (odd==1)
	{
	  Q = Q*H2 + (c[k-1]*I);
	}
      else
	{
	  P = P*H2 + (c[k-1] * I);
	}	   
      odd = 1 - odd;
    }
  if (odd==1)
    {
      Q = Q*U;
    }
  else
    {
      P = P*U;
    }
  
  Q -= P;
  
  // Invert Q: (store inverse Q in H2)
  if (n==2)
    {
      H2.set(0,0,Q(1,1));
      H2.set(0,1,-Q(0,1));
      H2.set(1,0,-Q(1,0));
      H2.set(1,1,Q(0,0));
      H2 *= 1./(Q(0,0)*Q(1,1)-Q(0,1)*Q(1,0)); // divide by det(H) 
    }
  else if (n==3)
    {
      H2.set(0,0,(Q(1,1)*Q(2,2)-Q(1,2)*Q(2,1)));
      H2.set(0,1,(Q(0,2)*Q(2,1)-Q(0,1)*Q(2,2)));
      H2.set(0,2,(Q(0,1)*Q(1,2)-Q(0,2)*Q(1,1)));
      H2.set(1,0,(Q(1,2)*Q(2,0)-Q(1,0)*Q(2,2)));
      H2.set(1,1,(Q(0,0)*Q(2,2)-Q(0,2)*Q(2,0)));
      H2.set(1,2,(Q(0,2)*Q(1,0)-Q(0,0)*Q(1,2)));
      H2.set(2,0,(Q(1,0)*Q(2,1)-Q(1,1)*Q(2,0)));
      H2.set(2,1,(Q(0,1)*Q(2,0)-Q(0,0)*Q(2,1)));
      H2.set(2,2,(Q(0,0)*Q(1,1)-Q(0,1)*Q(1,0)));
      H2 *= 1./(Q(0,0)*Q(1,1)*Q(2,2) + Q(0,1)*Q(1,2)*Q(2,0) + Q(0,2)*Q(1,0)*Q(2,1) 
		- Q(0,2)*Q(1,1)*Q(2,0) - Q(0,1)*Q(1,0)*Q(2,2) - Q(1,2)*Q(2,1)*Q(0,0)); 
      // divide by det(H)
    }
  else
    {
      cerr << "My matrix exponential works only for up to 3x3 matrices. If you need more, include MatrixExp.h and the boost library." << endl;
    }
  
  if (odd == 1)
    {
      U = -1.*((2.0 * H2*P) +I );
    }
  else
    {
      U = (2.0 * H2*P) + I;
    }
  
  //square result = U^(2^s)
  for(int i = 0; i < s; ++i)
    U = U*U;
  
  *this = U;
  
  return *this;
}

complex<double> Matrix::det()
{
  int n = this->getNDim();
  Matrix Q(n);
  Q = *this;
  complex<double> det;
 
  if (n==2)
    {
      det = Q(0,0)*Q(1,1)-Q(0,1)*Q(1,0);
    }
  else if (n==3)
    {
      det = Q(0,0)*Q(1,1)*Q(2,2) + Q(0,1)*Q(1,2)*Q(2,0) + Q(0,2)*Q(1,0)*Q(2,1) 
	- Q(0,2)*Q(1,1)*Q(2,0) - Q(0,1)*Q(1,0)*Q(2,2) - Q(1,2)*Q(2,1)*Q(0,0);
    }

  return det;
}

complex<double> Matrix::trace()
{
  int n = this->getNDim();
  Matrix Q(n);
  Q = *this;
  complex<double> trace;
 
  if (n==2)
    {
      trace = Q(0,0)+Q(1,1);
    }
  else if (n==3)
    {
      trace = Q(0,0)+Q(1,1)+Q(2,2);
    }

  return trace;
}

double Matrix::FrobeniusNorm()
{
  int n = this->getNDim();
  Matrix Q(n);
  Q = *this;
  double norm;
 
  for(int i=0; i<n; i++)
    { 
      for(int j=0; j<n; j++)
	{
	  norm = abs(Q(i,j))*abs(Q(i,j));
	}
    }
  
  norm = sqrt(norm);

  return norm;
}

double Matrix::OneNorm()
{
  int n = this->getNDim();
  Matrix Q(n);
  Q = *this;
  double norm[3];
  double onenorm;

  for(int j=0; j<n; j++)
    { 
      for(int i=0; i<n; i++)
	{ 
	  norm[j] = abs(Q(i,j));
	}
    }
  
  onenorm = max(norm[0], norm[1]);
  onenorm = max(onenorm, norm[2]);
  
  return onenorm;
}


Matrix& Matrix::inv()
{
  int n = this->getNDim();
  Matrix H2(n);
  Matrix Q(n);
  Q = *this;
 
  if (n==2)
    {
      H2.set(0,0,Q(1,1));
      H2.set(0,1,-Q(0,1));
      H2.set(1,0,-Q(1,0));
      H2.set(1,1,Q(0,0));
      H2 *= 1./(Q(0,0)*Q(1,1)-Q(0,1)*Q(1,0)); // divide by det(H) 
    }
  else if (n==3)
    {
      H2.set(0,0,(Q(1,1)*Q(2,2)-Q(1,2)*Q(2,1)));
      H2.set(0,1,(Q(0,2)*Q(2,1)-Q(0,1)*Q(2,2)));
      H2.set(0,2,(Q(0,1)*Q(1,2)-Q(0,2)*Q(1,1)));
      H2.set(1,0,(Q(1,2)*Q(2,0)-Q(1,0)*Q(2,2)));
      H2.set(1,1,(Q(0,0)*Q(2,2)-Q(0,2)*Q(2,0)));
      H2.set(1,2,(Q(0,2)*Q(1,0)-Q(0,0)*Q(1,2)));
      H2.set(2,0,(Q(1,0)*Q(2,1)-Q(1,1)*Q(2,0)));
      H2.set(2,1,(Q(0,1)*Q(2,0)-Q(0,0)*Q(2,1)));
      H2.set(2,2,(Q(0,0)*Q(1,1)-Q(0,1)*Q(1,0)));
      H2 *= 1./(Q(0,0)*Q(1,1)*Q(2,2) + Q(0,1)*Q(1,2)*Q(2,0) + Q(0,2)*Q(1,0)*Q(2,1) 
		- Q(0,2)*Q(1,1)*Q(2,0) - Q(0,1)*Q(1,0)*Q(2,2) - Q(1,2)*Q(2,1)*Q(0,0)); // divide by det(H)
    }

  *this = H2;
  return *this;
}


// Pade approximant of log(I+A) (I is unit matrix). good for A\sim I
Matrix& Matrix::logm_pade(const int m)
{
   const int n = this->getNDim();
   Matrix S(n,0.);
   Matrix A(n);
   A = *this;
   Matrix I(n,1.);
   Matrix D(n); // denominator
   Matrix invD(n); // denominator
   double xi;
   double wi;
   gsl_integration_glfixed_table *table;
   table = gsl_integration_glfixed_table_alloc(m);

   for (int i=0; i<m; i++)
     {
       gsl_integration_glfixed_point(0.,1.,i,&xi,&wi,table);
       D = I + xi*A;
       // compute inverse of D:
       invD = D;
       invD.inv();
       S = S + wi * ( A*invD );
     }
   
   *this = S;
   
   return *this;
}

// Matrix square root by product from Denman-Beavers (DB) iteration.
// computes principal square root X of the matrix A using the product form
// of the Denman-Beavers iteration. The matrix M tends to I. 
// scale specifies scaling: 0, no scaling. 1, determinant scaling (default)
// maxit is the number of iterations.
// Adabted from The Matrix Function Toolbox by Nick Higham (MATLAB code)
Matrix& Matrix::sqrtm(const int scale)
{
   const int n = this->getNDim();
   int sc = scale;
   double eps = 1e-2;
   double tol = sqrt(static_cast<double>(n))*1e-16/2.;
   double g;
   double Mres;
   double reldiff;
   Matrix X(n);
   Matrix Xold(n);
   Matrix M(n);
   Matrix invM(n);
   Matrix I(n,1.);
   Matrix Mr(n);
   Matrix XmXo(n);
   
   X = *this;
   M = *this;
   
   int maxit = 25; // maximal number of iterations

   for (int k = 0; k < maxit; k++)
     {
       if(sc == 1)
	 {
	   g = pow(abs(M.det()),-1./(2.*n));
	   X = g*X;
	   M = g*g*M;
	 }
   
       Xold = X;
       invM = M;
       invM.inv();
       
       X = X*(I + invM)/2.;
       M = 0.5*(I+ (M+invM)/2.);

       Mr = M-I;
       Mres = Mr.FrobeniusNorm();

       XmXo = X - Xold;
       
       reldiff = XmXo.FrobeniusNorm()/X.FrobeniusNorm();
       if(reldiff<eps)
	 sc = 0; // switch to no scaling

       if(Mres <= tol)
	 break;
     }

   *this = X;
   
   return *this;
}


// matrix logarithm using Pade approximant (inverse scaling and squaring)
// A.H. Al-Mohy and N.J. Higham, Improved Inverse Scaling and Squaring algorithms for the matrix logarithm, MIMS eprint 2011.83
// this is not the improved one, just standard. 
Matrix& Matrix::logm()
{
   const int n = this->getNDim();
   const Matrix I(n,1.);
   Matrix X(n);
   Matrix L(n);
   int k, p, itk;
   double normdiff;
   int j1, j2;
   Matrix M(n);
   int m;
   
   double xvals[16] =
     {
       1.586970738772063e-005,
       2.313807884242979e-003,
       1.938179313533253e-002,
       6.209171588994762e-002,
       1.276404810806775e-001,
       2.060962623452836e-001,
       2.879093714241194e-001,
       3.666532675959788e-001,
       4.389227326152340e-001,
       5.034050432047666e-001,
       5.600071293013720e-001,
       6.092525642521717e-001,
       6.519202543720032e-001,
       6.888477797186464e-001,
       7.208340678820352e-001,
       7.485977242539218e-001
     };

   X = *this;
   k = 0;
   p = 0;
   itk = 5;
   
   while(1)
     {
       M = X-I;
       normdiff = M.OneNorm();
       
       if(normdiff <= xvals[15])
	 {
	   p = p+1;
	   // set j1
	   for (int i=0; i<16; i++)
	     {
	       if (normdiff <= xvals[i])
		 {
		   j1 = i;
		   break;
		 }
	     }

	   // set j2
	   for (int i=0; i<16; i++)
	     {
	       if (normdiff/2. <= xvals[i])
		 {
		   j2 = i;
		   break;
		 }
	     }

	   if((2*static_cast<double>(j1-j2)/3. < itk) || (p==2))
	     {
	       m = j1;
	       break; // break while loop
	     }
	 }

       X.sqrtm(); // take the square root
              
       k = k+1;
       
     }// while(1) loop
   

   L = X-I;
   L.logm_pade(m);
   
   X = pow(2.,k)*L;
   
   *this  = X;

   return *this;
}


// Output matrix elements in text, used when outputting the Wilson lines
// H.M. 20160801
string Matrix::getElementsText()
{
    stringstream ss;
    for (int i=0; i<nn; i++)
    {
        ss << e[i].real() << " " << e[i].imag() << " ";
    }
    return ss.str();
}


