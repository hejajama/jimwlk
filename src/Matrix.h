#ifndef Matrix_h
#define Matrix_h

#include <complex>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <gsl/gsl_integration.h>  // include gsl for Gauss-Legendre nodes and weights for log Pade
#include <string>

using namespace std;

class Matrix
{
private:
    int ndim;
    int nn;
    //complex<double>* e;
    vector<complex<double> > e;

public:
    
    //constructor(s)
    Matrix(int n);
    Matrix(int n, double a);

    //destructor
    ~Matrix()
      {
          /*
        if (e != NULL)
        {
          delete[] e;
          e = NULL;
        }
        else
            cerr << "Warning, tried to delete memory twice!" << endl;
           */
      }

    void setRe(int i, double a) {e[i]=complex<double>(a,e[i].imag());};
    void setRe(int i, int j, double a) {e[j+ndim*i]=complex<double>(a,e[j+ndim*i].imag());};
    void setIm(int i, double a) {e[i]=complex<double>(e[i].real(),a);};    
    void setIm(int i, int j, double a) {e[j+ndim*i]=complex<double>(e[j+ndim*i].real(),a);};
  
    void set(int i, complex<double> a) {e[i]=a;};
    void set(int i, int j, complex<double> a) {e[j+ndim*i]=a;};
    
    complex<double> get(int i) {return e[i];};
    complex<double> get(int i, int j) {return e[j+ndim*i];};

    double getRe(int i) {return e[i].real();};
    double getIm(int i) {return e[i].imag();};    
    
    int getNDim() const {return ndim;}
    int getNN()  const {return nn;}

    Matrix& expm(double t = 1.0, const int p = 6);
    Matrix& logm();
    Matrix& logm_pade(const int m);
    Matrix& inv();
    Matrix& normAm(const int m);
    Matrix& sqrtm(const int scale = 1);
    complex<double> det();
    complex<double> trace();
    double FrobeniusNorm();
    double OneNorm();
    
    string getElementsText();   // Return string containing elements

    //    complex<double> trace() {if (ndim==3) return e[0]+e[4]+e[8]; else if(ndim==2) return e[0]+e[3]; };
    //operators:

    //()
    std::complex<double> operator () (const int i) const { return e[i]; }
    std::complex<double> operator () (const int i, const int j) const { return e[j+ndim*i]; }
    //std::complex<double> operator () (const int i, const int j) { return e[j+ndim*i]; }

    //=
    const Matrix& operator = (const Matrix& p) 
      {
	nn = p.getNN();
	if(&p != this ) 
	  {
	    for(int i=0; i<nn; i++) 
	      {
		e[i] = p.e[i];
	      }
	  }
	return *this;
      }

    //==
    bool operator == (const Matrix& p) const 
    {
      for(int i=0; i<nn; i++)
	if(e[i] != p.e[i]) return false;
      return true;
    }

    //!=
    bool operator != (const Matrix& p) const 
    {
      for(int i=0; i<nn; i++)
	if(e[i] != p.e[i]) return true;
      return false;
    }

    //+=
    Matrix& operator += (const Matrix& a) 
      {
	for(int i=0; i<nn; i++) e[i] += a.e[i];
	return *this;
      }

    //-=
    Matrix& operator -= (const Matrix& a) 
      {
	for(int i=0; i<nn; i++) e[i] -= a.e[i];
	return *this;
      }

    //*=
    Matrix& operator *= (const complex<double> a) 
      {
	for(int i=0; i<nn; i++) e[i] *= a;
	return *this;
      }

    // /=
    Matrix& operator /= (const complex<double> a) 
      {
	for(int i=0; i<nn; i++) e[i] /= a;
	return *this;
      }
    
    double square() const 
    {
      double tr = 0.0;
      for(int i=0; i<nn;i++) 
	{
	  tr += e[i].real()*e[i].real()+e[i].imag()*e[i].imag();
	}
      return 0.5*tr;
    }

    Matrix& imag();

    Matrix& conjg();

    
    //<<
    friend ostream& operator<<(ostream& os, const Matrix& p) 
    {
      for(int i=0; i<p.getNDim();i++) 
	{
	  for(int j=0; j<p.getNDim();j++) os << p(i,j);
	  if(i<p.getNDim()-1) os << std::endl;
	}
      return os;
    }


};

Matrix  operator + (const Matrix& a, const Matrix& b);
Matrix  operator - (const Matrix& a, const Matrix& b);
Matrix  operator - (const Matrix& a);
Matrix  operator / (const Matrix&  a, const Matrix& b);

Matrix operator * (const double a, const Matrix& b);
Matrix operator * (const std::complex<double> a,const Matrix& b);
Matrix operator * (const Matrix&  a, const double b);
Matrix operator * (const Matrix&  a, const Matrix& b);


#endif

