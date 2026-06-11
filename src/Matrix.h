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
    vector<complex<double> > e;

public:
    
    //constructor(s)
    Matrix(int n);
    Matrix(int n, double a);

    // OPT: explicit copy/move semantics. The user-declared copy assignment in the
    // original code suppressed the implicit move constructor/assignment, so every
    // temporary produced by the overloaded operators was deep-copied. With move
    // semantics, expression chains like A = B*C + D become cheap.
    Matrix(const Matrix&) = default;
    Matrix(Matrix&&) noexcept = default;
    Matrix& operator = (Matrix&&) noexcept = default;
    Matrix& operator = (const Matrix& p)
    {
        if (&p != this)
        {
            ndim = p.ndim;
            nn = p.nn;
            e = p.e;   // single vector assignment (memcpy-like), no per-element loop
        }
        return *this;
    }

    void setRe(int i, double a) {e[i]=complex<double>(a,e[i].imag());};
    void setRe(int i, int j, double a) {e[j+ndim*i]=complex<double>(a,e[j+ndim*i].imag());};
    void setIm(int i, double a) {e[i]=complex<double>(e[i].real(),a);};    
    void setIm(int i, int j, double a) {e[j+ndim*i]=complex<double>(e[j+ndim*i].real(),a);};
  
    void set(int i, complex<double> a) {e[i]=a;};
    void set(int i, int j, complex<double> a) {e[j+ndim*i]=a;};
    
    complex<double> get(int i) const {return e[i];};
    complex<double> get(int i, int j) const {return e[j+ndim*i];};

    double getRe(int i) const {return e[i].real();};
    double getIm(int i) const {return e[i].imag();};    
    
    int getNDim() const {return ndim;}
    int getNN()  const {return nn;}

    // OPT: direct access to storage for tight loops (FFT gather/scatter etc.)
    complex<double>* data() { return e.data(); }
    const complex<double>* data() const { return e.data(); }

    // OPT: set all elements to zero without allocating
    void setZero() { std::fill(e.begin(), e.end(), complex<double>(0.,0.)); }

    // OPT: *this += a*B  (axpy), avoids two temporaries of "*this = *this + a*B"
    void addMultiple(const complex<double> a, const Matrix& B)
    {
        for (int i=0; i<nn; i++) e[i] += a*B.e[i];
    }
    void addMultiple(const double a, const Matrix& B)
    {
        for (int i=0; i<nn; i++) e[i] += a*B.e[i];
    }

    // OPT: c = a*b written into existing storage (no allocation).
    // c must not alias a or b.
    static void mult(const Matrix& a, const Matrix& b, Matrix& c);

    Matrix& expm(double t = 1.0, const int p = 6);
    Matrix& logm();
    Matrix& logm_pade(const int m);
    Matrix& inv();
    Matrix& normAm(const int m);
    Matrix& sqrtm(const int scale = 1);
    complex<double> det() const;
    complex<double> trace() const;
    double FrobeniusNorm() const;
    double OneNorm() const;
    
    string getElementsText() const;   // Return string containing elements

    //operators:

    //()
    std::complex<double> operator () (const int i) const { return e[i]; }
    std::complex<double> operator () (const int i, const int j) const { return e[j+ndim*i]; }

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
      return !(*this == p);
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

    // /= by a real scalar (component-wise division, bitwise identical to
    // the original operator/(Matrix,double))
    Matrix& operator /= (const double a)
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
Matrix operator / (const Matrix&  a, const double b);

// OPT: rvalue overloads reuse the temporary's storage instead of allocating
inline Matrix operator + (Matrix&& a, const Matrix& b) { a += b; return std::move(a); }
inline Matrix operator + (const Matrix& a, Matrix&& b) { b += a; return std::move(b); }
inline Matrix operator + (Matrix&& a, Matrix&& b)      { a += b; return std::move(a); }
inline Matrix operator - (Matrix&& a, const Matrix& b) { a -= b; return std::move(a); }
inline Matrix operator * (const double s, Matrix&& a)  { a *= s; return std::move(a); }
inline Matrix operator * (Matrix&& a, const double s)  { a *= s; return std::move(a); }
inline Matrix operator * (const std::complex<double> s, Matrix&& a) { a *= s; return std::move(a); }
inline Matrix operator / (Matrix&& a, const double s)  { a /= s; return std::move(a); }

#endif

