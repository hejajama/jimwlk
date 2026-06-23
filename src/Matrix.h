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
    /// @brief Constructs a matrix with dimension n initialized to zero.
    /// @param n The dimension of the square matrix (n x n).
    Matrix(int n);
    
    /// @brief Constructs a matrix with dimension n where all elements are initialized to a.
    /// @param n The dimension of the square matrix (n x n).
    /// @param a The initial value for all matrix elements.
    Matrix(int n, double a);



    // Constructors that make expressions like A = B*C+D cheap 

    /// @brief Copy constructor (default).
    Matrix(const Matrix&) = default;
    
    /// @brief Move constructor.
    Matrix(Matrix&&) noexcept = default;
    
    /// @brief Move assignment operator.
    /// @param Other matrix to move from.
    /// @return Reference to this matrix.
    Matrix& operator = (Matrix&&) noexcept = default;
    
    /// @brief Copy assignment operator.
    /// @param p Matrix to copy from.
    /// @return Reference to this matrix.
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

    /// @brief Sets the real part of the i-th element (1D indexing).
    /// @param i Index in the flattened array.
    /// @param a New real value.
    void setRe(int i, double a) {e[i]=complex<double>(a,e[i].imag());}
    /// @brief Sets the real part of element (i,j) (2D indexing).
    /// @param i Row index.
    /// @param j Column index.
    /// @param a New real value.
    void setRe(int i, int j, double a) {e[j+ndim*i]=complex<double>(a,e[j+ndim*i].imag());}
    /// @brief Sets the imaginary part of the i-th element (1D indexing).
    /// @param i Index in the flattened array.
    /// @param a New imaginary value.
    void setIm(int i, double a) {e[i]=complex<double>(e[i].real(),a);}    
    /// @brief Sets the imaginary part of element (i,j) (2D indexing).
    /// @param i Row index.
    /// @param j Column index.
    /// @param a New imaginary value.
    void setIm(int i, int j, double a) {e[j+ndim*i]=complex<double>(e[j+ndim*i].real(),a);}
  
    /// @brief Sets the element at index i to value a (1D indexing).
    /// @param i Index in the flattened array.
    /// @param a New complex value.
    void set(int i, complex<double> a) {e[i]=a;}
    
    /// @brief Sets the element at position (i,j) to value a (2D indexing).
    /// @param i Row index.
    /// @param j Column index.
    /// @param a New complex value.
    void set(int i, int j, complex<double> a) {e[j+ndim*i]=a;}
    
    /// @brief Gets the element at index i (1D indexing).
    /// @param i Index in the flattened array.
    /// @return The complex value at index i.
    complex<double> get(int i) const {return e[i];}
    
    /// @brief Gets the element at position (i,j) (2D indexing).
    /// @param i Row index.
    /// @param j Column index.
    /// @return The complex value at position (i,j).
    complex<double> get(int i, int j) const {return e[j+ndim*i];}

    /// @brief Gets the real part of the element at index i.
    /// @param i Index in the flattened array.
    /// @return The real part of element at index i.
    double getRe(int i) const {return e[i].real();}
    
    /// @brief Gets the imaginary part of the element at index i.
    /// @param i Index in the flattened array.
    /// @return The imaginary part of element at index i.
    double getIm(int i) const {return e[i].imag();}    
    
    /// @brief Gets the dimension of the square matrix.
    /// @return The matrix dimension (n for an n x n matrix).
    int getNDim() const {return ndim;}
    
    /// @brief Gets the total number of elements in the matrix.
    /// @return The number of elements (ndim * ndim).
    int getNN()  const {return nn;}

    // OPT: direct access to storage for tight loops (FFT gather/scatter etc.)
    /// @brief Gets a pointer to the underlying data array.
    /// @return Pointer to the first element of the matrix data.
    complex<double>* data() { return e.data(); }
    
    /// @brief Gets a const pointer to the underlying data array.
    /// @return Const pointer to the first element of the matrix data.
    const complex<double>* data() const { return e.data(); }

    // OPT: set all elements to zero without allocating
    /// @brief Sets all matrix elements to zero without reallocating memory.
    void setZero() { std::fill(e.begin(), e.end(), complex<double>(0.,0.)); }

    /// @brief In-place addition of scaled matrix: this += a*B (complex scalar). Avoids temporary matrix allocation.
    /// @param a Complex scaling factor.
    /// @param B Matrix to add scaled.
    void addMultiple(const complex<double> a, const Matrix& B)
    {
        for (int i=0; i<nn; i++) e[i] += a*B.e[i];
    }
    
    /// @brief In-place addition of scaled matrix: this += a*B (real scalar).
    /// @param a Real scaling factor.
    /// @param B Matrix to add scaled.
    void addMultiple(const double a, const Matrix& B)
    {
        for (int i=0; i<nn; i++) e[i] += a*B.e[i];
    }

    // OPT: c = a*b written into existing storage (no allocation).
    // c must not alias a or b.
    /// @brief Static method for matrix multiplication c = a*b into pre-allocated storage.
    /// @param a First matrix multiplicand.
    /// @param b Second matrix multiplicand.
    /// @param c Result matrix (must not alias a or b).
    /// @note More efficient than c = a*b as it reuses existing storage.
    static void mult(const Matrix& a, const Matrix& b, Matrix& c);

    /// @brief Computes the matrix exponential exp(t*this).
    /// @param t Time/scaling parameter (default 1.0).
    /// @param p Order of the approximation (default 6).
    /// @return Reference to this matrix containing the result.
    Matrix& expm(double t = 1.0, const int p = 6);
    
    /// @brief Computes the matrix logarithm log(this).
    /// @return Reference to this matrix containing the result.
    Matrix& logm();
    
    /// @brief Computes the matrix logarithm using Padé approximation.
    /// @param m Order of the Padé approximation.
    /// @return Reference to this matrix containing the result.
    Matrix& logm_pade(const int m);
    
    /// @brief Computes the matrix inverse.
    /// @return Reference to this matrix containing the result.
    Matrix& inv();
    
    /// @brief Normalizes the matrix using a specific norm.
    /// @param m Norm type to use for normalization.
    /// @return Reference to this matrix after normalization.
    Matrix& normAm(const int m);
    
    /// @brief Computes the matrix square root.
    /// @param scale Scale parameter (default 1).
    /// @return Reference to this matrix containing the result.
    Matrix& sqrtm(const int scale = 1);
    
    /// @brief Computes the determinant of the matrix.
    /// @return The determinant as a complex number.
    complex<double> det() const;
    
    /// @brief Computes the trace (sum of diagonal elements).
    /// @return The trace as a complex number.
    complex<double> trace() const;
    
    /// @brief Computes the Frobenius norm of the matrix.
    /// @return The Frobenius norm (sqrt of sum of squared absolute values).
    double FrobeniusNorm() const;
    
    /// @brief Computes the one norm (maximum absolute column sum).
    /// @return The one norm.
    double OneNorm() const;
    
    /// @brief Returns a string representation of all matrix elements.
    /// @return String containing the matrix elements.
    string getElementsText() const;

    //operators:

    /// @brief Parenthesis operator for 1D access to elements.
    /// @param i Index in the flattened array.
    /// @return The element at index i.
    std::complex<double> operator () (const int i) const { return e[i]; }
    
    /// @brief Parenthesis operator for 2D access to elements.
    /// @param i Row index.
    /// @param j Column index.
    /// @return The element at position (i,j).
    std::complex<double> operator () (const int i, const int j) const { return e[j+ndim*i]; }

    /// @brief Equality comparison operator.
    /// @param p Matrix to compare with.
    /// @return True if all elements are equal, false otherwise.
    bool operator == (const Matrix& p) const 
    {
      if(ndim != p.ndim || nn != p.nn) return false;
      for(int i=0; i<nn; i++)
	if(e[i] != p.e[i]) return false;
      return true;
    }

    /// @brief Inequality comparison operator.
    /// @param p Matrix to compare with.
    /// @return True if any element differs, false otherwise.
    bool operator != (const Matrix& p) const 
    {
      return !(*this == p);
    }

    /// @brief In-place addition assignment operator.
    /// @param a Matrix to add to this matrix.
    /// @return Reference to this matrix.
    Matrix& operator += (const Matrix& a) 
      {
	for(int i=0; i<nn; i++) e[i] += a.e[i];
	return *this;
      }

    /// @brief In-place subtraction assignment operator.
    /// @param a Matrix to subtract from this matrix.
    /// @return Reference to this matrix.
    Matrix& operator -= (const Matrix& a) 
      {
	for(int i=0; i<nn; i++) e[i] -= a.e[i];
	return *this;
      }

    /// @brief In-place multiplication assignment by complex scalar.
    /// @param a Complex scalar to multiply this matrix by.
    /// @return Reference to this matrix.
    Matrix& operator *= (const complex<double> a) 
      {
	for(int i=0; i<nn; i++) e[i] *= a;
	return *this;
      }

    /// @brief In-place division assignment by complex scalar.
    /// @param a Complex scalar to divide this matrix by.
    /// @return Reference to this matrix.
    Matrix& operator /= (const complex<double> a) 
      {
	for(int i=0; i<nn; i++) e[i] /= a;
	return *this;
      }

    // /= by a real scalar (component-wise division, bitwise identical to
    // the original operator/(Matrix,double))
    /// @brief In-place division assignment by real scalar.
    /// @param a Real scalar to divide this matrix by.
    /// @return Reference to this matrix.
    Matrix& operator /= (const double a)
      {
	for(int i=0; i<nn; i++) e[i] /= a;
	return *this;
      }
    
    /// @brief Computes the squared Frobenius norm.
    /// @return Half the sum of squared absolute values of all elements.
    double square() const 
    {
      double tr = 0.0;
      for(int i=0; i<nn;i++) 
	{
	  tr += e[i].real()*e[i].real()+e[i].imag()*e[i].imag();
	}
      return 0.5*tr;
    }

    /// @brief Extracts the imaginary part of all elements.
    /// @return Reference to this matrix with imaginary parts as real values.
    Matrix& imag();

    /// @brief Computes the complex conjugate of all elements.
    /// @return Reference to this matrix containing the conjugate values.
    Matrix& conjg();

    
    /// @brief Stream insertion operator for output.
    /// @param os Output stream.
    /// @param p Matrix to output.
    /// @return Reference to the output stream.
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

/// @brief Addition operator for two matrices.
/// @param a First matrix.
/// @param b Second matrix.
/// @return Result of a + b.
Matrix  operator + (const Matrix& a, const Matrix& b);

/// @brief Subtraction operator for two matrices.
/// @param a First matrix.
/// @param b Second matrix to subtract.
/// @return Result of a - b.
Matrix  operator - (const Matrix& a, const Matrix& b);

/// @brief Unary negation operator.
/// @param a Matrix to negate.
/// @return Negative of the matrix.
Matrix  operator - (const Matrix& a);

/// @brief Division operator (matrix division, equivalent to a * inv(b)).
/// @param a Numerator matrix.
/// @param b Denominator matrix.
/// @return Result of a / b.
Matrix  operator / (const Matrix&  a, const Matrix& b);

/// @brief Multiplication operator: real scalar times matrix.
/// @param a Real scalar.
/// @param b Matrix to multiply.
/// @return Result of a * b.
Matrix operator * (const double a, const Matrix& b);

/// @brief Multiplication operator: complex scalar times matrix.
/// @param a Complex scalar.
/// @param b Matrix to multiply.
/// @return Result of a * b.
Matrix operator * (const std::complex<double> a,const Matrix& b);

/// @brief Multiplication operator: matrix times real scalar.
/// @param a Matrix to multiply.
/// @param b Real scalar.
/// @return Result of a * b.
Matrix operator * (const Matrix&  a, const double b);

/// @brief Multiplication operator: matrix times matrix.
/// @param a First matrix.
/// @param b Second matrix.
/// @return Result of a * b.
Matrix operator * (const Matrix&  a, const Matrix& b);

/// @brief Division operator: matrix divided by real scalar.
/// @param a Matrix to divide.
/// @param b Real scalar divisor.
/// @return Result of a / b.
Matrix operator / (const Matrix&  a, const double b);

// OPT: rvalue overloads reuse the temporary's storage instead of allocating
/// @brief Rvalue addition (reuses temporary storage).
/// @param a Rvalue reference to first matrix.
/// @param b Second matrix (const reference).
/// @return Sum reusing a's storage.
inline Matrix operator + (Matrix&& a, const Matrix& b) { a += b; return std::move(a); }

/// @brief Rvalue addition (reuses temporary storage).
/// @param a First matrix (const reference).
/// @param b Rvalue reference to second matrix.
/// @return Sum reusing b's storage.
inline Matrix operator + (const Matrix& a, Matrix&& b) { b += a; return std::move(b); }

/// @brief Rvalue addition (reuses temporary storage).
/// @param a First rvalue matrix.
/// @param b Second rvalue matrix.
/// @return Sum reusing a's storage.
inline Matrix operator + (Matrix&& a, Matrix&& b)      { a += b; return std::move(a); }

/// @brief Rvalue subtraction (reuses temporary storage).
/// @param a Rvalue reference to first matrix.
/// @param b Second matrix (const reference).
/// @return Difference reusing a's storage.
inline Matrix operator - (Matrix&& a, const Matrix& b) { a -= b; return std::move(a); }

/// @brief Rvalue multiplication by real scalar (reuses temporary storage).
/// @param s Real scalar.
/// @param a Rvalue reference to matrix.
/// @return Product reusing a's storage.
inline Matrix operator * (const double s, Matrix&& a)  { a *= s; return std::move(a); }

/// @brief Rvalue multiplication by real scalar (reuses temporary storage).
/// @param a Rvalue reference to matrix.
/// @param s Real scalar.
/// @return Product reusing a's storage.
inline Matrix operator * (Matrix&& a, const double s)  { a *= s; return std::move(a); }

/// @brief Rvalue multiplication by complex scalar (reuses temporary storage).
/// @param s Complex scalar.
/// @param a Rvalue reference to matrix.
/// @return Product reusing a's storage.
inline Matrix operator * (const std::complex<double> s, Matrix&& a) { a *= s; return std::move(a); }

/// @brief Rvalue division by real scalar (reuses temporary storage).
/// @param a Rvalue reference to matrix.
/// @param s Real scalar divisor.
/// @return Quotient reusing a's storage.
inline Matrix operator / (Matrix&& a, const double s)  { a /= s; return std::move(a); }

#endif

