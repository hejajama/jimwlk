// FFT.h is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.

#ifndef FFT_H
#define FFT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <cmath>
#include <complex>
#include <math.h>
#include "Matrix.h"
#include <vector>
#include <algorithm>
#include <functional>

#include <fftw3.h>

using namespace std;


template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::plus<T>());
    return result;
}

template <typename T>
std::vector<T> operator*(const std::vector<T>& a, const std::complex<double>& b)
{
  std::vector<T> result;
  result = a;
  int size = a.size();
  for (int i = 0; i<size; i++)
    {
      result[i]=b*a[i];
    }
  return result;
}

template <typename T>
std::vector<T> operator / (const std::vector<T>& a, const double b)
{
  std::vector<T> result;
  result = a;
  int size = a.size();
  for (int i = 0; i<size; i++)
    {
      result[i]=a[i]/b;
    }
  return result;
}


class FFT 
{
 private:
  fftw_complex *input, *output;
  fftw_plan p, pback;
  int n0, n1;

public:

  // Constructor.
  FFT(const int nn[]) 
    {
      n0 = nn[0];
      n1 = nn[1];
      input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nn[0] * nn[1]);
      output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nn[0] * nn[1]);
#ifdef JIMWLK_DETERMINISTIC_FFT
      // FFTW_MEASURE selects a plan by run-time benchmarking, which is not
      // deterministic between runs (and the resulting plans can differ in
      // floating point round-off). For validation / strict reproducibility
      // build with -DJIMWLK_DETERMINISTIC_FFT.
      p = fftw_plan_dft_2d(nn[0], nn[1], input, output, FFTW_FORWARD, FFTW_ESTIMATE);
      pback = fftw_plan_dft_2d(nn[0], nn[1], input, output, FFTW_BACKWARD, FFTW_ESTIMATE);
#else
      p = fftw_plan_dft_2d(nn[0], nn[1], input, output, FFTW_FORWARD, FFTW_MEASURE);
      pback = fftw_plan_dft_2d(nn[0], nn[1], input, output, FFTW_BACKWARD, FFTW_MEASURE);
#endif
    };
  // Destructor.
  ~FFT() 
    {
      fftw_destroy_plan(p);
      fftw_destroy_plan(pback); // BUGFIX: pback was leaked in the original
      fftw_free(input); fftw_free(output);
    };
  void fftnVector(vector<complex<double> > **data, vector<complex<double> > **outdata, const int nn[], const int ndim, const int isign);
  void fftnArray(complex<double> **data, complex<double> **outdata, const int nn[], const int ndim, const int isign, const int mDim);
  template<class T>
  void fftn(T **data, T **outdata, const int nn[], const int ndim, const int isign);

private:
  // OPT: shared driver. The four quadrant-copy loop nests of the original
  // are a single FFT-shift: newpos = ((i+n0/2)%n0)*n1 + ((j+n1/2)%n1).
  // The mapping (and therefore the result) is identical.
  void execute(const int isign)
  {
      if (isign==1)
        fftw_execute(p);
      else
        fftw_execute(pback);

      // if this is inverse transform, normalize.
      if(isign == -1)
        {
          const double ntot = static_cast<double>(n0)*static_cast<double>(n1);
          const int N = n0*n1;
          for(int i=0; i<N; i++)
            {
              output[i][0]/=ntot;
              output[i][1]/=ntot;
            }
        }
  }
};

#endif // FFT_H
