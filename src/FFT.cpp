// FFT.cpp is part of the JIMWLK solver.
// Copyright (C) 2011 Bjoern Schenke.
// This new version uses FFTW
#include "FFT.h"

//**************************************************************************
// FFT class.
//
// OPT: the original code copied the data into the FFTW buffer (and back)
// with four separate loop nests per quadrant; these are equivalent to a
// single FFT-shift index map, newpos = ((i+n0/2)%n0)*n1 + ((j+n1/2)%n1),
// which we precompute once per call and reuse for all mDim components.
// The data movement (and hence the numerical result) is identical to the
// original. The unused (and broken) fftnComplex routine was removed.
//**************************************************************************

void FFT::fftnVector(vector<complex<double> > **data, vector<complex<double> > **outdata, const int nn[], const int ndim, const int isign)
{
    (void) ndim; // always 2 here
    const int N0 = nn[0], N1 = nn[1];
    const int ntot = N0*N1;
    const int mDim = data[0]->size();
    // mDim is the size of the vector (how many rows)

    // precompute the shift map once
    vector<int> shift(ntot);
    for (int i=0; i<N0; i++)
      {
        const int is = ((i+N0/2)%N0)*N1;
        for (int j=0; j<N1; j++)
          shift[i*N1+j] = is + ((j+N1/2)%N1);
      }

    for (int k=0; k<mDim; k++)
      {
        for (int pos=0; pos<ntot; pos++)
          {
            const complex<double> v = (*data[pos])[k];
            const int newpos = shift[pos];
            input[newpos][0] = v.real();
            input[newpos][1] = v.imag();
          }

        execute(isign);

        for (int pos=0; pos<ntot; pos++)
          {
            const int newpos = shift[pos];
            (*outdata[pos])[k] = complex<double>(output[newpos][0],output[newpos][1]);
          }
      }
}

void FFT::fftnArray(complex<double> **data, complex<double>  **outdata, const int nn[], const int ndim, const int isign, const int mDim)
{
    (void) ndim; // always 2 here
    const int N0 = nn[0], N1 = nn[1];
    const int ntot = N0*N1;

    vector<int> shift(ntot);
    for (int i=0; i<N0; i++)
      {
        const int is = ((i+N0/2)%N0)*N1;
        for (int j=0; j<N1; j++)
          shift[i*N1+j] = is + ((j+N1/2)%N1);
      }

    for (int k=0; k<mDim; k++)
      {
        for (int pos=0; pos<ntot; pos++)
          {
            const complex<double> v = data[pos][k];
            const int newpos = shift[pos];
            input[newpos][0] = v.real();
            input[newpos][1] = v.imag();
          }

        execute(isign);

        for (int pos=0; pos<ntot; pos++)
          {
            const int newpos = shift[pos];
            outdata[pos][k] = complex<double>(output[newpos][0],output[newpos][1]);
          }
      }
}


// Performs Fast Fourier Transform of any object of class "T" (matrix or something else) using a wrapper for FFTW
// This routine takes data as a function of -x_max/2 to x_max/2 and returns it ordered similarly - no need to resort before or after!
template <class T>
void FFT::fftn(T **data, T **outdata, const int nn[], const int ndim, const int isign)
{
    (void) ndim; // always 2 here
    const int N0 = nn[0], N1 = nn[1];
    const int ntot = N0*N1;
    int mDim = data[0]->getNDim();
    mDim*=mDim;

    vector<int> shift(ntot);
    for (int i=0; i<N0; i++)
      {
        const int is = ((i+N0/2)%N0)*N1;
        for (int j=0; j<N1; j++)
          shift[i*N1+j] = is + ((j+N1/2)%N1);
      }

    for (int k=0; k<mDim; k++)
      {
        for (int pos=0; pos<ntot; pos++)
          {
            const int newpos = shift[pos];
            input[newpos][0] = data[pos]->getRe(k);
            input[newpos][1] = data[pos]->getIm(k);
          }

        execute(isign);

        for (int pos=0; pos<ntot; pos++)
          {
            const int newpos = shift[pos];
            outdata[pos]->set(k, complex<double>(output[newpos][0],output[newpos][1]));
          }
      }
}


// Define specializations of the template:
template void FFT::fftn(Matrix **data,Matrix **outdata, const int nn[], const int ndim, const int isign);
