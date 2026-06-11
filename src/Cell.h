#ifndef Cell_h
#define Cell_h

#include <complex>
#include <iostream>
#include <cstdlib>
#include <vector>

#include "Matrix.h"

using namespace std;

class Cell
{
private:

  int Nc;
  double xpos;  // the cells x position with respect to the lattice
  double ypos;  // the cells y position with respect to the lattice

  // OPT: matrices are stored by value (no extra pointer indirection, fewer
  // heap allocations, better locality). The Uy array (10 matrices, used only
  // for unequal-rapidity correlations) is allocated lazily on first use;
  // in the original it cost 10 heap-allocated matrices per cell whether or
  // not the feature was used.
  Matrix U;     // U is in the fundamental rep. (Nc*Nc matrix)
  Matrix Ui;    // Ui is the initial U in the fundamental rep. (Nc*Nc matrix)
  Matrix UA;    // UA is in the adjoint rep.
  vector<Matrix> Uy; // intermediate U's (for unequal y correlations), lazy

  void ensureUy()
  {
    if (Uy.empty())
      Uy.assign(10, Matrix(Nc,1.));
  }

public:
  Cell(int N);
  ~Cell();

  void setX(double in) { xpos = in; };
  void setY(double in) { ypos = in; };

  double getX() const { return xpos; };
  double getY() const { return ypos; };

  void setU(const Matrix& x) { U = x; };
  void setUi(const Matrix& x) { Ui = x; };
  void setUy(const int i, const Matrix& x) { ensureUy(); Uy[i] = x; };
  void setUA(const Matrix& x) { UA = x; };

  Matrix& getU() { return U; };
  Matrix& getUy(int i) { ensureUy(); return Uy[i]; };
  Matrix& getUi() { return Ui; };
  Matrix& getUA() { return UA; };

  void computeAdjointU();

};

#endif

