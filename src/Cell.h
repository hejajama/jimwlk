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
  // if cell is on the heap, all of the following will be too.
  Matrix* U;     // U is in the fundamental rep. (Nc*Nc matrix)
  Matrix* Ui;    // Ui is the initial U in the fundamental rep. (Nc*Nc matrix)
  Matrix** Uy;    // Uy is an intermediate U in the fundamental rep. (Nc*Nc matrix) (for unequal y correlations)
  Matrix* UA;    // UA is in the adjoint rep.
    
public:
  Cell(int N);
  ~Cell();

  void setX(double in) { xpos = in; };
  void setY(double in) { ypos = in; };

  double getX() { return xpos; };
  double getY() { return ypos; };

  void setU(const Matrix& x) { *U = x; };
  void setUi(const Matrix& x) { *Ui = x; };
  void setUy(const int i, const Matrix& x) { *Uy[i] = x; };
  void setUA(const Matrix& x) { *UA = x; };

  Matrix& getU() const { return *U; };
  Matrix& getUy(int i) const { return *Uy[i]; };
  Matrix& getUi() const { return *Ui; };
  Matrix& getUA() const { return *UA; };

  void computeAdjointU();

};

#endif

