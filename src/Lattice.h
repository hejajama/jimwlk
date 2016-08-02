#ifndef Lattice_h
#define Lattice_h

#include <complex>
#include <iostream>
#include <cstdlib>
#include <string>

#include "Matrix.h"
#include "Cell.h"
#include "Parameters.h"

// The Lattice class is a level higher than the Cell class
// It takes care of the overall structure of the lattice and the arrangement of individual cells
// "cells" is an array of pointers to individual cells of the lattice.
// The values of the quantities in a cell can be modified or retrieved by the public functions
// in both lattice and cells. The user has acces to the cell values via the lattice class only,
// because cells is a private object in Lattice.
// cells are initialized in the Lattice constructor.

using namespace std;

class Lattice
{
private:

  int size;             // the total number of cells (length*length)
  int Nc;               // the number of colors in SU(Nc): Determines the dimension of the used matrices
 
public:
  //constructor
  Lattice(Parameters *param,int N, int length);
  //destructor
  ~Lattice();

  //functions to access values within individual cells
  double getX(int pos) { return cells[pos]->getX(); };
  double getY(int pos) { return cells[pos]->getY(); };
  int getSize() { return size; };

  Cell** cells;         // the actual array of cells, the "lattice". cells is an array of pointers to cell objects
    
  void PrintWilsonLines(string filename);

};

#endif

