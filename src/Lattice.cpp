#include "Lattice.h"
#include <iostream>
#include <fstream>

//constructor
Lattice::Lattice(Parameters *param, int N, int length)
{
  Nc = N;
  size = length*length;

  cout << "Allocating square lattice of size " << length << "x" << length << " ...";

  // initialize the array of cells
  cells = new Cell*[size];
  for(int i=0; i<size; i++)
    {
      cells[i] = new Cell(Nc);
    }

  // x (i) is the outer loop, y (j) the inner
  for (int i=0; i<length; i++)
    for (int j=0; j<length; j++)
      {
	int pos = i*length+j;
	cells[pos]->setX((i-length/2.));
	cells[pos]->setY((j-length/2.));
	//cells[pos]->setX(i);
	//cells[pos]->setY(j);
      }
}

Lattice::~Lattice()
{
  for(int i=0; i<size; i++)
    delete cells[i];
  delete[] cells;
}

void Lattice::PrintWilsonLines(string filename)
{
    cout << "Saving Wilson lines into " << filename << endl;
    fstream output(filename.c_str(),ios::out);
    if (!output.good())
    {
	cerr << "Can't open file for saving, directory does not exist???" << endl;
	exit(1);
    }
    int length = sqrt(size);
    for (int yind=0; yind<length; yind++)
    {
        for (int xind=0; xind<length; xind++)
        {
            int pos = yind*length+xind;
            output << yind << " " << xind << " " << cells[pos]->getU().getElementsText() << endl;
        }
    }
    output.close();
}
