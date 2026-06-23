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
      }
}

Lattice::~Lattice()
{
  for(int i=0; i<size; i++)
    delete cells[i];
  delete[] cells;
}

void Lattice::PrintWilsonLines(string filename, Parameters *param)
{
    cout << "Saving Wilson lines into " << filename << endl;

    if (param->getInitMethod() == 11) //  binary
    {
        int N = param->getSize();
        int Nc = param->getNc();
        double L = param->getL();
        double a = L/N;
        double tmp = 0; // not implemented
        std::ofstream Outfile;
        Outfile.open(filename.c_str(), ios::out | ios::binary);
        Outfile.write((char *) &N ,sizeof(int));
        Outfile.write((char *) &Nc ,sizeof(int));
        Outfile.write((char *) &L ,sizeof(double));
        Outfile.write((char *) &a ,sizeof(double));
        Outfile.write((char *) &tmp ,sizeof(double));
        double val1[2];   // BUGFIX: was heap allocated and leaked
        for(int ix=0; ix<N; ix++)
        {
            for(int iy=0; iy<N; iy++)
            {
                // BUGFIX: the matrix loop bounds were hardcoded to 3
                // (wrong for Nc=2)
                for(int j=0; j<Nc; j++)
                {
                    for(int k=0; k<Nc; k++)
                    {
                        int indx = N*iy+ix;
                        val1[0] = (cells[indx]->getU()).getRe(j*Nc+k);
                        val1[1] = (cells[indx]->getU()).getIm(j*Nc+k);

                        Outfile.write((char *) val1 ,2*sizeof(double));
                    }
                }
            }
        }  
        
        if (Outfile.good()==false)
        {
            cerr << "Error when saving in file " << filename << endl;
            exit(1);
        } 
        Outfile.close();
    }
    else
    {
        // Text output.
        // BUGFIX: the original only wrote text output for initMethod==10,
        // so runs initialized with methods 1-4 (Gauss / MV / ...) silently
        // produced no Wilson line output at all. Text is now the default
        // for everything except the binary format (11).
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
                output << yind << " " << xind << " " << cells[pos]->getU().getElementsText() << "\n";
            }
        }
        output.close();
    }
}
