# JIMWLK

The code solves JIMWLK evolution equation using Langevin step method

Originally written by Björn Schenke, some modificatoins by Heikki Mäntysaari <heikki.mantysaari@jyu.fi>

The code can be built using cmake easily (needs GSL, FFTW and Boost libraries)
 mkdir build
 cd build
 cmake ..
 make

The binary file is generated in build/bin/jimwlk

The code reads one command line parameter which is the file name where all parameters are set. Tempate file is provided, see file input. The parameters in the input file are expained in detail in the source file src/Parameters.h

The initial condition can either be generated at the beginning (initMethod 1-3), 
or one can read in a set of Wilson lines generated using the IPGlasma code (https://github.com/schenke/ipglasma) 
(10: txt format, 11: binary format). The filename for the initial Wilson lines 
is the input_wline parameter. The evolved Wilson lines are saved in 
evolved_wilson_lines directory, and the filename contains the number 
of evolution steps s (filename: inputfile_steps_s). 

At fixed coupling, steps is related to the Bjorken x via 
steps =  alphas * log(x0/x)/ (pi*pi*ds), 
where ds is defined in input file and x0 can be taken to be 0.01 or so.
At running coupling, steps=log(x0/x)/(pi*pi*ds)


In the input file, if simpleLangevin is set to 1, the JIMWLK is solved using the method from 1212.4825 which does not require adjoint Wilson lines -> faster, and needs less memory

At running coupling, the parameter m that suppresses long distance Coulomb tails is given in GeV in the input file, see 1407.8458 for details
