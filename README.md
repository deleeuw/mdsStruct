# mdsStruct
This is code for a new version of the MDS smacof project (for now only the metric case). The C code is in the subdirectory ccode. smacofEngine.c has
a small main with hardcoded data, so "make" in the subdirectory builds an
executable. "make rlib" builds a shared library to be called from R.

There is R code for smacofR.R and smacofRC.R in the rcode subdirectory.
smacofR is the usual smacof, completely in R, using matrices throughout.
smacofRC calls the shared library and does all its computations in
compiled object code.

The number of iterations in smacof is determined by the natures of the
problem (how flat is my stress, how isolated is my local minimum) but
also by the parameters itmax (maximum number of iterations) and eps
(stop if the function decreases by less than eps in an iteration). The
subdirectory timing compares microbenchmark running times for various 
combinations of itmax and eps on six examples. 
