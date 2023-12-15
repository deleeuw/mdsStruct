CFLAGS = -O2 -Wall -Wextra -Werror -pedantic -std=c99
DFLAGS = -g3 -Wall -Wextra -Werror -pedantic -std=c99

sources = \
	smacofPrintRead.c \
	smacofRCTranslation.c \
	smacofUtils.c \
    smacofMatrixUtils.c \
	smacofLinearAlgebra.c \
	smacofAccelerate.c \
	smacofNorm.c \
	smacofInitial.c \
	smacofBases.c \
	smacofTransforms.c \
	smacofCCD.c \
	smacofIndices.c \
	smacofCore.c 

notyet = \
	smacofQP.c \
	smacofSort.c \
	smacofIsotone.c \
	smacofDerivatives.c 

objects = \
	smacofPrintRead.o \
	smacofRCTranslation.o \
	smacofUtils.o \
	smacofMatrixUtils.o \
	smacofLinearAlgebra.o \
	smacofAccelerate.o \
	smacofNorm.o \
	smacofInitial.o \
	smacofBases.o \
	smacofTransforms.o \
	smacofCCD.o \
	smacofIndices.o \
	smacofCore.o 
	
drivers = smacofSS.c \
	smacofSSEngine.c

includes = smacof.h

all: smacofSS

ofiles: $(includes) $(sources)
	$(CC) -c $(CFLAGS) $(sources)

clib: ofiles 
	ar rcs libsmacofSS.a $(objects)
	mv -f libsmacofSS.a /usr/local/lib
	rm -f *.o

smacofSS: clib
	$(CC) -o smacofSS $(CFLAGS) $(drivers) -lsmacofSS
	mv -f smacofSS /usr/local/bin
	
smacofSSDebug: $(includes) $(sources) $(drivers)
	$(CC) -o smacofSSDebug $(DFLAGS) $(drivers) $(sources)

rlib: ofiles
	R CMD SHLIB -o libsmacofSS.so smacofSSEngine.c $(LLIB)/libsmacofSS.a

format: $(includes) $(sources) $(drivers)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes) $(drivers)
	
clean:
	rm -f $(objects) ../*/*.o

pristine: clean
	rm -f smacofSS smacofSSDebug \
		libsmacofSS.a libsmacofSS.so libsmacofSS.0.dylib libsmacofSS.dylib

