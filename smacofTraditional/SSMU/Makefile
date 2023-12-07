CFLAGS = -O2 -Wall -Wextra -pedantic -std=c99
DFLAGS = -g3 -Wall -Wextra -pedantic -std=c99

sources = \
	../Common/smacofCommonPrintRead.c \
	../Common/smacofCommonRCTranslation.c \
	../Common/smacofCommonUtils.c \
    ../Common/smacofCommonMatrixUtils.c \
	../Common/smacofCommonLinearAlgebra.c \
	../Common/smacofCommonAccelerate.c \
	../Unweighted/smacofUnweightedNorm.c \
	../Unweighted/smacofUnweightedInitial.c \
	../Unweighted/smacofUnweightedCore.c 

objects = \
	smacofCommonPrintRead.o \
	smacofCommonRCTranslation.o \
	smacofCommonUtils.o \
	smacofCommonMatrixUtils.o \
	smacofCommonLinearAlgebra.o \
	smacofCommonAccelerate.o \
	smacofUnweightedNorm.o \
	smacofUnweightedInitial.o \
	smacofUnweightedCore.o 

drivers = smacofSSMU.c \
	smacofSSMUEngine.c

includes = ../Unweighted/smacofUnweighted.h ../Common/smacofCommon.h

all: smacofSSMU

ofiles: $(includes) $(sources)
	$(CC) -c $(CFLAGS) $(sources)

clib: ofiles 
	ar rcs libsmacofSSMU.a $(objects)
	mv -f libsmacofSSMU.a /usr/local/lib
	rm -f *.o

smacofSSMU: clib
	$(CC) -o smacofSSMU $(CFLAGS) $(drivers) -lsmacofSSMU
	mv -f smacofSSMU /usr/local/bin

smacofSSMUDebug: $(includes) $(sources) $(drivers)
	$(CC) -o smacofSSMUDebug $(DFLAGS) $(drivers) $(sources)

rlib: ofiles
	R CMD SHLIB -o libsmacofSSMU.so smacofSSMUEngine.c libsmacofSSMU.a

format: $(includes) $(sources)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes)
	
clean:
	rm -f $(objects) ../*/*.o

pristine: clean
	rm -f smacofSSMU smacofSSMUDebug \
		libsmacofSSMU.a libsmacofSSMU.so libsmacofSSMU.0.dylib libsmacofSSMU.dylib