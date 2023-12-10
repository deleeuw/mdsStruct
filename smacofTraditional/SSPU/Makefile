CFLAGS = -O2 -Wall -Wextra -Werror -pedantic -std=c99
DFLAGS = -g3 -Wall -Wextra -Werror -pedantic -std=c99

sources = \
	../Common/smacofCommonPrintRead.c \
	../Common/smacofCommonRCTranslation.c \
	../Common/smacofCommonUtils.c \
    ../Common/smacofCommonMatrixUtils.c \
	../Common/smacofCommonLinearAlgebra.c \
	../Common/smacofCommonAccelerate.c \
	../Common/smacofCommonBases.c \
	../Unweighted/smacofUnweightedNorm.c \
	../Unweighted/smacofUnweightedInitial.c \
	../Unweighted/smacofUnweightedCore.c \
	../Unweighted/smacofUnweightedCCD.c \
	../Unweighted/smacofUnweightedTransforms.c

objects = \
	smacofCommonPrintRead.o \
	smacofCommonRCTranslation.o \
	smacofCommonUtils.o \
	smacofCommonMatrixUtils.o \
	smacofCommonLinearAlgebra.o \
	smacofCommonAccelerate.o \
	smacofCommonBases.o \
	smacofUnweightedNorm.o \
	smacofUnweightedInitial.o \
	smacofUnweightedCore.o \
	smacofUnweightedCCD.o \
	smacofUnweightedTransforms.o

drivers = smacofSSPU.c smacofSSPUEngine.c

includes = ../Unweighted/smacofUnweighted.h ../Common/smacofCommon.h

all: smacofSSPU

ofiles: $(includes) $(sources)
	$(CC) -c $(CFLAGS) $(sources)

clib: ofiles 
	ar rcs libsmacofSSPU.a $(objects)
	mv -f libsmacofSSPU.a /usr/local/lib
	rm -f *.o

smacofSSPU: clib
	$(CC) -o smacofSSPU $(CFLAGS) $(drivers) -lsmacofSSPU
	mv -f smacofSSPU /usr/local/bin

smacofSSPUDebug: $(includes) $(sources) $(drivers)
	$(CC) -o smacofSSPUDebug $(DFLAGS) $(drivers) $(sources)

rlib: ofiles
	R CMD SHLIB -o libsmacofSSPU.so smacofSSPUEngine.c libsmacofSSPU.a

format: $(includes) $(sources)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes)
	
clean:
	rm -f $(objects) ../*/*.o

pristine: clean
	rm -f smacofSSPU smacofSSPUDebug libsmacofSSPU.* 