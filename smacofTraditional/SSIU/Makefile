CFLAGS = -O2 -Wall -Wextra -Werror -pedantic -std=c99
DFLAGS = -g3 -Wall -Wextra -Werror -pedantic -std=c99

sources = \
	../Common/smacofCommonPrintRead.c \
	../Common/smacofCommonRCTranslation.c \
	../Common/smacofCommonUtils.c \
    ../Common/smacofCommonMatrixUtils.c \
	../Common/smacofCommonLinearAlgebra.c \
	../Common/smacofCommonAccelerate.c \
	../Unweighted/smacofUnweightedNorm.c \
	../Unweighted/smacofUnweightedInitial.c \
	../Unweighted/smacofUnweightedCore.c \
	../Unweighted/smacofUnweightedTransforms.c

objects = \
	smacofCommonPrintRead.o \
	smacofCommonRCTranslation.o \
	smacofCommonUtils.o \
	smacofCommonMatrixUtils.o \
	smacofCommonLinearAlgebra.o \
	smacofCommonAccelerate.o \
	smacofUnweightedNorm.o \
	smacofUnweightedInitial.o \
	smacofUnweightedCore.o \
	smacofUnweightedTransforms.o

drivers = smacofSSIU.c smacofSSIUEngine.c

includes = ../Unweighted/smacofUnweighted.h ../Common/smacofCommon.h

all: smacofSSIU

ofiles: $(includes) $(sources)
	$(CC) -c $(CFLAGS) $(sources)

clib: ofiles 
	ar rcs libsmacofSSIU.a $(objects)
	mv -f libsmacofSSIU.a /usr/local/lib
	rm -f *.o

smacofSSIU: clib
	$(CC) -o smacofSSIU $(CFLAGS) $(drivers) -lsmacofSSIU
	mv -f smacofSSIU /usr/local/bin

smacofSSIUDebug: $(includes) $(sources) $(drivers)
	$(CC) -o smacofSSIUDebug $(DFLAGS) $(drivers) $(sources)

rlib: ofiles
	R CMD SHLIB -o libsmacofSSIU.so smacofSSIUEngine.c libsmacofSSIU.a

format: $(includes) $(sources)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes)
	
clean:
	rm -f $(objects) ../*/*.o

pristine: clean
	rm -f smacofSSIU smacofSSIUDebug libsmacofSSIU.* 