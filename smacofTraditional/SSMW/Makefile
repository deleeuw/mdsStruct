CFLAGS = -O2 -Wall -Wextra -Werror -pedantic -std=c99
DFLAGS = -g3 -Wall -Wextra -Werror -pedantic -std=c99

sources = \
	../Common/smacofCommonPrintRead.c \
	../Common/smacofCommonRCTranslation.c \
	../Common/smacofCommonUtils.c \
    ../Common/smacofCommonMatrixUtils.c \
	../Common/smacofCommonLinearAlgebra.c \
	../Common/smacofCommonAccelerate.c \
	../Weighted/smacofWeightedNorm.c \
	../Weighted/smacofWeightedInitial.c \
	../Weighted/smacofWeightedDerivatives.c \
	../Weighted/smacofWeightedCore.c 

objects = \
	smacofCommonPrintRead.o \
	smacofCommonRCTranslation.o \
	smacofCommonUtils.o \
	smacofCommonMatrixUtils.o \
	smacofCommonLinearAlgebra.o \
	smacofCommonAccelerate.o \
	smacofWeightedNorm.o \
	smacofWeightedInitial.o \
	smacofWeightedDerivatives.o \
	smacofWeightedCore.o 

drivers = smacofSSMW.c \
	smacofSSMWEngine.c

drovers = smacofSSMW.c \
	smacofSSMWEngineDebug.c

includes = ../Weighted/smacofWeighted.h ../Common/smacofCommon.h

all: smacofSSMW

ofiles: $(includes) $(sources)
	$(CC) -c $(CFLAGS) $(sources)

clib: ofiles 
	ar rcs libsmacofSSMW.a $(objects)
	mv -f libsmacofSSMW.a /usr/local/lib
	rm -f *.o

smacofSSMW: clib
	$(CC) -o smacofSSMW $(CFLAGS) $(drivers) -lsmacofSSMW
	mv -f smacofSSMW /usr/local/bin

smacofSSMWDebug: $(includes) $(sources) $(drivers)
	$(CC) -o smacofSSMW $(DFLAGS) $(drivers) $(sources)

rlib: ofiles
	R CMD SHLIB -o libsmacofSSMW.so smacofSSMWEngine.c libsmacofSSMW.a

format: $(includes) $(sources)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes)
	
clean:
	rm -f $(objects) ../*/*.o

pristine: clean
	rm -f smacofSSMW smacofSSMWDebug smacofSSMW.txt \
		libsmacofSSMW.a libsmacofSSMW.so libsmacofSSMW.0.dylib libsmacofSSMW.dylib