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
	smacofCore.c 

objects = \
	smacofPrintRead.o \
	smacofRCTranslation.o \
	smacofUtils.o \
	smacofMatrixUtils.o \
	smacofLinearAlgebra.o \
	smacofAccelerate.o \
	smacofNorm.o \
	smacofInitial.o \
	smacofCore.o 

drivers = smacofSSM.c \
	smacofSSMEngine.c

includes = smacof.h

all: smacofSSM

ofiles: $(includes) $(sources)
	$(CC) -c $(CFLAGS) $(sources)

clib: ofiles 
	ar rcs libsmacofSSM.a $(objects)
	mv -f libsmacofSSM.a /usr/local/lib
	rm -f *.o

smacofSSM: clib
	$(CC) -o smacofSSM $(CFLAGS) $(drivers) -lsmacofSSM
	mv -f smacofSSM /usr/local/bin
	
smacofSSMUDebug: $(includes) $(sources) $(drivers)
	$(CC) -o smacofSSMDebug $(DFLAGS) $(drivers) $(sources)

rlib: ofiles
	R CMD SHLIB -o libsmacofSSM.so smacofSSMEngine.c libsmacofSSM.a

format: $(includes) $(sources) $(drivers)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes) $(drivers)
	
clean:
	rm -f $(objects) ../*/*.o

pristine: clean
	rm -f smacofSSM smacofSSMDebug \
		libsmacofSSM.a libsmacofSSM.so libsmacofSSM.0.dylib libsmacofSSM.dylib