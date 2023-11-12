sources = Engines/smacofExperimentalEngine.c \
	Common/smacofCommonPrintRead.c \
	Common/smacofCommonLinearAlgebra.c \
	Common/smacofCommonUtils.c \
	Common/smacofCommonAccelerate.c \
	Common/smacofCommonRCTranslation.c \
	SSMU/smacofUnweightedNorm.c \
	SSMU/smacofUnweightedInitial.c \
	SSMU/smacofUnweightedCore.c 

objects = Engines/smacofExperimentalEngine.o \
	Common/smacofCommonPrintRead.o \
	Common/smacofCommonLinearAlgebra.o \
	Common/smacofCommonUtils.o \
	Common/smacofCommonAccelerate.o \
	Common/smacofCommonRCTranslation.o \
	SSMU/smacofUnweightedNorm.o \
	SSMU/smacofUnweightedInitial.o \
  	SSMU/smacofUnweightedCore.o 

includes = ./Include/smacof.h

all: clib rlib dylib

ofiles: $(includes) $(sources)
	$(CC) -c -O2 -g3 -Wall -Winline -Wextra -pedantic -std=c99 $(sources)

smacof: ofiles
	$(CC) -o smacof -O2 -g3 -Wall -Winline -Wextra -pedantic \
	-std=c99 $(objects)

clib: ofiles
	ar rcs libsmacof.a $(objects)

rlib: ofiles
	R CMD SHLIB -o libsmacof.so $(objects)

dylib: ofiles
	$(CC) -O2 -Wall -Wall -Winline -pedantic -Wextra -std=c99 -dynamiclib  \
		-current_version 0.0 -compatibility_version 0.1 \
		$(objects) -o libsmacof.0.dylib 
	ln -s libsmacof.0.dylib libsmacof.dylib

format: $(includes) $(sources)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes)
	
clean:
	rm -f $(objects)

pristine: clean
	rm -f libsmacof.a libsmacof.so libsmacof.0.dylib libsmacof.dylib