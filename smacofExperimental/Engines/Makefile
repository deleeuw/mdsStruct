sources = smacofExperimentalEngine.c \
	../Common/smacofCommonPrintRead.c \
	../Common/smacofCommonRCTranslation.c \
	../Common/smacofCommonUtils.c \
	../Common/smacofCommonLinearAlgebra.c \
	../Common/smacofCommonAccelerate.c \
	../SSMU/smacofUnweightedNorm.c \
	../SSMU/smacofUnweightedInitial.c \
	../SSMU/smacofUnweightedCore.c

objects = smacofExperimentalEngine.o \
	smacofCommonPrintRead.o \
	smacofCommonRCTranslation.o \
	smacofCommonUtils.o \
	smacofCommonLinearAlgebra.o \
	smacofCommonAccelerate.o \
	smacofUnweightedNorm.o \
	smacofUnweightedInitial.o \
	smacofUnweightedCore.o

includes = ../Include/smacof.h

all: smacof

ofiles: $(includes) $(sources)
	$(CC) -c -O2 -w -std=c99 $(sources)

smacof: ofiles
	$(CC) -o smacof -O2 -w \
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
	rm -f smacof libsmacof.a libsmacof.so libsmacof.0.dylib libsmacof.dylib