sources = ../smacofEngines/smacofSSMWEngine.c \
	../smacofEngines/smacofSSMUEngine.c \
	../smacofCommon/ccode/smacofCommonIndices.c \
	../smacofCommon/ccode/smacofCommonAccelerate.c \
	../smacofCommon/ccode/smacofCommonLinearAlgebra.c \
	../smacofCommon/ccode/smacofCommonPrintRead.c \
	../smacofCommon/ccode/smacofCommonUtils.c \
	../smacofCommon/ccode/smacofWeightedCore.c \
	../smacofCommon/ccode/smacofWeightedDerivatives.c \
	../smacofCommon/ccode/smacofWeightedInitial.c \
	../smacofCommon/ccode/smacofWeightedNorm.c \
	../smacofCommon/ccode/smacofUnweightedCore.c \
	../smacofCommon/ccode/smacofUnweightedInitial.c \
	../smacofCommon/ccode/smacofUnweightedNorm.c 


objects = smacofSSMWEngine.o \
	smacofSSMUEngine.o \
	smacofCommonIndices.o \
	smacofCommonAccelerate.o \
	smacofCommonLinearAlgebra.o \
	smacofCommonPrintRead.o \
	smacofCommonUtils.o \
	smacofWeightedCore.o \
	smacofWeightedDerivatives.o \
	smacofWeightedInitial.o \
	smacofWeightedNorm.o \
	smacofUnweightedCore.o \
	smacofUnweightedInitial.o \
	smacofUnweightedNorm.o

includes = ../smacofInclude/smacof.h

all: ofiles clib rlib dylib

ofiles: $(includes) $(sources)
	$(CC) -c -O2 -g -Wall -Winline -pedantic $(sources)

clib: ofiles
	ar rcs libsmacof.a $(objects)

rlib: ofiles
	R CMD SHLIB -o libsmacof.so $(objects)

dylib: ofiles
	$(CC) -O2 -Wall -Wall -Winline -pedantic -dynamiclib  \
		-current_version 0.0 -compatibility_version 0.1 \
		$(objects) -o libsmacof.0.dylib 
	ln -s libsmacof.0.dylib libsmacof.dylib

format: $(includes) $(sources)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" $(sources) $(includes)
	
clean:
	rm -f $(objects)

pristine: clean
	rm -rf libsmacof.a libsmacof.so *.dylib