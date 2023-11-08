sources = ../smacofEngines/smacofSSMWEngine.c \
	../smacofEngines/smacofSSMUEngine.c \
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

all: clib rlib

clib: $(includes) $(sources)
	@$(CC) -c -O2 -g -Wall -Winline -pedantic $(sources)
	@ar rcs libsmacof.a $(objects)

#dylib: $(includes) $(sources)
#	$(CC) -O2 -Wall -dynamiclib -current_version 1.0 -compatibility_version 1.0 \
#	$(sources) -o libsmacof.1.dylib
#	@ln -s libsmacof.1.dylib libsmacof.dylib
#	@rm -f *.o $(objects)

rlib: clib $(includes) $(sources)
	@R CMD SHLIB -o libsmacof.so $(sources)
	@rm -f *.o $(objects)

format: $(includes) $(sources)
	@clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" $(sources)
	
clean:
	@rm -f *.o $(objects) ../smacofCommon/ccode/*.o

pristine: clean
	rm -rf libsmacof.a libsmacof.so *.dylib