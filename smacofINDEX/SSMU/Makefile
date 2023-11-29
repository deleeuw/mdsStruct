sources = smacofSSMUEngine.c \
	../Common/smacofCommonPrint.c \
	../Common/smacofCommonUtils.c \
	../Common/smacofCommonLinearAlgebra.c \
	../Common/smacofCommonAccelerate.c \
	../Common/smacofCommonIndices.c \
	../Unweighted/smacofUnweightedNorm.c \
	../Unweighted/smacofUnweightedInitial.c \
	../Unweighted/smacofUnweightedCore.c \
	../Unweighted/smacofUnweightedRead.c 
 

objects = smacofSSMUEngine.o \
	smacofCommonPrint.o \
	smacofCommonUtils.o \
	smacofCommonLinearAlgebra.o \
	smacofCommonAccelerate.o \
	smacofCommonIndices.o \
	smacofUnweightedNorm.o \
	smacofUnweightedInitial.o \
	smacofUnweightedCore.o \
	smacofUnweightedRead.o 


includes = ../Unweighted/smacofUnweighted.h

all: smacof

ofiles: $(includes) $(sources)
	$(CC) -c -O2 -Wall -Wextra -pedantic -std=c99 $(sources)

debug: clean $(includes) $(sources)
	$(CC) -c -g3 -Wall -Wextra -pedantic -std=c99 $(sources)
	$(CC) -o smacofSSMU -g3  -Wall -Wextra -pedantic -std=c99 $(objects)

smacof: ofiles
	$(CC) -o smacofSSMU -O2  -Wall -Wextra -pedantic -std=c99 $(objects)

clib: ofiles 
	ar rcs libsmacofSSMU.a $(objects)

rlib: ofiles
	R CMD SHLIB -o libsmacofSSMU.so $(objects)

dylib: ofiles
	$(CC) -O2 -Wall -Wall -Winline -pedantic -Wextra -std=c99 -dynamiclib  \
		-current_version 0.0 -compatibility_version 0.1 \
		$(objects) -o libsmacofSSMU.0.dylib 
	ln -s libsmacofSSMU.0.dylib libsmacofSSMU.dylib

format: $(includes) $(sources)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes)

test: smacof
	./smacof > smacofSSMU.txt
	
clean:
	rm -f $(objects) ../../*.c

pristine: clean
	rm -f smacofSSMU smacofSSMU.txt \
		libsmacofSSMU.a libsmacofSSMU.so libsmacofSSMU.0.dylib libsmacofSSMU.dylib