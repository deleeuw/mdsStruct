sources = smacofSSMUEngine.c \
	../Common/smacofCommonPrintRead.c \
	../Common/smacofCommonRCTranslation.c \
	../Common/smacofCommonUtils.c \
	../Common/smacofCommonLinearAlgebra.c \
	../Common/smacofCommonAccelerate.c \
	../Unweighted/smacofUnweightedNorm.c \
	../Unweighted/smacofUnweightedInitial.c \
	../Unweighted/smacofUnweightedCore.c 

objects = smacofSSMUEngine.o \
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
	$(CC) -c -O2 -Wall -Wextra -pedantic -std=c99 $(sources)

debug: clean $(includes) $(sources)
	$(CC) -c -g3 -Wall -Wextra -pedantic -std=c99 $(sources)
	$(CC) -o smacofSSMUd -g3  -Wall -Wextra -pedantic -std=c99 $(objects)

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
	rm -f smacofSSMU smacofSSMUd smacofSSMU.txt \
		libsmacofSSMU.a libsmacofSSMU.so libsmacofSSMU.0.dylib libsmacofSSMU.dylib