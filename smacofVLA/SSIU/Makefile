sources = smacofSSIUEngine.c \
	../Common/smacofCommonPrintRead.c \
	../Common/smacofCommonRCTranslation.c \
	../Common/smacofCommonUtils.c \
	../Common/smacofCommonLinearAlgebra.c \
	../Common/smacofCommonAccelerate.c \
	../Unweighted/smacofUnweightedNorm.c \
	../Unweighted/smacofUnweightedInitial.c \
	../Unweighted/smacofUnweightedCore.c \
    ../Unweighted/smacofUnweightedTransforms.c

objects = smacofSSIUEngine.o \
	smacofCommonPrintRead.o \
	smacofCommonRCTranslation.o \
	smacofCommonUtils.o \
	smacofCommonLinearAlgebra.o \
	smacofCommonAccelerate.o \
	smacofUnweightedNorm.o \
	smacofUnweightedInitial.o \
	smacofUnweightedCore.o \
	smacofUnweightedTransforms.o

includes = ../Include/smacof.h

all: smacof

ofiles: $(includes) $(sources)
	$(CC) -c -O2 -Wall -Wextra -pedantic -std=c99 $(sources)

debug: clean $(includes) $(sources)
	$(CC) -c -g3 -Wall -Wextra -pedantic -std=c99 $(sources)
	$(CC) -o smacofSSIUd -g3  -Wall -Wextra -pedantic -std=c99 $(objects)

smacof: ofiles
	$(CC) -o smacofSSIU -O2  -Wall -Wextra -pedantic -std=c99 $(objects)

clib: ofiles 
	ar rcs libsmacofSSIU.a $(objects)

rlib: ofiles
	R CMD SHLIB -o libsmacofSSIU.so $(objects)

dylib: ofiles
	$(CC) -O2 -Wall -Wall -Winline -pedantic -Wextra -std=c99 -dynamiclib  \
		-current_version 0.0 -compatibility_version 0.1 \
		$(objects) -o libsmacofSSIU.0.dylib 
	ln -s libsmacofSSIU.0.dylib libsmacofSSIU.dylib

format: $(includes) $(sources)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes)

test: smacof
	./smacof > smacofSSIU.txt
	
clean:
	rm -f $(objects) ../../*.c

pristine: clean
	rm -f smacofSSIU smacofSSIUd smacofSSIU.txt \
		libsmacofSSIU.a libsmacofSSIU.so libsmacofSSIU.0.dylib libsmacofSSIU.dylib