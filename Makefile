sources = smacofEngines/smacofSSMWEngine.c \
	smacofCommon/ccode/smacofCommon*.c \
	smacofCommon/ccode/smacofWeighted*.c \
	smacofCommon/ccode/smacofUnweighted*.c 

objects = smacofEngines/smacofSSMWEngine.o \
	smacofCommon/ccode/smacofCommon*.o \
	smacofCommon/ccode/smacofWeighted*.o \
	smacofCommon/ccode/smacofUnweighted*.o
	
ofiles: smacofInclude/smacof.h $(sources)
	clang -c -O2 smacofEngines/*.c smacofCommon/ccode/*.c

clib: ofiles
	ar rcs smacofBinaries/libsmacof.a *.o

rlib: ofiles
	R CMD SHLIB -o smacofBinaries/smacofShareLib.so *.o

install: 
	cp smacofInclude/smacof.h /usr/local/include 
	cp smacofBinaries/libsmacof.a /usr/local/lib
	
clean:
	rm -rf *.o $(objects)

pristine: clean
	rm -rf smacofBinaries/* 

