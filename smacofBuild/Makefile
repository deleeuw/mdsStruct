sources = ../smacofEngines/smacof*Engine.c \
	../smacofCommon/ccode/smacofCommon*.c \
	../smacofCommon/ccode/smacofWeighted*.c \
	../smacofCommon/ccode/smacofUnweighted*.c 

objects = ../smacofEngines/smacof*Engine.o \
	../smacofCommon/ccode/smacofCommon*.o \
	../smacofCommon/ccode/smacofWeighted*.o \
	../smacofCommon/ccode/smacofUnweighted*.o

includes = ../smacofInclude/smacof.h

all: clib dylib rlib

clib: $(includes) $(sources)
	$(CC) -c -O2 -Wall $(sources)
	ar rcs libsmacof.a *.o
	rm -f *.o


dylib: $(includes) $(sources)
	$(CC) -c -g -Wall -fPIC -O2 ../smacofEngines/*.c ../smacofCommon/ccode/*.c 
	$(CC) -arch arm64 -dynamiclib -Wl,-headerpad_max_install_names \
		-undefined dynamic_lookup -o libsmacof.dylib *.o
	rm -f *.o

rlib: $(includes) $(sources)
	R CMD SHLIB -o smacofShareLib.so $(sources)
	rm -f *.o

install: 
	cp ../smacofInclude/smacof.h /usr/local/include 
	cp libsmacof.a /usr/local/lib
	cp libsmacof.dylib /usr/local/lib

format: $(includes) $(sources)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" $(sources)
	
clean:
	rm -f *.o

pristine: clean
	rm -rf ../smacofBinaries/* 

