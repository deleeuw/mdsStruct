smacof: Makefile smacofEngine.c smacofTorgerson.c smacofUtils.c smacofCore.c smacofNorm.c smacofPrint.c smacof.h
	gcc -o smacof -O2 smacofEngine.c smacofTorgerson.c smacofUtils.c smacofCore.c smacofNorm.c smacofPrint.c

rlib: 
	R CMD SHLIB smacofEngine.c smacofTorgerson.c smacofUtils.c smacofCore.c smacofNorm.c smacofPrint.c

clib:
	ar rcs libSmacof.a *.o

clean:
	rm -rf *.o *.o core a.out

pristine:
	rm -rf smacof smacofEngine.so libSmacof.a

format:
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" *.c *.h

