smacof: Makefile smacofEngine.c smacofTorgerson.c smacofUtils.c smacofCore.c smacofNorm.c smacof.h
	gcc -o smacof -O2 smacofEngine.c smacofTorgerson.c smacofUtils.c smacofCore.c smacofNorm.c

rlib: 
	R CMD SHLIB smacofEngine.c smacofTorgerson.c smacofUtils.c smacofCore.c smacofNorm.c

clib:
	ar rcs libSmacof.a *.o

clean:
	rm -rf *.o *.o core a.out

pristine:
	rm -rf smacof smacofEngine.so libSmacof.a

format:
	clang-format -i *.c *.h

