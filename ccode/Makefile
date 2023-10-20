smacofSSWR: Makefile smacofEngine.c smacofInitial.c smacofUtils.c smacofCore.c smacofNorm.c smacof.h
	clang -o smacof -O2 smacofEngine.c smacofInitial.c smacofUtils.c smacofCore.c smacofNorm.c

rlib: 
	R CMD SHLIB smacofEngine.c smacofInitial.c smacofUtils.c smacofCore.c smacofNorm.c

clean:
	rm -rf *.o *.o core a.out

pristine:
	rm -rf smacof smacof.so

format:
	clang-format -i *.c *.h

