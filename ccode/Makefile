smacof: Makefile smacofEngine.c smacofTorgerson.c smacofUtils.c \
smacofCore.c smacofNorm.c smacofPrint.c smacofLinearAlgebra.c smacof.h
	gcc -o smacof -O2 smacofEngine.c smacofTorgerson.c smacofUtils.c \
	smacofCore.c smacofNorm.c smacofLinearAlgebra.c smacofPrint.c

debug: Makefile smacofEngine.c smacofTorgerson.c smacofUtils.c \
smacofCore.c smacofNorm.c smacofPrint.c smacofLinearAlgebra.c smacof.h
	gcc -o smacofdebug -g smacofEngine.c smacofTorgerson.c smacofUtils.c \
	smacofCore.c smacofNorm.c smacofLinearAlgebra.c smacofPrint.c

rlib: 
	R CMD SHLIB smacofEngine.c smacofTorgerson.c smacofUtils.c \
	smacofCore.c smacofNorm.c smacofLinearAlgebra.c smacofPrint.c

clib:
	ar rcs libSmacof.a *.o

clean:
	rm -rf *.o *.dSYM core a.out

pristine: clean
	rm -rf smacof *.so *.a *.txt

format:
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" *.c *.h

