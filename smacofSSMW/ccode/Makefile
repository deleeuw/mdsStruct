smacofSSMW: smacof.h *.c
	$(CC) -o smacofSSMW -O2 *.c

debug: smacof.h *.c
	$(CC) -o smacofSSMW -g *.c

rlib: smacof.h *.c
	R CMD SHLIB -o smacofSSMW.so *.c
	
ofiles: smacof.h *.c
	clang -c -O2 *.c
		
clib: ofiles
	ar rcs libSmacof.a *.o

clean:
	rm -rf *.o *.dSYM core a.out

pristine: clean
	rm -rf smacofSSMW *.so *.a *.txt

format:
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" *.c *.h

