
ofiles: smacof.h *.c
	clang -c -O2 *.c

clib: ofiles
	ar rcs /usr/local/lib/libSmacof.a *.o
		
clean:
	rm -rf *.o *.dSYM core a.out

pristine: clean
	rm -rf *.so *.a *.txt

format:
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" *.c *.h

