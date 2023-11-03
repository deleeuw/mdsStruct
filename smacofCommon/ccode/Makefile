clib: ofiles
	ar rcs /usr/local/lib/libSmacof.a *.o
	
ofiles: include *.c
	clang -c -O2 *.c
		
include: smacof.h
	cp smacof.h /usr/local/include

clean:
	rm -rf *.o *.dSYM core a.out

pristine: clean
	rm -rf *.so *.a *.txt

format:
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" *.c *.h

