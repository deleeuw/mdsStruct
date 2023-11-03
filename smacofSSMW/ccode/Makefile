smacofSSMW: *.c \
../../smacofCommon/ccode/*.c \
../../smacofCommon/ccode/*.h \
/usr/local/lib/libSmacof.a
	$(CC) -o smacofSSMW -O2 *.c -lSmacof

rlib: ../../smacofCommon/ccode/*.h *.c
	R CMD SHLIB -o smacofSSMW.so *.c -lSmacof

clean: clean
	rm -rf smacofSSMW *.o *.dSym a.out core *.so *.a *.txt

format:
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" *.c

