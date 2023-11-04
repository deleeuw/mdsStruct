execs: Makefile \
*.c \
../../smacofCommon/ccode/*.c \
../../smacofCommon/ccode/*.h \
/usr/local/lib/libSmacof.a
	$(CC) -o smacofSSMW1 -O2 smacofSSMW1.c smacofSSMWEngine.c -lSmacof
	$(CC) -o smacofSSMW2 -O2 smacofSSMW2.c smacofSSMWEngine.c -lSmacof
	$(CC) -o smacofSSMW3 -O2 smacofSSMW3.c smacofSSMWEngine.c -lSmacof

rlib: ../../smacofCommon/ccode/*.h *.c
	R CMD SHLIB -o smacofSSMW.so smacofSSMWEngine.c -lSmacof

clean: 
	rm -rf *.o *.dSym 

pristine: clean
	rm -rf *.out smacofSSMW? *.so 

format:
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" *.c

