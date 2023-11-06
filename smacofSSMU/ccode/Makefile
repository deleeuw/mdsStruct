all: run

smacofSSMU1: smacofSSMU1.c 
	clang -o smacofSSMU1 -O2 smacofSSMU1.c -lsmacof

run: smacofSSMU1 
	./smacofSSMU1 > smacofSSMU1.txt

clean:
	rm -f smacofSSMU[1-9]