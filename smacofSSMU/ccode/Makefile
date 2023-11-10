all: run

smacofSSMU1: smacofSSMU1.c 
	clang -o smacofSSMU1 -O2 -Wall -pedantic smacofSSMU1.c -lsmacof -L/opt/homebrew/lib -lprofiler

run: smacofSSMU1 
	./smacofSSMU1 > smacofSSMU1.txt

clean:
	rm -f smacofSSMU[1-9]