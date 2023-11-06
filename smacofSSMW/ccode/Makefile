all: run

smacofSSMW1: smacofSSMW1.c /usr/local/lib/libsmacof.a
	clang -o smacofSSMW1 -O2 smacofSSMW1.c -lsmacof

smacofSSMW2: smacofSSMW2.c /usr/local/lib/libsmacof.a
	clang -o smacofSSMW2 -O2 smacofSSMW2.c -lsmacof

smacofSSMW3: smacofSSMW3.c /usr/local/lib/libsmacof.a
	clang -o smacofSSMW3 -O2 smacofSSMW3.c -lsmacof

run: smacofSSMW1 smacofSSMW2 smacofSSMW3
	./smacofSSMW1 > smacofSSMW1.txt
	./smacofSSMW2 > smacofSSMW2.txt
	./smacofSSMW3 > smacofSSMW3.txt

clean:
	rm -f smacofSSMW[1-9]