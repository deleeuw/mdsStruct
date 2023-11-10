
runRoskam: runRoskam.c 
	clang -o runRoskam -O2 -Wall -pedantic runRoskam.c -L/usr/local/lib -lsmacof

timer:
	hyperfine -m 100 -w 25 runRoskam

clean:
	rm -f runRoskam