
runEkman: runMorseU.c runMorseW.c 
	clang -o runMorseU -O2 -Wall -pedantic runMorseU.c -L/usr/local/lib -lsmacof
	clang -o runMorseW -O2 -Wall -pedantic runMorseW.c -L/usr/local/lib -lsmacof

timer:
	hyperfine -m 100 -w 25 runMorseU runMorseW

clean:
	rm -f runMorseU runMorseW *.o