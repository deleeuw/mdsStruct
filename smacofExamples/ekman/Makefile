
runEkman: runEkmanU.c runEkmanW.c 
	clang -o runEkmanU -O2 -Wall -pedantic runEkmanU.c -lsmacof
	clang -o runEkmanW -O2 -Wall -pedantic runEkmanW.c -lsmacof

timer:
	hyperfine -m 100 -w 25 -N runEkmanU runEkmanW

clean:
	rm -f runEkmanU runEkmanW *.o