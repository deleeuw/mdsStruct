
runEkman: runEqdistU.c runEqdistW.c 
	clang -o runEqdistU -O2 -Wall -pedantic runEqdistU.c -lsmacof
	clang -o runEqdistW -O2 -Wall -pedantic runEqdistW.c -lsmacof

timer:
	hyperfine -m 100 -w 25 runEqdistU runEqdistW

clean:
	rm -f runEqdistU runEqdistW *.o