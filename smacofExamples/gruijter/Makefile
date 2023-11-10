
runEkman: runGruijterU.c runGruijterW.c 
	clang -o runGruijterU -O2 -Wall -pedantic runGruijterU.c -lsmacof
	clang -o runGruijterW -O2 -Wall -pedantic runGruijterW.c -lsmacof

timer:
	hyperfine -m 100 -w 25 runGruijterU runGruijterW

clean:
	rm -f runGruijterU runGruijterW *.o