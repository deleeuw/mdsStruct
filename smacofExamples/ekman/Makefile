runEkman: runEkmanU.c runEkmanW.c 
	clang -o runEkmanU -O2 runEkmanU.c -L../../smacofBuild -lsmacof
	clang -o runEkmanW -O2 runEkmanW.c -L../../smacofBuild -lsmacof

timer: 
	/usr/bin/time -ph runEkmanU
	/usr/bin/time -ph runEkmanW


clean:
	rm -f runEkmanU runEkmanW *.o