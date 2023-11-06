runEkman: runEkmanU.c runEkmanW.c /usr/local/lib/libsmacof.a
	clang -o runEkmanU -O2 runEkmanU.c -lsmacof
	clang -o runEkmanW -O2 runEkmanW.c -lsmacof

timer: 
	/usr/bin/time -ph runEkmanU
	/usr/bin/time -ph runEkmanW


clean:
	rm -f runEkmanU runEkmanW *.o