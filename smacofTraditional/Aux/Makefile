sources = riddle.c \
	../Common/smacofCommonPrintRead.c \
	../Common/smacofCommonRCTranslation.c 

objects = riddle.o \
	smacofCommonPrintRead.o \
	smacofCommonRCTranslation.o 


includes = smacof.h

all: riddle

ofiles: $(includes) $(sources)
	$(CC) -c -O2 -Wall -Wextra -pedantic -std=c99 $(sources)

debug: clean $(includes) $(sources)
	$(CC) -c -g3 -Wall -Wextra -pedantic -std=c99 $(sources)
	$(CC) -o riddled -g3  -Wall -Wextra -pedantic -std=c99 $(objects)

smacof: ofiles
	$(CC) -o riddle -O2  -Wall -Wextra -pedantic -std=c99 $(objects)

format: $(includes) $(sources)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes)

test: smacof
	./riddle > riddle.txt
	
clean:
	rm -f $(objects)

pristine: clean
	rm -f riddle riddle.txt