CC=gcc
CFLAGS=-O3 -Wall

OBJECTS= ./hello.o ./tinystr.o ./tinyxml.o ./tinyxmlerror.o ./tinyxmlparser.o 

.PHONY: clean

all: a

lop: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o a

clean:
	rm -f ./*~ ./*.o a

