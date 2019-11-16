# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb -fopenmp
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff-openmp partdiff-openmp-zeilen partdiff-openmp-spalten partdiff-openmp-element
OBJS = askparams.o partdiff-openmp-zeilen.o partdiff-openmp-spalten.o partdiff-openmp-element.o

# Targets ...
all: partdiff-openmp partdiff-openmp-zeilen partdiff-openmp-spalten partdiff-openmp-element

partdiff: partdiff.o askparams.o
	$(CC) $(LFLAGS) -o $@ partdiff.o askparams.o $(LIBS)

partdiff-openmp: partdiff-openmp-zeilen.o askparams.o
	$(CC) $(LFLAGS) -o $@ partdiff-openmp-zeilen.o askparams.o $(LIBS)

partdiff-openmp-%: partdiff-openmp-%.o askparams.o
	$(CC) $(LFLAGS) -o $@ partdiff-openmp-$*.o askparams.o $(LIBS)

partdiff.o: partdiff.c
	$(CC) -c $(CFLAGS) -D NORMAL $< -o $@
partdiff-openmp-zeilen.o: partdiff.c
	$(CC) -c $(CFLAGS) -D ROW $< -o $@
partdiff-openmp-spalten.o: partdiff.c
	$(CC) -c $(CFLAGS) -D COLUMN $< -o $@
partdiff-openmp-element.o: partdiff.c
	$(CC) -c $(CFLAGS) -D SINGLE $< -o $@



askparams.o: askparams.c

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)
