CC=gcc
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb -fopenmp
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff partdiff_row partdiff_column partdiff_single
OBJS = partdiff.o partdiff_row.o partdiff_column.o partdiff_single.o askparams.o


all: askparams.o partdiff partdiff_row partdiff_column partdiff_single




partdiff: partdiff.o
	gcc $(CFLAGS) partdiff.o askparams.o $(LIBS)
partdiff_row: partdiff_row.o
	gcc $(CFLAGS) partdiff_row.o askparams.o $(LIBS)
partdiff_column: partdiff_column.o
	gcc $(CFLAGS) partdiff_column.o askparams.o $(LIBS)
partdiff_single: partdiff_single.o
	gcc $(CFLAGS) partdiff_single.o askparams.o $(LIBS)


askparams.o: askparams.c
	gcc $(CFLAGS) -c askparams.c
partdiff.o: partdiff.c
	gcc -D NORMAL $(CFLAGS) -c partdiff.c
partdiff_row.o: partdiff.c
	gcc -D ROW $(CFLAGS) -c partdiff.c -o partdiff_row.o
partdiff_column.o: partdiff.c
	gcc -D COLUMN $(CFLAGS) -c partdiff.c -o partdiff_column.o
partdiff_single.o: partdiff.c
	gcc -D SINGLE $(CFLAGS) -c partdiff.c -o partdiff_single.o

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)

