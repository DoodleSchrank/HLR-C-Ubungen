RM= rm -rf
CC = mpicc
CFLAGS = -Wpedantic -Wall -Wextra -ggdb -Og
INCLUDES =
LFLAGS =
LIBS = -lm

SRCS := timempi.c timempi2.c 
OBJS = $(SRCS:.c=.o)
MAIN = $(SRCS:.c=.x)

all: $(MAIN)
	@echo  Build ${MAIN}

$(MAIN):%.x: %.o
	@echo Building $@
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) $(LIBS) -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

.PHONY: clean cleanall

clean:
	$(RM) $(OBJS)
	$(RM) $(MAIN)
