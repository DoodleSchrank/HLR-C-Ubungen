%: %.o
	gcc -o $* $^.o
	
%.o: %.c
	gcc -c $^.c


.PHONY : clean

clean:
	rm -f *.o
	rm -f *~
