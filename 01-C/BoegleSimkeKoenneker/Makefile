all: map pointer

%: %.o
	gcc -o $^ $@
	
%.o: %.c
	gcc -c $<

clean:
	rm -f *.bin
	rm -f *.o
	rm -f *~
