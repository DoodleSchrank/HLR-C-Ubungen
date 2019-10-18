all: map.bin pointer.bin

%.bin: %.o
	gcc -o $@ $<
	
%.o: %.c
	gcc -c $<

clean:
	rm -f *.bin
	rm -f *.o
	rm -f *~
