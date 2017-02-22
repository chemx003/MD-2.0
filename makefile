main: main.c
	gcc -o c-main main.c auxillary_functions.c -lm -O3
	make clean

clean:
	rm *.dat
