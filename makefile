main: main.c auxillary_functions.c
	gcc -o c-main main.c auxillary_functions.c -lm -O3
	make clean

clean:
	rm *.dat 2>/dev/null
