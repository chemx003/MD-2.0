main: main.c auxillary_functions.c
	gcc -o c-main main.c auxillary_functions.c -lm 
	make clean

clean:
	rm *.dat
