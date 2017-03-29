main: main.c auxillary_functions.c
	make clean-data
	gcc -o c-main main.c auxillary_functions.c -llapack -lm -O3
clean:
	make clean-data clean-exe
clean-data:
	rm *.dat -f
clean-exe:
	rm c-main -f
