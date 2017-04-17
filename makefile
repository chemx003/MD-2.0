main: main.c auxillary_functions.c
	make clean-data
	gcc -o c-main main.c auxillary_functions.c -llapack -lm -O3

postp: postp.c
	gcc -o c-postp postp.c auxillary_functions.c -llapack -lm

clean:
	make clean-data clean-exe
clean-data:
	rm *.dat -f
clean-exe:
	rm c-main -f
	rm c-postp -f
