main: main.cpp
	g++ -o c-main main.cpp auxillary_functions.cpp 
	make clean

clean:
	rm *.dat
