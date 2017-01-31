/*-----------------------  Function Prototypes  ------------------------------*/
//  Calculate the total energy
	void calc_E();

//  Calculate the temperature of the system
	void calc_temp();

//  Checks if vectors are unit length
	int check_unit_length(double* x, double* y, double* z);

//  Generate a random double between dMin and dMax
	double dRand(double dMin, double dMax);

//  Print energies
	void print_energies();

// 	Print global variables
	void print_global_variables();

//  Print temperature
	void print_temp();

//  Write energies
	void write_energies(int time);

//  Write scalar order parameter
	void write_sop(double* ex, double* ey, double* ez, int time);

//  Write data for the pair correlation function
	void write_pcf(double* x, double* y, double* z, double histo[][2], int avg);

//  Write a single point to a file
	void write_point(double x, double y, double z,
					double ex, double ey, double ez);

//  Wrtie the temperature to a file
	void write_temp(double time);

//  Write the positions and orientations to a file
	void write_vectors(double* x, double* y, double* z,
					   double* ex, double* ey, double* ez);

//  Same as write vectors but with marker for color
	void write_vectors_colored(double* x, double* y, double* z,
							   double* ex, double* ey, double* ez,
							   double* color);
/*----------------------------------------------------------------------------*/
