/*-----------------------  Function Prototypes  ------------------------------*/
//  Calculate the director field at the current timestep
	void calc_dir_field(double* x, double* y, double* z,
				   		double* ex, double* ey, double* ez,
				   		double* x_dir, double* y_dir, double* z_dir,
				   		double* ex_dir, double* ey_dir, double* ez_dir,
				   		double* eigenval, double q[][3][4], int avg);

//  Calculate the total energy
	void calc_E();

//  Calculate the temperature of the system
	void calc_temp();

//  Checks if vectors are unit length
	int check_unit_length(double* x, double* y, double* z);

//  Generate a random double between dMin and dMax
	double dRand(double dMin, double dMax);

//  Returns the length of a vector
	double mag_vec(double x, double y, double z);

//  Set information = NAN if overlap with sphere
	void mark_particles(double* x, double* y, double* z,
						double* vx, double* vy, double* vz,
						double* ex, double* ey, double* ez,
						double* ux, double* uy, double* uz,
						double* fx, double* fy, double* fz,
						double* gx, double* gy, double* gz);

//  Newton Raphson method for finding roots
	double newton_raphson(double ex, double ey, double ez,
						  double ux, double uy, double uz, int i);

//  Print energies
	void print_energies();

// 	Print global variables
	void print_global_variables();

//  Print temperature
	void print_temp();

//  Function that returns the order parameter in the x
	double return_sopx(double* ex);

//  Function that returns the temperature
	double return_temp();

//  Reinitialize arrays to shorter length ... after marking NAN
	void resize(double* x, double* y, double* z,
				 double* vx, double* vy, double* vz,
				 double* ex, double* ey, double* ez,
				 double* ux, double* uy, double* uz,
				 double* fx, double* fy, double* fz,
				 double* gx, double* gy, double* gz);

//  Rescale velocities to temp_init
	void rescale(double* x, double* y, double* z,
				double* vx, double* vy, double* vz,
				double* ex, double* ey, double* ez,
				double* ux, double* uy, double* uz);

//  Write energies
	void write_energies(int time);

//  Write scalar order parameter
	void write_sop(double* ex, double* ey, double* ez, int time);

//  Write data for the pair correlation function
	void write_pcf(double* x, double* y, double* z, double histo[][2], int avg);

//  Write data for orientation correlation function
	void write_ocf(double* x, double* y, double* z, 
				   double* ex, double* ey, double* ez, 
				   double histo2[][2], int avg);

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
