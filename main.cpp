#include <iostream>

using namespace std;

/*-----------------------  Function Declarations  ----------------------------*/
//  Gay-Berne: Calulate the forces and torques
	void gb		(double x[], double y[], double z[],
				double ex[], double ey[], double ez[]);	

//  Initialize the simulation
	void init	(double* x, double* y, double* z,
			  	double* xOld, double* yOld, double* zOld,
			  	double* ex, double* ey, double* ez,
			  	double* exOld, double* eyOld, double* ezOld);

//  Calculate pair correlation function
	void pc		(double x[], double y[], double z[]);

// 	Print global variables
	void print_global_variables();

//  Numerically integrate Newton's Equations
	void verlet (double x[], double y[], double z[],
			  	double xOld[], double yOld[], double zOld[],
			  	double ex[], double ey[], double ez[],
			  	double exOld, double eyOld[], double ezOld[]);

//  These wil eventually get put into a header file
/*----------------------------------------------------------------------------*/


/*--------------------------  Global Variables  ------------------------------*/
//  Simulation Parameters
int 	N			= 25;			//  Number of particles

double 	num_steps 	= 10000, 		//  Number of timesteps
	   	dt 			= 0.0015, 		//  Length of time step
	   	temp_init 	= 1.5,			//  Initial temperature

	   	L			= 7.0,			//  Length of simulation box

	   	M			= 1.0,			//	Particle mass
	  	I			= 1.0;			//  Particle moment of inertia

//  Data
double	K, V, E, 					//Pot, kin, tot energies
		P;							//Pressure

/*----------------------------------------------------------------------------*/




//  Main function
int main(){
	
	//  Local variables
	double 	x[N], y[N], z[N], 				//  Particle coords at n
			xOld[N], yOld[N], zOld[N],		//	Particle coords at n-1

			ex[N], ey[N], ez[N],			//  Particle orient at n
			exOld[N], eyOld[N], ezOld[N],	//  Particle orient at n-1

			fx[N], fy[N], fz[N],			//  Forces
			gx[N], gy[N], gz[N];			//  Gorques

	//  Iteration
	print_global_variables();
	init(x, y, z, xOld, yOld, zOld,
		 ex, ey, ez, exOld, eyOld, ezOld);			//  Initialize

	/*for(){
		gb(); 		//  Calculate the forces and torques
		verlet(); 	//  Integrate the eqns of motion
	}
	
	//  Analysis & Post-Processing
	pc();			//  Calculate PCF*/
}

/*  Initialize particle positions, orientations, velocities, and angular 
 *  velocities 																  */
void init	(double* x, double* y, double* z,
			 double* xOld, double* yOld, double* zOld,
			 double* ex, double* ey, double* ez,
			 double* exOld, double* eyOld, double* ezOld){

	
}

void print_global_variables(){
	cout << "GLOBAL VARIABLES" << endl << endl;
	
	cout << "N = " << N << endl << endl;

	cout << "num_steps = " << num_steps << endl;
	cout << "dt = " << dt << endl;
	cout << "temp_init = " << temp_init << endl << endl;

	cout << "L = " << L << endl << endl;

	cout << "M = " << M << endl;
	cout << "I = " << I << endl;

	cout << "K = " << K << endl;
	cout << "V = " << V << endl;
	cout << "E = " << E << endl;
	cout << "P = " << P << endl;
}
