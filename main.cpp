#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <time.h>

using namespace std;

/*-----------------------  Function Declarations  ----------------------------*/
//  Calculate the total energy
	void calc_E();

//  Calculate the temperature of the system
	void calc_temp();

//  Checks if vectors are unit length
	int check_unit_length(double* x, double* y, double* z);

//  Generate a random double between dMin and dMax
	double dRand(double dMin, double dMax);

//  Gay-Berne: Calulate the forces and torques
	void gb		(double* x, double* y, double* z,
				double* ex, double* ey, double* ez);	

//  Initialize the simulation
	void init	(double* x, double* y, double* z,
			  	double* xOld, double* yOld, double* zOld,
			  	double* ex, double* ey, double* ez,
			  	double* exOld, double* eyOld, double* ezOld);

//  Calculate pair correlation function
	void pc		(double* x, double* y, double* z);

//  Print energies
	void print_energies();

// 	Print global variables
	void print_global_variables();

//  Print temperature
	void print_temp();

//  Numerically integrate Newton's Equations
	void verlet (double* x, double* y, double* z,
			  	double* xOld, double* yOld, double* zOld,
			  	double* ex, double* ey, double* ez,
			  	double* exOld, double* eyOld, double* ezOld);

//  Write the positions and orientations to a file
	void write_vectors(double* x, double* y, double* z,
					   double* ex, double* ey, double* ez);

//  These wil eventually get put into a header file
/*----------------------------------------------------------------------------*/


/*--------------------------  Global Variables  ------------------------------*/
//  Simulation Parameters
int 	N			= 25;			//  Number of particles

double 	num_steps 	= 10, 		//  Number of timesteps
	   	dt 			= 0.0015, 		//  Length of time step
	   	temp_init 	= 1.5,			//  Initial temperature

	   	L			= 7.0,			//  Length of simulation box

	   	M			= 1.0,			//	Particle mass
	  	I			= 1.0,			//  Particle moment of inertia

		KB			= 1.0;			//  Boltzmann Constant

//  Data
double	K, V, E, 					//  Pot, kin, tot energies
		P,							//  Pressure
		T;							//  Temperature

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
	srand(time(NULL));

	print_global_variables();

	init(x, y, z, xOld, yOld, zOld,
		 ex, ey, ez, exOld, eyOld, ezOld);			//  Initialize

	calc_E();
	calc_temp();

	print_energies();
	print_temp();

	write_vectors(x, y, z, ex, ey, ez);

	if(check_unit_length(ex, ey, ez) == -1){
		cout << "ORIENTATION VECTORS DO NOT HAVE UNIT LENGTH" << endl << endl;
	}

	if(check_unit_length(exOld, eyOld, ezOld) == -1){
		cout << "OLD ORIENTATION VECTORS DO NOT HAVE UNIT LENGTH" <<endl <<endl;
	}

	for(int i = 0; i < num_steps; i++) {
		/*gb(); 		//  Calculate the forces and torques
		verlet(); 	//  Integrate the eqns of motion	*/
	}
	
	//  Analysis & Post-Processing
	//pc();			//  Calculate PCF*/
}

/*  Initialize particle positions, orientations, velocities, and angular 
 *  velocities 																  */
void init	(double* x, double* y, double* z,
			 double* xOld, double* yOld, double* zOld,
			 double* ex, double* ey, double* ez,
			 double* exOld, double* eyOld, double* ezOld){

	double 	sumVx, sumVy, sumVz,			//  Set lin mtm = 0

			sumVx2, sumVy2, sumVz2, 		//  Set kinetic energy
			sumUx2, sumUy2, sumUz2,

			sfvx, sfvy, sfvz,				//  Scaling factor
			sfux, sfuy, sfuz,

			vx[N], vy[N], vz[N],			//  velocities
			ux[N], uy[N], uz[N],

			a,								//  Particle spacing

			mag;							//  Magnitude of orientation

	int		NUM_LINE,						//  Max amt of particles in a direc.
			p;								//  Particles placed

	/*  Initial setting of postions, orientations, velocities and angular 
	 *  velocities*/
	sumVx = sumVy = sumVz = 0.0;
	sumVx2 = sumVy2 = sumVz2 = 0.0;
	sumUx2 = sumUy2 = sumUz2 = 0.0;

	p=0; K=0;
	
	//  Max number of particles along a side of a cube length l with spacing a
	NUM_LINE = ceil(pow(N, 1.0 / 3.0));
	a = L / NUM_LINE;

	for(int i = 0; i < NUM_LINE; i++) {
		for(int j = 0; j < NUM_LINE; j++) {
			for(int k = 0; k< NUM_LINE; k++) {
				if(p<N){
					//  Assign lattice sites to particle
					x[p] = (i + 0.5 + dRand(-0.1, 0.1)) * a;
					y[p] = (j + 0.5 + dRand(-0.1, 0.1)) * a;
					z[p] = (k + 0.5 + dRand(-0.1, 0.1)) * a;

					//  Assign random orientations
					ex[p] = dRand(-1.0, 1.0);
					ey[p] = dRand(-1.0, 1.0);
					ez[p] = dRand(-1.0, 1.0); 

					mag = sqrt(ex[p]*ex[p] + ey[p]*ey[p] + ez[p]*ez[p]);

					ex[p] = ex[p] / mag;
					ey[p] = ey[p] / mag;
					ez[p] = ez[p] / mag;

					//  Assign random velocities and ang velocites
					vx[p] = dRand(-0.5, 0.5);
					vy[p] = dRand(-0.5, 0.5);
					vz[p] = dRand(-0.5, 0.5);

					ux[p] = dRand(-0.5, 0.5);
					uy[p] = dRand(-0.5, 0.5);
					uz[p] = dRand(-0.5, 0.5);

					//  Sum for corrections to energy and mtm
					sumVx = sumVx + vx[p];
					sumVy = sumVy + vy[p]; 
					sumVz = sumVz + vz[p];

					sumVx2 = sumVx2 + vx[p] * vx[p];
					sumVy2 = sumVy2 + vy[p] * vy[p];
					sumVz2 = sumVz2 + vz[p] * vz[p];

					sumUx2 = sumUx2 + ux[p] * ux[p];
					sumUy2 = sumUy2 + uy[p] * uy[p];
					sumUz2 = sumUz2 + uz[p] * uz[p];
				}

				//  Add to particle count
				p++;
			}
		}
	}

	//  Current CM velocity
	sumVx = sumVx / N; sumVy = sumVy / N; sumVz = sumVz / N;

	//  Current mean squared velocities
	sumVx2 = sumVx2 / N; sumVy2 = sumVy2 / N; sumVz2 = sumVz2 / N;
	sumUx2 = sumUx2 / N; sumUy2 = sumUy2 / N; sumUz2 = sumUz2 / N;

	//  Calculate scaling factors
	sfvx = sqrt(KB * temp_init / sumVx2);
	sfvy = sqrt(KB * temp_init / sumVy2);
	sfvz = sqrt(KB * temp_init / sumVz2);

	sfux = sqrt(KB * temp_init / sumUx2);
	sfuy = sqrt(KB * temp_init / sumUy2);
	sfuz = sqrt(KB * temp_init / sumUz2);

	//  Correct velocities and ang velocities
	for(int i = 0; i < N; i++) {
		//  Set cm mtm to 0
		vx[i] = vx[i] - sumVx;
		vy[i] = vy[i] - sumVy;
		vz[i] = vz[i] - sumVz;

		//  Scale velocities and ang velocites
		vx[i] = vx[i] * sfvx;
		vy[i] = vy[i] * sfvy;
		vz[i] = vz[i] * sfvz;

		ux[i] = ux[i] * sfux;
		uy[i] = uy[i] * sfuy;
		uz[i] = uz[i] * sfuz;

		//  Set old positions and orientations for verlet algo
		xOld[i] = x[i] - dt * vx[i];
		yOld[i] = y[i] - dt * vy[i];
		zOld[i] = z[i] - dt * vz[i];
		/*I don't think this is going to be unit length but correcting to 
		 * unit length might change the temperature/energy - TEST*/
		exOld[i] = ex[i] - dt * ux[i];
		eyOld[i] = ey[i] - dt * uy[i];
		ezOld[i] = ez[i] - dt * uz[i];

		//  Calculate the kinetic energy
		K = K + 0.5 * M * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])
			+ 0.5 * I * (ux[i]*ux[i] + uy[i]*uy[i] +vz[i]*vz[i]);
	}
}

//  Calculate the total energy
void calc_E(){
	E = K + V;
}

//  Calculate the temperature of the system
void calc_temp(){
	T = (2 * K) / (5 * N * KB);
}

/*  Check if vectors are unit length
 *  --->	returns -1 if not		*/
int check_unit_length(double* x, double* y, double* z){
	double mag;
	int unit = 1;

	for(int i = 0; i < N; i++) {
		mag = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];

		if(abs(mag - 1.0) > 0.001){
			unit = -1;
		}
	}

	return unit;
}

/*  Generate a random double between dMin and dMax							  */
double dRand(double dMin, double dMax){
	double d = (double) rand() / RAND_MAX;
	return dMin + d * (dMax - dMin);
}

void print_energies(){
	cout << "ENERGIES" << endl << endl;

	cout << "K = " << K << endl;
	cout << "V = " << V << endl;
	cout << "E = " << E << endl;

	cout << endl;
}

void print_global_variables(){
	cout << "GLOBAL VARIABLES" << endl << endl;
	
	cout << "N = " << N << endl << endl;

	cout << "num_steps = " << num_steps << endl;
	cout << "dt = " << dt << endl;
	cout << "temp_init = " << temp_init << endl << endl;

	cout << "L = " << L << endl << endl;

	cout << "M = " << M << endl;
	cout << "I = " << I << endl << endl;

	cout << "KB = " << KB << endl << endl;

	cout << "K = " << K << endl;
	cout << "V = " << V << endl;
	cout << "E = " << E << endl;
	cout << "P = " << P << endl;
	cout << "T = " << T << endl;

	cout << endl;
}

void print_temp(){
	cout << "TEMPERATURE" << endl << endl;

	cout << "T = " << T << endl;

	cout << endl;
}

//  Write the positions and orientations to a file for display puposes
void write_vectors(double* x, double* y, double* z,
				   double* ex, double* ey, double* ez){
	ofstream o;
	o.open("vector.dat", ios::app);
	for(int i = 0; i < N; i++){
		o << x[i] << "\t" << y[i] << "\t" << z[i] << "\t"
		  << ex[i] << "\t" << ey[i] << "\t" << ez[i] << endl;
	}
	//  Need extra lines for gnuplot to recognize blocks	
	o << endl << endl << endl;
}
