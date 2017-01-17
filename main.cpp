#include <iostream>

using namespace std;

/*-----------------------  Function Declarations  ----------------------------*/
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

// 	Print global variables
	void print_global_variables();

//  Numerically integrate Newton's Equations
	void verlet (double* x, double* y, double* z,
			  	double* xOld, double* yOld, double* zOld,
			  	double* ex, double* ey, double* ez,
			  	double* exOld, double* eyOld, double* ezOld);

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

	/*for() {
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

	double 	sumVx = sumVy = sumVz = 0.0,	//  Set lin mtm = 0

			sumVx2 = sumVy2 = sumVz2 = 0.0, //  Set kinetic energy
			sumUx2 = sumUy2 = sumUz2 = 0.0,

			sfvx, sfvy, sfvz,				//  Scaling factor
			sfux, sfuy, sfuz,

			vx[N], vy[N], vz[N],			//  velocities
			ux[N], uy[N], uz[N],

			mag;							//  Magnitude of orientation

	int		NUM_LINE,						//  Max amt of particles in a direc.
			p,								//  Particles placed
			a;								//  Particle spacing

	/*  Initial setting of postions, orientations, velocities and angular 
	 *  velocities															  */
	p=0;

	for(int i = 0; i < NUM_LINE; i++) {
		for(int j = 0; j < NUM_LINE; j++) {
			for(int k = 0; k< NUM_LINE; k++) {
				if(p<N){
					//  Assign lattice sites to particle
					x[p] = (i + 0.5 + dRand(-0.1, 0.1)) * a;
					y[p] = (j + 0.5 + dRand(-0.1, 0.1)) * a;
					z[p] = (k + 0.5 + dRand(-0.1, 0.1)) * a;

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
	sfvx = sqrt(kB * temp_init / sumVx2);
	sfvy = sqrt(kB * temp_init / sumVy2);
	sfvz = sqrt(kB * temp_init / sumVz2);

	sfux = sqrt(kB * temp_init / sumUx2);
	sfuy = sqrt(kB * temp_init / sumUy2);
	sfuz = sqrt(kB * temp_init / sumUz2);

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
	}
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
