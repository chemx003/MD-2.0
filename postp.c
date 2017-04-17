#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "auxillary_functions.h"

/*-----------------------  Function Declarations  ----------------------------*/

/*----------------------------------------------------------------------------*/

/*--------------------------  Global Variables  ------------------------------*/
//  Simulation Parameters
int 	N				= 100,			//  Number of particles
		pcf_bins		= 400,			//  Number of bins for pcf
		pcf_num_steps	= 10,			// 	Steps to avg pcf over
		num_bin_x 		= 6,			//  Director bins
		num_bin_y		= 6,			
		num_bin_z		= 6,
		num_steps 		= 100, 			//  Number of timesteps
		num_steps_eqbm	= 100;			//  Number of eqbm timesteps

double 	dt 				= 0.0015, 		//  Length of time step
	   	temp_init 		= 0.8,			//  Initial temperature
	   	xi = 0, eta = 0,				//  Thermostat variables

	   	L				= 48.5,//63.1,			//  Length of simulation box
	   	SL				= 16.5,//21.1,			//	Short length of the simulation box

	   	M				= 1.0,			//	Particle mass
	  	I				= 1.0,			//  Particle moment of inertia

	  	R				= 3.0,			//  Immersed sphere radius
	  	W				= 350000,		//  Anchoring coefficient

		KB				= 1.0,			//  Boltzmann Constant
		PI				= 3.14159265358979; //  Pi

//  Data
double	KT, KR, K, V, E, 				//  Pot, kin, tot energies
		P,								//  Pressure
		T;								//  Temperature


/*----------------------------------------------------------------------------*/

//  Main function
int main(){

	//  Local variables
	double 	x[N], y[N], z[N],
			ex[N], ey[N], ez[N],
			x_dir[num_bin_x*num_bin_y*num_bin_z],
			y_dir[num_bin_x*num_bin_y*num_bin_z],
			z_dir[num_bin_x*num_bin_y*num_bin_z],
			ex_dir[num_bin_x*num_bin_y*num_bin_z],
			ey_dir[num_bin_x*num_bin_y*num_bin_z],
			ez_dir[num_bin_x*num_bin_y*num_bin_z];

	double q[num_bin_x*num_bin_y*num_bin_z][3][4];
	double eigenval[num_bin_x*num_bin_y*num_bin_z][4];

	int step;

	//  Read the vector.dat file
	FILE* r;
	r = fopen("./vector.dat", "r");

	//  Write to another file for test
	FILE* o;
	o = fopen("./test_vector.dat", "a");


	while(feof(r)==0) {

		if(feof(r) == 0) {
			fscanf(r, "%i\t%i\n", &step, &N);
		}

		//  Read data in from vector.dat
		for(int i=0; i < N; i++) {
			if(feof(r) == 0){
				fscanf(r,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &x[i], &y[i], &z[i], 
						&ex[i], &ey[i], &ez[i]); 
			}
		}

		/*  Functions to operate on x y z ex ey ez ... maybe just print
		 *  after the sphere is removed and use that N? */
		if(step>=1000){
			calc_dir_field(x, y, z, ex, ey, ez, x_dir, y_dir, z_dir,
						ex_dir, ey_dir, ez_dir, eigenval, q, 0);
		}

		/*  This writes out the animation.dat file which should not have the
		 *  step number and number of particles at the top of each block */
		for(int i = 0; i < N; i++){
				fprintf(o, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x[i], y[i], z[i], 
						ex[i], ey[i], ez[i]); 
		}

			fprintf(o, "\n\n");
	}

	calc_dir_field(x, y, z, ex, ey, ez, x_dir, y_dir, z_dir,
					ex_dir, ey_dir, ez_dir, eigenval, q, 1);

	//  Close files
	fclose(r);
	fclose(o);
}
