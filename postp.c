#include <stdio.h>

/*-----------------------  Function Declarations  ----------------------------*/

/*----------------------------------------------------------------------------*/

/*--------------------------  Global Variables  ------------------------------*/
//  Simulation Parameters
int		N = 100;

/*----------------------------------------------------------------------------*/

//  Main function
int main(){

	//  Local variables
	double 	x[N], y[N], z[N],
			ex[N], ey[N], ez[N];

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

		/*  This writes out the animation.dat file which should not have the
		 *  step number and number of particles at the top of each block */
		for(int i = 0; i < N; i++){
				fprintf(o, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%i\n", x[i], y[i], z[i], 
						ex[i], ey[i], ez[i]); 
		}

			fprintf(o, "\n\n");
	}

	//  Close files
	fclose(r);
	fclose(o);
}
