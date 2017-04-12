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

	//  Read the vector.dat file
	FILE* r;
	r = fopen("./vector.dat", "r");

	//  Write to another file for test
	FILE* o;
	o = fopen("./test_vector.dat", "a");


	while(feof(r)==0) {

		//  Read data in from vector.dat
		for(int i=0; i < N; i++) {
			if(feof(r) == 0){
				fscanf(r,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &x[i], &y[i], &z[i], 
						&ex[i], &ey[i], &ez[i]); 
			}
		}

		//  Write data out to test_vector.dat
		for(int i = 0; i < N; i++){
				fprintf(o, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%i\n", x[i], y[i], z[i], 
						ex[i], ey[i], ez[i], i); 
		}

			fprintf(o, "\n\n");
	}

	//  Close files
	fclose(r);
	fclose(o);
}
