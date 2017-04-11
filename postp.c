#include <stdio.h>

/*-----------------------  Function Declarations  ----------------------------*/

/*----------------------------------------------------------------------------*/

/*--------------------------  Global Variables  ------------------------------*/
//  Simulation Parameters
int		N = 4096;

/*----------------------------------------------------------------------------*/

//  Main function
int main(){

	//  Local variables
	double 	x[N], y[N], z[N],
			ex[N], ey[N], ez[N];

	//  Read the vector.dat file
	FILE* r;
	r = fopen("./vector.dat", "r");

	for(int i=0; i < N; i++){
		if(feof(r) == 0){
			fscanf(r,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &x[i], &y[i], &z[i], 
					&ex[i], &ey[i], &ez[i]); 
		}
	}

	fclose(r);

	//  Write to another file for test
	FILE* o;
	o = fopen("./test_vector.dat", "a");

	for(int i = 0; i < N; i++){
		fprintf(o, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x[i], y[i], z[i], 
				ex[i], ey[i], ez[i]); 
	}

	fclose(o);
}
