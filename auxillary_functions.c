#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "auxillary_functions.h"

/*--------------------------  Global Variables  ------------------------------*/
//  Simulation Parameters
extern int 	N,								//  Number of particles
	   		pcf_bins,						//  Number of bins for pcf
	   		pcf_num_steps,					//  Steps to avg pcf over
	   		num_bin_x,						//  Director bins
	   		num_bin_y,
	   		num_bin_z;

extern double 	num_steps, 					//  Number of timesteps
				dt, 						//  Length of time step
				temp_init,					//  Initial temperature
				xi, eta,					//  Thermostat variables

				L, SL,						//  Length of simulation box

				M,							//	Particle mass
				I,							//  Particle moment of inertia
				R,							//  Sphere radius
				W,							//  Anchoring coefficient

				KB,							//  Boltzmann Constant
				PI;

//  Data
extern double	KT, KR, K, V, E, 			//  Pot, kin, tot energies
				P,							//  Pressure
				T;							//  Temperature

//  LAPACK FUNCTION TO DIAGONALIZE SYMMETRIC MATRIX
extern void dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA, 
				   double* W, double* WORK, int* LWORK, int* INFO);
/*----------------------------------------------------------------------------*/
//  Calculate the director field at the current timestep
void calc_dir_field(double* x, double* y, double* z,
				   double* ex, double* ey, double* ez,
				   double* x_dir, double* y_dir, double* z_dir,
				   double* ex_dir, double* ey_dir, double* ez_dir,
				   double* eigenval, double q[][3][4], int avg){
	/*  NEW GLOBAL VARIABLES
	 *  int		num_bin_x,
	 *  		num_bin_y,
	 *  		num_bin_z; */

	double 	len_bin_x, len_bin_y, len_bin_z, //  bin dimensions
			bin_x, bin_y, bin_z, //  center of current bin 
			dx, dy, dz;  //Distance comp of bin center to particle

	double //q[num_bin_x*num_bin_y*num_bin_z][3][4], 
		   qt[9], w[3], work[9]; //  q tensor and q tensor in form for lapack

	char jobz, uplo; //  See documentation for dysev

	int order, lda, lwork, info, bin_number,
		n_bin;	//  Number of particles in bin

	//  Initialize params 
	jobz = 'V';
	uplo = 'U';
	order = 3;
	lda = 3;
	for(int i = 0; i < 3; i++) {w[i] = 0.0;}
	for(int i = 0; i < 9; i++) {work[i] = 0.0;}
	lwork = 9;
	info = 0;
	n_bin = 0;
	bin_number = 0;

	//  Calculate bin lengths
	len_bin_x = L / num_bin_x;
	len_bin_y = SL / num_bin_y;
	len_bin_z = SL / num_bin_z;
	
	//  Loop over all bins
	for(int i = 0; i < num_bin_x; i++) {
		for(int j = 0; j < num_bin_y; j++) {
			for(int k = 0; k < num_bin_z; k++) {
				
				if(avg == 0){
					//  Center of the bin
					bin_x = len_bin_x*(0.5 + i);
					bin_y = len_bin_y*(0.5 + j);
					bin_z = len_bin_z*(0.5 + k);

					//  Store center of bin coords
					x_dir[bin_number] = bin_x;
					y_dir[bin_number] = bin_y;
					z_dir[bin_number] = bin_z;

					/*  For all bins check all particles to see
					 *  if they are in that bin */
					for(int  p = 0; p < N; p++) {

						// Calculate distance from bin center
						dx = fabs(x[p] - bin_x);
						dy = fabs(y[p] - bin_y);
						dz = fabs(z[p] - bin_z);
						
						//  Check to see if particle is in bin
						if(dx < len_bin_x/2 && dy < len_bin_y/2 && 
								dz < len_bin_z/2){
							//  Take total number of particles in bin for avg
							q[bin_number][0][3]++;

							//  Sum up Q-Tensor components ... only store upper
							//  later
							q[bin_number][0][0] += ex[p]*ex[p] - 1.0/3.0;
							q[bin_number][0][1] += ex[p]*ey[p];
							q[bin_number][0][2] += ex[p]*ez[p];

							q[bin_number][1][0] += ey[p]*ex[p];
							q[bin_number][1][1] += ey[p]*ey[p] - 1.0/3.0;
							q[bin_number][1][2] += ey[p]*ez[p];

							q[bin_number][2][0] += ez[p]*ex[p];
							q[bin_number][2][1] += ez[p]*ey[p];
							q[bin_number][2][2] += ez[p]*ez[p] - 1.0/3.0;
						}
					}
				}

				if(avg==1){
					//  Average Q-Tensor
					for(int a = 0; a < 3; a++) {
						for(int b = 0; b < 3; b++) {
							q[bin_number][a][b] = q[bin_number][a][b]
												/ q[bin_number][0][3];
						}
					}

					//  Convert Q-Tensor to proper format
					for(int a = 0; a < 3; a++) {
						for(int b = 0; b < 3; b++) {
							qt[b + 3*a] = q[bin_number][b][a];
						}
					}

					//  Call LAPACK
					dsyev_(&jobz, &uplo, &order, qt, &lda, w, work, &lwork, &info);

					//  Store eigenstuff and position of bin?
					eigenval[bin_number] = pow(w[0]*w[1]*w[2], 2)
						/pow(w[0]*w[0]+w[1]*w[1]+w[2]*w[2], 3) - 1.0/54;
					
					//double mag = sqrt(qt[6]*qt[6] + qt[7]*qt[7] + qt[8]*qt[8]);

					ex_dir[bin_number] = qt[6];///mag;
					ey_dir[bin_number] = qt[7];///mag;
					ez_dir[bin_number] = qt[8];///mag;

					ex_dir[bin_number] = ex_dir[bin_number]
						/mag_vec(ex_dir[bin_number], ey_dir[bin_number], ez_dir[bin_number]);
					ey_dir[bin_number] = ey_dir[bin_number]
						/mag_vec(ex_dir[bin_number], ey_dir[bin_number], ez_dir[bin_number]);
					ez_dir[bin_number] = ez_dir[bin_number]
						/mag_vec(ex_dir[bin_number], ey_dir[bin_number], ez_dir[bin_number]);
				}	

				//  Increment bin number
				bin_number++;
			}
		}
	}

	if(avg==1){
		/*  Here we'll want to write out the values to a .dat file maybe make
		 *  director stuff local to this function? */
		FILE *o, *p;
		o = fopen("director.dat", "a");
		p = fopen("biaxial.dat", "a");

		double mid_z = (floor(num_bin_z/2) + 0.5) * len_bin_z;

		printf("mid_z=%f\n\n", mid_z);

		for(int i = 0; i < bin_number; i++){
			if(fabs(z_dir[i] - mid_z) == 0 && isnan(ex_dir[i]) == 0){
				double mag = mag_vec(ex_dir[i], ey_dir[i], ez_dir[i]);
				ex_dir[i] = ex_dir[i]/mag;
				ey_dir[i] = ey_dir[i]/mag;
				ez_dir[i] = ez_dir[i]/mag;
				printf("length before = %f\n", mag);
				mag = mag_vec(ex_dir[i], ey_dir[i], ez_dir[i]);
				printf("length after = %f\n", mag);
				fprintf(o, "%f    %f    %f    %f\n", x_dir[i], y_dir[i], ex_dir[i], ey_dir[i]);
				fprintf(p, "%f    %f    %f\n", x_dir[i], y_dir[i], eigenval[i]);
			}
		}

		//  Need extra lines for gnuplot to recognize blocks	
		fprintf(o, "\n\n");
		fprintf(p, "\n\n");

		fclose(o);
		fclose(p);
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

		if(fabs(mag - 1.0) > 0.001){
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

//Returns the length of a vector
double mag_vec(double x, double y, double z){
	double magn;

	magn = sqrt(x*x + y*y + z*z);

	return magn;
}

//  Set information = NAN if overlap with sphere
void mark_particles(double* x, double* y, double* z,
			  		double* vx, double* vy, double* vz,
			  		double* ex, double* ey, double* ez,
			  		double* ux, double* uy, double* uz,
			  		double* fx, double* fy, double* fz,
					double* gx, double* gy, double* gz){

	for(int i = 0; i < N; i++){
		if(pow(x[i] - L/2.0, 2.0) + pow(y[i] - SL/2.0, 2.0) + 
				pow(z[i] - SL/2.0, 2.0) <= pow(R + 3.0, 2.0)){
			
			x[i] = NAN;
			y[i] = NAN; 
			z[i] = NAN;
			vx[i] = NAN;
			vy[i] = NAN;
			vz[i] = NAN;
			ex[i] = NAN;
			ey[i] = NAN;
			ez[i] = NAN; 
			ux[i] = NAN;
			uy[i] = NAN;
			uz[i] = NAN;
			fx[i] = NAN;
			fy[i] = NAN;
			fz[i] = NAN;
			gx[i] = NAN;
			gy[i] = NAN;
			gz[i] = NAN;
		}
	}
}

//  Newton raphson method for calculating the lagrange multiplier
double newton_raphson(double ex, double ey, double ez,
					  double ux, double uy, double uz, int i){
	double l, lOld, temp;
	double dot1, dot2;
	int count=0;

	lOld = ux*ux + uy*uy + uz*uz;

	dot1 = ex*ux + ey*uy + ez*uz;
	dot2 = ux*ux + uy*uy + uz*uz;
	

	if(isnan(ux)==0 && isnan(uy)==0 && isnan(uz)==0){
		while(fabs(l - lOld) > 0.01 && isnan(l)==0){
			temp = l;
			l = lOld - (lOld*lOld + 2*lOld*(dot1+1/dt) + dot2  + 2*dot1/dt)/
				(2*lOld + 2*dot1 + 1/dt);
			count++;
			lOld = temp; 
			if(count > 1000){
				printf("count = %i\n",count);
				printf("i = %i, l = %f\n",i,l);
			} 
		}
	}

	return l;
}

void print_energies(){
	printf("ENERGIES\n");

	printf("KT = %f\n", KT);
	printf("KR = %f\n", KR);
	printf("K = %f\n", K);
	printf("V = %f\n", V);
	printf("E = %f\n\n", E);
}

void print_global_variables(){
	printf("GLOBAL VARIABLES\n");
	
	printf("N = %i\n", N);

	printf("num_steps = %f\n", num_steps);
	printf("dt = %f\n", dt);
	printf("temp_init = %f\n\n", temp_init);

	printf("L = %f\n ", L);

	printf("M = %f\n", M);
	printf("I = %f\n\n", I);

	printf("KB = %f\n\n", KB);

	printf("KT = %f\n", KT);
	printf("KR = %f\n", KR);
	printf("K = %f\n", K);
	printf("V = %f\n", V);
	printf("E = %f\n", E);
	printf("P = %f\n", P);
	printf("T = %f\n\n", T);
}

void print_temp(){
	printf("TEMPERATURE\n");

	printf("T = %f\n\n",T);
}

//  Reinitialize arrays to shorter length ... after marking NAN
void resize(double* x, double* y, double* z,
			 double* vx, double* vy, double* vz,
			 double* ex, double* ey, double* ez,
			 double* ux, double* uy, double* uz,
			 double* fx, double* fy, double* fz,
			 double* gx, double* gy, double* gz){
	int	count, newN, p; //  Count up the number NANs			
	p=0; count = 0;

	printf("%i\n", N);

	for(int i = 0; i < N; i++) {
		if(isnan(x[i]) != 0){
			count++;
		}
	}

	printf("count = %i\n\n", count);

	//  Number of particles in the new simulation
	newN = N - count;

	//  Declare temporary arrays
	double tx[newN]; // = (double*) calloc(newN,sizeof(double));
	double ty[newN]; // = (double*) calloc(newN,sizeof(double));
	double tz[newN]; // = (double*) calloc(newN,sizeof(double));

	double tvx[newN]; // = (double*) calloc(newN,sizeof(double));
	double tvy[newN]; // = (double*) calloc(newN,sizeof(double));
	double tvz[newN]; // = (double*) calloc(newN,sizeof(double));

	double tex[newN]; //  = (double*) calloc(newN,sizeof(double));
	double tey[newN]; //  = (double*) calloc(newN,sizeof(double));
	double tez[newN]; //  = (double*) calloc(newN,sizeof(double));

	double tux[newN]; //  = (double*) calloc(newN,sizeof(double));
	double tuy[newN]; //  = (double*) calloc(newN,sizeof(double));
	double tuz[newN]; //  = (double*) calloc(newN,sizeof(double));

	double tfx[newN]; //  = (double*) calloc(newN,sizeof(double));
	double tfy[newN]; //  = (double*) calloc(newN,sizeof(double));
	double tfz[newN]; //  = (double*) calloc(newN,sizeof(double));

	double tgx[newN]; //  = (double*) calloc(newN,sizeof(double));
	double tgy[newN]; //  = (double*) calloc(newN,sizeof(double));
	double tgz[newN]; //  = (double*) calloc(newN,sizeof(double));
	
	//  Assign values to the temporary arrays
	for(int i = 0; i < N; i++) {
		if(isnan(x[i]) == 0){
			tx[p] = x[i];
			ty[p] = y[i];
			tz[p] = z[i];

			tvx[p] = vx[i];
			tvy[p] = vy[i];
			tvz[p] = vz[i];

			tex[p] = ex[i];
			tey[p] = ey[i];
			tez[p] = ez[i];

			tux[p] = ux[i];
			tuy[p] = uy[i];
			tuz[p] = uz[i];

			tfx[p] = fx[i];
			tfy[p] = fy[i];
			tfz[p] = fz[i];

			tgx[p] = gx[i];
			tgy[p] = gy[i];
			tgz[p] = gz[i];

			p++;
		}
	}

	//  Assign arrays to temp array pointer
	for(int i = 0; i < N; i++){
		if(i<newN){
			x[i] = tx[i];
			y[i] = ty[i];
			z[i] = tz[i];

			vx[i] = tvx[i];
			vy[i] = tvy[i];
			vz[i] = tvz[i];

			ex[i] = tex[i];
			ey[i] = tey[i];
			ez[i] = tez[i];

			ux[i] = tux[i];
			uy[i] = tuy[i];
			uz[i] = tuz[i];

			fx[i] = tfx[i];
			fy[i] = tfy[i];
			fz[i] = tfz[i];

			gx[i] = tgx[i];
			gy[i] = tgy[i];
			gz[i] = tgz[i];
		}
		else{
			x[i] = 0;
			y[i] = 0;
			z[i] = 0;

			vx[i] = 0;
			vy[i] = 0;
			vz[i] = 0;

			ex[i] = 0;
			ey[i] = 0;
			ez[i] = 0;

			ux[i] = 0;
			uy[i] = 0;
			uz[i] = 0;

			fx[i] = 0;
			fy[i] = 0;
			fz[i] = 0;

			gx[i] = 0;
			gy[i] = 0;
			gz[i] = 0;	
		}
	}

	
	for(int i = 0; i < N; i++){
		printf("x[p] = %f\n", x[i]);
	}

	//  Assign new number of particles
	N = newN;
}

//  Rescale the velocities to temp_init ... move this to auxilary functions 
void rescale(double* x, double* y, double* z,
			 double* vx, double* vy, double* vz,
			 double* ex, double* ey, double* ez,
			 double* ux, double* uy, double* uz){

	double 	sumVx2, sumVy2, sumVz2, 		//  Set kinetic energy
			sumUx2, sumUy2, sumUz2,

			sfvx, sfvy, sfvz,				//  Scaling factor
			sfux, sfuy, sfuz;
	
	//  Set this equal to zero
	sumVx2 = sumVy2 = sumVz2 = 0.0;
	sumUx2 = sumUy2 = sumUz2 = 0.0;

	for(int p = 0; p < N; p++) {
		//  Sum for corrections to energy and mtm
		sumVx2 = sumVx2 + vx[p] * vx[p];
		sumVy2 = sumVy2 + vy[p] * vy[p];
		sumVz2 = sumVz2 + vz[p] * vz[p];

		sumUx2 = sumUx2 + ux[p] * ux[p];
		sumUy2 = sumUy2 + uy[p] * uy[p];
		sumUz2 = sumUz2 + uz[p] * uz[p];
	}

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
		//  Scale velocities and ang velocites
		vx[i] = vx[i] * sfvx;
		vy[i] = vy[i] * sfvy;
		vz[i] = vz[i] * sfvz;

		ux[i] = ux[i] * sfux;
		uy[i] = uy[i] * sfuy;
		uz[i] = uz[i] * sfuz; 

		//  Calculate the kinetic energy
		KT = KT + 0.5 * M * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
		KR = KR + 0.5 * I * (ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i]);
		K = KT + KR;
	}
}

//  Funtion that returns the order parameter in x
double return_sopx(double* ex){
	double sopx = 0;

	for(int i = 0; i < N; i++){
		sopx = sopx + (3.0*ex[i]*ex[i] - 1)/2.0;
	}

	sopx = sopx / N;

	return sopx;
}

//  Function that returns the temperature
double return_temp(){
	double temp;
	calc_temp();

	temp = T;

	return temp;
}

//  Write the energies to a file 'time' refers to the current timestep
void write_energies(int time){
//	ofstream o;
//	o.fopen("energy.dat", ios::app);
	FILE* o;
	o = fopen("energy.dat", "a");

	fprintf(o, "%i\t%f\t%f\t%f\n", time, V, K, E);

	fclose(o);
}

//  Writes the scalar order parameter to a file
void write_sop(double* ex, double* ey, double* ez, int time){
	double sopx, sopy, sopz; //order param wrt axes
	double exi, eyi, ezi, mag;

	sopx = 0.0; sopy = 0.0; sopz = 0.0;

	for(int i = 0; i < N; i++) {

		//wrt x
		sopx = sopx + (3.0*ex[i]*ex[i] - 1)/2.0;
		sopy = sopy + (3.0*ey[i]*ey[i] - 1)/2.0;
		sopz = sopz + (3.0*ez[i]*ez[i] - 1)/2.0;
	}

	sopx = sopx / N; sopy = sopy /N; sopz = sopz / N;

	FILE* o;
	o = fopen("sop.dat", "a");

	if(fabs(sopx)<3.0){
		fprintf(o,"%i\t%f\t%f\t%f\n", time, sopx, sopy, sopz);
	}

	fclose(o);
}

//  Writes orientation pair correlation
void write_ocf(double* x, double* y, double* z,
			   double* ex, double* ey, double* ez,
			   double histo2[][2], int avg){
	double dR = 0.025,		//  Width of bin
		   R = dR,			//  Bin label
		   dx, dy, dz,		//  Components of seperation vector
		   r;				//  Seperation distance

	for(int b = 0; b < pcf_bins; b++) {
		for(int i = 0; i < N-1; i++) {
			for(int j = i+1; j < N; j++) {
				dx = x[i] - x[j];
				dy = y[i] - y[j];
				dz = z[i] - z[j];

				//  Apply minimum image critereon (trying a different method)
				if(fabs(dx) > 0.5*L){
					dx = dx - L*(dx / fabs(dx));
				}

				if(fabs(dy) > 0.5*SL){
					dy = dy - SL*(dy / fabs(dy));
				}
				
				if(fabs(dz) > 0.5*SL){
					dz = dz - SL*(dz / fabs(dz));
				}
				
				r = sqrt(dx*dx + dy*dy + dz*dz);

				if(r < R+dR && r > R){
					double op = (3*pow(ex[i]*ex[j] + ey[i]*ey[j] 
									+ ez[i]*ez[j],2) - 1)/2;
					histo2[b][1] = histo2[b][1] + op;
					histo2[b][0]++;

					//printf("%f\n", op);
				}
			}
		}
		R = R + dR;
	}

	//  Write to file
	if(avg == 1){
		FILE* p;
		p = fopen("ocf.dat","a");

		R = dR;

		for(int b = 0; b < pcf_bins; b++) {
			if(histo2[b][0] != 0){histo2[b][1] = histo2[b][1] / histo2[b][0];}
			histo2[b][0] = R;
			fprintf(p, "%f\t%f\n", histo2[b][0], histo2[b][1]);	
			R = R + dR;
		}
		
		fclose(p);
	}
}

/*  Writes data for pair correlation function to a file
 *  if(avg == 0) store data
 *  if(avg == 1) average data and print 									  */
void write_pcf(double* x, double* y, double* z,
			   double histo[][2], int avg){
	double 	dR = 0.025,		//  Width of bin
			R = dR,			//  Bin label -incremented in loop
			dx, dy, dz,		//  Components of seperation vector
			r;				//  Seperation distance

	for(int b = 0; b < pcf_bins; b++) {
		histo[b][0] = R;

		for(int i = 0; i < N-1; i++) {
			for(int j = i+1; j < N; j++) {
				dx = x[i] - x[j];
				dy = y[i] - y[j];
				dz = z[i] - z[j];

				//  Apply minimum image critereon (trying a different method)
				if(fabs(dx) > 0.5*L){
					dx = dx - L*(dx / fabs(dx));
				}

				if(fabs(dy) > 0.5*SL){
					dy = dy - SL*(dy / fabs(dy));
				}
				
				if(fabs(dz) > 0.5*SL){
					dz = dz - SL*(dz / fabs(dz));
				}
	
				r = sqrt(dx*dx + dy*dy + dz*dz);

				if(r > R && r < R+dR){
					histo[b][1]++;
				}
			}
		}
		R = R + dR; //  Increment bin size
	}
	
	//  Write to file
	if(avg == 1){
		FILE* p;
		p = fopen("pcf.dat", "a");

		R = dR;
		for(int b = 1; b < pcf_bins; b++) {
			histo[b][1] = histo[b][1] / pcf_num_steps;
			histo[b][1] = 2*histo[b][1] / (4*PI*R*R*dR*N*N/(L*SL*SL));
			fprintf(p,"%f\t%f\n", histo[b][0], histo[b][1]);
			R = R + dR;
		}
		
		fclose(p);
	}
}

//  Write one point to a file
void write_point(double x, double y, double z,
				 double ex, double ey, double ez){
	FILE* o;
	o = fopen("points.dat", "a");

	fprintf(o, "%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, ex, ey, ez);

	fclose(o);
}

//  Writes the temperature to a file
void write_temp(double time){
	FILE* o; 
	o = fopen("temp.dat", "a");

	calc_temp();

	fprintf(o,"%f\t%f\n", time, T);

	fclose(o);
}

//  Write the positions and orientations to a file for display puposes
void write_vectors(double* x, double* y, double* z,
				   double* ex, double* ey, double* ez){
	FILE* o;
	o = fopen("vector.dat", "a");

	for(int i = 0; i < N; i++){
		fprintf(o, "%f\t%f\t%f\t%f\t%f\t%f\n", x[i], y[i], z[i], 
				ex[i], ey[i], ez[i]);
	}

	//  Need extra lines for gnuplot to recognize blocks	
	fprintf(o, "\n\n");

	fclose(o);
}

/*  same as write, but extra integer to mark color
void write_vectors_colored(double* x, double* y, double* z,
				   double* ex, double* ey, double* ez,
				   double* color){
	FILE* o;
	o = fopen("vector.dat", "a");

	for(int i = 0; i < N; i++){
		o, x[i],"\t",y[i],"\t",z[i],"\t"
		 ,ex[i],"\t",ey[i],"\t",ez[i],"\t" 
		 ,(int) color[i]);
	}

	//  Need extra lines for gnuplot to recognize blocks	
	o,endl,endl);

	o.fclose();
}*/
