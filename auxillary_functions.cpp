#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

/*--------------------------  Global Variables  ------------------------------*/
//  Simulation Parameters
extern int 	N,								//  Number of particles
	   		pcf_bins,						//  Number of bins for pcf
	   		pcf_num_steps;					//  Steps to avg pcf over

extern double 	num_steps, 					//  Number of timesteps
	   	dt, 								//  Length of time step
	   	temp_init,							//  Initial temperature
	   	xi, eta,							//  Thermostat variables

	   	L, SL,									//  Length of simulation box

	   	M,									//	Particle mass
	  	I,									//  Particle moment of inertia

		KB,									//  Boltzmann Constant
		PI;

//  Data
extern double	KT, KR, K, V, E, 					//  Pot, kin, tot energies
				P,							//  Pressure
				T;							//  Temperature

/*----------------------------------------------------------------------------*/

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
		while(abs(l - lOld) > 0.01 && isnan(l)==0){
			temp = l;
			l = lOld - (lOld*lOld + 2*lOld*(dot1+1/dt) + dot2  + 2*dot1/dt)/
				(2*lOld + 2*dot1 + 1/dt);
			count++;
			lOld = temp; 
			if(count > 1000){
				cout << count << endl;
				cout << "i = " << i <<",l = " <<  l << endl << endl;
			} 
		}
	}

	return l;
}

void print_energies(){
	cout << "ENERGIES" << endl << endl;

	cout << "KT = " << KT << endl;
	cout << "KR = " << KR << endl;
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

	cout << "KT = " << KT << endl;
	cout << "KR = " << KR << endl;
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

//  Rescale the velocities to temp_init .... move this to auxilary functions 

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
	ofstream o;
	o.open("energy.dat", ios::app);
	
	o << time << "\t" << V << "\t" << K << "\t" << E << endl;

	o.close();
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

	ofstream o;
	o.open("sop.dat", ios::app);

	if(abs(sopx)<3.0){
		o << time << " " << sopx << " " << sopy << " " << sopz << endl;
	}

	o.close();
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
		histo2[b][0] = R;

		for(int i = 0; i < N-1; i++) {
			for(int j = i+1; j < N; j++) {
				dx = x[i] - x[j];
				dy = y[i] - y[j];
				dz = z[i] - z[j];
				
				r = sqrt(dx*dx + dy*dy + dz*dz);

				if(r > R and r < R+dR){
					histo2[b][1] = histo2[b][1] 
						+ (3*(ex[i]*ex[j] + ey[i]*ey[j] + ez[i]*ez[j]) - 1)/2;
				}
			}
		}
		R = R + dR;
	}

	//  Write to file
	if(avg == 1){
		ofstream p;
		p.open("ocf.dat");

		R = dR;

		for(int b = 0; b < pcf_bins; b++) {
			histo2[b][1] = histo2[b][1] / pcf_num_steps;
			histo2[b][1] = 2*histo2[b][1] / (4*PI*R*R*dR*N*N/(L*SL*SL));
			p << histo2[b][0] << " " << histo2[b][1] << endl;	
			R = R + dR;
		}
		
		p.close();
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
				
				r = sqrt(dx*dx + dy*dy + dz*dz);

				if(r > R and r < R+dR){
					histo[b][1]++;
				}
			}
		}
		R = R + dR; //  Increment bin size
	}
	
	//  Write to file
	if(avg == 1){
		ofstream p;
		p.open("pcf.dat");

		R = dR;
		for(int b = 1; b < pcf_bins; b++) {
			histo[b][1] = histo[b][1] / pcf_num_steps;
			histo[b][1] = 2*histo[b][1] / (4*PI*R*R*dR*N*N/(L*SL*SL));
			p << histo[b][0] << " " << histo[b][1] << endl;
			R = R + dR;
		}
		
		p.close();
	}
}

//  Write one point to a file
void write_point(double x, double y, double z,
				 double ex, double ey, double ez){
	ofstream o;
	o.open("points.dat", ios::app);

	o << x << "    " << y << "     " << z << "     " 
	  << ex << "    " << ey << "    " << ez << "     " << endl;

	o.close();
}

//  Writes the temperature to a file
void write_temp(double time){
	ofstream o; 
	o.open("temp.dat", ios::app);

	calc_temp();

	o << time << "     " << T << endl;

	o.close();
}

//  Write the positions and orientations to a file for display puposes
void write_vectors(double* x, double* y, double* z,
				   double* ex, double* ey, double* ez){
	ofstream o;
	o.open("vector.dat", ios::app);

	for(int i = 0; i < N; i++){
		o << x[i] << "\t" << y[i] << "\t" << z[i] << "       "
		  << ex[i] << "\t" << ey[i] << "\t" << ez[i] << endl;
	}

	//  Need extra lines for gnuplot to recognize blocks	
	o << endl << endl << endl;

	o.close();
}

//  same as write, but extra integer to mark color
void write_vectors_colored(double* x, double* y, double* z,
				   double* ex, double* ey, double* ez,
				   double* color){
	ofstream o;
	o.open("vector.dat", ios::app);

	for(int i = 0; i < N; i++){
		o << x[i] << "\t" << y[i] << "\t" << z[i] << "\t"
		  << ex[i] << "\t" << ey[i] << "\t" << ez[i] << "\t" 
		  << (int) color[i] << endl;
	}

	//  Need extra lines for gnuplot to recognize blocks	
	o << endl << endl << endl;

	o.close();
}



