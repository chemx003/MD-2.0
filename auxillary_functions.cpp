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

	   	L,									//  Length of simulation box

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

	for(int i = 0; i < N; i++) {

		//wrt x
		sopx = sopx + (3.0*ex[i]*ex[i] - 1)/2.0;
		sopy = sopy + (3.0*ey[i]*ey[i] - 1)/2.0;
		sopz = sopz + (3.0*ez[i]*ez[i] - 1)/2.0;
	}

	sopx = sopx / N; sopy = sopy /N; sopz = sopz / N;

	ofstream o;
	o.open("sop.dat", ios::app);

	o << time << " " << sopx << " " << sopy << " " << sopz << endl;

	o.close();
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
			histo[b][1] = 2*histo[b][1] / (4*PI*R*R*dR*N*N/(L*L*L));
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



