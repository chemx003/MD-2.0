#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

/*--------------------------  Global Variables  ------------------------------*/
//  Simulation Parameters
extern int 	N;								//  Number of particles

extern double 	num_steps, 					//  Number of timesteps
	   	dt, 								//  Length of time step
	   	temp_init,							//  Initial temperature

	   	L,									//  Length of simulation box

	   	M,									//	Particle mass
	  	I,									//  Particle moment of inertia

		KB;									//  Boltzmann Constant

//  Data
extern double	K, V, E, 					//  Pot, kin, tot energies
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

//  Write the energies to a file 'time' refers to the current timestep
void write_energies(int time){
	ofstream o;
	o.open("energy.dat", ios::app);
	
	o << time << "\t" << V << "\t" << K << "\t" << E << endl;

	o.close();
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

	o.close();
}
