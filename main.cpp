#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "auxillary_functions.h"
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
using namespace std;

/*-----------------------  Function Declarations  ----------------------------*/
//  Gay-Berne: Calulate the forces and torques
	void gb		(double* x, double* y, double* z,
				double* ex, double* ey, double* ez,
				double* fx, double* fy, double* fz,
				double* gx, double* gy, double* gz);	

//  Initialize the simulation
	void init	(double* x, double* y, double* z,
			  	double* vx, double* vy, double* vz,
			  	double* ex, double* ey, double* ez,
			  	double* ux, double* uy, double* uz);

//  Calculate pair correlation function
	void pc		(double* x, double* y, double* z);



//  Numerically integrate Newton's Equations
	void iterate (double* x, double* y, double* z,
			  	double* vx, double* vy, double* vz,
			  	double* ex, double* ey, double* ez,
			  	double* ux, double* uy, double* uz,
			  	double* fx, double* fy, double* fz,
			  	double* gx, double* gy, double* gz);

//  Newton Raphson
	double newton_raphson(double ex, double ey, double ez,
							double ux, double uy, double uz, int i);

//  rescale velocities
	void rescale(double* x, double* y, double* z,
			  	double* vx, double* vy, double* vz,
			  	double* ex, double* ey, double* ez,
			  	double* ux, double* uy, double* uz);

/*----------------------------------------------------------------------------*/


/*--------------------------  Global Variables  ------------------------------*/
//  Simulation Parameters
int 	N				= 216,			//  Number of particles
		pcf_bins		= 200,			//  Number of bins for pcf
		pcf_num_steps	= 200;			// 	Steps to avg pcf over

double 	num_steps 		= 100000, 		//  Number of timesteps
	   	dt 				= 0.0015, 		//  Length of time step
	   	temp_init 		= 0.9,			//  Initial temperature
	   	xi = 0, eta = 0,				//  Thermostat variables

	   	L				= 18.1,			//  Length of simulation box
	   	SL				= 6.1,			//	Short length of the simulation box

	   	M				= 1.0,			//	Particle mass
	  	I				= 1.0,			//  Particle moment of inertia

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
	double 	x[N], y[N], z[N], 				//  Particle coords at n
			vx[N], vy[N], vz[N],			//	Particle coords at n-1

			ex[N], ey[N], ez[N],			//  Particle orient at n
			ux[N], uy[N], uz[N],			//  Particle orient at n-1

			fx[N], fy[N], fz[N],			//  Forces
			gx[N], gy[N], gz[N],			//  Gorques
			histo[pcf_bins][2];				//  Histogram for pcf

	//  Iteration
	srand(1);

	print_global_variables();

	init(x, y, z, vx, vy, vz,ex, ey, ez, ux, uy, uz);	//  Initialize

	//  Calculate the forces and torques
	gb(x, y, z, ex, ey, ez, fx, fy, fz, gx, gy, gz);

	write_vectors(x, y, z, ex, ey, ez);

	calc_E(); print_energies();
	calc_temp(); print_temp();

	for(int i = 0; i < num_steps; i++) {

		iterate(x, y, z, vx, vy, vz,
			   ex, ey, ez, ux, uy, uz,
			   fx, fy, fz,
			   gx, gy, gz); 	//  Integrate the eqns of motion		

		calc_E(); write_energies(i);

		if(num_steps-i <= pcf_num_steps){
			write_pcf(x, y, z, histo, 0);
		}

		if(i%100 == 0) {
			cout << i << endl;

			calc_E(); print_energies();
			calc_temp(); print_temp();
		}

		/*if(i < 20000) {
			L = L - 0.0005;
			if(i%100 == 0 && i>10000) {
				rescale(x,y,z,xOld,yOld,zOld,ex,ey,ez,exOld,eyOld,ezOld);
			}
		}*/

		write_sop(ex, ey, ez, i);
		write_temp(i);

		if(i%10 == 0) {write_vectors(x, y, z, ex, ey, ez);}
	}

	calc_E(); print_energies();
	calc_temp(); print_temp();

	//  Analysis & Post-Processing
	write_pcf(x, y, z, histo, 1);
}

//  Gay-Berne: Calulate the forces and torques
void gb		(double* x, double* y, double* z,
			double* ex, double* ey, double* ez,
			double* fx, double* fy, double* fz,
			double* gx, double* gy, double* gz){
	
	double	mu, nu,					//  Exponents -adustable parameters
			kappa, xappa, 			//  Ratio of length and energy parameters
			chi, xhi, 				// 	Shape anisotropy parameters

			rc, cutterm, dcutterm,	//  Cuttoff distance

			dx, dy, dz,				//	Components of rij
			rij, rij2,				// 	Magnitude of rij
			hx, hy, hz,				//	Components of rij unit vector

			si, sj, sij,			//  Scalar products with e_i, e_j, and r_ij
			sp, sm,					//  Sums and differences of si and sj
			spchi, smchi,			// 	sp and sm divided by the term with chi
			spxhi, smxhi,			//  " with xhi

			sigma,					//  Distance parameter

			//  The well depth and its derivatives
			epsilon, eps1, eps2, deps_dsi, deps_dsj, deps_dsij,

			//  Rho and its multiples and derivatives
			rho, rho6, rho12, rhoterm, drhoterm,
			rhoc, rho6c, rho12c, 

			//  Derivatives of sigma
			dsig_dsi, dsig_dsj, dsig_dsij, prefac,

			//  Potential bewtween pair and derivatives
			pot, dpot_drij, dpot_dsi, dpot_dsj, dpot_dsij,

			//  Forces and torques between pairs
			fxi, fyi, fzi, gx1, gy1, gz1, gx2, gy2, gz2;

	/*----------------------Set Parameters-----------------------*/
	mu = 2.0; nu = 1.0;
	kappa = 3.0; xappa = 5.0; //  xappa is a touch suspicious

	chi = (kappa*kappa - 1.0) / (kappa*kappa + 1.0);
	xhi = (pow(xappa, 1.0/mu) - 1.0) / (pow(xappa, 1.0/mu) + 1.0);

	rc = 6.0;

	//  Example chooses eps0 and sigma_s equal 1 by units...
	/*-----------------------------------------------------------*/

	ofstream p;
	p.open("force.dat", ios::app);

	//  Resetting quantiies
	for(int i = 0; i < N; i++) {
		fx[i] = 0; fy[i] = 0; fz[i] = 0;
		gx[i] = 0; gy[i] = 0; gz[i] = 0;
	}
	V = 0; P = 0;

	//  Calculations
	for(int i = 0; i < N-1 ; i++) {
		for(int j = i+1; j<N; j++) {
			//  Components of rij
			dx = x[i] - x[j];
			dy = y[i] - y[j];
			dz = z[i] - z[j];

			//  Apply minimum image critereon (trying a different method)
			if(abs(dx) > 0.5*L){
				dx = dx - L*(dx / abs(dx));
			}

			if(abs(dy) > 0.5*SL){
				dy = dy - SL*(dy / abs(dy));
			}
			
			if(abs(dz) > 0.5*SL){
				dz = dz - SL*(dz / abs(dz));
			}

			//  magnitude of rj
			rij2 = dx*dx + dy*dy + dz*dz;
			rij = sqrt(rij2);

			//  if(rij <= 0.95){cout << "rij = " << rij << endl;}

			//  Calculate the unit vector of rij
			hx = dx / rij;
			hy = dy / rij;
			hz = dz / rij;

			if(rij < rc){
				//  Dot products, their sum/difference, and the fractions
				si = hx*ex[i] + hy*ey[i] + hz*ez[i]; 
				sj = hx*ex[j] + hy*ey[j] + hz*ez[j];
				sp = si + sj; sm = si - sj;
				sij = ex[i]*ex[j] + ey[i]*ey[j] + ez[i]*ez[j];
				spchi = sp / (1.0 + chi*sij);
				smchi = sm / (1.0 - chi*sij);
				spxhi = sp / (1.0 + xhi*sij);
				smxhi = sm / (1.0 - xhi*sij);

				//  Distance parameter (sigma)
				sigma = 1.0 / sqrt(1.0 - 0.5*chi*(sp*spchi + sm*smchi));

				//  Well depth (epsilon) !!Watch the notation!!
				eps1 = 1.0 / sqrt(1.0 - (chi*chi*sij*sij));
				eps2 = 1.0 - 0.5*xhi*(sp*spxhi + sm*smxhi);
				epsilon = pow(eps1, nu) * pow(eps2, mu);

				/*  Potential at rij										  */
				rho = rij - sigma + 1.0;
				rho6 = 1.0 / pow(rho, 6);
				rho12 = rho6 * rho6;
				rhoterm = 4.0 * (rho12 - rho6);
				drhoterm = -24.0 * (2.0*rho12 - rho6) / rho;
				pot = epsilon * rhoterm;

				//  Potential at rc
				rhoc = rc - sigma + 1.0;
				rho6c = 1.0 / pow(rhoc, 6);
				rho12c = rho6c * rho6c;
				cutterm = 4.0 * (rho12c - rho6c);
				dcutterm = -24.0 * (2.0*rho12c - rho6c) / rhoc;
				pot = pot - epsilon*cutterm;

				//  Derivatives of sigma
				prefac = 0.5 * chi * pow(sigma, 3);
				dsig_dsi = prefac * (spchi + smchi);
				dsig_dsj = prefac * (spchi - smchi);
				prefac = prefac * (0.5*chi);
				dsig_dsij = -prefac * (spchi*spchi - smchi*smchi);

				//  Derivatives of epsilon
				prefac = -mu * xhi * pow(eps1, nu) * pow(eps2, mu-1);
				deps_dsi = prefac * (spxhi + smxhi);
				deps_dsj = prefac * (spxhi - smxhi);
				prefac = prefac * (0.5*xhi);
				deps_dsij = -prefac * (spxhi*spxhi - smxhi*smxhi);
				deps_dsij = deps_dsij 
							+ nu*chi*chi*pow(eps1, nu+2)*pow(eps2, mu)*sij;

				//  Derivatives of the potential at rij
				dpot_drij = epsilon*drhoterm;
				dpot_dsi = rhoterm*deps_dsi - epsilon*drhoterm*dsig_dsi;
				dpot_dsj = rhoterm*deps_dsj - epsilon*drhoterm*dsig_dsj;
				dpot_dsij = rhoterm*deps_dsij - epsilon*drhoterm*dsig_dsij;

			/*	p<<i<< "     " << j << "     " << x[i] << "     " << x[j] 
					<< "     " << -rij << "     " << sigma
					<< "     " << rho << "     " << -dpot_drij*hx << endl;  */

				//  Forces at rij
				fxi = -dpot_drij*hx - dpot_dsi*(ex[i] - si*hx)/rij
					  - dpot_dsj*(ex[j] - sj*hx)/rij;
				fyi = -dpot_drij*hy - dpot_dsi*(ey[i] - si*hy)/rij
					  - dpot_dsj*(ey[j] - sj*hy)/rij;
				fzi = -dpot_drij*hz - dpot_dsi*(ez[i] - si*hz)/rij
					  - dpot_dsj*(ez[j] - sj*hz)/rij;

				//  Torques at rij
				gx1 = dpot_dsi*hx + dpot_dsij*ex[j];
				gy1 = dpot_dsi*hy + dpot_dsij*ey[j];
				gz1 = dpot_dsi*hz + dpot_dsij*ez[j];

				gx2 = dpot_dsj*hx + dpot_dsij*ex[i];
				gy2 = dpot_dsj*hy + dpot_dsij*ey[i];
				gz2 = dpot_dsj*hz + dpot_dsij*ez[i]; 

			/*	if(isnan(gx1)!=0 && i== 302){
					cout << "BEFORE CUTOFF" << endl;
					cout << "fxi = " << fxi << endl;
					cout << "rij = " << rij << endl;
					cout << "rhoterm = " << rhoterm << endl;
					cout << "rho = " << rho << endl;
					cout << "rho6 = " << rho6 << endl;
					cout << "rho12 = " << rho12 << endl;
					cout << "sigma = " << sigma << endl;
					cout << "sigma root term = " << 
						1.0 - 0.5*chi*(sp*spchi + sm*smchi) << endl;
					cout << "si = " << si << endl;
					cout << "sj = " << sj << endl;
					cout << "sp = " << sp << endl;
					cout << "sm = " << sm << endl;
					cout << "spchi, smshi = " << spchi << ", " << smchi 
						<< endl;
					cout << "(xi, yi, zi) = (" << x[i] << ", " << y[i] 
						<< ", " << z[i] << ")" << endl;
					cout << "(xj, yj, zj) = (" << x[j] << ", " << y[j] 
						<< ", " << z[j] << ")" << endl;
					cout << "(hx, hy, hz) = (" << hx << ", " << hy << ", " 
						<< hz << ")" << endl;
					cout << "(exi, eyi, ezi) = (" << ex[i] << ", " 
						<< ey[i] << ", " << ez[i] << ")" << endl;
					cout << "(exj, eyj, ezj) = (" << ex[j] << ", "
						<< ey[j] << ", " << ez[j] << ")" << endl;

					double eimag = sqrt(ex[i]*ex[i] + ey[i]*ey[i] + ez[i]*ez[i]);
					double ejmag = sqrt(ex[j]*ex[j] + ey[j]*ey[j] + ez[j]*ez[j]);

					cout << "eimag = " << eimag << endl;
					cout << "ejmag = " << ejmag << endl;
					cout << "(gxi, gyi, gzi) = (" << gx[i] << "," 
						<< gy[i] << "," << gz[i] << ")" << endl;
					cout << "spchi = " << spchi << endl; 
					cout << "smchi = " << smchi << endl;
					cout << "chi = " << chi <<endl;
					cout << "dpot_drij = " << dpot_drij << endl;
					cout << "dpot_dsi = " << dpot_dsi << endl;
					cout << "dpot_dsj = " << dpot_dsj << endl;
					cout << "dpot_dsij = " << dpot_dsij << endl;
					cout << "i = " << i << endl;
					cout << "j = " << j << endl << endl;
				}*/

				//  Derivatives of the potential at the cuttoff
				dpot_drij = epsilon * dcutterm;
				dpot_dsi = cutterm*deps_dsi - epsilon*dcutterm*dsig_dsi;
				dpot_dsj = cutterm*deps_dsj - epsilon*dcutterm*dsig_dsj;
				dpot_dsij = cutterm*deps_dsij - epsilon*dcutterm*dsig_dsij;

				//  Forces at cuttoff
				fxi = fxi + dpot_drij*hx + dpot_dsi*(ex[i] - si*hx)/rij
					  + dpot_dsj*(ex[j] - sj*hx)/rij;
				fyi = fyi + dpot_drij*hy + dpot_dsi*(ey[i] - si*hy)/rij
					  + dpot_dsj*(ey[j] - sj*hy)/rij;
				fzi = fzi + dpot_drij*hz + dpot_dsi*(ez[i] - si*hz)/rij
					  + dpot_dsj*(ez[j] - sj*hz)/rij;

				//  Torques at cuttoff
				gx1 = gx1 - dpot_dsi*hx - dpot_dsij*ex[j];
				gy1 = gy1 - dpot_dsi*hy - dpot_dsij*ey[j];
				gz1 = gz1 - dpot_dsi*hz - dpot_dsij*ez[j];

				gx2 = gx2 - dpot_dsj*hx - dpot_dsij*ex[i];
				gy2 = gy2 - dpot_dsj*hy - dpot_dsij*ey[i];
				gz2 = gz2 - dpot_dsj*hz - dpot_dsij*ez[i];

				/*if(j == 255 && fxi < -10){
					cout << "AFTER CUTOFF" << endl;
					cout << "fxi = " << fxi << endl;
					cout << "rij = " << rij << endl;
					cout << "rhoterm = " << rhoterm << endl;
					cout << "rhoc = " << rhoc << endl;
					cout << "rho6c = " << rho6c << endl;
					cout << "rho12c = " << rho12c << endl;
					cout << "sigma = " << sigma << endl;
					cout << "si = " << si << endl;
					cout << "sj = " << sj << endl;
					cout << "sp = " << sp << endl;
					cout << "sm = " << sm << endl;
					cout << "(xi, yi, zi) = (" << x[i] << ", " << y[i] 
						<< ", " << z[i] << ")" << endl;
					cout << "(xj, yj, zj) = (" << x[j] << ", " << y[j] 
						<< ", " << z[j] << ")" << endl;
					cout << "(hx, hy, hz) = (" << hx << ", " << hy << ", " 
						 << hz << ")" << endl;
					cout << "(exi, eyi, ezi) = (" << ex[i] << ", " 
						<< ey[i] << ", " << ez[i] << ")" << endl;
					cout << "(exj, eyj, ezj) = (" << ex[j] << ", "
						<< ey[j] << ", " << ez[j] << ")" << endl;
					cout << "spchi = " << spchi << endl; 
					cout << "smchi = " << smchi << endl;
					cout << "chi = " << chi <<endl;
					cout << "dpot_drij = " << dpot_drij << endl;
					cout << "dpot_dsi = " << dpot_dsi << endl;
					cout << "dpot_dsj = " << dpot_dsj << endl;
					cout << "i = " << i << endl;
					cout << "j = " << j << endl << endl;
				}*/

				//  Write forces and torques
				fx[i] = fx[i] + fxi;
				fy[i] = fy[i] + fyi;
				fz[i] = fz[i] + fzi;

				fx[j] = fx[j] - fxi;
				fy[j] = fy[j] - fyi;
				fz[j] = fz[j] - fzi;		

				gx[i] = gx[i] - gx1;
				gy[i] = gy[i] - gy1;
				gz[i] = gz[i] - gz1;

				gx[j] = gx[j] - gx2;
				gy[j] = gy[j] - gy2;
				gz[j] = gz[j] - gz2;

				// Calculate potential
				V = V + pot;
			}
		}
	}
}

/*  Initialize particle positions, orientations, velocities, and angular 
 *  velocities 																  */
void init	(double* x, double* y, double* z,
			 double* vx, double* vy, double* vz,
			 double* ex, double* ey, double* ez,
			 double* ux, double* uy, double* uz){

	double 	sumVx, sumVy, sumVz,			//  Set lin mtm = 0

			sumVx2, sumVy2, sumVz2, 		//  Set kinetic energy
			sumUx2, sumUy2, sumUz2,

			sfvx, sfvy, sfvz,				//  Scaling factor
			sfux, sfuy, sfuz,

			a,b,								//  Particle spacing

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
	b = SL / NUM_LINE;

	cout << "SPACING = " << a << endl << endl;

	for(int i = 0; i < NUM_LINE; i++) {
		for(int j = 0; j < NUM_LINE; j++) {
			for(int k = 0; k< NUM_LINE; k++) {
				if(p<N){
					//  Assign lattice sites to particle
					x[p] = (i + 0.5) * a;
					y[p] = (j +	0.5) * b;
					z[p] = (k + 0.5) * b;

					//  Assign random orientations
					ex[p] = 1; //dRand(-1.0, 1.0);
					ey[p] = 0; //dRand(-1.0, 1.0);
					ez[p] = 0; //dRand(-1.0, 1.0); 

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

		//  Calculate the kinetic energy
		KT = KT + 0.5 * M * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
		KR = KR + 0.5 * I * (ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i]);
		K = KT + KR;
	}

	//write_vectors(x,y,z,fx,fy,fz);
}

/*  Calculate the new positions and orientations. Update the kinetic energy.
 *  Applies periodic boundary conditions.*/
void iterate	(double* x, double* y, double* z,
			  	double* vx, double* vy, double* vz,
			  	double* ex, double* ey, double* ez,
			  	double* ux, double* uy, double* uz,
			  	double* fx, double* fy, double* fz,
			  	double* gx, double* gy, double* gz){

	double	xNew, yNew, zNew,			//  Place holders for position
			exNew, eyNew, ezNew,		//  Place holders for orientation

			vxi, vyi, vzi,				//  Velocities at current timestep
			uxi, uyi, uzi,				//  Ang. velocities at current timestep

			d1, d2, b, lm;					//  Dot products and the lagrange mult.

	double  color[N];

	//  Reset velocities and kinetic energies
	vxi = vyi = vzi = 0.0;
	uxi = uyi = uzi = 0.0;
	KT = KR = K = 0.0;

	//  Calculations
	for(int i = 0; i < N; i++) {
		vx[i] = vx[i]*(1 - xi*dt/2) + 0.5*fx[i]*dt/M;
		vy[i] = vy[i]*(1 - xi*dt/2) + 0.5*fy[i]*dt/M;
		vz[i] = vz[i]*(1 - xi*dt/2) + 0.5*fz[i]*dt/M;

		ux[i] = ux[i]*(1 - xi*dt/2) + 0.5*gx[i]*dt/I;
		uy[i] = uy[i]*(1 - xi*dt/2) + 0.5*gy[i]*dt/I;
		uz[i] = uz[i]*(1 - xi*dt/2) + 0.5*gz[i]*dt/I;

		d1 = ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i];
		d2 = ex[i]*ux[i] + ey[i]*uy[i] + ez[i]*uz[i];
		b = d2 + 1/dt;

		double lm1 = -b + sqrt(b*b - d1 - 2*d2/dt);
		double lm2 = -b - sqrt(b*b - d1 - 2*d2/dt);
		lm = lm1;
		if(abs(lm2) < abs(lm1)){ lm = lm2; }

		//cout << "lm = " << lm << endl;

		double mag = ex[i]*ex[i] + ey[i]*ey[i] + ez[i]*ez[i];
		if(abs(1-mag) > 0.1){
			cout << "i = " << i << ",mag = " << mag << endl;
		}

		if(isnan(lm)!=0){
			cout << "LM IS NAN for i = " << i << endl;
			double mag = ex[i]*ex[i] + ey[i]*ey[i] + ez[i]*ez[i];
			cout << "mag = " << mag << endl;
			cout << "b = " << b << endl;
			cout << "d2 = " << d2 << endl;
			cout << "(ux,uy,uz) = (" << ux[i]  << "," << uy[i]
				<< "," << uz[i] << ")" << endl;
			cout << "(gx,gy,gz) = (" << gx[i] << "," << gy[i]
				<< "," << gz[i] << ")";
			cout << endl << endl;
			kill(getpid(), SIGKILL);
		}

		x[i] = x[i] + vx[i]*dt;
		y[i] = y[i] + vy[i]*dt;
		z[i] = z[i] + vz[i]*dt;
 
		ex[i] = ex[i] + (ux[i] + lm*ex[i])*dt;
		ey[i] = ey[i] + (uy[i] + lm*ey[i])*dt;
		ez[i] = ez[i] + (uz[i] + lm*ez[i])*dt;

		//  Calculate the kinetic energy
		KT = KT + 0.5 * M * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
		KR = KR + 0.5 * I * (ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i]);
		K = KT + KR;
		//if(mag != 1.0){ cout << mag << endl; }
	}

	xi = xi + (2*K - (5*N-3)*temp_init)*dt;

	KT = KR = K = 0;

	gb(x, y, z, ex, ey, ez, fx, fy, fz, gx, gy, gz);

	for(int i = 0; i < N; i++) {
		vx[i] = (vx[i] + 0.5*fx[i]*dt/M)/(1 + xi*dt/2);
		vy[i] = (vy[i] + 0.5*fy[i]*dt/M)/(1 + xi*dt/2);
		vz[i] = (vz[i] + 0.5*fz[i]*dt/M)/(1 + xi*dt/2);

		d1 = ux[i]*ex[i] + uy[i]*ey[i] + uz[i]*ez[i];
		lm = 2*d1/dt + gx[i]*ex[i] + gy[i]*ey[i] + gz[i]*ez[i];

		ux[i] = (ux[i] + 0.5*(gx[i])*dt/I - 0.5*lm*ex[i]*dt/I)/(1 + xi*dt/2);
		uy[i] = (uy[i] + 0.5*(gy[i])*dt/I - 0.5*lm*ey[i]*dt/I)/(1 + xi*dt/2);
		uz[i] = (uz[i] + 0.5*(gz[i])*dt/I - 0.5*lm*ez[i]*dt/I)/(1 + xi*dt/2);

		//  Calculate the kinetic energy
		KT = KT + 0.5 * M * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
		KR = KR + 0.5 * I * (ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i]);
		K = KT + KR;
	}

	//write_vectors(x,y,z,fx,fy,fz);

	//  Apply periodic boundary conditions
	for(int i = 0; i < N; i++) {
		if(x[i] < 0.0){
			x[i] = x[i] + L;
		}
		else if(x[i] > L){
			x[i] = x[i] - L;
		}

		if(y[i] < 0.0){
			y[i] = y[i] + SL;
		}
		else if(y[i] > SL){
			y[i] = y[i] - SL;
		}

		if(z[i] < 0.0){
			z[i] = z[i] + SL;
		}
		else if(z[i] > SL){
			z[i] = z[i] - SL;
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

