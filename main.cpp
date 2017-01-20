#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "auxillary_functions.h"

using namespace std;

/*-----------------------  Function Declarations  ----------------------------*/
//  Gay-Berne: Calulate the forces and torques
	void gb		(double* x, double* y, double* z,
				double* ex, double* ey, double* ez,
				double* fx, double* fy, double* fz,
				double* gx, double* gy, double* gz);	

//  Initialize the simulation
	void init	(double* x, double* y, double* z,
			  	double* xOld, double* yOld, double* zOld,
			  	double* ex, double* ey, double* ez,
			  	double* exOld, double* eyOld, double* ezOld);

//  Calculate pair correlation function
	void pc		(double* x, double* y, double* z);



//  Numerically integrate Newton's Equations
	void verlet (double* x, double* y, double* z,
			  	double* xOld, double* yOld, double* zOld,
			  	double* ex, double* ey, double* ez,
			  	double* exOld, double* eyOld, double* ezOld,
			  	double* fx, double* fy, double* fz,
			  	double* gx, double* gy, double* gz);
/*----------------------------------------------------------------------------*/


/*--------------------------  Global Variables  ------------------------------*/
//  Simulation Parameters
int 	N			= 100;			//  Number of particles

double 	num_steps 	= 2000, 		//  Number of timesteps
	   	dt 			= 0.0015, 		//  Length of time step
	   	temp_init 	= 1.5,			//  Initial temperature

	   	L			= 6.0,			//  Length of simulation box

	   	M			= 1.0,			//	Particle mass
	  	I			= 1.0,			//  Particle moment of inertia

		KB			= 1.0;			//  Boltzmann Constant

//  Data
double	K, V, E, 					//  Pot, kin, tot energies
		P,							//  Pressure
		T;							//  Temperature

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
	srand(time(NULL));

	print_global_variables();

	init(x, y, z, xOld, yOld, zOld,
		 ex, ey, ez, exOld, eyOld, ezOld);			//  Initialize

	write_vectors(x, y, z, ex, ey, ez);

	calc_E(); print_energies();
	calc_temp(); print_temp();

	for(int i = 0; i < num_steps; i++) {
		gb(x, y, z,
		   ex, ey, ez,
		   fx, fy, fz,
		   gx, gy, gz); 		//  Calculate the forces and torques

		verlet(x, y, z, xOld, yOld, zOld,
			   ex, ey, ez, exOld, eyOld, ezOld,
			   fx, fy, fz,
			   gx, gy, gz); 	//  Integrate the eqns of motion	

		write_vectors(x, y, z, ex, ey, ez);
	}

	calc_E(); print_energies();
	calc_temp(); print_temp();

	//  Analysis & Post-Processing
	//pc();			//  Calculate PCF*/
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
	xhi = (pow(xappa, 1.0/mu) - 1.0) / (pow(xappa, 1.0/mu) - 1.0);

	rc = 3.0;

	//  Example chooses eps0 and sigma_s equal 1 by units...
	/*-----------------------------------------------------------*/

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

			if(abs(dy) > 0.5*L){
				dy = dy - L*(dy / abs(dy));
			}
			
			if(abs(dz) > 0.5*L){
				dz = dz - L*(dz / abs(dz));
			}

			//  magnitude of rij
			rij2 = dx*dx + dy*dy + dz*dz;
			rij = sqrt(rij2);

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
				sigma = 1.0 / sqrt(1.0 - 0.5*xhi*(sp*spchi + sm*smchi));

				//  Well depth (epsilon) !!Watch the notation!!
				eps1 = 1.0 / sqrt(1.0 - (chi*chi*sij*sij));
				eps2 = 1.0 - 0.5*xhi*(sp*spxhi + sm*smxhi);
				epsilon = pow(eps1, nu) * pow(eps2, mu);

				/*  Potential at rij - !!!something is off here !!!
				 *  some of the rij are < 1									*/
				rho = rij - sigma + 1.0;
				rho6 = 1.0 / pow(rho, 6);
				rho12 = rho6 * rho6;
				rhoterm = 4.0 * (rho12 - rho6);
				drhoterm = -24.0 * (2.0*rho12 - rho6) / rho;
				pot = epsilon * rhoterm;

				//  Potential at rc
				rho = rc - sigma + 1.0;
				rho6 = 1.0 / pow(rho, 6);
				rho12 = rho6 * rho6;
				cutterm = 4.0 * (rho12 - rho6);
				dcutterm = -24.0 * (2.0*rho12 - rho6) / rho;
				pot = pot - epsilon*cutterm;

				//  Derivatives of sigma
				prefac = 0.5 * chi * pow(sigma, 3);
				dsig_dsi = prefac * (spchi + smchi);
				dsig_dsj = prefac * (spchi + smchi);
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
				dpot_drij = epsilon * drhoterm;
				dpot_dsi = rhoterm*deps_dsi - epsilon*drhoterm*dsig_dsi;
				dpot_dsj = rhoterm*deps_dsj - epsilon*drhoterm*dsig_dsj;
				dpot_dsij = rhoterm*deps_dsij - epsilon*drhoterm*dsig_dsij;

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

				//  Derivatives of the potential at the cuttoff
				dpot_drij = epsilon * dcutterm;
				dpot_dsi = cutterm*deps_dsi - epsilon*dcutterm*dsig_dsi;
				dpot_dsj = cutterm*deps_dsj - epsilon*dcutterm*dsig_dsj;
				dpot_dsij = cutterm*deps_dsij - epsilon*dcutterm*dsig_dsij;

				//  Forces at cuttoff
				fxi = fxi + dpot_dsi*(ex[i] - si*hx)/rij
					  + dpot_dsj*(ex[j] - sj*hx)/rij;
				fyi = fyi + dpot_dsi*(ey[i] - si*hy)/rij
					  + dpot_dsj*(ey[j] - sj*hy)/rij;
				fzi = fzi + dpot_dsi*(ez[i] - si*hz)/rij
					  + dpot_dsj*(ez[j] - sj*hz)/rij;

				//  Torques at cuttoff
				gx1 = gx1 - dpot_dsi*hx - dpot_dsij*ex[j];
				gy1 = gy1 - dpot_dsi*hy - dpot_dsij*ey[j];
				gz1 = gz1 - dpot_dsi*hz - dpot_dsij*ez[j];

				gx2 = gx2 - dpot_dsj*hx - dpot_dsij*ex[i];
				gy2 = gy2 - dpot_dsj*hy - dpot_dsij*ey[i];
				gz2 = gz2 - dpot_dsj*hz - dpot_dsij*ez[i];

				// Calculate potential
				V = V + pot;
			}
		}
	}
}

/*  Initialize particle positions, orientations, velocities, and angular 
 *  velocities 																  */
void init	(double* x, double* y, double* z,
			 double* xOld, double* yOld, double* zOld,
			 double* ex, double* ey, double* ez,
			 double* exOld, double* eyOld, double* ezOld){

	double 	sumVx, sumVy, sumVz,			//  Set lin mtm = 0

			sumVx2, sumVy2, sumVz2, 		//  Set kinetic energy
			sumUx2, sumUy2, sumUz2,

			sfvx, sfvy, sfvz,				//  Scaling factor
			sfux, sfuy, sfuz,

			vx[N], vy[N], vz[N],			//  velocities
			ux[N], uy[N], uz[N],

			a,								//  Particle spacing

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

	for(int i = 0; i < NUM_LINE; i++) {
		for(int j = 0; j < NUM_LINE; j++) {
			for(int k = 0; k< NUM_LINE; k++) {
				if(p<N){
					//  Assign lattice sites to particle
					x[p] = (i + 0.5) * a;
					y[p] = (j + 0.5) * a;
					z[p] = (k + 0.5) * a;

					//  Assign random orientations
					ex[p] = dRand(-1.0, 1.0);
					ey[p] = dRand(-1.0, 1.0);
					ez[p] = dRand(-1.0, 1.0); 

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

		//  Set old positions and orientations for verlet algo
		xOld[i] = x[i] - dt * vx[i];
		yOld[i] = y[i] - dt * vy[i];
		zOld[i] = z[i] - dt * vz[i];
		/*  I don't think this is going to be unit length but correcting to 
		 * unit length might change the temperature/energy - TEST*/
		exOld[i] = ex[i] - dt * ux[i];
		eyOld[i] = ey[i] - dt * uy[i];
		ezOld[i] = ez[i] - dt * uz[i];

		//  Calculate the kinetic energy
		K = K + 0.5 * M * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])
			+ 0.5 * I * (ux[i]*ux[i] + uy[i]*uy[i] +vz[i]*vz[i]);
	}
}

/*  Calculate the new positions and orientations. Update the kinetic energy.
 *  Applies periodic boundary conditions.*/
void verlet	(double* x, double* y, double* z,
			  	double* xOld, double* yOld, double* zOld,
			  	double* ex, double* ey, double* ez,
			  	double* exOld, double* eyOld, double* ezOld,
			  	double* fx, double* fy, double* fz,
			  	double* gx, double* gy, double* gz){

	double	xNew, yNew, zNew,			//  Place holders for position
			exNew, eyNew, ezNew,		//  Place holders for orientation

			vxi, vyi, vzi,				//  Velocities at current timestep
			uxi, uyi, uzi,				//  Ang. velocities at current timestep

			d1, d2, lm;					//  Dot products and the lagrange mult.

	//  Reset kinetic energy
	K = 0;

	//  Calculations
	for(int i = 0; i < N; i++) {
		//  New postions
		xNew = 2*x[i] - xOld[i] + fx[i]*dt*dt/M;
		yNew = 2*y[i] - yOld[i] + fy[i]*dt*dt/M;
		zNew = 2*z[i] - zOld[i] + fz[i]*dt*dt/M;

		//  Store old positions and update current positions
		xOld[i] = x[i]; yOld[i] = y[i]; zOld[i] = z[i];
		x[i] = xNew; y[i] = yNew; z[i] = zNew;

		//  Velocities at the current timestep
		vxi = (x[i] - xOld[i]) / dt;
		vyi = (y[i] - yOld[i]) / dt;
		vzi = (z[i] - zOld[i]) / dt; 

		//  New Orientations
		exNew = 2*ex[i] - exOld[i] + gx[i]*dt*dt/I;
		eyNew = 2*ey[i] - eyOld[i] + gy[i]*dt*dt/I;
		ezNew = 2*ez[i] - ezOld[i] + gz[i]*dt*dt/I;

		//  Calculate the lagrange multiplier
		d1 = ex[i]*exNew + ey[i]*eyNew + ez[i]*ezNew;
		d2 = exNew*exNew + eyNew*eyNew + ezNew*ezNew;
		lm = -d1 + sqrt(d1*d1 - d2 + 1.0);

		//  Adjust new orientations
		exNew = exNew + ex[i]*lm; 
		eyNew = eyNew + ey[i]*lm;
		ezNew = ezNew + ez[i]*lm;

		//  Store old orientations and update current orientations
		exOld[i] = ex[i]; eyOld[i] = ey[i]; ezOld[i] = ez[i];
		ex[i] = exNew; ey[i] = eyNew; ez[i] = ezNew;

		//  Ang. velocities at the current timestep
		uxi = (ex[i] - exOld[i]) / dt;
		uyi = (ey[i] - eyOld[i]) / dt;
		uzi = (ez[i] - ezOld[i]) / dt;

		//  Update Kinetic Energy
		K = K + 0.5 * M * (vxi*vxi + vyi*vyi + vzi*vzi)
			+ 0.5 * I * (uxi*uxi + uyi*uyi +vzi*vzi);
	}

	//  Apply periodic boundary conditions
	for(int i = 0; i < N; i++) {
		if(x[i] < 0.0){
			x[i] = x[i] + L;
		}
		else if(x[i] > L){
			x[i] = x[i] - L;
		}

		if(y[i] < 0.0){
			y[i] = y[i] + L;
		}
		else if(y[i] > L){
			y[i] = y[i] - L;
		}

		if(z[i] < 0.0){
			z[i] = z[i] + L;
		}
		else if(z[i] > L){
			z[i] = z[i] - L;
		}
	}
}
