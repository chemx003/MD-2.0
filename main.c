#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "auxillary_functions.h"
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>

/*-----------------------  Function Declarations -----------------------------*/
//  Gay-Berne: calulate the forces and torques
	void gb		(double* x, double* y, double* z,
				double* ex, double* ey, double* ez,
				double* fx, double* fy, double* fz,
				double* gx, double* gy, double* gz,
				int sphere);	

//  Initialize the simulation
	void init	(double* x, double* y, double* z,
			  	double* vx, double* vy, double* vz,
			  	double* ex, double* ey, double* ez,
			  	double* ux, double* uy, double* uz);

//  Numerically integrate equations of motion
	void iterate (double* x, double* y, double* z,
			  	double* vx, double* vy, double* vz,
			  	double* ex, double* ey, double* ez,
			  	double* ux, double* uy, double* uz,
			  	double* fx, double* fy, double* fz,
			  	double* gx, double* gy, double* gz,
			  	int sphere);

/*----------------------------------------------------------------------------*/


/*--------------------------  Global Variables  ------------------------------*/
//  Simulation Parameters
int 	N				= 4096,				//  Number of particles
		pcf_bins		= 400,				//  Number of bins for pcf
		pcf_num_steps	= 10,				// 	Steps to avg pcf over
		num_bin_x 		= 11,				//  Director bins
		num_bin_y		= 11,			
		num_bin_z		= 11,
		num_steps 		= 75000, 			//  Number of timesteps
		num_steps_eqbm	= 25000;			//  Number of eqbm timesteps

double 	dt 				= 0.0015, 			//  Length of time step
	   	temp_init 		= 0.8,				//  Initial temperature
	   	xi = 0, eta = 0,					//  Thermostat variables

	   	L				= 48.5,				//  Length of sim box
	   	SL				= 16.5,				//	Short length of sim box

	   	M				= 1.0,				//	Particle mass
	  	I				= 1.0,				//  Particle moment of inertia

	  	R				= 3.0,				//  Immersed sphere radius
	  	W				= 35000,			//  Anchoring coefficient

		KB				= 1.0,				//  Boltzmann Constant
		PI				= 3.14159265358979; //  Pi

//  Data
double	KT, KR, K, V, E, 					//  Pot, kin, tot energies
		P,									//  Pressure
		T,									//  Temperature
		density;							//  Density

/*----------------------------------------------------------------------------*/




//  Main function
int main(){
	
	//  Local variables
	double 	x[N], y[N], z[N], 				//  Particle coords
			vx[N], vy[N], vz[N],			//	Particle vel.
			ex[N], ey[N], ez[N],			//  Particle orient
			ux[N], uy[N], uz[N],			//  Particle ang. vel.

			fx[N], fy[N], fz[N],			//  Forces
			gx[N], gy[N], gz[N],			//  Gorques

			//  Bin coordinates
			x_dir[num_bin_x*num_bin_y*num_bin_z],
			y_dir[num_bin_x*num_bin_y*num_bin_z],
			z_dir[num_bin_x*num_bin_y*num_bin_z],
			
			//  Local director
			ex_dir[num_bin_x*num_bin_y*num_bin_z],
			ey_dir[num_bin_x*num_bin_y*num_bin_z],
			ez_dir[num_bin_x*num_bin_y*num_bin_z],	

			histo[pcf_bins][2],				//  Histogram for pcf
			histo2[pcf_bins][2],			//  Histogram for ocf
			avg_temp, avg_sop;				//  Holders for avgs

	double q[num_bin_x*num_bin_y*num_bin_z][3][4];		//  Q-tensor
	double eigenval[num_bin_x*num_bin_y*num_bin_z][4];	//  Eigenvalues

	int		c_temp, c_sop;					//  Counters for avgs

	int 	num_bins = num_bin_x*num_bin_y*num_bin_z;	//  Total bins

	clock_t start = clock(), diff;

	//  Initialize variables
	avg_temp = 0.0; avg_sop = 0.0; c_temp = 0.0; c_sop = 0.0;

	for(int i = 0; i < num_bins; i++) {
		for(int j = 0; j < 3; j++) {
			for(int k = 0; k < 4; k++) {
				q[i][j][k] = 0.0;
			}
		}
	}

	//  Set random seed
	srand(1);

	//  Check for correct initialization
	print_global_variables();

	//  Initialize simulation and calc first forces and torques
	init(x, y, z, vx, vy, vz,ex, ey, ez, ux, uy, uz);
	gb(x, y, z, ex, ey, ez, fx, fy, fz, gx, gy, gz, 0);
	calc_E(); print_energies();
	calc_temp(); print_temp();
	
	//  Equilibriation of nematic liquid crystal
	for(int i = 0; i < num_steps_eqbm; i++) {

		iterate(x, y, z, vx, vy, vz,
			   ex, ey, ez, ux, uy, uz,
			   fx, fy, fz,
			   gx, gy, gz, 0); 		

		calc_E(); write_energies(i);

		if(i%100 == 0) {
			printf("EQBM: %i\n", i);

			calc_E(); print_energies();
			calc_temp(); print_temp();

			density = N/(L*SL*SL);
			printf("density = %f\n\n", density);

			double sop = return_sopx(ex);
			printf("sop = %f\n\n", sop);

			//  Print current time
			diff = clock() - start;
			int msec = diff*1000 / CLOCKS_PER_SEC;
			printf("Time taken %d seconds %d milliseconds\n\n", 
					msec/1000, msec%1000);

			write_vectors(x, y, z, ex, ey, ez, i);
		}

		write_sop(ex, ey, ez, i);
		write_temp(i);
	}
	
	//  Insert sphere
	mark_particles(x, y, z, vx, vy, vz, 
					ex, ey, ez, ux, uy, uz, 
					fx, fy, fz, gx, gy, gz);

	resize(x, y, z, vx, vy, vz, 
			ex, ey, ez, ux, uy, uz, 
			fx, fy, fz, gx, gy, gz);

	//  Calculate the forces and torques w/ sphere
	gb(x, y, z, ex, ey, ez, fx, fy, fz, gx, gy, gz, 1);

	//  Equilibriate with sphere
	for(int i = num_steps_eqbm; i < num_steps + num_steps_eqbm; i++) {

		iterate(x, y, z, vx, vy, vz,
			   ex, ey, ez, ux, uy, uz,
			   fx, fy, fz,
			   gx, gy, gz, 1);		

		calc_E(); write_energies(i);

		if(i%100 == 0) {
			printf("%i\n", i);

			calc_E(); print_energies();
			calc_temp(); print_temp();

			density = N/(L*SL*SL - 4.0*PI*R*R*R/3.0);
			printf("density = %f\n", density);

			double sop = return_sopx(ex);
			printf("sop = %f\n\n", sop);

			//  Print current time
			diff = clock() - start;
			int msec = diff*1000 / CLOCKS_PER_SEC;
			printf("Time taken %d seconds %d milliseconds\n\n", 
					msec/1000, msec%1000);

			write_vectors(x, y, z, ex, ey, ez, i);
		}

		write_sop(ex, ey, ez, i);
		write_temp(i);

		//  Calculate average temp and sopx
		if(i>30000){
			avg_temp = avg_temp + return_temp();
			c_temp++;

			if(fabs(return_sopx(ex))<3.0){
				avg_sop = avg_sop + return_sopx(ex);
				c_sop++;
			}
		}

		//  Director field 
		if(i >= num_steps + num_steps_eqbm - 50){
			calc_dir_field(x, y, z, ex, ey, ez, x_dir, y_dir, z_dir,
				ex_dir, ey_dir, ez_dir, eigenval, q, 0);
		} 

		//  PCF & OCF
		if(num_steps + num_steps_eqbm - i <= pcf_num_steps){
			write_pcf(x, y, z, histo, 0);
			write_ocf(x, y, z, ex, ey, ez, histo2, 0);
		}

	}	
	calc_E(); print_energies();
	calc_temp(); print_temp();

	//  Calc average temps
	avg_temp = avg_temp/c_temp;
	avg_sop = avg_sop/c_sop;
	printf("AVG_TEMP = %f\n", avg_temp);
	printf("AVG_SOPX = %f\n\n", avg_sop);

	//  Average PCF & OCF data
	write_pcf(x, y, z, histo, 1);
	write_ocf(x, y, z, ex, ey, ez, histo2, 1); 

	//  Average director/q-tensor data
	calc_dir_field(x, y, z, ex, ey, ez, x_dir, y_dir, z_dir,
					ex_dir, ey_dir, ez_dir, eigenval, q, 1); 

	//  Print total time
	diff = clock() - start;
	int msec = diff*1000 / CLOCKS_PER_SEC;
	printf("Time taken %d seconds %d milliseconds\n\n", msec/1000, msec%1000);
}

//  Gay-Berne: Calulate the forces and torques
void gb		(double* x, double* y, double* z,
			double* ex, double* ey, double* ez,
			double* fx, double* fy, double* fz,
			double* gx, double* gy, double* gz, 
			int sphere){
	
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
			fxi, fyi, fzi, gx1, gy1, gz1, gx2, gy2, gz2,

			//  Sphere quantities
			rhoS, rhoS18, rhoS19, sigmaS,
			r, r2, r6, r7, re, re5, re6, root, root3,
			chiS;

	/*----------------------Set Parameters-----------------------*/
	mu = 2.0; nu = 1.0;
	kappa = 3.0; xappa = 5.0; //  xappa is a touch suspicious

	chi = (kappa*kappa - 1.0) / (kappa*kappa + 1.0);
	xhi = (pow(xappa, 1.0/mu) - 1.0) / (pow(xappa, 1.0/mu) + 1.0);

	rc = 6.0;
	/*-----------------------------------------------------------*/

	//  Initialize to zero
	for(int i = 0; i < N; i++) {
		fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0;
		gx[i] = 0.0; gy[i] = 0.0; gz[i] = 0.0;
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
			if(fabs(dx) > 0.5*L){
				dx = dx - L*(dx / fabs(dx));
			}

			if(fabs(dy) > 0.5*SL){
				dy = dy - SL*(dy / fabs(dy));
			}
			
			if(fabs(dz) > 0.5*SL){
				dz = dz - SL*(dz / fabs(dz));
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
	
	//  Calculate sphere interactions
	if(sphere == 1){
		for(int i = 0; i < N; i++){
			//  Components of seperation vector
			dx = x[i] - L/2.0;
			dy = y[i] - SL/2.0;
			dz = z[i] - SL/2.0;

			//  Magnitude of seperation
			r2 = dx*dx + dy*dy + dz*dz;
			r = sqrt(r2);

			if(r < R + 6.0){
				// printf("in sphere terms\n\n");

				r6 = r*r*r*r*r*r;
				r7 = r6*r;

				//  Components of the unit vector
				hx = dx / r; 
				hy = dy / r;
				hz = dz / r;

				//  Scalar product of r and e
				re = hx*ex[i] + hy*ey[i] + hz*ez[i];
				re5 = re*re*re*re*re;
				re6 = re5*re;

				chiS = 8.0 / (9.0 + 4*R*R);
				sigmaS = sqrt((1 + 4*R*R)/2);
				root = sqrt(1.0 - chiS*re*re);
				root3 = root*root*root;
				double sigSbyR = sigmaS/root;
				
				//  Potential at R
				rhoS = r - sigSbyR + 1.0;
				rhoS18 = pow(1.0/rhoS, 18);
				rhoS19 = rhoS18/rhoS;
				
				//  Force due to sphere
				fxi = 72 * rhoS19 * (hx - chiS/(2*root3)*(ex[i]/r - re*hx/r));
				fyi = 72 * rhoS19 * (hy - chiS/(2*root3)*(ey[i]/r - re*hy/r));
				fzi = 72 * rhoS19 * (hz - chiS/(2*root3)*(ez[i]/r - re*hz/r));

				//printf("%f\t%f\t%f\t%f\t\n", fx[i],fxi,r,rhoS);

				//  Torque due to sphere
				gx1 = 72 * rhoS19 * chiS*sigmaS/(2*root3) * hx;
				gy1 = 72 * rhoS19 * chiS*sigmaS/(2*root3) * hy;
				gz1 = 72 * rhoS19 * chiS*sigmaS/(2*root3) * hz;

				//  Force due to surface anchoring
				fxi = fxi - 6*W*(re6/r7*hx - re5/r6*(ex[i]/r - re*hx/r));
				fyi = fyi - 6*W*(re6/r7*hy - re5/r6*(ey[i]/r - re*hy/r));
				fzi = fzi - 6*W*(re6/r7*hz - re5/r6*(ez[i]/r - re*hz/r));

				//  Torqe due to surface anchoring
				gx1 = gx1 - 6*W*re5/r6*hx;
				gy1 = gy1 - 6*W*re5/r6*hy;
				gz1 = gz1 - 6*W*re5/r6*hz;

				fx[i] = fx[i] + fxi;
				fy[i] = fy[i] + fyi;
				fz[i] = fz[i] + fzi;

				gx[i] = gx[i] - gx1;
				gy[i] = gy[i] - gy1;
				gz[i] = gz[i] - gz1;

				V = V + 4*rhoS18 - W*re/6;
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

			sfvx, sfvy, sfvz,				//  Scaling factors
			sfux, sfuy, sfuz,

			a,b,							//  Particle spacing

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

	printf("SPACING = %f, %f\n\n", a,b);

	for(int i = 0; i < NUM_LINE; i++) {
		for(int j = 0; j < NUM_LINE; j++) {
			for(int k = 0; k< NUM_LINE; k++) {

				//  Generate lattice site
				double q = (i + 0.5) * a;
				double r = (j +	0.5) * b;
				double s = (k + 0.5) * b;

				if(p<N) {
					//  Assign lattice site to particle 
					x[p] = q;
					y[p] = r;
					z[p] = s;

					//  Assign orientations
					ex[p] = 1;
					ey[p] = 0; 
					ez[p] = 0;

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
					
					//  Add to particle count
					p++;
				}
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
}

/*  Calculate the new positions and orientations. Update the kinetic energy.
 *  Applies periodic boundary conditions.*/
void iterate	(double* x, double* y, double* z,
			  	double* vx, double* vy, double* vz,
			  	double* ex, double* ey, double* ez,
			  	double* ux, double* uy, double* uz,
			  	double* fx, double* fy, double* fz,
			  	double* gx, double* gy, double* gz,
			  	int sphere){

	double	xNew, yNew, zNew,			//  Place holders for position
			exNew, eyNew, ezNew,		//  Place holders for orientation

			vxi, vyi, vzi,				//  Velocities at current timestep
			uxi, uyi, uzi,				//  Ang. velocities at current timestep

			d1, d2, b, lm;				//  Dot products and the lagrange mult.

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

		//  Calculate lagrange multipliers
		d1 = ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i];
		d2 = ex[i]*ux[i] + ey[i]*uy[i] + ez[i]*uz[i];
		b = d2 + 1/dt;

		double lm1 = -b + sqrt(b*b - d1 - 2*d2/dt);
		double lm2 = -b - sqrt(b*b - d1 - 2*d2/dt);

		lm = lm1;
		if(fabs(lm2) < fabs(lm1)){ lm = lm2; }

		double mag = ex[i]*ex[i] + ey[i]*ey[i] + ez[i]*ez[i];

		//  Check for unit length
		if(fabs(1 - mag) > 0.1){
			printf("i = %i\t mag = %f\n", i, mag);
		}

		//  Check that lagrange multiplier is real
		if(isnan(lm)!=0){
			printf("LM IS NAN for i = %i\n", i);
			double mag = ex[i]*ex[i] + ey[i]*ey[i] + ez[i]*ez[i];
			printf("mag = %f\n", mag);
			printf("b = %f\n", b);
			printf("d2 = %f\n", d2);
			printf("(ux,uy,uz) = (%f,%f,%f)\n", ux[i], uy[i], uz[i]);
			printf("(gx,gy,gz) =  (%f,%f,%f)\n\n\n)", gx[i], gy[i], gz[i]);
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
	}
	
	//  Thermostat term
	xi = xi + (2*K - (5*N-3)*temp_init)*dt;

	KT = KR = K = 0;

	//  Calculate forces at the next timestep
	gb(x, y, z, ex, ey, ez, fx, fy, fz, gx, gy, gz, sphere);

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
