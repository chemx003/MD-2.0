void sphere-terms();

void sphere-terms(double* x, double* y, double* z,
				  double* ex, double* ey, double* ez){
	double 	rhoS, rhoS18, rhoS19, sigmaS,
			r, r2, re,
			chiS, R,
			dx, dy, dz, hx, hy, hz;

	//  Components of seperation vector
	dx = x[i] - L/2.0;
	dy = y[i] - SL/2.0;
	dz = z[i] - SL/2.0;

	//  Magnitude of seperation
	r2 = dx*dx + dy*dy + dz*dz;
	r = sqrt(r2);

	//  Components of the unit vector
	hx = dx / r; 
	hy = dy / r;
	hz = dz / r;

	//  Scalar product of r and e
	re = hx*ex[i] + hy*ey[i] + hz*ez[i];

	chiR = 8.0 / (9.0 + 4*R);
	sigmaS = sqrt((1 + 4*R)/2);
	
	//  Potential at R
	rhoS = r - sigmaS*sqrt(1 - chiS*re) + 1.0;
	rhoS18 = pow(1.0/rhoS, 18);
	rhoS19 = rhoS18/rhoS;
	
}
