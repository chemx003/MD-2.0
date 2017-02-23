void sphere-terms();

void sphere-terms(double* x, double* y, double* z,
				  double* ex, double* ey, double* ez){

	double 	rhoS, rhoS18, rhoS19, sigmaS,
			r, r2, re,
			chiS, R, W,
			dx, dy, dz, hx, hy, hz
			fx, fy, fz,
			gx, gy, gz;

	//  Components of seperation vector
	dx = x[i] - L/2.0;
	dy = y[i] - SL/2.0;
	dz = z[i] - SL/2.0;

	//  Magnitude of seperation
	r2 = dx*dx + dy*dy + dz*dz;
	r = sqrt(r2);
	r6 = r*r*r*r*r*r;
	r7 = r6*r;

	//  Components of the unit vector
	hx = dx / r; 
	hy = dy / r;
	hz = dz / r;

	//  Scalar product of r and e
	re = hx*ex[i] + hy*ey[i] + hz*ez[i];
	re5 = re*re*re*re*re
	re6 = re5*re;

	chiS = 8.0 / (9.0 + 4*R);
	sigmaS = sqrt((1 + 4*R)/2);
	root = sqrt(1 - chiS*re);
	root3 = root*root*root;
	
	//  Potential at R
	rhoS = r - sigmaS*sqrt(1 - chiS*re) + 1.0;
	rhoS18 = pow(1.0/rhoS, 18);
	rhoS19 = rhoS18/rhoS;
	
	//  Force due to sphere
	fx = -72 * rhoS19 * (hx - chiS/(2*root3)*(ex[i]/r - re*hx/r));
	fy = -72 * rhoS19 * (hy - chiS/(2*root3)*(ey[i]/r - re*hx/r));
	fz = -72 * rhoS19 * (hz - chiS/(2*root3)*(ez[i]/r - re*hz/r));

	//  Torque due to sphere
	gx = 72 * rhoS19 * chiS/(2*root3) * hx;
	gy = 72 * rhoS19 * chiS/(2*root3) * hy;
	gz = 72 * rhoS19 * chiS/(2*root3) * hz;

	//  Force due to surface anchoring
	fx = fx - 6*W*(re6/r7*hx - re5/r6*(ex[i]/r - re*hx/r));
	fy = fy - 6*W*(re6/r7*hy - re5/r6*(ey[i]/r - re*hy/r));
	fz = fz - 6*W*(re6/r7*hz - re5/r6*(ez[i]/r - re*hz/r));

	//  Torqe due to surface anchoring
	gx = gx - 6*W*re5/r6*hx;
	gy = gy - 6*W*re5/r6*hy;
	gz = gz - 6*W*re5/r6*hz;
}
