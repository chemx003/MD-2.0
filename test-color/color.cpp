#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

//  Function declarations
double drand(double dmin, double dmax);

int main(){
	ofstream o;
	o.open("plot.dat");

	for(int i = 0; i < 75; i++) {
		o << drand(0.0, 10.0) << "   " << drand(0.0, 10.0) 
			<< "    " << drand (0.0, 10.0) << "   "  << i%2 << endl;
	}

	o.close();
}

double drand(double dmin, double dmax){
	double d = (double) rand() / RAND_MAX;
	return dmin + d * (dmax - dmin);
}
