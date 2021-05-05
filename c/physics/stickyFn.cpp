#include <math.h>
#include <stdlib.h>
#include <iostream>
namespace bd{

void stickyF(double E, double rho, double beta, double k0, double& f, double& fprime) {
	//evaluate the sticky function and derivative for root finding
	double S = beta*2*rho*rho*E; double Sq = sqrt(S);
	f = exp(beta*E)/Sq-k0;
	fprime = beta*exp(E*beta)/S*(Sq-rho*rho/Sq);
}

double stickyNewton(double E, double rho, double k0, double beta) {
	//newton method to solve for E
	int cutoff = 2000; double tol = 1e-10; double x0 = E; double x1;
	double f = 0; double fprime = 0;
	for (int step = 0; step < cutoff; step++) {
		stickyF(x0, rho, beta, k0, f, fprime);
		x1 = x0 - f/fprime;
		if (abs(x1-x0)<tol) {
			return x1;
		}
		x0 = x1;
		//std::cout << x0 << "\n";
	}
	return x1;
}

}