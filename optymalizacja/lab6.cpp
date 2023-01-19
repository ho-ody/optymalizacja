#include "opt_alg.h"
#include "lab.h"
namespace l6 {
	
	matrix f6(matrix x, matrix ud1, matrix ud2) {
		matrix y(1, 1);
		y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;
		return y;
	}

	void testowa_f_celu_a() {
		ofstream output("_output.csv");
		int N = 2, mi = 20, lambda = 40, Nmax = 10000;
		double epsilon = 1e-3;
		double sigma_v = 0.1;
		matrix limits(N, 2), sigma(N, 1, sigma_v);
		limits(0, 0) = limits(1, 0) = -5;
		limits(0, 1) = limits(1, 1) = 5;

		for (int i = 0; i < 100; i++) {
			solution opt = EA(f6, N, limits, mi, lambda, sigma, epsilon, Nmax);
			output << t(opt.x(0)) << t(opt.x(1)) << t(opt.y(0)) << solution::f_calls << endl;
			solution::clear_calls();
		}
	}




}