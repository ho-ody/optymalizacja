#include "opt_alg.h"
#include "lab.h"
namespace l6 {
	
	matrix f6(matrix x, matrix ud1, matrix ud2) {
		matrix y;
		y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;
		return y;
	}
	matrix df(double t, matrix Y, matrix ud1, matrix ud2) {
		double m1 = 5, m2 = 5, k1 = 1, k2 = 1, F = 1;
		double b1 = (ud1)(0), b2 = (ud1)(1);
		matrix dY(4, 1);
		dY(0) = Y(1);
		dY(1) = (-b1 * Y(1) - b2 * (Y(1) - Y(3)) - k1 * Y(0) - k2 * (Y(0) - Y(2))) / m1;
		dY(2) = Y(3);
		dY(3) = (F + b2 * (Y(1) - Y(3)) + k2 * (Y(0) - Y(2))) / m2;
		return dY;
	}
	matrix fR(matrix x, matrix ud1, matrix ud2) {
		matrix y;
		//x=[b1 / b2] y=roznica pomiedzy polozeniami
		int N = 1001;
		static matrix X(N, 2);
		if (solution::f_calls == 1)
		{
			ifstream S("polozenia.txt");
			S >> X;
			S.close();
			//DUZE X TO POLOZENIA Z WYKRESU
		}
		matrix Y0(4, new double[4]{ 0,0,0,0 });
		matrix temp=x[0];
		// i wtedy zamiast &x[0] &temp
		matrix* Y = solve_ode(df, 0, 0.1, 100, Y0, temp);
		for (int i = 0; i < N; i++)
		{
			y = y + abs(X(i, 0) - Y[1](i, 0))
				+ abs(X(i, 1) - Y[1](i, 2)); //ity punkt na wykresie 
		}
		y = y / (2 * N);
		return y;
	}


	void testowa_f_celu_a() {
		ofstream output("_output.csv");
		int N = 2, mi = 20, lambda = 40, Nmax = 10000;
		double epsilon = 1e-3;
		double sigma_v = 1;
		matrix limits(N, 2), sigma(N, 1, sigma_v);
		limits(0, 0) = limits(1, 0) = -5;
		limits(0, 1) = limits(1, 1) = 5;

		for (int i = 0; i < 100; i++) {
			solution opt = EA(f6, N, limits, mi, lambda, sigma, epsilon, Nmax);
			output << t(opt.x(0)) << t(opt.x(1)) << t(opt.y(0)) << solution::f_calls << endl;
			solution::clear_calls();
		}
	}

	void problem_rzeczywi() {
		ofstream output("_output.csv");
		int N = 2, mi = 20, lambda = 40, Nmax = 10000;
		double epsilon = 1e-3;
		double sigma_v = 1;
		matrix limits(N, 2), sigma(N, 1, sigma_v);
		limits(0, 0) = limits(1, 0) = -0.1;
		limits(0, 1) = limits(1, 1) = 3;

		//for (int i = 0; i < 100; i++) {
		//solution opt = EA(fR, N, limits, mi, lambda, sigma, epsilon, Nmax);
		//cerr << t(opt.x(0)) << t(opt.x(1)) << t(opt.y(0)) << solution::f_calls << endl;
		//solution::clear_calls();
		//}
		solution opt = EA(f6, N, limits, mi, lambda, sigma, epsilon, Nmax);
		opt.x(0) = 1.19871;
		opt.x(1) = 2.13469;

		matrix Y0(4, new double[4]{ 0,0,0,0 });
		matrix* s_X = solve_ode(df, 0, 0.1, 100, Y0, opt.x);
		cerr << endl;
		cerr << endl;

		//cerr << Y[0][0](1000) << endl;
		//cerr << Y[0][1](1000);
		//cerr << Y[1][0](1000) << endl;
		//cerr << Y[1][1](1000) << endl;
		for (int i = 0; i < 1001; i++) {
			output << t(s_X[1][0](i));
			output << t(s_X[1][1](i));
			//output << t(s_X2[1][0](i));
			//output << t(s_X2[1][1](i));
			output << endl;
		}

	}


}