#include"opt_alg.h"
#include "lab.h"
namespace l2 {
	matrix f2(matrix x, matrix ud1, matrix ud2) {
		double y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;
		return matrix(y);
	}
	matrix df(double t, matrix Y, matrix ud1, matrix ud2) {
		double mr = 1, mc = 9, l = 0.5, b = 0.5, a_ref = M_PI, w_ref =0;
		double I = mr * l * l / 3 + mc * l * l;
		double k1 = ud2(0), k2 = ud2(1);
		double M = k1 * (a_ref - Y(0)) + k2 * (w_ref - Y(1));
		matrix dY(2, 1);
		dY(0) = Y(1);
		dY(1) = (M - b * Y(1)) / I;
		return dY;
	}
	matrix fR(matrix x, matrix ud1, matrix ud2) {
		matrix y;
		matrix Y0(2, 1);
		matrix* Y = solve_ode(df, 0, 0.1, 100, Y0, ud1, x);
		int n = get_len(Y[0]);
		double a_ref = M_PI, w_ref = 0;
		y = 0;
		for (int i = 0; i < n; i++)
			y = y + 10 * pow(a_ref - Y[1](i, 0), 2) + pow(w_ref - Y[1](i, 1), 2) + pow(x(0) * (a_ref - Y[1](i, 0)) + x(1) * (w_ref - Y[1](i, 1)), 2);
		y = y * 0.1;
		return y;
	}
	void test_zbiez_metod() {
		double alphaHJ = 0.5, epsilon = 1e-3, s = 0.5;
		double alphaR = 2, beta = 0.5;
		//matrix s01(2, 1, s);
		int Nmax = 1000;

		matrix x0 = matrix(2, 1, 0.5);
		matrix Xs_HJ1 = trans(x0);
		matrix Xs_Rosen1 = trans(x0);

		//przyklad1
		solution optHJ1 = HJ(f2, x0, s, alphaHJ, epsilon, Nmax, Xs_HJ1);
		cerr << "przykad1==============================================================\nMetoda Hooke'Jeevsa\n" << optHJ1; solution::clear_calls();
		solution optR1 = Rosen(f2, x0, matrix(2, 1, s), alphaR, beta, epsilon, Nmax, Xs_Rosen1);
		cerr << "\nMetoda Rosenbrocka\n" << optR1 << endl; solution::clear_calls();

		//przyklad2
		s = 0.5;
		alphaHJ = 0.5;
		alphaR = 2;
		beta = 0.5;
		epsilon = 1e-3;
		Nmax = 1000;
		x0 = matrix(2, 1, 0.75);
		matrix Xs_HJ2 = trans(x0);
		matrix Xs_Rosen2 = trans(x0);

		solution optHJ2 = HJ(f2, x0, s, alphaHJ, epsilon, Nmax, Xs_HJ2);
		cerr << "przykad2==============================================================\nMetoda Hooke'Jeevsa\n" << optHJ2; solution::clear_calls();
		solution optR2 = Rosen(f2, x0, matrix(2, 1, s), alphaR, beta, epsilon, Nmax, Xs_Rosen2);
		cerr << "\nMetoda Rosenbrocka\n" << optR2 << endl; solution::clear_calls();
		//przyklad3
		s = 0.25;
		alphaHJ = 0.5;
		alphaR = 2;
		beta = 0.5;
		epsilon = 1e-3;
		Nmax = 1000;
		x0 = matrix(2, 1, 0.75);
		matrix Xs_HJ3 = trans(x0);
		matrix Xs_Rosen3 = trans(x0);

		solution optHJ3 = HJ(f2, x0, s, alphaHJ, epsilon, Nmax, Xs_HJ3);
		cerr << "przykad3==============================================================\nMetoda Hooke'Jeevsa\n" << optHJ3; solution::clear_calls();
		solution optR3 = Rosen(f2, x0, matrix(2, 1, s), alphaR, beta, epsilon, Nmax, Xs_Rosen3);
		cerr << "\nMetoda Rosenbrocka\n" << optR3 << endl; solution::clear_calls();
		//przyklad4
		s = 0.25;
		alphaHJ = 0.25;
		alphaR = 1.5;
		beta = 0.125;
		epsilon = 1e-8;
		Nmax = 1000;
		x0 = matrix(2, 1, 1.0);
		matrix Xs_HJ4 = trans(x0);
		matrix Xs_Rosen4 = trans(x0);

		solution optHJ4 = HJ(f2, x0, s, alphaHJ, epsilon, Nmax, Xs_HJ4);
		cerr << "przykad4==============================================================\nMetoda Hooke'Jeevsa\n" << optHJ4; solution::clear_calls();
		solution optR4 = Rosen(f2, x0, matrix(2, 1, s), alphaR, beta, epsilon, Nmax, Xs_Rosen4);
		cerr << "\nMetoda Rosenbrocka\n" << optR4 << endl; solution::clear_calls();
		//przyklad5
		s = 0.25;
		alphaHJ = 0.25;
		alphaR = 1.5;
		beta = 0.125;
		epsilon = 1e-8;
		Nmax = 1000;
		x0 = matrix(2, 1, -1.0);
		matrix Xs_HJ5 = trans(x0);
		matrix Xs_Rosen5 = trans(x0);

		solution optHJ5 = HJ(f2, x0, s, alphaHJ, epsilon, Nmax, Xs_HJ5);
		cerr << "przykad5==============================================================\nMetoda Hooke'Jeevsa\n" << optHJ5; solution::clear_calls();
		solution optR5 = Rosen(f2, x0, matrix(2, 1, s), alphaR, beta, epsilon, Nmax, Xs_Rosen5);
		cerr << "\nMetoda Rosenbrocka\n" << optR5 << endl; solution::clear_calls();
	}


	void testowa_f_celu_a() {
		ofstream output("_output.csv");
		double s = 0.5; //d³ugoœæ kroku:  0.5   0.25   0.125
		double epsilon = 1e-25, alphaHJ = 0.5, alphaR = 2, beta = 0.5;
		int Nmax = 1000;

		matrix x0(2, 1);
		matrix s0(2, 1, s);

		for (int j = 0; j < 3; j++) {
			switch (j) {
			case 0: s = 0.5; break;
			case 1: s = 0.25; break;
			case 2: s = 0.125; break;
			}
			matrix s0(2, 1, s);
			for (int i = 0; i < 100; i++) {
				x0 = 2 * rand_mat(2, 1) - 1;
				output << t(x0(0)) << t(x0(1));

				solution optHJ = HJ(f2, x0, s, alphaHJ, epsilon, Nmax);
				output << t(optHJ.x(0)) << t(optHJ.x(1)) << t(optHJ.y(0)) << t(optHJ.f_calls) << ";";
				solution::clear_calls();

				solution optR = Rosen(f2, x0, s0, alphaR, beta, epsilon, Nmax);
				output << t(optR.x(0)) << t(optR.x(1)) << t(optR.y(0)) << t(optHJ.f_calls) << t(s) << endl;
				solution::clear_calls();
			}
		}
	}
	void testowa_f_celu_b() {
		double s = 0.125, epsilon = 1e-5, alphaHJ = 0.5, alphaR = 2, beta = 0.5;
		int Nmax = 1000;

		matrix x0(2, 1);
		matrix s0(2, 1, s);
		x0(0) = -0.0161099; //wybrany przypadek
		x0(1) = 0.21992;

		solution optHJ = HJ(f2, x0, s, alphaHJ, epsilon, Nmax); 
		solution optR = Rosen(f2, x0, s0, alphaR, beta, epsilon, Nmax);
	}
	void problem_rzeczywi() {
		double alphaHJ = 0.5, epsilon = 1e-4, s = 0.1, alphaR = 2, beta = 0.5; int Nmax = 1000;
		matrix s0(2, 1, s);

		matrix x0 = 11 * rand_mat(2, 1);
		matrix Xs_HJ = trans(x0);
		matrix Xs_Rosen = trans(x0);

		solution opt_HJ = HJ(fR, x0, s, alphaHJ, epsilon, Nmax);
		cout << "Optymalne k1, k2 (HJ):" << endl << opt_HJ << endl;
		solution::clear_calls();
		solution optRos = Rosen(fR, x0, s0, alphaR, beta, epsilon, Nmax);
		cout << "Optymalne k1, k2 (Ros):" << endl << optRos << endl;
		solution::clear_calls();

		matrix Y0(2, 1);
		matrix* Y_HJ = solve_ode(df, 0, 0.1, 100, Y0, matrix(0), opt_HJ.x);
		matrix* Y_Ros = solve_ode(df, 0, 0.1, 100, Y0, matrix(0), optRos.x);


		ofstream output("_symul.csv");

		output << t("HJ-a") << t("HJ-w") << t("R-a") << t("R-w") << endl;
		for (int i = 0; i < 1000; i++) {
			output << t(Y_HJ[1][0](i));
			output << t(Y_HJ[1][1](i));
			output << t(Y_Ros[1][0](i));
			output << t(Y_Ros[1][1](i)) << endl;
		}
	}
}