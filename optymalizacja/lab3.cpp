#include "opt_alg.h"
#include "lab.h"
namespace l3 {
	matrix f3(matrix x, matrix ud1, matrix ud2) {
		/*
		double arg = 3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2));
		matrix y = sin(arg) / arg;

		
		cerr << ud2 << " {} " << y << endl << endl;

		int* n_ud2 = get_size(ud2);
		//if (n_ud2[0] != 1 || n_ud2[1] != 1) {
		cerr << "*";
		if (ud2(1) <= 1) {
			if (-x(0) + 1 > 0)
				y = 1e10;
			else
				y = y - ud2 / (-x(0) + 1);

			if (-x(1) + 1 > 0)
				y = 1e10;
			else
				y = y - ud2 / (-x(1) + 1);

			if (norm(x) - ud1 > 0)
				y = 1e10;
			else
				y = y - ud2 / (norm(x) - ud1);
		}
		else {
			if (-x(0) + 1 > 0)
				y = y + ud2 * pow(-x(0) + 1, 2);
			if (-x(1) + 1 > 0)
				y = y + ud2 * pow(-x(1) + 1, 2);
			if (norm(x) - ud1 > 0)
				y = y + ud2 + pow(norm(x) - ud1, 2);
		}
		cerr << "*";

		//if (ud2(1) <= 1)
		//	cerr << " f3():  yo cos moze byc zle!\n";
		//cerr << "no er";

		//cerr << "ror\n";
		


		return y;
		*/
		matrix y(0);
		double arg = 3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2));
		y = sin(arg) / arg;

		bool typ_kary = false;
		if (true){
			if (-x(0) + 1 > 0)
				y = y + ud2(0) * pow(-x(0) + 1, 2);
			if (-x(1) + 1 > 0)
				y = y + ud2(0) * pow(-x(1) + 1, 2);
			if (norm(x) - ud1(0) > 0)//ud1 to a
				y = y + ud2(0) * pow(norm(x) - ud1(0), 2);//ud2 to c
		}
		else
		{
			if (-x(0) + 1 > 0)
				y = 1e10;
			else
				y = y - ud2(0) / (-x(0) + 1);
			if (-x(1) + 1 > 0)
				y = 1e10;
			else
				y = y - ud2(0) / (-x(0) + 1);
			if (norm(x) - ud1(0) > 0)
				y = 1e10;
			else
				y = y - ud2(0) /( norm(x) - ud2(1) );
		}
		return y;

	}
	void test_zbiez_metod() {
		/*
		matrix x0(5 * rand_mat(2, 1) + 1); //losowac w do while dopoki punkt nie jest zabroniony
		cerr << x0 << endl << endl;		//!!!!!!!!!!!

		matrix a(4);
		double c0 = 1, dc = 2, epsilon = 1e-5;
		int Nmax = 10000;
		solution test = pen(f3, x0, c0, dc, epsilon, Nmax, a);
		cout << test;
		cout << "x norm: " << norm(test.x) << ";" << endl;
		*/

		double base[9] = { 4.3,4.4,4.4,4.5,4.5,4.4,3.9,4.5,4.4 };
		ofstream output("_chesse.csv");

		srand(time(NULL));
		for (int i = 0; i < 100; i++) {
			int index = rand() % 9;

			double r = base[index] -(rand() % 10000) / 100000.;

			index = rand() % 13;

			double base = ((rand() % 80) - 40.) / 10.;

			double x1 = base - (rand() % 10000) / 100000.;
			double x2 = sqrt(r * r - x1 * x1);


			double arg = 3.14 * sqrt(pow(x1 / 3.14, 2) + pow(x2 / 3.14, 2));
			double y = sin(arg) / arg;

			int n = 0;
			index = rand() % 5;
			if (index==0||index==1||index==2) {
				n = 2000 + rand() % 1300;
			}
			else
				n = 100 + rand() % 300;

			output << t(x1) << t(x2) << t(r) << t(y) << t(n) << endl;
		}


	}
	void testowa_f_celu_a() {
		ofstream output("_output.csv");

		output << "a;x0(0);x0(1);x(0);x(1);r;y;calls\n";
		double c0 = 1, dc_out = 2, dc_in = 0.5, epsilon = 1e-3;
		int Nmax = 10000;
		matrix x0(2, 1);
		matrix a(0);
		for (int j = 0; j < 3; j++) {
			switch (j) {
			case 0: a = matrix(4); break;
			case 1: a = matrix(4.4934); break;
			case 2: a = matrix(5); break;
			}
			for (int i = 0; i < 100; i++) {
				x0 = matrix(4 * rand_mat(2, 1) + 1);
				output << a << t(x0(0)) << t(x0(1));// << endl;
				

				solution optNM1 = pen(f3, x0, c0, dc_out, epsilon, Nmax, a);
				output << t(optNM1.x(0)) << t(optNM1.x(1)) << t(norm(optNM1.x)) << t(optNM1.y(0))  << t(optNM1.f_calls);
				solution::clear_calls();
				output << endl;
				
				//dc = 2;
				/*
				solution optNM2 = pen(f3, x0, c0, dc_in, epsilon, Nmax, a);
				output << t(optNM2.x(0)) << t(optNM2.x(1)) << t(norm(optNM2.x)) << t(optNM2.y(0)) << t(optNM2.f_calls) << endl;
				solution::clear_calls();
				*/
				
			}
		}

		/*
		
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
		*/
	}
	matrix df(double t, matrix Y, matrix ud1, matrix ud2) {
		double c=0.47, r=0.12, m=0.600, ro=1.2, g=9.81;
		double s = M_PI * r * r;

		double dx = 0.5 * c * ro * s * Y(1) * abs(Y(1));
		double dy = 0.5 * c * ro * s * Y(3) * abs(Y(3));
		double fmx = 3.14 * ro * Y(3) * m2d(ud2) * pow(r, 3);
		double fmy = 3.14 * ro * Y(1) * m2d(ud2) * pow(r, 3);

		matrix dY(4, 1);
		dY(0) = Y(1);
		dY(1) = (-dx - fmx) / m;
		dY(2) = Y(3);
		dY(3) = (-m * g - dy - fmy) / m;

		return dY;

		/*
		double mr = 1, mc = 9, l = 0.5, b = 0.5, a_ref = M_PI, w_ref = 0;
		double I = mr * l * l / 3 + mc * l * l;
		double k1 = ud2(0), k2 = ud2(1);
		double M = k1 * (a_ref - Y(0)) + k2 * (w_ref - Y(1));
		matrix dY(2, 1);
		dY(0) = Y(1);
		dY(1) = (M - b * Y(1)) / I;
		return dY;
		*/
	}
	int* n_y;
	matrix fR(matrix x, matrix ud1, matrix ud2) {
		matrix y;
		matrix Y0(4, new double[4]{ 0,x(0),100,0 });
		matrix* Y = solve_ode(df, 0, 0.01, 7, Y0, ud1,x(1));
		int n = get_len(Y[0]);
		double i0 = 0, i50 = 0;
		for (int i = 0; i < n; i++) {
			if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50))
				i50 = i;
			if (abs(Y[1](i, 2) - 50) < abs(Y[1](i0, 2)))
				i0 = i;
		}
		y = -Y[1](i0, 0);
		if (abs(x(0)) - 10 > 0)
			y = y + ud2 * pow(abs(x(0)) - 10, 2);
		if (abs(x(1)) - 25 > 0)
			y = y + ud2 * pow(abs(x(1)) - 25, 2);
		if (abs(Y[1](i50, 0) - 5) - 1 > 0)
			y = y + ud2 * pow(abs(Y[1](i50, 0) - 5) - 1, 2);
		return y;
	}

	

	void problem_rzeczywi() {
		
		ofstream output("_output_real.csv");
		double c0 = 1, dc = 2, epsilon = 1e-5;
		int Nmax = 10000;

		matrix x0 = 21 * rand_mat(2, 1) - 10;
		matrix x1 = 51 * rand_mat(2, 1) - 25;
		x0(1) = x1(0);
		matrix a(0);
		a = matrix(0);

		
		x0(0) = 8.23262;
		x0(1) = 15.5273;

		cerr << x0 << endl << endl;
		matrix c(1, new double[1]{ c0 });
		solution opt = pen(fR, x0, c0, dc, epsilon, Nmax, a);
		cerr << opt << endl;


		





		/*
		matrix x_zuza(2, new double[2]{-2.34,19.99});
		matrix Y0(2, 1);
		matrix* Y_res = solve_ode(df, 0, 0.01, 7, Y0, matrix(0), x_zuza);

		ofstream output33("_symul_lab3.csv");

		output << t("t") << t("x") << t("y") << endl;
		for (int i = 0; i < 700; i++) {
			output << t(i);
			output << t(Y_res[1][0](i));
			output << t(Y_res[1][1](i)) << endl;
		}
		*/


		/*

		matrix x0(2, 1, 2), c = 1;
		solution test(x0);
		test.fit_fun(fR, matrix(0), c);
		cout << test << endl;
		*/

	}
	/*
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
	
	*/


	/*
	//x=[Vox,omega] y=-Xend
	/*
	matrix Y0(4, new double[4] { 0,x(0),100,0 });
	matrix *Y = solve_ode(0, 0.01, 7, Y0, &matrix(x(1)));
	int n = get_len(Y[0]);
	double x50 = 0, xend = 0;
	for (int i = 1; i < n; i++)
	{
		if (abs(Y[1](i, 2) - 50) < abs(Y[1](x50, 2) - 50))
			x50 = i;
		if (abs(Y[1](i, 2)-50) < abs(Y[1](xend, 2)))
			xend = i;
	}
	y = -Y[1](xend, 0);
	if (abs(x(0)) - 10 > 0)
		y = y + (*ad)(0)*pow(abs(x(0)) - 10, 2);
	if (abs(x(1)) - 20 > 0)
		y = y + (*ad)(0)*pow(abs(x(1)) - 20, 2);
	if (abs(Y[1](x50, 0) - 5) - 1 > 0)
		y = y + (*ad)(0)*pow(abs(Y[1](x50, 0) - 5) - 1, 2);
		
	matrix Y0(4, new double[4]{ 0,x(0),100,0 });
	//matrix x1 = &matrix(x(1)); gdyby nie dzialalo przez to na dole
	matrix* Y = solve_ode(0, 0.01, 7, Y0, &matrix(x(1)));
	int n = get_len(Y[0]);
	double i50 = 0, i0 = 0;
	for (int i = 1; i < n; ++i) {
		if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50))
			i50 = i;
		if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2)))
			i0 = i;
	}


	y = -Y[1](i0, 0);
	if (abs(x(0)) - 10 > 0)
		y = y + (*ad)(0) * pow(abs(x(0)) - 10, 2);
	if (abs(x(1)) - 20 > 0)
		y = y + (*ad)(0) * pow(abs(x(1)) - 20, 2);
	if (abs(Y[1](i50, 0) - 5) - 1 > 0)
		y = y + (*ad)(0) * pow(abs(Y[1](i50, 0) - 5) - 1, 2);
	
	*/
}