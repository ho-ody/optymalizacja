#include "opt_alg.h"
#include "lab.h"
namespace l4 {
	matrix f4(matrix x, matrix ud1, matrix ud2) {
		matrix y(0);
		if (isnan(ud2(0, 0))) {
			y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
		}		
		else {
			y = f4(ud2[0] + x * ud2[1], ud1, NAN);
		}		
		return y;
	}
	matrix gf4(matrix x, matrix ud1, matrix ud2) {
		matrix g(2, 1);
		g(0) = 10 * x(0) + 8 * x(1) - 34;
		g(1) = 8 * x(0) + 10 * x(1) - 38;
		return g;
	}
	matrix hf4(matrix x, matrix ud1, matrix ud2) {
		matrix H(2, 2);
		H(0, 0) = 10;
		H(0, 1) = 8;
		H(1, 0) = 8;
		H(1, 1) = 10;
		return H;
	}
	matrix fR(matrix x, matrix ud1, matrix ud2) {
		int m = 100, n = get_len(x);
		static matrix X(n, m), Y(1, m);
		if (solution::f_calls == 1) {
			ifstream S("XData.txt");
			S >> X;
			S.close();
			S.open("YData.txt");
			S >> Y;
			S.close();
		}
		double h;
		matrix y(0);
		for (int i = 0; i < m; i++)
		{
			h = (trans(x) * X[i])();
			h = 1 / (1 + exp(-h));
			y = y - Y(0, i) * log(h) - (1 - Y(0, i)) * log(1 - h);
		}
		y = y / m;
		return y;
	}
	matrix gfR(matrix x, matrix ud1, matrix ud2) {
		int m = 100;
		int n = get_len(x);
		static matrix X(n, m), Y(1, m);
		if (solution::g_calls == 1)
		{
			ifstream S("XData.txt");
			S >> X;
			S.close();
			S.open("YData.txt");
			S >> Y;
			S.close();
		}
		double h;
		matrix g(n, 1);
		for (int j = 0; j < n; ++j)
		{
			for (int i = 0; i < m; ++i)
			{
				h = (trans(x) * X[i])();
				h = 1 / (1 + exp(-h));
				g(j) = g(j) + X(j, i) * (h - Y(0, i));
			}
			g(j) = g(j) / m;
		}
		return g;
	}
	void test_zbiez_metod() {
		double epsilon = 1e-5, h0 = 0.05;
		int Nmax = 10000;

		//for (int i = 0; i < 100; i++)
		matrix x0 = 21 * rand_mat(2, 1) - 10;
		x0(0) = 2.05956;
		x0(1) = 7.97644;

		cerr << x0 << endl << endl;

		h0 = 0.05;
		cerr << "  @:. h0 stalokrokwe 0.05\n";//h0 stalokrokwe 0.05
		solution SDopt0 = SD(f4, gf4, x0, h0, epsilon, Nmax);
		cerr << SDopt0 << endl; solution::clear_calls();

		solution CGopt0 = CG(f4, gf4, x0, h0, epsilon, Nmax);
		cerr << CGopt0 << endl; solution::clear_calls();

		solution NTopt0 = Newton(f4, gf4, hf4, x0, h0, epsilon, Nmax);
		cerr << NTopt0 << endl; solution::clear_calls();

		h0 = 0.12;
		cerr << "  @:. h0 stalokrokwe 0.12\n";//h0 stalokrokwe 0.12
		solution SDopt1 = SD(f4, gf4, x0, h0, epsilon, Nmax);
		cerr << SDopt1 << endl; solution::clear_calls();

		solution CGopt1 = CG(f4, gf4, x0, h0, epsilon, Nmax);
		cerr << CGopt1 << endl; solution::clear_calls();

		solution NTopt1 = Newton(f4, gf4, hf4, x0, h0, epsilon, Nmax);
		cerr << NTopt1 << endl; solution::clear_calls();

		h0 = -0.05;
		//system("cls");
		cerr << "  @:. h0 zmiennokrokowe -0.05, obliczane zlotym podzialem\n";//h0 zmiennokrokowe -0.05, obliczane z³otym podzia³em
		solution SDopt2 = SD(f4, gf4, x0, h0, epsilon, Nmax);
		cerr << SDopt2 << endl; solution::clear_calls();

		solution CGopt2 = CG(f4, gf4, x0, h0, epsilon, Nmax);
		cerr << CGopt2 << endl; solution::clear_calls();

		solution NTopt2 = Newton(f4, gf4, hf4, x0, h0, epsilon, Nmax);
		cerr << NTopt2 << endl; solution::clear_calls();
		
	}

	void testowa_f_celu_a() {
		ofstream output("_output.csv");
		double h0, epsilon = 1e-5;
		int Nmax = 10000;
		for (int i = 0; i < 3; i++) {
			switch (i) {
			case 0: h0 = 0.05; break;
			case 1: h0 = 0.12; break;
			case 2: h0 = -1.0; break;
			}
			for (int j = 0; j < 100; j++) {
				matrix x0 = 20 * rand_mat(2, 1) - 10;
				output << t(x0(0)) << t(x0(1));

				solution SDopt = SD(f4, gf4, x0, h0, epsilon, Nmax);
				output << t(SDopt.x(0)) << t(SDopt.x(1)) << t(SDopt.y(0)) << t(SDopt.f_calls) << t(SDopt.g_calls);
				solution::clear_calls();

				solution CGopt = CG(f4, gf4, x0, h0, epsilon, Nmax);
				output << t(CGopt.x(0)) << t(CGopt.x(1)) << t(CGopt.y(0)) << t(CGopt.f_calls) << t(CGopt.g_calls);
				solution::clear_calls();

				solution NTopt = Newton(f4, gf4, hf4, x0, h0, epsilon, Nmax);
				output << t(NTopt.x(0)) << t(NTopt.x(1)) << t(NTopt.y(0)) << t(NTopt.f_calls) << t(NTopt.g_calls) << t(NTopt.H_calls);
				solution::clear_calls();

				output << endl;
			}
		}
	}
	void testowa_f_celu_b() {
		double epsilon = 1e-5;
		int Nmax = 10000;
		matrix x0 = 21 * rand_mat(2, 1) - 10;
		//metoda najszybszego spadku
		solution SDopt0 = SD(f4, gf4, x0, 0.05, epsilon, Nmax); solution::clear_calls();
		solution SDopt1 = SD(f4, gf4, x0, 0.12, epsilon, Nmax); solution::clear_calls();
		solution SDopt2 = SD(f4, gf4, x0, -1.0, epsilon, Nmax); solution::clear_calls();
		//metoda gradientów sprzê¿onych
		solution CGopt0 = CG(f4, gf4, x0, 0.05, epsilon, Nmax); solution::clear_calls();
		solution CGopt1 = CG(f4, gf4, x0, 0.12, epsilon, Nmax); solution::clear_calls();
		solution CGopt2 = CG(f4, gf4, x0, -1.0, epsilon, Nmax); solution::clear_calls();
		//metoda newtona
		solution NTopt0 = Newton(f4, gf4, hf4, x0, 0.05, epsilon, Nmax); solution::clear_calls();
		solution NTopt1 = Newton(f4, gf4, hf4, x0, 0.12, epsilon, Nmax); solution::clear_calls();
		solution NTopt2 = Newton(f4, gf4, hf4, x0, -1.0, epsilon, Nmax); solution::clear_calls();
	}

	ofstream output_real("_output_real.csv");
	void single_call(double h0) {
		double epsilon = 1e-5;
		int Nmax = 10000;

		matrix x0(3, new double[3]{ 0,0,0 });
		solution solR = CG(fR, gfR, x0, h0, epsilon, Nmax);
	
		output_real << t(h0) << t(solR.x(0)) << t(solR.x(1)) << t(solR.x(2)) << t(solR.y(0));

		matrix X(3, 100), Y(1, 100);
		ifstream S("Xdata.txt");
		S >> X;
		S.close();
		S.open("Ydata.txt");
		S >> Y;
		S.close();
		double p = 0.0;
		for (int i = 0; i < 100; i++) {
			h0 = (trans(solR.x) * X[i])();
			h0 = 1 / (1 + exp(-h0));
			if ((h0 >= 0.5 && Y[i] == 1) || (h0 < 0.5 && Y[i] == 0)) {
				p++;
			}
		}
		p = p / 100;

		output_real << t(p) << t(solR.g_calls) << endl;
		solution::clear_calls();
	}
	void problem_rzeczywi() {
		single_call(0.01);
		single_call(0.001);
		single_call(0.0001);
	}




}