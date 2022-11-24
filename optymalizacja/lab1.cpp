#include"opt_alg.h"
namespace l1 {

	matrix f1(matrix x, matrix ud1, matrix ud2) {
		double y = -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * 3.14, 2)) + 0.002 * pow(0.1 * x(), 2);
		return matrix(y);
	}

	void test_zbiez_metod() {
		//input data
		double x0 = 10;
		double d = 1;
		double alpha = 1.5;
		double epsilon = 1e-5;
		double gamma = 1e-5;
		double Nmax = 1000;

		//calculations
		double* ab_range = expansion(f1, x0, d, alpha, Nmax);
		double a = *ab_range, b = *++ab_range;
		cerr << "expansion -> [" << a << "," << b << "]\n";

		a = -10; b = 1;
		epsilon = 0.0001;
		gamma = 1e-7;
		solution::clear_calls();

		solution fib_ = fib(f1, a, b, epsilon);
		cerr << "fibonacci -> \n" << fib_;
		solution::clear_calls();

		solution lag_ = lag(f1, a, b, epsilon, gamma, Nmax);
		cerr << "lagrange  -> \n" << lag_;
	}
	matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
		double a = 0.98, b = 0.63, g = 9.81, PA = 0.75, TA = 90, PB = 1, DB = 0.00365665, Fin = 0.01, Tin = 10, DA = ud2();
		matrix dY(3, 1);
		double FAout = a * b * DA * sqrt(2 * g * Y(0) / PA);
		if (Y(0) <= 0) FAout = 0;
		double FBout = a * b * DB * sqrt(2 * g * Y(1) / PB);
		if (Y(1) <= 0) FBout = 0;

		dY(0) = -FAout;
		dY(1) = FAout + Fin - FBout;
		dY(2) = Fin / Y(1) * (Tin - Y(2)) + FAout / Y(1) * (TA - Y(2));

		return dY;
	}
	matrix fR(matrix x, matrix ud1, matrix ud2) {
		matrix y;
		double t[3] = { 5,1,10 };
		matrix Y0 = matrix(3, t);
		matrix* Y = solve_ode(df1, 0, 1, 1000, Y0, ud1, x);

		//zapis VA, VB, TB do pliku -> symulacja
		//ofstream file;file.open("sym.csv");file<<Y[1];file.close();

		double max_ = Y[1](0.2);
		for (int i = 1; i < get_len(Y[0]); i++)
			if (max_ < Y[1](i, 2))
				max_ = Y[1](i, 2);
		y = abs(max_ - 50);
		return y;
	}

	template <typename T>
	string to_string_with_precision(const T a_value, const int n = 22)
	{
		std::ostringstream out;
		out.precision(n);
		out << std::fixed << a_value;
		return out.str();
	}
	string stringify(double in) {
		string s = to_string_with_precision(in);
		std::replace(s.begin(), s.end(), '.', ',');
		return s;
	}

	void testowa_f_celu_a() {
		ofstream part1("_part1.txt");
		ofstream part2("_part2.txt");

		double x0, d = 1, epsilon = 1e-5, gamma = 1e-11, Nmax = 1000, alpha;
		random_device R;
		alpha = 2.2; //1.5, 2.2, 2.9
		for (int i = 0; i < 1; i++) {
			//exp
			x0 = 200.0 * R() / R.max() - 100; //x0 z przedzialu -100 ; 100
			double* ab_range = expansion(f1, x0, d, alpha, Nmax);
			double a = *ab_range, b = *++ab_range;
			part1 << stringify(x0) << "\t" << stringify(a) << "\t" << stringify(b) << "\t" << solution::f_calls << "\t"; //zapis do pliku part1
			solution::clear_calls();
			//przedzia³ na sztywno
			//a = -100; b = 100;
			//fib
			solution s_fib = fib(f1, a, b, epsilon);
			part1 << stringify(m2d(s_fib.x)) << "\t" << stringify(m2d(s_fib.y)) << "\t" << solution::f_calls << endl; //zapis do pliku part1
			solution::clear_calls();
			//lag
			solution s_lag = lag(f1, a, b, epsilon, gamma, Nmax);
			part2 << stringify(m2d(s_lag.x)) << "\t" << stringify(m2d(s_lag.y)) << "\t" << solution::f_calls << endl;	//zapis do pliku part2
			solution::clear_calls();
		}
	}
	void problem_rzeczy_b() {
		double x0, d = 1, epsilon = 1e-4, gamma = 1e-6, Nmax = 1000, alpha = 2;
		//fib
		solution o_fib = fib(fR, 0.0001, 0.01, epsilon);
		cerr << "fib->\n" << o_fib << endl << endl;
		solution::clear_calls();
		//lag
		solution o_lag = lag(fR, 0.0001, 0.01, epsilon, gamma, Nmax);
		cerr << "lag->\n" << o_lag << endl << endl;
	}

}