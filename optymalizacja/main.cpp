/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
*********************************************/

#include"opt_alg.h"

void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

#define M_PI 3.14159265358979323846
#define M_E 2.71828182845904523536

using namespace std;

matrix f1(matrix x, matrix ud1, matrix ud2) {
	double y = -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * 3.14, 2)) + 0.002 * pow(0.1 * x(), 2);
	return matrix(y);
}
matrix f1_old(matrix x, matrix ud1, matrix ud2) {
	double y = -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * M_PI, 2)) + 0.002 * pow(0.1 * x(), 2);
	return matrix(y);
}
void lab1_rzeczywiste();
void testowa_f_celu_a();
void print_sol(solution in) {
	cerr << "\tx =       " << m2d(in.x) << endl;
	cerr << "\ty =       " << m2d(in.y) << endl;
	cerr << "\tf_calls = " << solution::f_calls << endl;
}
void reset_calls() {
	solution::f_calls = 0;
	solution::g_calls = 0;
	solution::H_calls = 0;
}

void lab1() {
	//lab1_rzeczywiste();
	testowa_f_celu_a();
	return;
	//input data
	double x0		= 10;
	double d		= 1;
	double alpha	= 1.5;
	double epsilon	= 1e-5;
	double gamma	= 1e-5;
	double Nmax		= 1000;

	//calculations
	double* ab_range = expansion(f1, x0, d, alpha, Nmax);
	double a = *ab_range, b = *++ab_range;
	cerr << "expansion -> [" << a << "," << b << "]\n";
	
	
	a = -10; b = 1;
	epsilon = 0.0001;
	gamma = 1e-7;
	reset_calls();
	solution temp = fib(f1, a, b, epsilon);

	cerr << "fibonacci -> " << endl;
	print_sol(temp);

	reset_calls();
	temp =lag(f1, a, b, epsilon, gamma, Nmax);
	cerr << "lagrange  -> " << endl;
	print_sol(temp);
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
	
	//double a = 0.98, b = 0.63, g = 9.81, VA = 5, PA = 0.75, DB = 36.5665, VB = 1, PB = 1, Fin = 0.01, Tin = 10, TA = 90;
	//double a = 0.98, b = 0.63, g = 9.81, PA = 1, TA = 90, PB = 1, DB = 0.00365665, Fin = 0.01, Tin = 10, DA = ud2();
	double a = 0.98, b = 0.63, g = 9.81, PA = 0.75, TA = 90, PB = 1, DB = 36.5665, Fin = 0.01, Tin = 10, DA = ud2();


	matrix dY(3, 1);
	
	double FAout = a * b * DA * sqrt(2 * g * Y(0) / PA);
	if (Y(0) <= 0)
		FAout = 0;
	double FBout = a * b * DB * sqrt(2 * g * Y(1) / PB);
	if (Y(1) <= 0)
		FBout = 0;

	dY(0) = -FAout;
	dY(1) = FAout + Fin - FBout;
	dY(2) = Fin / Y(1) * (Tin - Y(2)) + FAout / Y(1) * (TA-Y(2));
	
	return dY;
}
matrix fR(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	double t[3] = {5,1,10};
	matrix Y0 = matrix(3, t);
	matrix* Y = solve_ode(df1,0,1,1000,Y0,ud1,x);

	double max_ = Y[1](0.2);

	for (int i = 1; i < get_len(Y[0]); i++) {
		if (max_ < Y[1](i, 2))
			max_ = Y[1](i, 2);
		/*
		if (Y[1](i,2) > max_) {
			max_ = Y[1](i, 2);
		}
		*/
	}
	y = abs(max_ - 50);

	return y;
}

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 22)
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
		cerr << i << " ";
		x0 = 200.0 * R() / R.max() - 100; //x0 z przedzialu -100 ; 100

		part1 << stringify(x0) << "\t";
		double* ab_range = expansion(f1, x0, d, alpha, Nmax);
		double a = *ab_range, b = *++ab_range;
		a = -100; b = 100;
		part1 << stringify(a) << "\t";
		part1 << stringify(b) << "\t";
		part1 << solution::f_calls << "\t";

		solution::clear_calls();
		//fib
		cerr << "b";
		solution s_fib = fib(f1, a, b, epsilon);
		part1 << stringify(m2d(s_fib.x)) << "\t";
		part1 << stringify(m2d(s_fib.y)) << "\t";
		part1 << solution::f_calls << endl;
		solution::clear_calls();

		//cerr << s_fib.y << " vs " << f1(s_fib.x, matrix(0), matrix(0)) << endl;

		//lag
		cerr << "\n\n\n\n\nc";
		solution s_lag = lag(f1, a, b, epsilon, gamma, Nmax);
		part2 << stringify(m2d(s_lag.x)) << "\t";
		part2 << stringify(m2d(s_lag.y)) << "\t";
		part2 << solution::f_calls << endl;
		solution::clear_calls();
		cerr << "d\n";

	}



}

void lab1_rzeczywiste() {

	double x0, d = 1, epsilon = 1e-5, gamma = 1e-200, Nmax = 100, alpha = 2;

	random_device R;
	x0 = 99.0 * R() / R.max() + 1; //random x0
	cerr << "x0=" << x0 << endl;

	double* ab_range = expansion(f1, x0, d, alpha, Nmax);
	double a = *ab_range, b = *++ab_range;
	cerr << "expansion -> [" << a << "," << b << "]\n";

	//a = 11, b = 21;

	//fib(fR, a, b, epsilon);
	solution s_fib = fib(fR, a, b, epsilon);
	cerr << "fibonacci -> " << endl;
	cerr << s_fib << endl;
	solution::clear_calls();


	solution s_lag = lag(fR, a, b, epsilon, gamma, Nmax);
	cerr << "lagrange  -> " << endl;
	cerr << s_lag << endl;
	solution::clear_calls();

}


void lab2() {

}

void lab3() {

}

void lab4() {

}

void lab5(){

}

void lab6() {

}