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
	double y = -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * M_PI, 2)) + 0.002 * pow(0.1 * x(), 2);
	return matrix(y);
}


double fit_fun(double x) {
	return -cos(0.1 * x) * pow(M_E, -(0.1 * x - 2 * M_PI) * (0.1 * x - 2 * M_PI)) + 0.002 * (0.1 * x) * (0.1 * x);
}

#include <iomanip>
void lab1()
{
	std::pair<double, double> result = metoda_ekspansji(fit_fun, 100, 2);

	cerr << setprecision(16);

	cerr << "a = " << result.first << " b = " << result.second << endl;

	double* rest = expansion(f1,100,2,1.1,1000);

	cerr << "a = " << rest[0] << " b = " << rest[1] << endl;




	cerr << "fib_min = " << metoda_fibonacci(fit_fun, result.first, result.second) << endl;

	double fib_vv = m2d(fib(f1, rest[0], rest[1], 1e-5).x);

	cerr << "fib_new = " << fib_vv << endl;

	double lag_vv = m2d(lag(f1, rest[0], rest[1],1e-5,1.1,1000).x);


	cerr << "lag_min = " << metoda_lagrangea(fit_fun, result.first, result.second, (result.first+ result.second)/2) << endl;
	cerr << "lag_new = " << lag_vv << endl;

	double x0 = -100;
	double d = 1;
	double alpha = 2;
	double epsilon = 1e-5;
	double gamma = 1e-5;
	double Nmax = 10009;

	double* ab_range = expansion(f1, x0, d, alpha, Nmax);
	double a = *ab_range;
	double b = *++ab_range;
	cerr << "expansion -> [" << a << "," << b << "]\n";
	double fib_v = m2d(fib(f1, a, b, epsilon).x);
	cerr << "fibonacci -> " << fib_v << endl;
	double lag_v = m2d(lag(f1, a, b, epsilon, gamma, Nmax).x);
	cerr << "lagrange  -> " << lag_v << endl;
	




}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
