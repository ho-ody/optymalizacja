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
void lab1() {
	//input data
	double x0		= -100;
	double d		= 1;
	double alpha	= 2;
	double epsilon	= 1e-5;
	double gamma	= 1e-5;
	double Nmax		= 10009;

	//calculations
	double* ab_range = expansion(f1, x0, d, alpha, Nmax);
	double a = *ab_range, b = *++ab_range;
	cerr << "expansion -> [" << a << "," << b << "]\n";
	double fib_v = m2d(fib(f1, a, b, epsilon).x);
	cerr << "fibonacci -> " << fib_v << endl;
	double lag_v = m2d(lag(f1, a, b, epsilon, gamma, Nmax).x);
	cerr << "lagrange  -> " << lag_v << endl;
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