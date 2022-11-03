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
//return;
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
	
	double a = 0.98, b = 0.63, g = 9.81, VA = 5, PA = 0.75, DB = 36.5665, VB = 1, PB = 1, Fin = 0.01, Tin = 10, TA = 90;
	matrix dY(3, 1);
	
	double FAout = a * b * m2d(ud2) * sqrt(2 * g * Y(1) / PA);
	double FBout = a * b * DB * sqrt(2 * g * Y(2) / PB);

	dY(0) = -1. * FAout;
	dY(1) = FAout + Fin - FBout;
	dY(2) = Fin / Y(1) * (Tin - Y(2)) + FAout / Y(1) * (TA-Y(2));
	
	return dY;
}
matrix fR(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	matrix Y0 = matrix(3, new double[3]{ 5,1,10 });
	matrix* Y = solve_ode(df1,0,1,1000,Y0,ud1,x);

	double max_ = -999999999999;
	for (int i = 0; i < get_len(Y[1])-1; i++) {
		if (Y[1](i,2) > max_) {
			max_ = Y[1](i, 2);
		}
	}
	y = abs(max_ - 50);

	return y;
}

void lab1_rzeczywiste() {

	double a = 1, b = 100, epsilon = 1e-5;

	//fib(fR, a, b, epsilon);
	double fib_v = m2d(fib(fR, a, b, epsilon).x);
	cerr << fib_v << endl;



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