#pragma once
#include"ode_solver.h"

std::pair<double, double> metoda_ekspansji(double f(double), double x0, double d);
double metoda_fibonacci(double f(double), double a, double b);
double metoda_lagrangea(double f(double), double a, double b, double c);

double cacl_nth_element_fib(int n);

template <typename T>
string transform(T i) {
	stringstream ss;
	ss << i;
	string s = ss.str();
	replace(s.begin(), s.end(), '.', ',');
	s += ";";
	return s;
}
template <typename T>
string t(T i) {
	return transform(i);
}