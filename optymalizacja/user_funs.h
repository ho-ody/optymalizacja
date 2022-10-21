#pragma once

#include"ode_solver.h"

std::pair<double, double> metoda_ekspansji(double f(double), double x0, double d);
double metoda_fibonacci(double f(double), double a, double b);
double metoda_lagrangea(double f(double), double a, double b, double c);

double cacl_nth_element_fib(int n);