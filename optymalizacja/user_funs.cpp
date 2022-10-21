#include"user_funs.h"

double fib_c(int n)
{
	return round(pow((1 + sqrt(5)) / 2, n) / sqrt(5));
}
double cacl_nth_element_fib(int n)
{
	return round(pow((1 + sqrt(5)) / 2, n) / sqrt(5));
}

std::pair<double,double> metoda_ekspansji(double f(double), double x0, double d) {
	int nmax = 1000;
	int i = 0;
	double epsilon = 1e-10;
	double alfa = 1.1;

	double x1 = x0 + d;
	if (abs(f(x0) - f(x1)) < epsilon) {
		return std::make_pair(x0, x1);
	}
	if (f(x1) > f(x0)) {
		d = -d;
		x1 = x0 + d;
		if (f(x1) >= f(x0)) {
			return std::make_pair(x1, x0 - d);
		}
	}
	double x_i = x0, x_iplus = x1, x_iminus;

	i = 1;
	for (int n = 0; n < nmax, f(x_i) > f(x_iplus); n++, i++) {
		x_iminus = x_i;
		x_i = x_iplus;
		x_iplus = x0 + pow(alfa, i) * d;
	}
	if (d > 0) {
		return std::make_pair(x_iminus, x_iplus);
	}
	return std::make_pair(x_iplus, x_iminus);
}






double metoda_fibonacci(double f(double), double a, double b) {
	double epsilon = 1e-5;
	int imax = 10000000, i;
	for (i = 0; i < imax; i++)
		if (fib_c(i) > (b - a) / epsilon)
			break;
	double k = i;
	
	double a_i = a, b_i = b;
	double c_i = b_i - (fib_c(k - 1) / fib_c(k)) * (b_i - a_i);
	double d_i = a_i + b_i - c_i;
	double a_iplus, b_iplus, c_iplus, d_iplus;

	for (int i = 0; i <= k - 3; i++) {
		
		if (f(c_i) < f(d_i)) {
			a_iplus = a_i;
			b_iplus = d_i;
		}
		else {
			b_iplus = b_i;
			a_iplus = c_i;
		}

		c_iplus = b_iplus - (fib_c(k - i - 2) / fib_c(k - i - 1)) * (b_iplus - a_iplus);
		d_iplus = a_iplus + b_iplus - c_iplus;

		a_i = a_iplus;
		b_i = b_iplus;
		c_i = c_iplus;
		d_i = d_iplus;
	}
	return c_iplus;
}

double metoda_lagrangea(double f(double), double a, double b, double c) {
	/*
	double epsilon = 1e-5;
	int i = 0, imax = 100000;
	double a_i = a, b_i = b, c_i = c, a_iplus, b_iplus, c_iplus, d_i, d_iplus, d_iminus;
	
	double gamma = 0.1;
	do {
		double l = f(a_i) * (pow(b_i, 2) - pow(c_i, 2)) + f(b_i) * (pow(c_i, 2) - pow(a_i, 2)) + f(c_i) * (pow(a_i, 2) - pow(b_i, 2));
		double m = f(a_i) * (b_i - c_i) + f(b_i) * (c_i - a_i) + f(c_i) * (a_i - b_i);
		if (m <= 0)
			cerr << " @ metoda_lagrangea: error0\n";
		d_i = 0.5 * l / m;
		if (a_i < d_i && d_i < c_i) {
			if (f(d_i) < f(c_i)) {
				a_iplus = a_i;
				c_iplus = d_i;
				b_iplus = c_i;
			}
			else {
				a_iplus = d_i;
				c_iplus = c_i;
				b_iplus = b_i;
			}
		}
		else {
			if (a_i < d_i && d_i < c_i) {
				if (f(d_i) < f(c_i)) {
					a_iplus = c_i;
					c_iplus = d_i;
					b_iplus = b_i;
				}
				else {
					a_iplus = a_i;
					c_iplus = c_i;
					b_iplus = d_i;
				}
			}
			else {
				cerr << " @ metoda_lagrangea: error1\n";
			}
		}
		i++;
		if (i > imax) {
			cerr << " @ metoda_lagrangea: error2\n";
			break;
		}

		a_i = a_iplus;
		b_i = b_iplus;
		c_i = c_iplus;
		d_iminus = d_i;
		//d_iminus = d_i;
		//d_i = d_iplus;

	} while (b_i - a_i < epsilon || abs(d_i- d_iminus) < gamma);
	
	return d_i;
	*/
	return -1;
}