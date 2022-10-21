#include"opt_alg.h"

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };
		solution X0(x0);
		solution X1(matrix(x0 + d));

		if (X1.fit_fun(ff) == X0.fit_fun(ff)) {
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}
		if (X1.y > X0.y) {
			d = -d;
			X1.x = x0 + d;
			if (X1.fit_fun(ff) >= X0.fit_fun(ff)) {
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x - d);
				return p;
			}
		}
		solution X2;
		for (int i = 1; !(X2.fit_fun(ff) > X1.fit_fun(ff) || solution::f_calls >= Nmax);i++) {
			X0 = X1;
			X1 = X2;
			X2.x = x0 + pow(alpha, i) * d;
		}
		if (d > 0) {
			p[0] = X0.x(0);
			p[1] = X2.x(0);
		}
		else {
			p[0] = X2.x(0);
			p[1] = X0.x(0);
		}
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		double k = 0;
		for (;;k++)
			if (cacl_nth_element_fib(k) > (b - a) / epsilon)
				break;
		solution A(a), B(b), C, D;
		C.x = B.x - (cacl_nth_element_fib(k - 1) / cacl_nth_element_fib(k)) * (B.x - A.x);
		D.x = A.x + B.x - C.x;
		for (int i = 0; i <= k - 3; i++) {
			if (C.fit_fun(ff) < D.fit_fun(ff))
				B.x = D.x;
			else
				A.x = C.x;
			C.x = B.x - (cacl_nth_element_fib(k - i - 2) / cacl_nth_element_fib(k - i - 1)) * (B.x - A.x);
			D.x = A.x + B.x - C.x;
		}
		return C;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		
		solution A(a), B(b), C((a + b) / 2), D(0.);
		matrix l, m;
		while (1) {
			A.fit_fun(ff); //fit_fun() oblicza funkcje w celu
			B.fit_fun(ff);
			C.fit_fun(ff);
			l = A.y * (pow(B.x, 2) - pow(C.x,2)) + B.y * (pow(C.x,2)-pow(A.x,2)) + C.y*(pow(A.x,2)-pow(B.x,2));
			m = A.y * (B.x - C.x) + B.y * (C.x - A.x) + C.y * (A.x - B.x);
			if (m <= 0) {
				C.x = NAN;
				C.y = NAN;
				return C;
			}
			D.x = 0.5 * l / m;
			D.fit_fun(ff);
			if (A.x < C.x && C.x < D.x) {
				if (D.y < C.y) {
					A = C;
					C = D;
				}
				else {
					B = D;
				}
			}
			else if (A.x < D.x && D.x < C.x) {
				if (D.y > C.y) {
					B = C;
					C = D;
				}
				else {
					A = D;
				}
			}
			else {
				C.x = NAN;
				C.y = NAN;
				return C;
			}
			if (Nmax <= D.f_calls - 3 || A.x - B.x < epsilon || (D.x - C.x > -gamma && D.x - C.x < gamma)) {
				C.fit_fun(ff);
				return C;
			}
		}


		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
