#include"opt_alg.h"

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };
		solution X0(x0);
		solution X1(matrix(x0 + d));
		X0.fit_fun(ff);
		X1.fit_fun(ff);
		if (X0.y == X1.y) {
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}
		if (X1.y > X0.y) {
			d = -d;
			X1.x = X0.x + d;
			X1.fit_fun(ff);
			X0.fit_fun(ff);
			if ( X1.y >= X0.y) {
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x - d);
				return p;
			}
		}
		solution X2;
		for (int i = 1;;i++) {
			X2.x = x0 + pow(alpha, i) * d;
			X2.fit_fun(ff);
			if (X2.y >= X1.y || solution::f_calls >= Nmax)
				break;
			X0 = X1;
			X1 = X2;
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
		solution Xopt;
		Xopt.ud = b - a;
		int n = static_cast<int>(ceil(log2(sqrt(5) * (b - a) / epsilon) / log2((1 + sqrt(5)) / 2)));
		int* F = new int[n] {1, 1};
		for (int i = 2; i < n; ++i)
			F[i] = F[i - 2] + F[i - 1];
		solution A(a), B(b), C, D;
		C.x = B.x - 1.0 * F[n - 2] / F[n - 1] * (B.x - A.x);
		D.x = A.x + B.x - C.x;
		C.fit_fun(ff, ud1, ud2);
		D.fit_fun(ff, ud1, ud2);
		for (int i = 0; i <= n - 3; ++i)
		{
			if (C.y < D.y)
				B = D;
			else
				A = C;
			C.x = B.x - 1.0 * F[n - i - 2] / F[n - i - 1] * (B.x - A.x);
			D.x = A.x + B.x - C.x;
			C.fit_fun(ff, ud1, ud2);
			D.fit_fun(ff, ud1, ud2);
		}
		Xopt = C;
		Xopt.flag = 0;
		return Xopt;
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
	Xopt.ud = b - a;
	solution A(a), B(b), C, D, D_old(a);
	C.x = (a + b) / 2;
	A.fit_fun(ff, ud1, ud2);
	B.fit_fun(ff, ud1, ud2);
	C.fit_fun(ff, ud1, ud2);
	double l, m;
	while (true)
	{
		l = m2d(A.y * (pow(B.x) - pow(C.x)) + B.y * (pow(C.x) - pow(A.x)) + C.y * (pow(A.x) - pow(B.x)));
		m = m2d(A.y * (B.x - C.x) + B.y * (C.x - A.x) + C.y * (A.x - B.x));
		if (m <= 0)
		{
			Xopt = D_old;
			Xopt.flag = 2;
			return Xopt;
		}
		D.x = 0.5 * l / m;
		D.fit_fun(ff, ud1, ud2);
		if (A.x <= D.x && D.x <= C.x)
		{
			if (D.y < C.y)
			{
				B = C;
				C = D;
			}
			else
				A = D;
		}
		else if (C.x <= D.x && D.x <= B.x)
		{
			if (D.y < C.y)
			{
				A = C;
				C = D;
			}
			else
				B = D;
		}
		else
		{
			Xopt = D_old;
			Xopt.flag = 2;
			return Xopt;
		}
		Xopt.ud.add_row((B.x - A.x)());
		if (B.x - A.x < epsilon || abs(D.x() - D_old.x()) < gamma)
		{
			Xopt = D;
			Xopt.flag = 0;
			break;
		}
		if (solution::f_calls > Nmax)
		{
			Xopt = D;
			Xopt.flag = 1;
			break;
		}
		D_old = D;
	}
	return Xopt;

	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}


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

//ofstream output_HJ("_outputHJ.csv");	//zapis iteracje
//ofstream output_R("_outputR.csv");	//zapis iteracje
solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		
		solution XB(x0), XB_old, X;
		XB.fit_fun(ff, ud1, ud2);
		//output_HJ << t(XB.x(0)) << t(XB.x(1)) << endl;	//zapis iteracje
		while (true)
		{
			X = HJ_trial(ff, XB, s, ud1, ud2); //odpalenie etapu probnego
			if (X.y < XB.y) //sprawdzamy czy etap probny przyniosl poprawe
			{
				while (true) //etap roboczy wykonywany co chwile
				{
					XB_old = XB;
					XB = X;
					X.x = 2 * XB.x - XB_old.x;
					X.fit_fun(ff, ud1, ud2);
					X = HJ_trial(ff, X, s, ud2, ud2);
					if (X.y >= XB.y)
					break; //przerwanie etapu roboczego
					if (solution::f_calls > Nmax)
						return XB;
				}
			}
			else //zmniejszamy dlugosc kroku
				s *= alpha;
			if (s<epsilon || solution::f_calls>Nmax) { //warunki stopu
				XB.flag = 0;
				return XB;
			}
			//output_HJ << t(XB.x(0)) << t(XB.x(1)) << endl;	//zapis iteracje
		}
		Xopt.flag = 1;
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
		int n = get_dim(XB); //dlugosc wektora X
		matrix D = ident_mat(n); //kierunki macierz jednostkowa
		//Etap probny konczy sie porazka gdy zostajemy na podstawowym punkcie, nastyepnie zmniejszamy kroki az krok bedzie mniejszy
			// od epsilon
			//Po zakonczeniu etapu probnego, jesli znalezlismy nowy punkt wykonujemy etap roboczy, odbicie lustrzane starej bazy wzgledem nowej bazy
			// z punkty X (odbicia) odpalamy etap probny, to co zwroci  porownujemy z punktem symetrii, jesli jest lepszy wykonujemy kolejny raz etap roboczy
			// jesli jest gorszy anulujemy etap roboczy , wracamy do bazy i rozpoczynamy iteracje 
		solution X;
		for (int i = 0; i < n; ++i)
		{
			X.x = XB.x + s * D[i];
			X.fit_fun(ff, ud1, ud2);
			if (X.y < XB.y)
				XB = X;
			else
			{
				X.x = XB.x - s * D[i];
				X.fit_fun(ff, ud1, ud2);
				if (X.y < XB.y)
					XB = X;
			}
		}
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
		solution X(x0), Xt; //xt - xtemporary
		int n = get_dim(X);
		matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n); //l-wzgledne przesuniecie, p - macierz przechowujaca ilosc porazek, s - dlugosci krokow, d - macierz kierunkow
		X.fit_fun(ff, ud1, ud2);
		//output_R << t(X.x(0)) << t(X.x(1)) << endl;	//zapis iteracje
		while (true)
		{
			for (int i = 0; i < n; ++i)
			{
				Xt.x = X.x + s(i) * D[i];									// [] -> () jd
				Xt.fit_fun(ff, ud1, ud2);
				if (Xt.y < X.y)
				{
					X = Xt;
					l(i) += s(i); //dodajemy dlugosc ktora przebylismy		// =  ->  +=
					s(i) *= alpha; // przypadek gdy poczatkowy krok jest udany
				}
				else
				{
					++p(i); //zwiekszamy licznik porazek
					s(i) *= -beta; //zmniejszamy krok
				}
			}
			bool change = true;
			for (int i = 0; i < n; ++i)
				if (p(i) == 0 || l(i) == 0)					
				{
					change = false; //przypadek w ktorym nie trzeba zmieniac bazy
						break;
				}
			if (change) //przypadek obrotu
			{
				//pierwszy kierunek to pierwsza kolumna z macierzy Q podzielona przez dlugosc tej kolumny
				matrix Q(n, n), v(n, 1);
				for (int i = 0; i < n; ++i)
					for (int j = 0; j <= i; ++j)
						Q(i, j) = l(i);
				Q = D * Q;
				v = Q[0] / norm(Q[0]);
				D.set_col(v, 0);
				for (int i = 1; i < n; ++i)
				{
					matrix temp(n, 1);
					for (int j = 0; j < i; ++j)
						temp = temp + trans(Q[i]) * D[j] * D[j];
					v = (Q[i] - temp) / norm(Q[i] - temp);
					D.set_col(v, i);
				}
				s = s0;
				l = matrix(n, 1);
				p = matrix(n, 1);
			}
			double max_s = abs(s(0));
			for (int i = 1; i < n; ++i)
				if (max_s < abs(s(i)))
					max_s = abs(s(i));
			if (max_s<epsilon || solution::f_calls>Nmax) { // warunek stopu
				//X.fit_fun(ff, ud1, ud2);
				X.flag = 0;
				return X;
			}
			//output_R << t(X.x(0)) << t(X.x(1)) << endl;	//zapis iteracje
		}
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
