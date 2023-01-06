#include"opt_alg.h"

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };
		solution X0(x0);
		solution X1(matrix(x0 + d));
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);
		if (X0.y == X1.y) {
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}
		if (X1.y > X0.y) {
			d = -d;
			X1.x = X0.x + d;
			X1.fit_fun(ff, ud1, ud2);
			X0.fit_fun(ff, ud1, ud2);
			if ( X1.y >= X0.y) {
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x - d);
				return p;
			}
		}
		solution X2;
		for (int i = 1;;i++) {
			X2.x = x0 + pow(alpha, i) * d;
			X2.fit_fun(ff, ud1, ud2);
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

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c0, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
		solution X(x0), X1;
		matrix c(2, new double[2]{ c0, dc }); //wysylane do  fitfun c to waga przy karze, dc do zorientowania jaka kara
		//matrix c(1, new double[1]{ c0 });
		while (true)
		{
			X1 = sym_NM(ff, X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c);
			if (norm(X.x - X1.x) < epsilon) {
				X1.flag = 0;
				break;
			}
			if (solution::f_calls > Nmax) {
				X1.flag = 1;
				break;
			}

			/*
			if (solution::f_calls > Nmax || norm(X.x - X1.x) < epsilon) {
				X1.flag = 0;
				return X1;
			}		
			*/
			c(0) *= dc; //obliczenie kary
			X = X1;
		}
		// X1.y _l;
		return X1;
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
		int n = get_len(x0);
		matrix D = ident_mat(n); //do generowania punktow macierz jednostkowa 
		int N = n + 1; //liczba wierzcholkow sympleksow
		solution* S = new solution[N]; //sympleks
		S[0].x = x0;
		S[0].fit_fun(ff, ud1, ud2);
		for (int i = 1; i < N; ++i)
		{
			S[i].x = S[0].x + s * D[i - 1];
			S[i].fit_fun(ff, ud1, ud2);

			int* n_y = get_size(S[i].y);
			//cerr << n_y[0] << " " << n_y[1] << endl;
			//cerr << S[i].y << endl;

			/*
			if (n_y[0] != 1 || n_y[1] != 1)
				cerr << "mamy problem\n";
			else
				cerr << "my nie lol\n";
			cerr << S[i] << endl;
			*/

		}
		

		solution PR, PE, PN; //odbicie, ekspansja ,zawezenie
		matrix pc; //srodek ciezkosci sympleksu
		int i_min, i_max; //najlepszy/najgorszy wierzcholek
		while (true)
		{
			i_min = i_max = 0;
			//cerr << "no er";
			for (int i = 1; i < N; ++i)
			{
				/*
				cerr << S[i_min] << endl;
				cerr << S[i] << endl;
				cerr << endl;
				*/
				if (S[i_min].y > S[i].y)
					i_min = i;
				if (S[i_max].y < S[i].y)
					i_max = i;
			}
			//cerr << "ror\n";
			pc = matrix(n, 1);
			for (int i = 0; i < N; ++i)
				if (i != i_max)
					pc = pc + S[i].x;
			//pc = pc / (n); //mamy pc
			pc = pc / (N - 1);
			PR.x = pc + alpha * (pc - S[i_max].x);

			PR.fit_fun(ff, ud1, ud2);
			if (S[i_min].y <= PR.y && PR.y < S[i_max].y)
				S[i_max] = PR;
			else if (PR.y < S[i_min].y)
			{
				PE.x = pc + gamma * (PR.x - pc);
				PE.fit_fun(ff, ud1, ud2);
				if (PR.y <= PE.y) //Jesli PE nie jest lepsze to bierzemy PR
					S[i_max] = PR;
				else
					S[i_max] = PE;
			}
			else
			{
				PN.x = pc + beta * (S[i_max].x - pc);
				PN.fit_fun(ff, ud1, ud2);
				if (PN.y < S[i_max].y)
					S[i_max] = PN;
				else
				{
					for (int i = 0; i < N; ++i)
						if (i != i_min)
						{
							S[i].x = delta * (S[i].x + S[i_min].x);
							S[i].fit_fun(ff, ud1, ud2);
						}
				}
			}
			double max_s = norm(S[0].x - S[i_min].x);
			for (int i = 1; i < N; ++i)
				if (max_s < norm(S[i].x - S[i_min].x))
					max_s = norm(S[i].x - S[i_min].x);

			if (solution::f_calls > Nmax || max_s < epsilon) { //!!!!!!!!!!!!!!!!!!!!
				S[i_min].flag = 0;
				S[i_min].fit_fun(ff, ud1, ud2);
//				cerr << "returning: " << S[i_min].y;
				return S[i_min];
			}		

			////cerr << S[i_min] << endl << endl;
		}
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

/*
ofstream output_SD_0_05("_output_SD_0_05.csv");
ofstream output_SD_0_12("_output_SD_0_12.csv");
ofstream output_SD_m_zk("_output_SD_m_zk.csv");
ofstream output_CG_0_05("_output_CG_0_05.csv");
ofstream output_CG_0_12("_output_CG_0_12.csv");
ofstream output_CG_m_zk("_output_CG_m_zk.csv");
ofstream output_Newton_0_05("_output_Newton_0_05.csv");
ofstream output_Newton_0_12("_output_Newton_0_12.csv");
ofstream output_Newton_m_zk("_output_Newton_m_zk.csv");
void write(string x1, string x2, double h0, string method) {
	if (method == "SD") {
		if (h0 == 0.05)
			output_SD_0_05 << x1 << x2 << endl;
		if (h0 == 0.12)
			output_SD_0_12 << x1 << x2 << endl;
		if (h0 < 0)
			output_SD_m_zk << x1 << x2 << endl;
	}
	else if (method == "CG") {
		if (h0 == 0.05)
			output_CG_0_05 << x1 << x2 << endl;
		if (h0 == 0.12)
			output_CG_0_12 << x1 << x2 << endl;
		if (h0 < 0)
			output_CG_m_zk << x1 << x2 << endl;
	}
	else if (method == "Newton") {
		if (h0 == 0.05)
			output_Newton_0_05 << x1 << x2 << endl;
		if (h0 == 0.12)
			output_Newton_0_12 << x1 << x2 << endl;
		if (h0 < 0)
			output_Newton_m_zk << x1 << x2 << endl;
	}
}
*/

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		int n = get_len(x0);
		solution X, X1;
		X.x = x0;
		matrix d(n, 1);
		matrix P(2, 2);
		solution h;
		double* ab; //przedzial ab
		//write(t(X.x(0)), t(X.x(1)), h0, "SD");
		while (true)
		{
			X.grad(gf); //gradient
			d = -X.g;
			if (h0 < 0) //wersja zmiennokrokowa
			{
				P(0, 0) = X.x(0);
				P(1, 0) = X.x(1);
				P(0, 1) = d(0);
				P(1, 1) = d(1);

				ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);			
				h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);	

				X1.x = X.x + h.x * d;
			}
			else
				X1.x = X.x + h0 * d;
			//write(t(X1.x(0)), t(X1.x(1)), h0, "SD");
			if (solution::f_calls > Nmax || solution::g_calls > Nmax || norm(X1.x-X.x) < epsilon)
			{
				X1.fit_fun(ff, ud1, ud2);
				return X1;
			}
			X = X1;
		}
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
		int n = get_len(x0);
		solution X, X1;
		X.x = x0;
		matrix d(n, 1);
		matrix P(2, 2);
		solution h;
		double* ab, beta;
		X.grad(gf);
		d = -X.g;
		//write(t(X.x(0)), t(X.x(1)), h0, "CG");
		while (true)
		{
			if (h0 < 0)
			{
				P(0, 0) = X.x(0);
				P(1, 0) = X.x(1);
				P(0, 1) = d(0);
				P(1, 1) = d(1);
				ab = expansion(ff,0, 1, 1.2, Nmax, ud1, P);
				h = golden(ff,ab[0], ab[1], epsilon, Nmax, ud1, P);
				X1.x = X.x + h.x * d;
			}
			else
				X1.x = X.x + h0 * d;
			//write(t(X1.x(0)), t(X1.x(1)), h0, "CG");
			if (solution::f_calls > Nmax || solution::g_calls > Nmax || norm(X1.x - X.x) < epsilon)
			{
				X1.fit_fun(ff, ud1);
				return X1;
			}
			X1.grad(gf);
			beta = pow(norm(X1.g), 2) / pow(norm(X.g), 2);
			d = -X1.g + beta * d;
			X = X1;
		}
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
		int n = get_len(x0);
		solution X, X1;
		X.x = x0;
		matrix d(n, 1);
		matrix P(2, 2);
		solution h;
		double* ab;
		//write(t(X.x(0)), t(X.x(1)), h0, "Newton");
		while (true)
		{
			X.grad(gf);
			X.hess(Hf);
			d = -inv(X.H) * X.g;
			if (h0 < 0)
			{
				P(0, 0) = X.x(0);
				P(1, 0) = X.x(1);
				P(0, 1) = d(0);
				P(1, 1) = d(1);
				ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
				h = golden(ff,ab[0], ab[1], epsilon, Nmax, ud1, P);
				X1.x = X.x + h.x * d;
			}
			else
				X1.x = X.x + h0 * d;	
			//write(t(X1.x(0)), t(X1.x(1)), h0, "Newton");
			if (solution::f_calls > Nmax || solution::g_calls > Nmax || norm(X1.x - X.x) < epsilon)
			{
				X1.fit_fun(ff,ud1);
				return X1;
			}
			X = X1;
		}
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
		double alfa = (sqrt(5) - 1) / 2;
		solution A, B, C, D;
		A.x = a;
		B.x = b;
		C.x = B.x - alfa * (B.x - A.x);
		C.fit_fun(ff, ud1, ud2);
		D.x = A.x + alfa * (B.x - A.x);
		D.fit_fun(ff, ud1, ud2);
		while (true)
		{
			if (C.y < D.y)
			{
				B = D;
				D = C;
				C.x = B.x - alfa * (B.x - A.x);
				C.fit_fun(ff, ud1, ud2);
			}
			else
			{
				A = C;
				C = D;
				D.x = A.x + alfa * (B.x - A.x);
				D.fit_fun(ff, ud1, ud2);
			}
			if (B.x - A.x<epsilon || solution::f_calls>Nmax)
			{
				A.x = (A.x + B.x) / 2;
				A.fit_fun(ff, ud1, ud2);
				return A;
			}
		}
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
		int n = get_len(x0);
		matrix D = ident_mat(n);
		matrix A(2,2);
		solution X, P, h;
		X.x = x0;
		double* ab;
		while (true)
		{	
			P = X;
			for (int i = 0; i < n; ++i)
			{
				A(0, 0) = P.x(0);
				A(1, 0) = P.x(1);
				A(0, 1) = D[i](0);
				A(1, 1) = D[i](1);

				ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, A);

				h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, A);
				P.x = P.x + h.x * D[i];
			}
			if (norm(X.x - P.x) < epsilon || solution::f_calls > Nmax)
			{
				P.fit_fun(ff, ud1);
				return P;
			}
			for (int i = 0; i < n - 1; ++i)
				D.set_col(D[i + 1], i); 
			D.set_col(P.x - X.x, n - 1); 

			A(0, 0) = P.x(0);				
			A(1, 0) = P.x(1);				
			A(0, 1) = D[n-1](0);		
			A(1, 1) = D[n-1](1);	

			ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, A);
			h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, A);
			X.x = P.x + h.x * D[n - 1];
		}
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
