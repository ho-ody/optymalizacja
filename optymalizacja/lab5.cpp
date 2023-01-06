#include "opt_alg.h"
#include "lab.h"
namespace l5 {
	matrix f5(matrix x, matrix ud1, matrix ud2) {
		matrix y(2, 1);
		if (isnan(ud2(0, 0))) {
			y = matrix(2, 1); //bo mamy dwie funkcje celu;
			/*
			x(1);
			y(0) = ud1[0]() * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2));
			y(1) = 1 / ud1[0]() * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2));
			*/
			double a = ud1[0]();
			y(0) = a * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2));
			y(1) = 1 / a * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2));
			//cerr << "Y:" << y << endl << endl;
		}
		else {
			//cerr << ud2[0] << endl << endl;
			//cerr << ud2[1] << endl << endl;
			//cerr << ud1 << endl << endl;
			//cerr << x << endl << endl;

			y = f5(ud2[0] + x * ud2[1], ud1, NAN);

//			cerr << ud1[0] << endl << endl;
//			cerr << ud1[1] << endl << endl;

			double w = ud1[1]();
			//y = ud1[1]() * y(0) + (1 - ud1[1]()) * y(1);
			y = w * y(0) + (1 - w) * y(1);

			//cerr << "Y:" << y << endl << endl;

			/*
			solution temp;
			temp.x = ad[0] + x * ad[1]; //x + alfa*d //tak jak tydzien temu
			temp.fit_fun(ud); //tu obliczamy poprostu y1,y2 z tego ifa powyzej
			y = ud[1]() * temp.y(0) + (1 - ud[1]()) * temp.y(1); //w * y1 + ...
			--f_calls;
			*/
		}


		return y;
	}

	void test_zbiez_metod() {
		double w = 0.0;
		double epsilon = 1e-5;
		int Nmax = 5000;
		int a = 1;

		for (int i = 0; i < 1; i++) {
			matrix x0 = 20 * rand_mat(2, 1) - 10;
			x0(0) = 0.; x0(1) = 0.;
			matrix ud1(1,2);
		//	ud1(0) = a;		//!!
		//	ud1(1) = w;		//!!

			ud1(0, 0) = a;
			//ud1(1, 0) = a;
			//ud1(0, 1) = w;
			ud1(0, 1) = w;

//!!			cerr << ud1[0] << endl << endl;
//!!			cerr << ud1[1] << endl << endl;

			solution powel = Powell(f5, x0, epsilon, Nmax, ud1);
			cerr << x0(0) << ";" << x0(1) << ";" << powel.x(0) << ";" << powel.x(1) << ";" << powel.y(0) << ";" << powel.y(1) << ";" << solution::f_calls << endl;
			solution::clear_calls();
			w += 0.01;
		}




	}















}