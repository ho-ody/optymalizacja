#include "opt_alg.h"
#include "lab.h"
namespace l5 {
	matrix f5(matrix x, matrix ud1, matrix ud2) {
		matrix y(2, 1);
		if (isnan(ud2(0, 0))) {
			y = matrix(2, 1);

			double a = ud1[0]();
			y(0) = a * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2));
			y(1) = 1 / a * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2));
		}
		else {
			y = f5(ud2[0] + x * ud2[1], ud1, NAN);

			double w = ud1[1]();
			y = w * y(0) + (1 - w) * y(1);
		}
		return y;
	}

	matrix fR(matrix x, matrix ud1, matrix ud2) {
		matrix y;
		if (isnan(ud2(0, 0))) {
			y = matrix(3, 1);
			double ro = 7800, P = 1e3, E = 207e9;
			y(0) = ro * x(0) * M_PI * pow(x(1), 2) / 4;
			y(1) = 64 * P * pow(x(0), 3) / (3 * E * M_PI * pow(x(1), 4));
			y(2) = 32 * P * x(0) / (M_PI * pow(x(1), 3));
		}
		else {
			solution T;
			T = ud2[0] + x * ud2[1];
			T.y = fR(T.x, ud1, NAN);

			//normalizacja
			double m_min = 0.061261056745001,	   m_max = 15.3152641862502;
			double o_min = 0.00000524878137604189, o_max = 0.00328048836002618;
			matrix yn(2, 1);
			yn(0) = (T.y(0) - m_min) / (m_max - m_min);
			yn(1) = (T.y(1) - o_min) / (o_max - o_min);

			y = ud1() * yn(0) + (1 - ud1()) * yn(1);

			double c = 1e10;
			if (T.x(0) < 0.1)	//l < 100mm
				y = y + c * pow(0.1 - T.x(0), 2);
			if (T.x(0) > 1)		//l > 1000mm
				y = y + c * pow(T.x(0) - 1, 2);
			if (T.x(1) < 0.01)	//d < 10mm
				y = y + c * pow(0.01 - T.x(1), 2);
			if (T.x(1) > 0.05)	//d > 50mm
				y = y + c * pow(T.x(1) - 0.05, 2);
			if (T.y(1) > 0.005)	//ugiecie > 5mm
				y = y + c * pow(T.y(1) - 0.005, 2);
			if (T.y(2) > 300e6)	//naprezenie > 300MPA
				y = y + c * pow(T.y(2) - 300e6, 2);
		}
		return y;
	}

	void testowa_f_celu_a() {
		ofstream output("_output.csv");
		double w = 0.0;
		double epsilon = 1e-5;
		int Nmax = 1000;
		int a = 1;

		for (int i = 0; i < 101; i++) {
			matrix x0 = 20 * rand_mat(2, 1) - 10;
			output << t(w) << t(x0(0)) << t(x0(1));
			for (int j = 0; j < 3; j++) {
				switch (j) {
				case 0: a = 1; break;
				case 1: a = 10; break;
				case 2: a = 100; break;
				}
				matrix ud1(1, 2);
				ud1(0, 0) = a;
				ud1(0, 1) = w;

				solution powel = Powell(f5, x0, epsilon, Nmax, ud1);
				output << t(powel.x(0)) << t(powel.x(1)) << t(powel.y(0)) << t(powel.y(1)) << t(solution::f_calls);
				solution::clear_calls();		
			}
			w += 0.01;
			output << endl;
		}
	}

	void problem_rzeczywi() {
		ofstream output_real("_output_real.csv");
		matrix ld0(2, 1);
		double w = 0.0;
		double epsilon = 1e-5;
		int Nmax = 1000;
		
		for (int i = 0; i < 101; i++) {
			matrix ud1(1, 1);
			ud1(0, 0) = w;
			
			ld0 = rand_mat(2, 1);
			ld0(0) = ld0(0) * (1. - 0.1) + 0.1;
			ld0(1) = ld0(1) * (0.05 - 0.01) + 0.01;
			
			solution solR = Powell(fR, ld0, epsilon, Nmax, ud1);

			output_real << t(w) << t(ld0(0)) << t(ld0(1)) << t(solR.x(0)) << t(solR.x(1)) << t(solR.y(0)) << t(solR.y(1)) << t(solution::f_calls) << endl;
			solution::clear_calls();
			w += 0.01;
		}
	}

}