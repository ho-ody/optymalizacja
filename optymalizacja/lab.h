#pragma once
#define M_PI 3.14159265358979323846
#define M_E 2.71828182845904523536

#include "user_funs.h"
namespace l1 {
	void test_zbiez_metod();	//spradzanie poprawnosci obliczen metody fibbonaciego i lagrangea
	void testowa_f_celu_a();	//rozwi¹zanie dla funkcji testowej, punkt 5.a. z pdf'a
	void problem_rzeczy_b();		//rozwi¹zanie dla przypadku rzeczywistego, punkt 5.b. z pdf'a
}
namespace l2 {
	void test_zbiez_metod();
	void testowa_f_celu_a();
	void testowa_f_celu_b();
	void problem_rzeczywi();
}
namespace l3 {
	void test_zbiez_metod();
	void testowa_f_celu_a();
	void problem_rzeczywi();
}
namespace l4 {
	void test_zbiez_metod();
	void testowa_f_celu_a();
	void testowa_f_celu_b();
	void problem_rzeczywi();
}
namespace l5 {
	void testowa_f_celu_a();
	//void testowa_f_celu_a();
	//void testowa_f_celu_b();
	void problem_rzeczywi();
}