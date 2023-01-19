/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
*********************************************/

#include "opt_alg.h"
#include "lab.h"

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
		lab6();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void testowa_f_celu_a();
void problem_rzeczy_b();
void test_zbiez_metod();
void lab1() {
	//l1::test_zbiez_metod();	//spradzanie poprawnosci obliczen metody fibbonaciego i lagrangea
	//l1::testowa_f_celu_a();	//rozwi¹zanie dla funkcji testowej, punkt 5.a. z pdf'a
	l1::problem_rzeczy_b();		//rozwi¹zanie dla przypadku rzeczywistego, punkt 5.b. z pdf'a
}

void lab2() {
	//l2::test_zbiez_metod();
	//l2::testowa_f_celu_a();
	//l2::testowa_f_celu_b();
	l2::problem_rzeczywi();
}

void lab3() {
	//l3::test_zbiez_metod();
	l3::testowa_f_celu_a();
	//l3::problem_rzeczywi();
}

void lab4() {
	l4::test_zbiez_metod();
	//l4::testowa_f_celu_a();
	//l4::testowa_f_celu_b();
	//l4::problem_rzeczywi();
}

void lab5(){
	//l5::testowa_f_celu_a();
	l5::problem_rzeczywi();
}

void lab6() {
	l6::testowa_f_celu_a();
}
