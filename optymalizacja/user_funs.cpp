#include"user_funs.h"

double cacl_nth_element_fib(int n)
{
	return round(pow((1 + sqrt(5)) / 2, n) / sqrt(5));
}