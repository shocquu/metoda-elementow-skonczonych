#pragma once

typedef double (*vFunctionCall)(double arg);
typedef double (*vFunctionCall2)(double arg1, double arg2);

// f(x) = 5x^2 + 3x + 6
double f(double x);

// f(x, y) = 5x^2*y^2 + 3xy + 6
double f(double x, double y);
