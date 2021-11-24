#pragma once
#include <iostream>
#include <iomanip>
using namespace std;

//#define N 2
//#define M 2

void transpose(double* src, double* dst, const int N, const int M);
void inverseMatrix(double invM[2][2], double M[2][2], double detM);
double detOfMatrix(double matrix[2][2]);
void printMatrix(double M[4][4]);
void printMatrix(double M[2][2]);
void printMatrix(double M[4]);
//void printMatrix(double** matrix);