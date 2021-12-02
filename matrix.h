#pragma once
#include <iostream>
#include <iomanip>
using namespace std;

//#define N 4
//#define M 2

void transpose(double* src, double* dst, const int M, const int N);
void inverseMatrix(double invM[2][2], double M[2][2], double detM);
double detOfMatrix(double matrix[2][2]);
void printMatrix(double M[16][16]);
void printMatrix(double M[4][4]);
void printMatrix(double M[2][2]);
void printMatrix(double M[4]);
void printMatrix(double** matrix, const int M, const int N);
void printMatrix4x4(double** matrix);