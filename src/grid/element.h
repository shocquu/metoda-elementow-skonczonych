#pragma once

/**
 * Element uniwersalny.
 *
 * @param id[4] - identyfikatory wêz³ów elementu
 * @param H[4][4] - macierz ..
 * @param Hbc[4][4] - macierz ..
 * @param C[4][4] - macierz ciep³a w³aœciwego
 * @param P[16] - wektor obci¹¿eñ / warunków brzegowych
 */
struct Element {
	int index = 0;
	int id[4] = { 0, 0, 0, 0 };
	double J[2][2] = { 0 }, invJ[2][2] = { 0 };
	double H[4][4] = { 0 }, Hbc[4][4] = { 0 }, C[4][4] = { 0 };
	double P[4] = { 0 };
	double detJ = 0;
	double temperature = 0;
};