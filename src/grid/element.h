#pragma once

/**
 * Element uniwersalny.
 *
 * @param id[4] - identyfikatory w�z��w elementu
 * @param H[4][4] - macierz ..
 * @param Hbc[4][4] - macierz ..
 * @param C[4][4] - macierz ciep�a w�a�ciwego
 * @param P[16] - wektor obci��e� / warunk�w brzegowych
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