#include "element4_2d.h"

Element4_2D::Element4_2D(int integralPoints) {
	this->p = integralPoints * integralPoints;
	Gauss g;

	for (int i = 0; i < integralPoints; i++) {
		for (int j = i; j < integralPoints + i; j++) {
			ksi.push_back(g.x[integralPoints - 1][j % integralPoints]); // warto�ci ksi dla n-tego punktu ca�kowania
			eta.push_back(g.x[integralPoints - 1][i % integralPoints]); // warto�ci eta dla n-tego punktu ca�kowania
		}
	}

	auto k = ksi.begin();
	auto e = eta.begin();
	while (k != ksi.end() and e != eta.end()) {
		nMatrix.push_back(nKsiEta(*k, *e));
		dNdKsiMatrix.push_back(dNdKsi(*e));
		dNdEtaMatrix.push_back(dNdEta(*k));
		*k++, *e++;
	}
}

/**
 * Zwraca tablic� funkcji kszta�tu.
 *
 * @param ksi - warto�ci ksi dla funkcji
 * @param eta - warto�ci eta dla funkcji
 */

std::array<double, 4> Element4_2D::nKsiEta(double ksi, double eta) {
	std::array<double, 4> N = { 0 };
	N[0] = 0.25 * (1 - ksi) * (1 - eta);
	N[1] = 0.25 * (1 + ksi) * (1 - eta);
	N[2] = 0.25 * (1 + ksi) * (1 + eta);
	N[3] = 0.25 * (1 - ksi) * (1 + eta);
	return N;
}

/**
 * Zwraca tablic� pochodnych funkcji kszta�tu po eta.
 *
 * @param ksi - warto�ci ksi dla funkcji
 */
std::array<double, 4> Element4_2D::dNdEta(double ksi) {
	std::array<double, 4> dEta = { 0 };
	dEta[0] = -0.25 * (1 - ksi);
	dEta[1] = -0.25 * (1 + ksi);
	dEta[2] = 0.25 * (1 + ksi);
	dEta[3] = 0.25 * (1 - ksi);
	return dEta;
}

/**
 * Zwraca tablic� pochodnych funkcji kszta�tu po ksi.
 *
 * @param eta - warto�ci eta dla funkcji
 */
std::array<double, 4> Element4_2D::dNdKsi(double eta) {
	std::array<double, 4> dKsi = { 0 };
	dKsi[0] = -0.25 * (1 - eta);
	dKsi[1] = 0.25 * (1 - eta);
	dKsi[2] = 0.25 * (1 + eta);
	dKsi[3] = -0.25 * (1 + eta);
	return dKsi;
}
