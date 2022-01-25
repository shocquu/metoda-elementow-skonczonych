#include "element4_2d.h"

Element4_2D::Element4_2D(int integralPoints) {
	this->p = integralPoints * integralPoints;
	Gauss g;

	for (int i = 0; i < integralPoints; i++) {
		for (int j = i; j < integralPoints + i; j++) {
			ksi.push_back(g.x[integralPoints - 1][j % integralPoints]); // wartoœci ksi dla n-tego punktu ca³kowania
			eta.push_back(g.x[integralPoints - 1][i % integralPoints]); // wartoœci eta dla n-tego punktu ca³kowania
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
 * Zwraca tablicê funkcji kszta³tu.
 *
 * @param ksi - wartoœci ksi dla funkcji
 * @param eta - wartoœci eta dla funkcji
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
 * Zwraca tablicê pochodnych funkcji kszta³tu po eta.
 *
 * @param ksi - wartoœci ksi dla funkcji
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
 * Zwraca tablicê pochodnych funkcji kszta³tu po ksi.
 *
 * @param eta - wartoœci eta dla funkcji
 */
std::array<double, 4> Element4_2D::dNdKsi(double eta) {
	std::array<double, 4> dKsi = { 0 };
	dKsi[0] = -0.25 * (1 - eta);
	dKsi[1] = 0.25 * (1 - eta);
	dKsi[2] = 0.25 * (1 + eta);
	dKsi[3] = -0.25 * (1 + eta);
	return dKsi;
}
