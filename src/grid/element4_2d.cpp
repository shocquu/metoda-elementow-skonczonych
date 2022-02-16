#include "element4_2d.h"

Element4_2D::Element4_2D(int integralPoints) {
	this->n = integralPoints;
	this->p = n * n;
	Gauss g;

	for (int i = 0; i < n; i++) {
		w.push_back(g.w[n - 1][i]); // wagi

		for (int j = i; j < n + i; j++) {
			ksi.push_back(g.x[n - 1][j % n]); // wartoœci ksi dla n-tego punktu ca³kowania
			eta.push_back(g.x[n - 1][i % n]); // wartoœci eta dla n-tego punktu ca³kowania
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

	// @TO_FIX Na sztywno
	double points[8][2] = {
		{ ksi[0], -1 }, { ksi[0], -1 }, // pc1
		{ 1, eta[1] }, { 1, eta[1] },   // pc2
		{ ksi[2], 1 }, { ksi[2], 1 },   // pc3
		{ -1, eta[3] }, { -1, eta[3] }, // pc4
	};
	Side sides[4] = {
		{ points[0], points[1] }, // dó³
		{ points[2], points[3] }, // prawo
		{ points[4], points[5] }, // góra
		{ points[6], points[7] }, // lewo
	};

	// Inicjalizacja macierzy funkcji N dla nowych wartoœci ksi i eta
	for (int i = 0; i < p; i++) {
		std::array<double, 4> Nrow1 = nKsiEta(sides[i].pc1[0], sides[i].pc1[1]);
		std::array<double, 4> Nrow2 = nKsiEta(sides[i].pc2[0], sides[i].pc2[1]);

		for (size_t j = 0; j < p; j++) {
			Npc1[i][j] = Nrow1[j];
			Npc2[i][j] = Nrow2[j];
		}
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
