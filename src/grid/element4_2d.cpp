#include "element4_2d.h"

Element4_2D::Element4_2D(int integralPoints) {
	this->p = integralPoints * integralPoints;
	this->ksi = new double[p];
	this->eta = new double[p];
	Gauss g;

	int it = 0;
	for (int i = 0; i < integralPoints; i++) {
		for (int j = i; j < integralPoints + i; j++) {
			ksi[it] = g.x[integralPoints - 1][j % integralPoints]; // wartoœci ksi dla n-tego punktu ca³kowania
			eta[it] = g.x[integralPoints - 1][i % integralPoints]; // wartoœci eta dla n-tego punktu ca³kowania
			it++;
		}		
	}

	this->nMatrix = new double* [p];
	this->ksiMatrix = new double* [p];
	this->etaMatrix = new double* [p];

	for (int i = 0; i < p; i++) {
		nMatrix[i] = new double[p];
		ksiMatrix[i] = new double[p];
		etaMatrix[i] = new double[p];
	}

	fillMatrices();
}

Element4_2D::~Element4_2D() {
	for (int i = 0; i < p; i++) {
		delete etaMatrix[i];
		delete ksiMatrix[i];
	}
	delete[] etaMatrix, ksiMatrix;
}

/**
 * Zwraca odleg³oœæ miêdzy dwoma punktami.
 *
 * @param x1, y1 - wspó³rzêdna x i y pierwszego punktu
 * @param x2, y2 - wspó³rzêdna x i y drugiego punktu
 */
double Element4_2D::distance(double x1, double y1, double x2, double y2) {
	return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

/**
 * Zwraca tablicê funkcji kszta³tu.
 *
 * @param ksi - wartoœci ksi dla funkcji
 * @param eta - wartoœci eta dla funkcji
 */
double* Element4_2D::nKsiEta(double ksi, double eta) {
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
double* Element4_2D::dNdEta(double ksi) {
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
double* Element4_2D::dNdKsi(double eta) {
	dKsi[0] = -0.25 * (1 - eta);
	dKsi[1] = 0.25 * (1 - eta);
	dKsi[2] = 0.25 * (1 + eta);
	dKsi[3] = -0.25 * (1 + eta);

	return dKsi;
}

/**
 * Zwraca Jakobian przekszta³cenia po eta.
 *
 * @param pc - punkt ca³kowania
 * @param xy - wymno¿enie po "x" lub po "y"
 */
double Element4_2D::dXYdEta(int pc, double xy) {
	return dEta[pc] * xy;
}

/**
 * Zwraca sumê Jakobianu przekszta³cenia po eta.
 *
 * @param xy - wymno¿enie po "x" lub po "y"
 */
double Element4_2D::dXYdEta(double* xy) {
	double sum = 0;
	for (size_t i = 0; i < 4; i++) sum += dEta[i] * xy[i];
	return sum;
}

/**
 * Zwraca Jakobian przekszta³cenia po ksi.
 *
 * @param pc - punkt ca³kowania
 * @param xy - wymno¿enie po "x" lub po "y"
 */
double Element4_2D::dXYdKsi(int i, double xy) {
	return dKsi[i] * xy;
}

/**
 * Zwraca sumê Jakobianu przekszta³cenia po ksi.
 *
 * @param xy - wymno¿enie po "x" lub po "y"
 */
double Element4_2D::dXYdKsi(double* xy) {
	double sum = 0;
	for (size_t i = 0; i < 4; i++) sum += dKsi[i] * xy[i];
	return sum;
}

/**
 * Wype³nia macierz pochodnych funkcji N po eta i ksi.
 *
 * @param xy - wymno¿enie po "x" lub po "y"
 */
void Element4_2D::fillMatrices() {
	for (int x = 0; x < p; x++)	
		for (int y = 0; y < 4; y++)	{
			nMatrix[x][y] = nKsiEta(ksi[x], eta[x])[y];
			etaMatrix[x][y] = dNdEta(ksi[x])[y];
			ksiMatrix[x][y] = dNdKsi(eta[x])[y];
		}
	
}
