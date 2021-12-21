#include <math.h>
#include "../utils/matrix.h"
#include "element4_2d.h"
#include "grid.h"

Element4_2D::Element4_2D(double* ksi, double* eta, int n = 2) {
	this->p = n * n;
	this->ksi = ksi;
	this->eta = eta;
	nMatrix = new double* [p];
	ksiMatrix = new double* [p];
	etaMatrix = new double* [p];

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
 * Ustawia w przekazanej tablicy sumê macierzy C dla ka¿dego punktu ca³kowania.
 *
 * @param globalH - tablica do zapisu zagregowanej macierzy H
 * @param globalC - tablica do zapisu zagregowanej macierzy C
 * @param globalP - tablica do zapisu zagregowanej macierzy P
 * @param currEl - element, na którym metoda ma operowaæ
 */
void Element4_2D::aggregate(double**& globalH, double**& globalC, double*& globalP, Element currEl) {
	int nodesId[4] = { currEl.id[0] - 1, currEl.id[1] - 1, currEl.id[2] - 1, currEl.id[3] - 1 };
	int sideNodes[4][4][2] = {									// pomocnicza macierz ID wêz³ów do rozmieszczenia elementów
		{ { nodesId[0], nodesId[0] }, { nodesId[0], nodesId[1] }, { nodesId[0], nodesId[2] }, { nodesId[0], nodesId[3] } },
		{ { nodesId[1], nodesId[0] }, { nodesId[1], nodesId[1] }, { nodesId[1], nodesId[2] }, { nodesId[1], nodesId[3] } },
		{ { nodesId[2], nodesId[0] }, { nodesId[2], nodesId[1] }, { nodesId[2], nodesId[2] }, { nodesId[2], nodesId[3] } },
		{ { nodesId[3], nodesId[0] }, { nodesId[3], nodesId[1] }, { nodesId[3], nodesId[2] }, { nodesId[3], nodesId[3] } },
	};

	for (int i = 0; i < this->p; i++) {
		for (int j = 0; j < this->p; j++) {
			int rowIndex = sideNodes[i][j][0];
			int colIndex = sideNodes[i][j][1];
			globalH[rowIndex][colIndex] += currEl.H[i][j] + currEl.Hbc[i][j];
			globalC[rowIndex][colIndex] += currEl.C[i][j];
		}

		int rowIndex = nodesId[i];
		globalP[rowIndex] += currEl.P[i];		
	}
}

/**
 * Ustawia w przekazanej tablicy sumê macierzy Hbc oraz wektor P dla ka¿dego punktu ca³kowania.
 *
 * @param Hbc - tablica, do której zapisaæ sumê macierzy Hbc z ka¿dego punktu ca³kowania
 * @param P - tablica, do której zapisaæ wektor P
 * @param alpha - ...
 * @param ambientTemp - temperatura otoczenia
 */
void Element4_2D::calcHbc(double Hbc[4][4], double P[4], Grid grid, Element currEl, int elIndex, double alpha = 25, double ambientTemp = 1200) {
	struct Side { double* pc1, * pc2; };
	double Npc1[4][4] = { 0 }, Npc2[4][4] = { 0 };
	double w[2] = { 1, 1 };
	double ksi, eta;
	ksi = eta = 1 / sqrt(3);

	double points[8][2] = {
		{ -ksi, -1 }, { ksi, -1 }, // pc1
		{ 1, -eta }, { 1, eta },   // pc2
		{ ksi, 1 }, { -ksi, 1 },   // pc3
		{ -1, eta }, { -1, -eta }, // pc4
	};
	Side sides[4] = {
		{ points[0], points[1] }, // dó³
		{ points[2], points[3] }, // prawo
		{ points[4], points[5] }, // góra
		{ points[6], points[7] }, // lewo
	};

	/// @TODO - mo¿na policzyæ raz dla wszystkich œcian
	// Inicjalizacja macierzy funkcji N dla nowych wartoœci ksi i eta
	for (int i = 0; i < this->p; i++) {
		double* Nrow = this->nKsiEta(sides[i].pc1[0], sides[i].pc1[1]);
		Npc1[i][0] = Nrow[0];
		Npc1[i][1] = Nrow[1];
		Npc1[i][2] = Nrow[2];
		Npc1[i][3] = Nrow[3];
		Nrow = this->nKsiEta(sides[i].pc2[0], sides[i].pc2[1]);
		Npc2[i][0] = Nrow[0];
		Npc2[i][1] = Nrow[1];
		Npc2[i][2] = Nrow[2];
		Npc2[i][3] = Nrow[3];
	}

	// Wype³nianie macierzy
	for (int i = 0; i < this->p; i++) { // i - œciana elementu
		int n1Index = currEl.id[i] - 1;
		int n2Index = currEl.id[(i + 1) % 4] - 1;

		// Sprawdzanie warunku brzegowego na œcianach elementu
		if (grid.nodes[n1Index].bc > 0 && grid.nodes[n2Index].bc > 0) {
			double L = distance(grid.nodes[n1Index].x, grid.nodes[n1Index].y, grid.nodes[n2Index].x, grid.nodes[n2Index].y);
			double locDetJ = L / 2;

			for (int j = 0; j < this->p; j++) {
				for (int k = 0; k < this->p; k++) {
					Hbc[j][k] += w[0] * (Npc1[i][j] * Npc1[i][k]) * alpha * locDetJ;
					Hbc[j][k] += w[1] * (Npc2[i][j] * Npc2[i][k]) * alpha * locDetJ;
				}

				P[j] += (w[0] * Npc1[i][j] + w[1] * Npc2[i][j]) * alpha * ambientTemp * locDetJ;
			}
		}
	}
}

/**
 * Tworzenie macierzy H dla ka¿dego punktu ca³kowania. Wymaga poprzedzenia funkcj¹
 * `fillMatrices`, by operowaæ na uzupe³nionych macierzach.
 *
 * @param H - macierz, do której zapisaæ wynik
 * @param detJ - wyznacznik macierzy J
 * @param k - wspó³czynnik przewodzenia
 */
void Element4_2D::calcH(double H[4][4], double k = 30) {
	for (int i = 0; i < this->p; i++) {
		double rowsX[4] = { dNdX[i][0],	dNdX[i][1], dNdX[i][2],	dNdX[i][3] };
		double rowsY[4] = { dNdY[i][0], dNdY[i][1], dNdY[i][2], dNdY[i][3] };

		for (int j = 0; j < this->p; j++)
			for (int l = 0; l < this->p; l++)
				H[j][l] += k * (rowsX[j] * rowsX[l] + rowsY[j] * rowsY[l]) * detJ;
	}
}

/**
 * Ustawia w przekazanej tablicy sumê macierzy C (pojemnoœci cieplnej) dla ka¿dego punktu ca³kowania.
 *
 * @param c - ciep³o w³aœciwe
 * @param ro - gêstoœæ materia³u
 */
void Element4_2D::calcC(double C[4][4], double c, double ro, Element currEl, int currElId, Grid grid) {	
	for (int i = 0; i < this->p; i++) {
		double* N = nMatrix[i];

		for (int j = 0; j < this->p; j++)
			for (int k = 0; k < this->p; k++)
				C[j][k] += c * ro * (N[j] * N[k]) * detJ;
	}
}

/**
 * Oblicza macierz Jakobiego.
 *
 * @param J - macierz Ÿród³owa do zapisania danych
 * @param el4 - struktura z metodami do liczenia pochodnych
 * @param pc - punkt ca³kowania
 * @param x, y - wspó³rzêdna x i y
 */
void Element4_2D::calcJ(double J[2][2], int pc, double x, double y) {
	J[0][0] += this->dXYdKsi(pc, x);
	J[0][1] += this->dXYdEta(pc, x);
	J[1][0] += this->dXYdKsi(pc, y);
	J[1][1] += this->dXYdEta(pc, y);
}

/**
 * Uzupe³nia Jakobian dla ka¿dego punktu ca³kowania i-tego elementu.
 *
 * @param i - punkt ca³kowania
 * @param el4 - struktura z danymi o elemencie
 * @param grid - siatka, na której wykonywana jest procedura liczenia Jakobianu
 */
void Element4_2D::jacobian(Grid grid, Element &currEl) {
	// Obliczanie macierzy Jacobiego
	for (int i = 0; i < this->p; i++)	{
		int nodeId = currEl.id[i] - 1;
		double x = grid.nodes[nodeId].x;
		double y = grid.nodes[nodeId].y;

		this->calcJ(currEl.J, i, x, y);
	}

	detJ = currEl.detJ = detOfMatrix(currEl.J);
	inverseMatrix(currEl.invJ, currEl.J, currEl.detJ);

	// Wype³nij tablicê pochodnych funkcji N po X i Y - jakobian
	for (int j = 0; j < this->p; j++) {
		for (int k = 0; k < this->p; k++) {
			dNdX[j][k] = currEl.invJ[0][0] * this->ksiMatrix[j][k] + currEl.invJ[0][1] * this->etaMatrix[j][k];
			dNdY[j][k] = currEl.invJ[1][0] * this->ksiMatrix[j][k] + currEl.invJ[1][1] * this->etaMatrix[j][k];
		}
	}
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
