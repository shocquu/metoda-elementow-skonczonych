#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include "functions.h"
#include "matrix.h"
#include "gauss.h"

using namespace std;

#define LAB_NO 10

//#define SHOW_MATRIX
//#define SHOW_DETAILS
//#define TEST

/**
 * Wêze³ bêd¹cy czêœci¹ elementu.
 * 
 * @param x, y - wspó³rzêdne wêz³ów
 * @param bc - informacja o obecnoœci warunku brzegowego
 */
struct Node {
	int id = -1;
	short int bc = 0;
	double x = 0, y = 0;
};

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
	int id[4] = { 0, 0, 0, 0 };
	double H[4][4] = { 0 }, C[4][4] = { 0 }, Hbc[4][4] = { 0 };
	double P[16] = { 0 };	
};

/**
 * Siatka elementów.
 *
 * @param width, height - wymiary siatki
 * @param nW, nH - liczba wêz³ów odpowiednio na osi X i osi Y
 * @param nN - ³¹czna liczba wêz³ów w siatce
 * @param nE - liczba elementów, z których sk³ada siê siatka
 */
struct Grid {
	double width, height;
	int nW, nH, nN, nE;
	Node* nodes;
	Element* elements;

	Grid() {
		this->width = 0;
		this->height = 0;
		this->nW = 0;
		this->nH = 0;
		this->nN = 0;
		this->nE = 0;
		this->nodes = NULL;
		this->elements = NULL;
	}
	Grid(double height, double width, int nH, int nW) {
		this->height = height;
		this->width = width;
		this->nH = nH;
		this->nW = nW;
		this->nN = nH * nW;
		this->nE = (nW - 1) * (nH - 1);
		this->nodes = new Node[nN];
		this->elements = new Element[nE];

		setNodesCoords();
		setElementsNodes();
		setNodesBC();
	}

	~Grid() {
		// ERROR
		//delete[] this->nodes;
		//delete[] this->elements;
	}

	void setNodesCoords() {
		int i = 0;

		for (int x = 0; x < nW; x++) {
			for (int y = 0; y < nH; y++) {
				nodes[i].id = i;
				++i;
				nodes[i].x = nodes[i - 1].x;
				nodes[i].y = nodes[i - 1].y + height / (static_cast<double>(nH) - 1);
			}

			nodes[i].x = nodes[i - 1].x + width / (static_cast<double>(nW) - 1);
			nodes[i].y = 0;
		}
	}

	void printNodesCoords() {
		for (int i = 0; i < nN; i++)
			std::cout << "Node #" << i + 1 << ":\t[" << this->nodes[i].x << ", " << this->nodes[i].y << "]\n";
	}

	void setElementsNodes() {
		int n = 1, nhSum = nH;

		for (int i = 0; i < nE; i++) {
			elements[i].id[0] = n;
			elements[i].id[1] = elements[i].id[0] + nH;
			elements[i].id[2] = elements[i].id[1] + 1;
			elements[i].id[3] = elements[i].id[0] + 1;	
			n++;

			if (n == nhSum) {
				nhSum += nH;
				n++;
			}
		}
	}

	void printElementsNodes(string type = "inline") {
		for (int i = 0; i < nE; i++) {
			if (type == "inline") {
				std::cout << "El #" << i + 1 << "\t{" << elements[i].id[0] << ", " << elements[i].id[1] << ", " << elements[i].id[2] << ", " << elements[i].id[3] << "}\n";
			} else {
				std::cout << "(" << elements[i].id[3] << ")---(" << elements[i].id[2] << ")\n";
				std::cout << " |\t|\n";
				std::cout << "(" << elements[i].id[0] << ")---(" << elements[i].id[1] << ")\n\n";
			}			
		}
	}

	void setNodesBC() {
		for (int i = 0; i < nN; i++) {
			if (nodes[i].x == 0 || nodes[i].x == width)  nodes[i].bc = 1;
			if (nodes[i].y == 0 || nodes[i].y == height) nodes[i].bc = 1;
		}
	}

	void printBCNodes() {
		std::cout << "\nBC nodes:\n";
		for (int i = 0; i < nN; i++)
			if(nodes[i].bc == 1)
				std::cout << i + 1 << "  ";
		
		std::cout << "\n";
	}
};

/**
 * Dwuwymiarowy, czteropunktowy element.
 * 
 * @param N[4] - tablica funkcji kszta³tu
 * @param dKsi[4], dEta[4] - tablica pochodnych funkcji N kolejno po Ksi i Eta
 * @param ksiMatrix, etaMatrix - tablice pochodnych funkcji dla wszystkich punktów ca³kowania
 * @param ksi, eta - schematy ca³kowania
 * @param p - liczba punktów ca³kowania
 */
struct Element4_2D {
	double dNdX[4][4] = { 0 }, dNdY[4][4] = { 0 };
	double dNdX_T[1][4] = { 0 }, dNdY_T[1][4] = { 0 };
	double J[2][2] = { 0 }, invJ[2][2] = { 0 };
	double N[4], dKsi[4], dEta[4];
	double **ksiMatrix, **etaMatrix;
	double *ksi, *eta;
	double detJ = 0;
	int p = 4;

	Element4_2D(double* ksi, double* eta, int n = 2) {
		this->p = n * n;
		this->ksi = ksi;
		this->eta = eta;
		ksiMatrix = new double* [p];
		etaMatrix = new double* [p];

		for (int i = 0; i < p; i++) {
			ksiMatrix[i] = new double[p];
			etaMatrix[i] = new double[p];
		}

		fillMatrices();
	}

	~Element4_2D() {
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
	double distance(double x1, double y1, double x2, double y2)	{
		return sqrt( pow(x2 - x1, 2) + pow(y2 - y1, 2) );
	}

	/**
	 * Ustawia w przekazanej tablicy sumê macierzy C dla ka¿dego punktu ca³kowania.
	 *
	 * @param globalH - tablica do zapisu zagregowanej macierzy H
	 * @param globalC - tablica do zapisu zagregowanej macierzy C
	 * @param globalP - tablica do zapisu zagregowanej macierzy P
	 * @param currEl - element, na którym metoda ma operowaæ
	 */
	void aggregate(double**& globalH, double**& globalC, double*& globalP, Element currEl) {
		int nodes[4] = { currEl.id[0], currEl.id[1], currEl.id[2], currEl.id[3] };
		int tempMatrix[4][4][2] = {							 // pomocnicza macierz do rozmieszczenia elementów
			{ { nodes[0], nodes[0] }, { nodes[0], nodes[1] }, { nodes[0], nodes[2] }, { nodes[0], nodes[3] } },
			{ { nodes[1], nodes[0] }, { nodes[1], nodes[1] }, { nodes[1], nodes[2] }, { nodes[1], nodes[3] } },
			{ { nodes[2], nodes[0] }, { nodes[2], nodes[1] }, { nodes[2], nodes[2] }, { nodes[2], nodes[3] } },
			{ { nodes[3], nodes[0] }, { nodes[3], nodes[1] }, { nodes[3], nodes[2] }, { nodes[3], nodes[3] } },
		};

		for (int i = 0; i < this->p; i++) {
			for (int j = 0; j < this->p; j++) {
				int rowIndex = tempMatrix[i][j][0] - 1;
				int colIndex = tempMatrix[i][j][1] - 1;
				globalH[rowIndex][colIndex] += currEl.H[i][j]; // H czy Hbc?
				globalC[rowIndex][colIndex] += currEl.C[i][j];
			}

			int rowIndex = nodes[i] - 1;
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
	void calcHbc(double Hbc[4][4], double P[4], Grid grid, Element currEl, int elIndex, double alpha = 25, double ambientTemp = 1200) {
		struct Side { double *pc1, *pc2; };
		double Npc1[4][4] = { 0 }, Npc2[4][4] = {0};
		double w[2] = { 1, 1 };
		double ksi, eta;
		ksi = eta = 1/sqrt(3);
		
		double points[8][2] = {
			{ -ksi, -1 }, { ksi, -1 }, //pc1
			{ 1, -eta }, { 1, eta },   //pc2
			{ ksi, 1 }, { -ksi, 1 },   //pc3
			{ -1, eta }, { -1, -eta }, //pc4
		};
		Side sides[4] = {
			{ points[0], points[1] }, //dol
			{ points[2], points[3] }, //prawo
			{ points[4], points[5] }, //gora
			{ points[6], points[7] }, //lewo
		};

		#ifdef TEST
			double x1 = sides[3].pc1[0];
			double y1 = sides[3].pc1[1];
			double x2 = sides[3].pc2[0];
			double y2 = sides[3].pc2[1];
			//double L = (y2 - 0.5773) / (y2 - y1);
			double L = (1 - y2) / 2;
			double localDetJ = L / 2;

			std::cout << localDetJ << std::endl;
		#endif

		// DO OPTYMALIZACJI - mo¿na policzyæ raz dla wszystkich œcian
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
		for (int i = 0; i < this->p; i++) { // i - punkt ca³kowania
			int n1Index = currEl.id[i] - 1;
			int n2Index = currEl.id[(i + 1) % 4] - 1; 
			double L = distance(grid.nodes[n1Index].x, grid.nodes[n1Index].y, grid.nodes[n2Index].x, grid.nodes[n2Index].y);
			double locDetJ = L / 2;

			// Sprawdzanie warunku brzegowego na œcianach elementu
			if ( grid.nodes[n1Index].bc > 0 && grid.nodes[n2Index].bc > 0 ) {
				//cout << "Element " << elIndex + 1 << " | Sciana " << i + 1 << " | Punkty ==> ";
				//cout << grid.nodes[n1Index].id + 1 << ", " << grid.nodes[n2Index].id + 1 << "\n";

				for (int j = 0; j < this->p; j++) {
					for (int k = 0; k < this->p; k++) {
						Hbc[j][k] = w[0] * (Npc1[i][j] * Npc1[i][k]);
						Hbc[j][k] += w[1] * (Npc2[i][j] * Npc2[i][k]);
						Hbc[j][k] *= alpha * locDetJ;
					}

					P[j] = w[0] * Npc1[i][j];
					P[j] += w[1] * Npc2[i][j];
					P[j] *= alpha * ambientTemp * locDetJ;
				}

				#ifdef SHOW_MATRIX
					std::cout << "Pow_" << i + 1 << "\n";
					cout << "DetJ: " << locDetJ << "\n";
					printMatrix(Hbc);
				#endif
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
	void calcH(double H[4][4], int k = 30) {
		for (int i = 0; i < this->p; i++) {
			double rowsX[4] = { dNdX[i][0],	dNdX[i][1], dNdX[i][2],	dNdX[i][3] };
			double rowsY[4] = { dNdY[i][0], dNdY[i][1], dNdY[i][2], dNdY[i][3] };

#ifdef SHOW_DETAILS
			std::cout << k << " * (" << rowsX[0] << " * " << dNdX_T[0][0] << " + " << rowsY[0] << " * " << dNdY_T[0][0] << ") * " << detJ << " = ";
			std::cout << k << " * (" << rowsX[0] << " * " << rowsX[0] << " + " << rowsY[0] << " * " << rowsY[0] << ") * " << detJ << " = ";
			std::cout << k * (rowsX[0] * rowsX[0] + rowsY[0] * rowsY[0]) * detJ << "\n";
#endif

			for (int j = 0; j < this->p; j++)
				for (int l = 0; l < this->p; l++)
					H[j][l] += k * (rowsX[j] * rowsX[l] + rowsY[j] * rowsY[l]) * detJ;
		}
#ifdef SHOW_DETAILS
		std::cout << "\n";
#endif
#ifdef SHOW_MATRIX
		printMatrix(H);
#endif
	}

	/**
	 * Ustawia w przekazanej tablicy sumê macierzy C (pojemnoœci cieplnej) dla ka¿dego punktu ca³kowania.
	 *
	 * @param c - ciep³o w³aœciwe
	 * @param ro - gêstoœæ materia³u
	 */
	void calcC(double C[4][4], double c = 700, double ro = 7800) {
		for (int i = 0; i < this->p; i++) {
			double* N = nKsiEta(ksi[i], eta[i]);

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
	void calcJ(double J[2][2], int pc, double x, double y) {
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
	void jacobian(Grid grid, Element currEl) {

#ifdef TEST
		/// Dane do testów
		double x[4] = { 0, 0.025, 0.025, 0 };
		double y[4] = { 0, 0, 0.025, 0.025 };
		J[0][0] += el4.dXYdKsiSum(x);
		J[0][1] += el4.dXYdEtaSum(x);
		J[1][0] += el4.dXYdKsiSum(y);
		J[1][1] += el4.dXYdEtaSum(y);
#else
		/// W³aœciwe dane
		// Stwórz macierz Jakobiego ze wszystkich punktów ca³kowania
		for (int j = 0; j < this->p; j++)
		{
			int id = currEl.id[j] - 1;
			double x = grid.nodes[id].x;
			double y = grid.nodes[id].y;

			this->calcJ(J, j, x, y);
		}
#endif

		this->detJ = detOfMatrix(J);
		inverseMatrix(invJ, J, detJ);

		// Wype³nij tablicê pochodnych funkcji N po X i Y - Jakobian?
		for (int j = 0; j < this->p; j++) {
			for (int k = 0; k < this->p; k++) {
				dNdX[j][k] = invJ[0][0] * this->ksiMatrix[j][k] + invJ[0][1] * this->etaMatrix[j][k];
				dNdY[j][k] = invJ[1][0] * this->ksiMatrix[j][k] + invJ[1][1] * this->etaMatrix[j][k];
			}
		}
	}

	/**
	 * Zwraca wyznacznik macierzy Jakobiego.
	 *
	 * @returns Wyznacznik macierzy J
	 */
	double getDetJ() {
		return this->detJ;
	}

	/**
	 * Zwraca tablicê funkcji kszta³tu.
	 * 
	 * @param ksi - wartoœci ksi dla funkcji
	 * @param eta - wartoœci eta dla funkcji
	 */
	double* nKsiEta(double ksi, double eta) {
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
	double* dNdEta(double ksi) {
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
	double* dNdKsi(double eta) {
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
	double dXYdEta(int pc, double xy) {
		return dEta[pc] * xy;
	}

	/**
	 * Zwraca Jakobian przekszta³cenia po ksi.
	 *
	 * @param pc - punkt ca³kowania
	 * @param xy - wymno¿enie po "x" lub po "y"
	 */
	double dXYdKsi(int i, double xy) {
		return dKsi[i] * xy;
	}

	/**
	 * Zwraca sumê Jakobianu przekszta³cenia po eta.
	 *
	 * @param xy - wymno¿enie po "x" lub po "y"
	 */
	double dXYdEtaSum(double *xy) {
		double sum = 0;
		for (size_t i = 0; i < 4; i++) sum += dEta[i] * xy[i];	
		return sum;
	}

	/**
	 * Zwraca sumê Jakobianu przekszta³cenia po ksi.
	 *
	 * @param xy - wymno¿enie po "x" lub po "y"
	 */
	double dXYdKsiSum(double* xy) {
		double sum = 0;
		for (size_t i = 0; i < 4; i++) sum += dKsi[i] * xy[i];
		return sum;
	}

	/**
	 * Wype³nia macierz pochodnych funkcji N po eta i ksi.
	 *
	 * @param xy - wymno¿enie po "x" lub po "y"
	 */
	void fillMatrices() {
		for (int x = 0; x < p; x++)
		{
			for (int y = 0; y < 4; y++)
			{
				etaMatrix[x][y] = dNdEta(ksi[x])[y];
				ksiMatrix[x][y] = dNdKsi(eta[x])[y];
			}
		}
	}	
};
 
double* globalP;
double** globalH, **globalC;

/**
 * Funkcja inicjuj¹ca globalne macierze do przechowywania zagregowanych tablic/wektorów.
 *
 * @param N - wielkoœæ macierzy
 */
void initGlobalMatrices(const int N) {
	globalP = new double  [N];
	globalH = new double* [N];
	globalC = new double* [N];

	for (int i = 0; i < N; i++) {
		globalH[i] = new double[N];
		globalC[i] = new double[N];
		globalP[i] = 0;

		for (int j = 0; j < N; j++) {
			globalH[i][j] = 0;
			globalC[i][j] = 0;
		}
	}
}

/**
 * Funkcja zwalniaj¹ca z pamiêci zinicjowane funkcj¹ initGlobalMatrices() macierze.
 *
 * @param N - wielkoœæ macierzy
 */
void destroyGlobalMatrices(const int N) {
	for (int i = 0; i < N; i++) {
		delete globalH[i];
		delete globalC[i];
	}
	delete[] globalH;
}

/**
 * Analizuje plik pod wskazan¹ œcie¿k¹ i zapisuje dane do przekazanych referencji.
 * Plik musi byæ w formacie takim jak na zajêciach, tj.
 * zmienne i ich wartoœci oddzielone spacj¹, informacje o elementach, wêz³ach i warunkach brzegu
 * poprzedzone asteriksem (*) z wartoœciami w nowej linii oddzielonymi spacj¹.
 * 
 * Przyk³ad:
 *  SimulationTime 500
 *  SimulationStepTime 50
 *  *Node
 *     1,  0.100000001, 0.00499999989      
 *  *Element, type=DC2D4
 *	 1,  1,  2,  6,  5
 *	*BC
 *	1, 2, 3, 4, 5, 8, 9, 12, 13, 14, 15, 16
 * 
 * @param path - œcie¿ka do pliku tekstowego
 * @param grid - siatka do zapisu danych
 * @param simulationTime, simulationStepTime, ... - referencje do zapisu wartoœci
 */
void parseTextFile(string path, Grid& grid, double& simulationTime, double& simulationStepTime, double& conductivity, double& alpha, double& density, double& specificHeat, double& initialTemp, double& ambientTemp) {
	int nodesNo = 0, elementsNo = 0;
	Element *elements = NULL;
	Node* nodes = NULL;

	string line, temp;
	ifstream file;

	file.open(path);

	if (!file.is_open()) return;
	while (getline(file, line)) {
		if (regex_match(line, std::regex("SimulationTime \\d+"))) {
			temp = regex_replace(line, std::regex("SimulationTime "), "$1");			
			simulationTime = stod(temp);
		}
		else if (regex_match(line, std::regex("SimulationStepTime \\d+"))) {
			temp = regex_replace(line, std::regex("SimulationStepTime "), "$1");
			simulationStepTime = stod(temp);
		}
		else if (regex_match(line, std::regex("Conductivity \\d+"))) {
			temp = regex_replace(line, std::regex("Conductivity "), "$1");
			conductivity = stod(temp);
		}
		else if (regex_match(line, std::regex("Alfa \\d+"))) {
			temp = regex_replace(line, std::regex("Alfa "), "$1");
			alpha = stod(temp);
		}
		else if (regex_match(line, std::regex("Tot \\d+"))) {
			temp = regex_replace(line, std::regex("Tot "), "$1");
			ambientTemp = stod(temp);
		}
		else if (regex_match(line, std::regex("InitialTemp \\d+"))) {
			temp = regex_replace(line, std::regex("InitialTemp "), "$1");
			initialTemp = stod(temp);
		}
		else if (regex_match(line, std::regex("Density \\d+"))) {
			temp = regex_replace(line, std::regex("Density "), "$1");
			density = stod(temp);
		}
		else if (regex_match(line, std::regex("SpecificHeat \\d+"))) {
			temp = regex_replace(line, std::regex("SpecificHeat "), "$1");
			specificHeat = stod(temp);
		}
		else if (regex_match(line, std::regex("Nodes number \\d+"))) {
			temp = regex_replace(line, std::regex("Nodes number "), "$1");
			nodesNo = stoi(temp);
		}
		else if (regex_match(line, std::regex("Elements number \\d+"))) {
			temp = regex_replace(line, std::regex("Elements number "), "$1");
			elementsNo = stoi(temp);
		}

		// Odczyt wêz³ów
		if (regex_match(line, std::regex("\\*Node"))) {			
			nodes = new Node[nodesNo];

			while (getline(file, line)) {
				if (regex_match(line, std::regex("\\s+\\d+,\\s+[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+),\\s+[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)"))) {
					istringstream ss(line);
					string token;
					string results[3];
					int i = 0;

					while (getline(ss, token, ',')) {
						results[i] = token;
						i = i > 3 ? 0 : i + 1;
					}

					int id = stoi(results[0]) - 1;
					double x = stod(results[1]);
					double y = stod(results[2]);

					nodes[id].id = id;
					nodes[id].x = x;
					nodes[id].y = y;
				}					
				else
					break;
			}
		}

		// Odczyt elementów
		if (regex_match(line, regex("\\*Element, type=DC2D4"))) { // uwzglêdniæ type !!!
			elements = new Element[elementsNo];

			while (getline(file, line)) {
				if (regex_match(line, regex("\\s*\\d+,\\s+\\d+,\\s+\\d+,\\s+\\d+,\\s+\\d+"))) {
					istringstream ss(line);
					string token;
					string results[5];
					int i = 0;

					while (getline(ss, token, ',')) {
						results[i] = token;
						i = i > 5 ? 0 : i + 1;						
					}

					int index = stoi(results[0]) - 1;
					int id1 = stoi(results[1]);
					int id2 = stoi(results[2]);
					int id3 = stoi(results[3]);
					int id4 = stoi(results[4]);

					elements[index].id[0] = id1;
					elements[index].id[1] = id2;
					elements[index].id[2] = id3;
					elements[index].id[3] = id4;
				}					
				else
					break;
			}
		}

		// Odczyt warunków brzegowych
		if (regex_match(line, std::regex("\\*BC"))) {
			while (getline(file, line)) {
				if (nodes == NULL)
					return;

				string token;
				istringstream ss(line);
				while (getline(ss, token, ',')) {
					int nodeId = stoi(token) - 1;
					nodes[nodeId].bc = 1;
				}
			}
		}		
	}

	file.close();

	grid.nN = nodesNo;
	grid.nE = elementsNo;
	grid.nodes = nodes;
	grid.elements = elements;

	// delete nodes, elements;
}

int main() {
	double a = 1/sqrt(3);
	double ksiSchema[4] = { -a, a, a, -a };
	double etaSchema[4] = { -a, -a, a, a };
	Element4_2D el4(ksiSchema, etaSchema, 2);
	Gauss gauss;

	#if LAB_NO >= 10
		Grid grid;
		double simTime, simStepTime, conductivity, alpha, density, specificHeat, initialTemp, ambientTemp;
		parseTextFile("samples/MES_31_31_v2.txt", grid, simTime, simStepTime, conductivity, alpha, density, specificHeat, initialTemp, ambientTemp);

		cout << "SimulationTime " << simTime << "\n";
		cout << "SimulationStepTime  " << simStepTime << "\n";
		cout << "Conductivity " << conductivity << "\n";
		cout << "Alfa " << alpha << "\n";
		cout << "Tot " << ambientTemp << "\n";
		cout << "InitialTemp " << initialTemp << "\n";
		cout << "Density " << density << "\n";
		cout << "SpecificHeat " << specificHeat << "\n";
		cout << "Nodes number " << grid.nN << "\n";
		cout << "Elements number " << grid.nE << "\n";
		grid.printNodesCoords();
		grid.printElementsNodes();
		grid.printBCNodes();

	#elif LAB_NO >= 6 && LAB_NO < 10
		Grid grid(0.1, 0.1, 4, 4);
		const int N = grid.nN;

		initGlobalMatrices(N);

		for (int i = 0; i < grid.nE; i++) {
			Element currEl = grid.elements[i];

			el4.jacobian(grid, currEl);
			el4.calcH(currEl.H, 25);
			el4.calcC(currEl.C, 700, 7800);
			el4.calcHbc(currEl.Hbc, currEl.P, grid, currEl, i, 300, 1200);
			el4.aggregate(globalH, globalC, globalP, currEl);
		}		
		
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// Oczekiwane => 12000, 12000, 12000, 12000, 12000, 0, 0, 12000, 12000, 0, 0, 12000, 12000, 12000, 12000, 12000
		// Otrzymane  => 6000, 12000, 12000, 6000, 6000, 0, 0, 6000, 6000, 0, 0, 12000, 6000, 12000, 6000, 6000
		double* testGlobalP = { new double[16] { 12000, 12000, 12000, 12000, 12000, 0, 0, 12000, 12000, 0, 0, 12000, 12000, 12000, 12000, 12000 } };
		// printMatrix(globalP, N);
				
		double** Ccaret;
		double** Hcaret;
		double* Pcaret, *vectorC, *t1;
		double* t0 = new double[N];
		for (int i = 0; i < N; i++) t0[i] = 100;

		cout << "Time[s]   MinTemp[s]   MaxTemp[s]\n";

		for (size_t dT = 50, i = 0; dT <= 500; dT += 50, i++) {
			Ccaret = dT > 0 ? divide(globalC, dT, N, N) : divide(globalC, 1, N, N);
			Hcaret = add(globalH, Ccaret, N, N);
			Hcaret = add(globalH, Ccaret, N, N);
			vectorC = multiply(Ccaret, t0, N);
			Pcaret = add(globalP, vectorC, N, N);
			t1 = gauss.elimination(Hcaret, Pcaret, N);
			t0 = t1;

			/*cout << std::string(89, '_') << " Iteration " << i << " " << std::string(89, '_') << "\n";
			cout << std::string(85, '_') << " Matrix ([H]+[C]/dT) " << std::string(85, '_') << "\n";
			printMatrix(Hcaret, N, N);
			cout << std::string(80, '_') << " Vector ({P}+{[C]/dT}*{T0}) " << std::string(81, '_') << "\n";
			printMatrix(Pcaret, N);
			cout << "\n";*/

			double minTemp = 9999, maxTemp = 0;
			//getMinAndMaxElement(t1, minTemp, maxTemp, N);
			for (size_t k = 0; k < N; k++)	{
				if (t1[k] < minTemp) minTemp = t1[k];
				if (t1[k] > maxTemp) maxTemp = t1[k];
			}
			cout << setw(7) << dT << "   "  << setw(10) << minTemp << "   " << setw(10) << maxTemp << "\n";
		}

		delete[] t0;
		destroyGlobalMatrices(N);		

	#elif LAB_NO == 4 || LAB_NO == 5
		Grid grid(0.025f, 0.025f, 4, 4);

		for (size_t i = 0; i < grid.nE; i++) {
			el4.jacobian(grid, grid.elements[i]);

			#if LAB_NO == 4
				el4.calcH(grid.elements[i].H);
			#elif LAB_NO == 5
				el4.calcHbc(grid.elements[i].Hbc, grid, 25);
			#endif
		}

	#elif LAB_NO == 3
		std::cout << std::setw(16) << " " << ">>> Funkcje ksztaltu dN/dKsi <<<\n";
		printMatrix4x4(el4.ksiMatrix);
		std::cout << std::setw(16) << " " << ">>> Funkcje ksztaltu dN/dEta <<<\n";
		printMatrix4x4(el4.etaMatrix);
	#elif LAB_NO == 2
		Gauss gauss;
		int n = 3;
		double interval1d = gauss.quadrature1d(-1, 1, (vFunctionCall)f, n);		// 15.3333 OK | https://www.wolframalpha.com/input/?i2d=true&i=Integrate%5B5Power%5Bx%2C2%5D%2B3x%2B6%2C%7Bx%2C-1%2C1%7D%5D
		double interval2d = gauss.quadrature2d((vFunctionCall2)f, n);		    // 26.2222 OK | https://www.wolframalpha.com/input/?i2d=true&i=Integrate%5BIntegrate%5B5Power%5Bx%2C2%5DPower%5By%2C2%5D+%2B+3xy+%2B+6%2C%7Bx%2C-1%2C1%7D%5D%2C%7By%2C-1%2C1%7D%5D
		std::cout << "Kwadratura 1D dla " << n << "p: " << interval1d << std::endl;
		std::cout << "Kwadratura 2D dla " << n << "p: " << interval2d << std::endl;
	#elif LAB_NO == 1
		Grid grid(0.2f, 0.1f, 5, 4);
		grid.printNodesCoords();
		grid.printElementsNodes();
	#endif
	system("pause");
	return 0;
}