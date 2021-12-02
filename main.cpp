#include <iostream>
#include <iomanip>
#include "functions.h"
#include "matrix.h"
#include "gauss.h"

using namespace std;

#define LAB_NO 6

//#define SHOW_MATRIX
//#define SHOW_DETAILS
//#define TEST

/**
 * Wêze³ bêd¹cy czêœci¹ elementu.
 * 
 * @attrib x, y - wspó³rzêdne wêz³ów
 * @attrib bc - informacja o obecnoœci warunku brzegowego
 */
struct Node {
	float x = 0, y = 0;
	short int bc;
};

/**
 * Element uniwersalny.
 * 
 * @attrib id[4] - identyfikatory wêz³ów elementu 
 */
struct Element {
	int id[4] = { 0, 0, 0, 0 };
	double H[4][4] = { 0 }, Hbc[4][4] = { 0 };

	int sideDown[2] = { id[0], id[1] };
	int sideRight[2] = { id[1], id[2] };
	int sideUp[2] = { id[2], id[3] };
	int sideLeft[2] = { id[3], id[0] };
};

/**
 * Siatka elementów.
 *
 * @attrib width, height - wymiary siatki
 * @attrib nW, nH - liczba wêz³ów odpowiednio na osi X i osi Y
 * @attrib nN - ³¹czna liczba wêz³ów w siatce
 * @attrib nE - liczba elementów, z których sk³ada siê siatka
 */
struct Grid {
	float width, height;
	int nW, nH, nN, nE;
	Node* nodes;
	Element* elements;

	Grid(float height, float width, int nH, int nW) {
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
	}

	~Grid() {
		// ERROR
		//delete[] this->nodes;
		//delete[] this->elements;
	}

	void setNodesCoords() {
		int i = 0;

		for (size_t x = 0; x < nW; x++)
		{
			for (size_t y = 0; y < nH; y++)
			{
				++i;
				nodes[i].y = nodes[i - 1].y + height / (nH - 1);
				nodes[i].x = nodes[i - 1].x;
			}

			nodes[i].x = nodes[i - 1].x + width / (nW - 1);
			nodes[i].y = 0;
		}
	}

	void printNodesCoords() {
		for (size_t i = 0; i < nN; i++)
		{
			std::cout << "Node #" << i + 1 << ":\t[" << this->nodes[i].x << ", " << this->nodes[i].y << "]\n";
		}
	}

	void setElementsNodes() {
		int n = 1, nhSum = nH;

		for (size_t i = 0; i < nE; i++)
		{
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

	void printElementsNodes() {
		for (size_t i = 0; i < nE; i++)
		{
			std::cout << "El #" << i + 1 << "\t{" << elements[i].id[0] << ", " << elements[i].id[1] << ", " << elements[i].id[2] << ", " << elements[i].id[3] << "}\n";

			/*std::cout << "(" << elements[i].id[3] << ")---(" << elements[i].id[2] << ")\n";
			std::cout << " |\t|\n";
			std::cout << "(" <<  elements[i].id[0] << ")---(" << elements[i].id[1] << ")\n\n";*/
		}
	}
};

/**
 * Dwuwymiarowy, czteropunktowy element.
 * 
 * @attrib N[4] - tablica funkcji kszta³tu
 * @attrib dKsi[4], dEta[4] - tablica pochodnych funkcji N kolejno po Ksi i Eta
 * @attrib etaMatrix, ksiMatrix - tablice pochodnych funkcji dla wszystkich punktów ca³kowania
 * @attrib ksi, eta - schematy ca³kowania
 * @attrib p - liczba punktów ca³kowania
 */
struct Element4_2D {
	double N[4], dKsi[4], dEta[4];
	double **ksiMatrix, **etaMatrix;
	double *ksi, *eta;	
	double dNdX[4][4] = { 0 }, dNdY[4][4] = { 0 };
	double dNdX_T[1][4] = { 0 }, dNdY_T[1][4] = { 0 };
	double J[2][2] = { 0 }, invJ[2][2] = { 0 };
	double detJ = 0;
	int p;

	Element4_2D(double* ksi, double* eta, int n = 2) {
		this->p = n * n;
		this->ksi = ksi;
		this->eta = eta;
		ksiMatrix = new double* [p];
		etaMatrix = new double* [p];

		for (size_t i = 0; i < p; i++) {
			ksiMatrix[i] = new double[4];
			etaMatrix[i] = new double[4];
		}

		fillMatrices();
	}

	~Element4_2D() {
		for (size_t i = 0; i < p; i++)
		{
			// Exception thrown
			//delete etaMatrix[i];
			//delete ksiMatrix[i];
		}
		//delete[] etaMatrix, ksiMatrix;
	}

	/// <summary>
	/// Zwraca odleg³oœæ miêdzy dwoma punktami.
	/// </summary>
	/// <param name="x1">Wspó³rzêdna X pierwszego punktu</param>
	/// <param name="y1">Wspó³rzêdna Y pierwszego punktu</param>
	/// <param name="x2">Wspó³rzêdna X drugiego punktu</param>
	/// <param name="y2">Wspó³rzêdna Y drugiego punktu</param>
	double distance(double x1, double y1, double x2, double y2)	{
		return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) * 1.0);
	}

	/// <summary>
	/// Agregacja macierzy H.
	/// </summary>
	/// <param name="globalMatrix">Macierz do zapisu danych</param>
	/// <param name="currElIndex">Element siatki, na którym metoda ma operowaæ</param>
	void aggregate(double**& globalMatrix, Element currEl) {
		double nodes[4] = { currEl.id[0], currEl.id[1], currEl.id[2], currEl.id[3] };
		double tempMatrix[4][4][2] = {
			{ { nodes[0], nodes[0] }, { nodes[0], nodes[1] }, { nodes[0], nodes[2] }, { nodes[0], nodes[3] } },
			{ { nodes[1], nodes[0] }, { nodes[1], nodes[1] }, { nodes[1], nodes[2] }, { nodes[1], nodes[3] } },
			{ { nodes[2], nodes[0] }, { nodes[2], nodes[1] }, { nodes[2], nodes[2] }, { nodes[2], nodes[3] } },
			{ { nodes[3], nodes[0] }, { nodes[3], nodes[1] }, { nodes[3], nodes[2] }, { nodes[3], nodes[3] } },
		};
		for (size_t j = 0; j < 4; j++)
		{
			for (size_t k = 0; k < 4; k++)
			{
				int rowIndex = tempMatrix[j][k][0] - 1;
				int colIndex = tempMatrix[j][k][1] - 1;
				globalMatrix[rowIndex][colIndex] += currEl.H[j][k];
			}
		}
	}

	/// <summary>
	/// Ustawia w przekazanej tablicy sumê macierzy Hbc dla ka¿dego punktu ca³kowania.
	/// </summary>
	/// <param name="Hbc">Tablica, do której zapisaæ wynik</param>
	/// <param name="grid">Siatka, na któej metoda ma operowaæ</param>
	/// <param name="currEl">Index obecnego elementu siatki</param>
	/// <param name="alpha">wspó³czynnik ...</param>
	void calcHbc(double Hbc[4][4], Grid grid, int currEl, double alpha = 25) {
		struct Side { double *pc1, *pc2; };
		double Npc1[4][4], Npc2[4][4];
		double w[2] = { 1, 1 };

		double ksi = 0.577350; // ???
		double eta = 0.577350; 
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

		// Inicjalizacja macierzy funkcji N dla nowych wartoœci ksi i eta
		for (size_t i = 0; i < this->p; i++)
		{
			/*int id = grid.elements[currEl].id[i] - 1;
			double x = grid.nodes[id].x;
			double y = grid.nodes[id].y;*/			

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
		for (size_t i = 0; i < this->p; i++) { // i - punkt ca³kowania
			double* N1 = this->nKsiEta(sides[i].pc1[0], sides[i].pc1[1]);
			double* N2 = this->nKsiEta(sides[i].pc2[0], sides[i].pc2[1]);

			for (size_t j = 0; j < 4; j++)
				for (size_t k = 0; k < 4; k++) {
					Hbc[j][k] = alpha * w[0] * (Npc1[i][j] * Npc1[i][k]);
					Hbc[j][k] += alpha * w[1] * (Npc2[i][j] * Npc2[i][k]);
					Hbc[j][k] *= 0.0125; // localDetJ
				}

			#ifdef SHOW_MATRIX
				std::cout << "Pow_" << i + 1 << "\n";
				printMatrix(Hbc);
			#endif	
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
		for (size_t i = 0; i < 4; i++) {
			double rowsX[4] = { dNdX[i][0],	dNdX[i][1], dNdX[i][2],	dNdX[i][3] };
			double rowsY[4] = { dNdY[i][0], dNdY[i][1], dNdY[i][2], dNdY[i][3] };
			//transpose(dNdX, *dNdX_T, 1, 4); transpose(dNdY, *dNdY_T, 1, 4);

#ifdef SHOW_DETAILS
			std::cout << k << " * (" << rowsX[0] << " * " << dNdX_T[0][0] << " + " << rowsY[0] << " * " << dNdY_T[0][0] << ") * " << detJ << " = ";
			std::cout << k << " * (" << rowsX[0] << " * " << rowsX[0] << " + " << rowsY[0] << " * " << rowsY[0] << ") * " << detJ << " = ";
			std::cout << k * (rowsX[0] * rowsX[0] + rowsY[0] * rowsY[0]) * detJ << "\n";
#endif

			for (size_t j = 0; j < 4; j++)
				for (size_t l = 0; l < 4; l++)
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
	void jacobian(int i, Grid grid) {

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
		for (size_t j = 0; j < this->p; j++)
		{
			int id = grid.elements[i].id[j] - 1;
			double x = grid.nodes[id].x;
			double y = grid.nodes[id].y;

			this->calcJ(J, j, x, y);
		}
#endif

		this->detJ = detOfMatrix(J);
		inverseMatrix(invJ, J, detJ);

		// Wype³nij tablicê pochodnych funkcji N po X i Y - Jakobian?
		for (size_t j = 0; j < this->p; j++) {
			for (size_t k = 0; k < this->p; k++) {
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

	// Funkcje kszta³tu
	double* nKsiEta(double ksi, double eta) {
		N[0] = 0.25 * (1 - ksi) * (1 - eta);
		N[1] = 0.25 * (1 + ksi) * (1 - eta);
		N[2] = 0.25 * (1 + ksi) * (1 + eta);
		N[3] = 0.25 * (1 - ksi) * (1 + eta);
		
		return N;
	}

	// Pochodna funkcji N wzglêdem eta
	double* dNdEta(double ksi) {
		dEta[0] = -0.25 * (1 - ksi);
		dEta[1] = -0.25 * (1 + ksi);
		dEta[2] = 0.25 * (1 + ksi);
		dEta[3] = 0.25 * (1 - ksi);

		return dEta;
	}

	// Pochodna funkcji N wzglêdem ksi
	double* dNdKsi(double eta) {
		dKsi[0] = -0.25 * (1 - eta);
		dKsi[1] = 0.25 * (1 - eta);
		dKsi[2] = 0.25 * (1 + eta);
		dKsi[3] = -0.25 * (1 + eta);

		return dKsi;
	}

	// Jakobian przekszta³cenia po eta
	double dXYdEta(int i, double xy) {
		return dEta[i] * xy;
	}

	// Jakobian przekszta³cenia po ksi
	double dXYdKsi(int i, double xy) {
		return dKsi[i] * xy;
	}

	// Suma jakobianu przekszta³cenia po eta
	double dXYdEtaSum(double *x) {
		double sum = 0;
		for (size_t i = 0; i < 4; i++) sum += dEta[i] * x[i];	
		return sum;
	}

	// Suma jakobianu przekszta³cenia po ksi
	double dXYdKsiSum(double* x) {
		double sum = 0;
		for (size_t i = 0; i < 4; i++) sum += dKsi[i] * x[i];
		return sum;
	}

	// Wype³nianie macierzy pochodnych funkcji N po eta i ksi
	void fillMatrices() {
		for (size_t x = 0; x < p; x++)
		{
			for (size_t y = 0; y < 4; y++)
			{
				etaMatrix[x][y] = dNdEta(ksi[x])[y];
				ksiMatrix[x][y] = dNdKsi(eta[x])[y];
			}
		}
	}	
};

int main() {
	double a = 0.577350;		// 1/sqrt(3)	
	double ksiSchema[4] = { -a, a, a, -a };
	double etaSchema[4] = { -a, -a, a, a };
	Element4_2D el4(ksiSchema, etaSchema, 2);
	Grid grid(0.100f, 0.100f, 4, 4);

	#if LAB_NO == 6
		// Init
		const int N = grid.nN;
		double **notSoGlobalH = new double *[N];
		for (int i = 0; i < N; ++i)
			notSoGlobalH[i] = new double[N];
		for (size_t i = 0; i < N; i++)
			for (size_t j = 0; j < N; j++)
				notSoGlobalH[i][j] = 0;

		for (size_t i = 0; i < grid.nE; i++) {
			el4.jacobian(i, grid);
			el4.calcH(grid.elements[i].H, 25);
			el4.aggregate(notSoGlobalH, grid.elements[i]);		

			cout << std::string(73, '_') << " Iteration " << i << std::string(73, '_') << "\n";
			printMatrix(notSoGlobalH, N, N);
		}

		// Destroy
		for (size_t i = 0; i < N; i++)
			delete notSoGlobalH[i];
		delete[] notSoGlobalH;

	#elif LAB_NO == 4 || LAB_NO == 5
		for (size_t i = 0; i < grid.nE; i++) {
			double H[4][4] = { 0 }, Hbc[4][4] = { 0 };
			jac.jacobian(i, el4, grid);

			#if LAB_NO == 4
				jac.calcH(H, jac.getDetJ());
			#elif LAB_NO == 5
				el4.calcHbc(grid.elements[i].Hbc, grid, 25);
			#endif

		}

	#elif LAB_NO == 3
		double a = -0.577350;
		double ksiSchema[4] = { a, -a, -a, a };
		double etaSchema[4] = { a, a, -a, -a };
		Element4_2D el4(ksiSchema, etaSchema, 2);
		std::cout << std::setw(16) << " " << ">>> Funkcje ksztaltu dN/dKsi <<<\n";
		el4.printMatrix(el4.ksiMatrix);
		std::cout << std::setw(16) << " " << ">>> Funkcje ksztaltu dN/dEta <<<\n";
		el4.printMatrix(el4.etaMatrix);
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