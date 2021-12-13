#include <iostream>
#include <iomanip>
#include "functions.h"
#include "matrix.h"
#include "gauss.h"

using namespace std;

#define LAB_NO 9

//#define SHOW_MATRIX
//#define SHOW_DETAILS
//#define TEST

/**
 * W�ze� b�d�cy cz�ci� elementu.
 * 
 * @attrib x, y - wsp�rz�dne w�z��w
 * @attrib bc - informacja o obecno�ci warunku brzegowego
 */
struct Node {
	float x = 0, y = 0;
	short int bc = -1;
};

/**
 * Element uniwersalny.
 * 
 * @attrib id[4] - identyfikatory w�z��w elementu
 * @attrib H[4][4] - macierz ..
 * @attrib Hbc[4][4] - macierz ..
 * @attrib C[4][4] - macierz ciep�a w�a�ciwego
 * @attrib P[16] - wektor obci��e� / warunk�w brzegowych
 */
struct Element {
	int id[4] = { 0, 0, 0, 0 };
	double H[4][4] = { 0 }, C[4][4] = { 0 }, Hbc[4][4] = { 0 };
	double P[16] = { 0 };	
};

/**
 * Siatka element�w.
 *
 * @attrib width, height - wymiary siatki
 * @attrib nW, nH - liczba w�z��w odpowiednio na osi X i osi Y
 * @attrib nN - ��czna liczba w�z��w w siatce
 * @attrib nE - liczba element�w, z kt�rych sk�ada si� siatka
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

		for (int x = 0; x < nW; x++) {
			for (int y = 0; y < nH; y++) {
				++i;
				nodes[i].x = nodes[i - 1].x;
				nodes[i].y = nodes[i - 1].y + height / (nH - 1);				
			}

			nodes[i].x = nodes[i - 1].x + width / (nW - 1);
			nodes[i].y = 0;
		}
	}

	void printNodesCoords() {
		for (int i = 0; i < nN; i++)
			std::cout << "Node #" << i + 1 << ":\t[" << this->nodes[i].x << ", " << this->nodes[i].y << "]\n";
	}

	void setElementsNodes() {
		int n = 1, nhSum = nH;

		// G�upie, ale dzia�a
		for (size_t i = 0; i < nH; i++)	{
			nodes[i].bc = 1;
			nodes[nN - i - 1].bc = 1;
		}

		for (int i = 0; i < nE; i++) {
			elements[i].id[0] = n;
			elements[i].id[1] = elements[i].id[0] + nH;
			elements[i].id[2] = elements[i].id[1] + 1;
			elements[i].id[3] = elements[i].id[0] + 1;			
			n++;

			if (n == nhSum) {
				nhSum += nH;
				nodes[n].bc = 1;
				nodes[n + nH - 1].bc = 1;
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

	void printBCNodes() {
		std::cout << "\nBC nodes:\n";
		for (size_t i = 0; i < nN; i++)
			if(nodes[i].bc == 1)
				std::cout << i + 1 << "  ";
		
		std::cout << "\n";
	}
};

/**
 * Dwuwymiarowy, czteropunktowy element.
 * 
 * @attrib N[4] - tablica funkcji kszta�tu
 * @attrib dKsi[4], dEta[4] - tablica pochodnych funkcji N kolejno po Ksi i Eta
 * @attrib etaMatrix, ksiMatrix - tablice pochodnych funkcji dla wszystkich punkt�w ca�kowania
 * @attrib ksi, eta - schematy ca�kowania
 * @attrib p - liczba punkt�w ca�kowania
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
	 * Zwraca odleg�o�� mi�dzy dwoma punktami.
	 *
	 * @param x1, y1 - wsp�rz�dna x i y pierwszego punktu
	 * @param x2, y2 - wsp�rz�dna x i y drugiego punktu
	 */
	double distance(double x1, double y1, double x2, double y2)	{
		return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) * 1.0);
	}

	/**
	 * Ustawia w przekazanej tablicy sum� macierzy C dla ka�dego punktu ca�kowania.
	 *
	 * @param globalH - tablica do zapisu zagregowanej macierzy H
	 * @param globalC - tablica do zapisu zagregowanej macierzy C
	 * @param globalP - tablica do zapisu zagregowanej macierzy P
	 * @param currEl - element, na kt�rym metoda ma operowa�
	 */
	void aggregate(double**& globalH, double**& globalC, double*& globalP, Element currEl) {
		int nodes[4] = { currEl.id[0], currEl.id[1], currEl.id[2], currEl.id[3] };
		int tempMatrix[4][4][2] = {							 // pomocnicza macierz do rozmieszczenia element�w
			{ { nodes[0], nodes[0] }, { nodes[0], nodes[1] }, { nodes[0], nodes[2] }, { nodes[0], nodes[3] } },
			{ { nodes[1], nodes[0] }, { nodes[1], nodes[1] }, { nodes[1], nodes[2] }, { nodes[1], nodes[3] } },
			{ { nodes[2], nodes[0] }, { nodes[2], nodes[1] }, { nodes[2], nodes[2] }, { nodes[2], nodes[3] } },
			{ { nodes[3], nodes[0] }, { nodes[3], nodes[1] }, { nodes[3], nodes[2] }, { nodes[3], nodes[3] } },
		};

		for (size_t j = 0; j < this->p; j++) {
			for (size_t k = 0; k < this->p; k++) {
				int rowIndex = tempMatrix[j][k][0] - 1;
				int colIndex = tempMatrix[j][k][1] - 1;
				globalH[rowIndex][colIndex] += currEl.H[j][k];
				globalC[rowIndex][colIndex] += currEl.C[j][k];
			}

			int rowIndex = nodes[j] - 1;
			globalP[rowIndex] += currEl.P[j];
		}
	}

	/**
	 * Ustawia w przekazanej tablicy sum� macierzy Hbc oraz wektor P dla ka�dego punktu ca�kowania.
	 *
	 * @param Hbc - tablica, do kt�rej zapisa� sum� macierzy Hbc z ka�dego punktu ca�kowania
	 * @param P - tablica, do kt�rej zapisa� wektor P
	 * @param alpha - ...
	 * @param ambientTemp - temperatura otoczenia
	 */
	void calcHbc(double Hbc[4][4], double P[4], Grid grid, double alpha = 25, double ambientTemp = 1200) {
		struct Side { double *pc1, *pc2; };
		double Npc1[4][4] = { 0 }, Npc2[4][4] = {0};
		double w[2] = { 1, 1 };
		double ksi, eta;

		ksi = eta = 1/sqrt(3);

		// TEST
		int elSides[4][4][2] = {
			{ grid.elements[0].id[0], grid.elements[0].id[1] },
			{ grid.elements[0].id[1], grid.elements[0].id[2] },
			{ grid.elements[0].id[2], grid.elements[0].id[3] },
			{ grid.elements[0].id[3], grid.elements[0].id[0] }
		};

		int inde = elSides[3][2][1];
		int inde2 = elSides[3][1][1];
		//cout << inde2 << "\n";
		

		double locDetJ = 0.025 / 2; // node.width & node.height		
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

		// Inicjalizacja macierzy funkcji N dla nowych warto�ci ksi i eta
		for (int i = 0; i < this->p; i++) {
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

		// Wype�nianie macierzy
		for (int i = 0; i < this->p; i++) { // i - punkt ca�kowania
			/* nt n1Index = currEl.id[i];
			int n2Index = currEl.id[(i + 1) % 4];
			if( nodes[n1Index].bc == 1 && nodes[n2Index].bc == 1 ) { */

			//cout << "__________________ Sciana " << i + 1 << " __________________\n";

			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					Hbc[j][k] = w[0] * (Npc1[i][j] * Npc1[i][k]);
					Hbc[j][k] += w[1] * (Npc2[i][j] * Npc2[i][k]);
					Hbc[j][k] *= alpha * locDetJ;
				}

				

				P[i] = w[0] * Npc1[i][j];	//P_pc1
				P[i] += w[1] * Npc2[i][j];	//P_pc2
				P[i] *= alpha * locDetJ * ambientTemp;
				
								
				//cout << P[i] << " = (" << w[0] << " * " << Npc1[i][j] << " + " << w[1] << " * " << Npc2[i][j] << ") * " << alpha << " * " << ambientTemp << " * " << locDetJ << "\n";
			}

			/*double Nrow[4] = { Npc1[3][0], Npc1[3][1], Npc1[3][2], Npc1[3][3] };
			double Nrow2[4] = { Npc2[3][0], Npc2[3][1], Npc2[3][2], Npc2[3][3] };
			cout << "pc1\nksi: " << sides[3].pc1[0] << "  eta: "<< sides[3].pc1[1] << "\t";
			printMatrix(Nrow);
			cout << "pc2\nksi: " << sides[3].pc2[0] << "  eta: " << sides[3].pc2[1] << "\t";
			printMatrix(Nrow2);*/

			//cout << "\n";
			
			#ifdef SHOW_MATRIX
				std::cout << "Pow_" << i + 1 << "\n";
				printMatrix(Hbc);
			#endif

			// }

			//cout << P[i] << "\t";
		}

		//cout << P[0] << "\t" << P[1] << "\t" << P[2] << "\t" << P[3] << "\n";
		//printMatrix(P);
		//cout << "\n";
	}
	
	/**
	 * Tworzenie macierzy H dla ka�dego punktu ca�kowania. Wymaga poprzedzenia funkcj�
	 * `fillMatrices`, by operowa� na uzupe�nionych macierzach.
	 *
	 * @param H - macierz, do kt�rej zapisa� wynik
	 * @param detJ - wyznacznik macierzy J
	 * @param k - wsp�czynnik przewodzenia
	 */
	void calcH(double H[4][4], int k = 30) {
		for (size_t i = 0; i < this->p; i++) {
			double rowsX[4] = { dNdX[i][0],	dNdX[i][1], dNdX[i][2],	dNdX[i][3] };
			double rowsY[4] = { dNdY[i][0], dNdY[i][1], dNdY[i][2], dNdY[i][3] };

#ifdef SHOW_DETAILS
			std::cout << k << " * (" << rowsX[0] << " * " << dNdX_T[0][0] << " + " << rowsY[0] << " * " << dNdY_T[0][0] << ") * " << detJ << " = ";
			std::cout << k << " * (" << rowsX[0] << " * " << rowsX[0] << " + " << rowsY[0] << " * " << rowsY[0] << ") * " << detJ << " = ";
			std::cout << k * (rowsX[0] * rowsX[0] + rowsY[0] * rowsY[0]) * detJ << "\n";
#endif

			for (size_t j = 0; j < this->p; j++)
				for (size_t l = 0; l < this->p; l++)
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
	 * Ustawia w przekazanej tablicy sum� macierzy C (pojemno�ci cieplnej) dla ka�dego punktu ca�kowania.
	 *
	 * @param c - ciep�o w�a�ciwe
	 * @param ro - g�sto�� materia�u
	 */
	void calcC(double C[4][4], double c = 700, double ro = 7800) {
		for (size_t i = 0; i < this->p; i++) {
			double* N = nKsiEta(ksi[i], eta[i]);

			for (size_t j = 0; j < this->p; j++)
				for (size_t k = 0; k < this->p; k++)
					C[j][k] += c * ro * (N[j] * N[k]) * detJ;
		}
	}

	/**
	 * Oblicza macierz Jakobiego.
	 *
	 * @param J - macierz �r�d�owa do zapisania danych
	 * @param el4 - struktura z metodami do liczenia pochodnych
	 * @param pc - punkt ca�kowania
	 * @param x, y - wsp�rz�dna x i y
	 */
	void calcJ(double J[2][2], int pc, double x, double y) {
		J[0][0] += this->dXYdKsi(pc, x);
		J[0][1] += this->dXYdEta(pc, x);
		J[1][0] += this->dXYdKsi(pc, y);
		J[1][1] += this->dXYdEta(pc, y);
	}
	
	/**
	 * Uzupe�nia Jakobian dla ka�dego punktu ca�kowania i-tego elementu.
	 *
	 * @param i - punkt ca�kowania
	 * @param el4 - struktura z danymi o elemencie
	 * @param grid - siatka, na kt�rej wykonywana jest procedura liczenia Jakobianu
	 */
	void jacobian(Grid grid, Element currEl) {

#ifdef TEST
		/// Dane do test�w
		double x[4] = { 0, 0.025, 0.025, 0 };
		double y[4] = { 0, 0, 0.025, 0.025 };
		J[0][0] += el4.dXYdKsiSum(x);
		J[0][1] += el4.dXYdEtaSum(x);
		J[1][0] += el4.dXYdKsiSum(y);
		J[1][1] += el4.dXYdEtaSum(y);
#else
		/// W�a�ciwe dane
		// Stw�rz macierz Jakobiego ze wszystkich punkt�w ca�kowania
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

		// Wype�nij tablic� pochodnych funkcji N po X i Y - Jakobian?
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
	 * Zwraca tablic� funkcji kszta�tu.
	 * 
	 * @param ksi - warto�ci ksi dla funkcji
	 * @param eta - warto�ci eta dla funkcji
	 */
	double* nKsiEta(double ksi, double eta) {
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
	double* dNdEta(double ksi) {
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
	double* dNdKsi(double eta) {
		dKsi[0] = -0.25 * (1 - eta);
		dKsi[1] = 0.25 * (1 - eta);
		dKsi[2] = 0.25 * (1 + eta);
		dKsi[3] = -0.25 * (1 + eta);

		return dKsi;
	}

	/**
	 * Zwraca Jakobian przekszta�cenia po eta.
	 *
	 * @param pc - punkt ca�kowania
	 * @param xy - wymno�enie po "x" lub po "y"
	 */
	double dXYdEta(int pc, double xy) {
		return dEta[pc] * xy;
	}

	/**
	 * Zwraca Jakobian przekszta�cenia po ksi.
	 *
	 * @param pc - punkt ca�kowania
	 * @param xy - wymno�enie po "x" lub po "y"
	 */
	double dXYdKsi(int i, double xy) {
		return dKsi[i] * xy;
	}

	/**
	 * Zwraca sum� Jakobianu przekszta�cenia po eta.
	 *
	 * @param xy - wymno�enie po "x" lub po "y"
	 */
	double dXYdEtaSum(double *xy) {
		double sum = 0;
		for (size_t i = 0; i < 4; i++) sum += dEta[i] * xy[i];	
		return sum;
	}

	/**
	 * Zwraca sum� Jakobianu przekszta�cenia po ksi.
	 *
	 * @param xy - wymno�enie po "x" lub po "y"
	 */
	double dXYdKsiSum(double* xy) {
		double sum = 0;
		for (size_t i = 0; i < 4; i++) sum += dKsi[i] * xy[i];
		return sum;
	}

	/**
	 * Wype�nia macierz pochodnych funkcji N po eta i ksi.
	 *
	 * @param xy - wymno�enie po "x" lub po "y"
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

void initGlobalMatrices(const int N) {
	globalP = new double[N];
	globalH = new double* [N];
	globalC = new double* [N];
	for (int i = 0; i < N; i++) {
		globalP[i] = 0;
		globalH[i] = new double[N];
		globalC[i] = new double[N];
		for (int j = 0; j < N; j++) {
			globalH[i][j] = 0;
			globalC[i][j] = 0;
		}
	}
}

void destroyGlobalMatrices(const int N) {
	for (int i = 0; i < N; i++) {
		delete globalH[i];
		delete globalC[i];
	}
	delete[] globalH;
}

int main() {
	double a = 1/sqrt(3);
	double ksiSchema[4] = { -a, a, a, -a };
	double etaSchema[4] = { -a, -a, a, a };
	Element4_2D el4(ksiSchema, etaSchema, 2);
	Gauss gauss;

	#if LAB_NO >= 6		
		Grid grid(0.100f, 0.100f, 4, 4);		
		const int N = grid.nN;

		initGlobalMatrices(N);

		// Chwilowe - DO ZMIANY
		double t0[16] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100 };
		double* testGlobalP = { new double[16]{ 120000, 120000, 120000, 120000, 120000, 0, 0, 120000, 120000, 0, 0, 120000, 120000, 120000, 120000, 120000 } };

		for (int i = 0; i < grid.nE; i++) {
			Element currEl = grid.elements[i];

			el4.jacobian(grid, currEl);
			el4.calcH(currEl.H, 25);
			el4.calcC(currEl.C, 700, 7800);
			el4.calcHbc(currEl.Hbc, currEl.P, grid, 300, 1200);
			el4.aggregate(globalH, globalC, globalP, currEl);
		}

		

		for (size_t dT = 1, i = 1; dT <= 500; dT += 50, i++) {
			cout << std::string(73, '_') << " Iteration " << i << " " << std::string(73, '_') << "\n";

			double** CdT = divide(globalC, dT, N, N);
			double** HplusCdT = add(globalH, CdT, N, N);
			double* PplusCdT = add(globalP, CdT, N, N);

			for (size_t i = 0; i < 16; i++)
				PplusCdT[i] = testGlobalP[i] * t0[i]; // TEST

			//printMatrix(HplusCdT, N, N);
			//printMatrix(PplusCdT, N);
		}

		grid.printBCNodes();
		//grid.printElementsNodes();

		// H*t + P = 0
		//double* result = gauss.elimination(globalH, globalP, N); // t-> ???

		/*for (size_t i = 0; i < N; i++)
			cout << globalP[i] << "  ";					
		cout << "\n";*/

		destroyGlobalMatrices(grid.nN);		

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