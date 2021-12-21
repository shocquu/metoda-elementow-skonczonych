#include "grid.h"

using namespace std;

Grid::Grid() {
	this->width = 0;
	this->height = 0;
	this->nW = 0;
	this->nH = 0;
	this->nN = 0;
	this->nE = 0;
	this->nodes = nullptr;
	this->elements = nullptr;
}
Grid::Grid(double height, double width, int nH, int nW) {
	this->height = height;
	this->width = width;
	this->nH = nH;
	this->nW = nW;
	this->nN = nH * nW;
	this->nE = (nW - 1) * (nH - 1);
	this->nodes = new Node[nN];
	this->elements = new Element[nE];

	initMatrices();
	setNodesCoords();
	setElementsNodes();
	setNodesBC();
}

Grid::~Grid() {
	// ERROR
	//delete[] this->nodes;
	//delete[] this->elements;
}

void Grid::initMatrices() {
	const int N = nN;
	globalP = new double[N];
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
 * Odczytuje plik pod wskazan¹ œcie¿k¹ i zapisuje dane do przekazanych referencji.
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
void Grid::readFromFile(string path, double& simulationTime, double& simulationStepTime, double& conductivity, double& alpha, double& density, double& specificHeat, double& initialTemp, double& ambientTemp) {
	int nodesNo = 0, elementsNo = 0;
	Element* elements = NULL;
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

	this->nN = nodesNo;
	this->nE = elementsNo;
	this->nodes = nodes;
	this->elements = elements;

	// @TODO
	// delete nodes, elements;
}

/**
 * Ustawia w przekazanej tablicy sumê macierzy C dla ka¿dego punktu ca³kowania.
 *
 * @param globalH - tablica do zapisu zagregowanej macierzy H
 * @param globalC - tablica do zapisu zagregowanej macierzy C
 * @param globalP - tablica do zapisu zagregowanej macierzy P
 * @param currEl - element, na którym metoda ma operowaæ
 */
void Grid::aggregate(Element currEl) {
	int p = 4; //TEMP


	int nodesId[4] = { currEl.id[0] - 1, currEl.id[1] - 1, currEl.id[2] - 1, currEl.id[3] - 1 };
	int sideNodes[4][4][2] = {									// pomocnicza macierz ID wêz³ów do rozmieszczenia elementów
		{ { nodesId[0], nodesId[0] }, { nodesId[0], nodesId[1] }, { nodesId[0], nodesId[2] }, { nodesId[0], nodesId[3] } },
		{ { nodesId[1], nodesId[0] }, { nodesId[1], nodesId[1] }, { nodesId[1], nodesId[2] }, { nodesId[1], nodesId[3] } },
		{ { nodesId[2], nodesId[0] }, { nodesId[2], nodesId[1] }, { nodesId[2], nodesId[2] }, { nodesId[2], nodesId[3] } },
		{ { nodesId[3], nodesId[0] }, { nodesId[3], nodesId[1] }, { nodesId[3], nodesId[2] }, { nodesId[3], nodesId[3] } },
	};

	for (int i = 0; i < p; i++) {
		for (int j = 0; j < p; j++) {
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
void Grid::calcHbc(Element4_2D &el4, Element &currEl, int elIndex, double alpha = 25, double ambientTemp = 1200) {
	int p = 4; // !!!!
	
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
	for (int i = 0; i < p; i++) {
		double* Nrow = el4.nKsiEta(sides[i].pc1[0], sides[i].pc1[1]);
		Npc1[i][0] = Nrow[0];
		Npc1[i][1] = Nrow[1];
		Npc1[i][2] = Nrow[2];
		Npc1[i][3] = Nrow[3];
		Nrow = el4.nKsiEta(sides[i].pc2[0], sides[i].pc2[1]);
		Npc2[i][0] = Nrow[0];
		Npc2[i][1] = Nrow[1];
		Npc2[i][2] = Nrow[2];
		Npc2[i][3] = Nrow[3];
	}

	// Wype³nianie macierzy
	for (int i = 0; i < p; i++) { // i - œciana elementu
		int n1Index = currEl.id[i] - 1;
		int n2Index = currEl.id[(i + 1) % 4] - 1;

		// Sprawdzanie warunku brzegowego na œcianach elementu
		if (nodes[n1Index].bc > 0 && nodes[n2Index].bc > 0) {
			double L = distance(nodes[n1Index].x, nodes[n1Index].y, nodes[n2Index].x, nodes[n2Index].y);
			double locDetJ = L / 2;

			for (int j = 0; j < p; j++) {
				for (int k = 0; k < p; k++) {
					currEl.Hbc[j][k] += w[0] * (Npc1[i][j] * Npc1[i][k]) * alpha * locDetJ;
					currEl.Hbc[j][k] += w[1] * (Npc2[i][j] * Npc2[i][k]) * alpha * locDetJ;
				}

				currEl.P[j] += (w[0] * Npc1[i][j] + w[1] * Npc2[i][j]) * alpha * ambientTemp * locDetJ;
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
void Grid::calcH(Element4_2D &el4, Element &currEl, double k = 30) {
	int p = 4; // !!!!

	for (int i = 0; i < p; i++) {
		double rowsX[4] = { el4.dNdX[i][0],	el4.dNdX[i][1], el4.dNdX[i][2],	el4.dNdX[i][3] };
		double rowsY[4] = { el4.dNdY[i][0], el4.dNdY[i][1], el4.dNdY[i][2], el4.dNdY[i][3] };

		for (int j = 0; j < p; j++)
			for (int l = 0; l < p; l++)
				currEl.H[j][l] += k * (rowsX[j] * rowsX[l] + rowsY[j] * rowsY[l]) * currEl.detJ;
	}
}

/**
 * Ustawia w przekazanej tablicy sumê macierzy C (pojemnoœci cieplnej) dla ka¿dego punktu ca³kowania.
 *
 * @param c - ciep³o w³aœciwe
 * @param ro - gêstoœæ materia³u
 */
void Grid::calcC(Element4_2D &el4, Element &currEl, double c = 700, double ro = 7800) {
	int p = 4; // !!!!
	
	for (int i = 0; i < p; i++) {
		double* N = el4.nMatrix[i];

		for (int j = 0; j < p; j++)
			for (int k = 0; k < p; k++)
				currEl.C[j][k] += c * ro * (N[j] * N[k]) * currEl.detJ;
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
void Grid::calcJ(Element4_2D &el4, double J[2][2], int pc, double x, double y) {
	J[0][0] += el4.dXYdKsi(pc, x);
	J[0][1] += el4.dXYdEta(pc, x);
	J[1][0] += el4.dXYdKsi(pc, y);
	J[1][1] += el4.dXYdEta(pc, y);
}

/**
 * Uzupe³nia Jakobian dla ka¿dego punktu ca³kowania i-tego elementu.
 *
 * @param i - punkt ca³kowania
 * @param el4 - struktura z danymi o elemencie
 * @param grid - siatka, na której wykonywana jest procedura liczenia Jakobianu
 */
void Grid::jacobian(Element4_2D &el4, Element &currEl) {
	int p = 4; // !!!!

	// Obliczanie macierzy Jacobiego
	for (int i = 0; i < p; i++) {
		int nodeId = currEl.id[i] - 1;
		double x = nodes[nodeId].x;
		double y = nodes[nodeId].y;

		//this->calcJ(el4, currEl.J, i, x, y);
		currEl.J[0][0] += el4.dXYdKsi(i, x);
		currEl.J[0][1] += el4.dXYdEta(i, x);
		currEl.J[1][0] += el4.dXYdKsi(i, y);
		currEl.J[1][1] += el4.dXYdEta(i, y);
	}

	currEl.detJ = detOfMatrix(currEl.J);
	inverseMatrix(currEl.invJ, currEl.J, currEl.detJ);

	// Wype³nij tablicê pochodnych funkcji N po X i Y - jakobian
	for (int j = 0; j < p; j++) {
		for (int k = 0; k < p; k++) {
			el4.dNdX[j][k] = currEl.invJ[0][0] * el4.ksiMatrix[j][k] + currEl.invJ[0][1] * el4.etaMatrix[j][k];
			el4.dNdY[j][k] = currEl.invJ[1][0] * el4.ksiMatrix[j][k] + currEl.invJ[1][1] * el4.etaMatrix[j][k];
		}
	}
}

/**
 * Zwraca odleg³oœæ miêdzy dwoma punktami.
 *
 * @param x1, y1 - wspó³rzêdna x i y pierwszego punktu
 * @param x2, y2 - wspó³rzêdna x i y drugiego punktu
 */
double Grid::distance(double x1, double y1, double x2, double y2) {
	return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

/**
 * Ustawia wspó³rzêdne wierzcho³ków siatki. 
 */
void Grid::setNodesCoords() {
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

/**
 * Wyœwietla wspó³rzêdne wierzcho³ków siatki.
 */
void Grid::printNodesCoords() {
	for (int i = 0; i < nN; i++)
		std::cout << "Node #" << i + 1 << ":\t[" << this->nodes[i].x << ", " << this->nodes[i].y << "]\n";
}

/**
 * Ustawia wêz³y nale¿¹ce do ka¿dego z elementów siatki.
 */
void Grid::setElementsNodes() {
	int n = 1, nhSum = nH;

	for (int i = 0; i < nE; i++) {
		elements[i].index = i;
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

/**
 * Wyœwietla wêz³y nale¿¹ce do ka¿dego z elementów siatki.
 */
void Grid::printElementsNodes() {
	for (int i = 0; i < nE; i++) {
		//if (type == "inline") {
			std::cout << "El #" << i + 1 << "\t{" << elements[i].id[0] << ", " << elements[i].id[1] << ", " << elements[i].id[2] << ", " << elements[i].id[3] << "}\n";
		//}
		//else {
		//	std::cout << "(" << elements[i].id[3] << ")---(" << elements[i].id[2] << ")\n";
		//	std::cout << " |\t|\n";
		//	std::cout << "(" << elements[i].id[0] << ")---(" << elements[i].id[1] << ")\n\n";
		//}
	}
}

/**
 * Ustawia warunki brzegowe dla wêz³ów granicznych siatki.
 */
void Grid::setNodesBC() {
	for (int i = 0; i < nN; i++) {
		if (nodes[i].x == 0 || nodes[i].x == width)  nodes[i].bc = 1;
		if (nodes[i].y == 0 || nodes[i].y == height) nodes[i].bc = 1;
	}
}

/**
 * Wyœwietla wêz³y z warunka brzegowymi.
 */
void Grid::printBCNodes() {
	std::cout << "\nBC nodes:\n";
	for (int i = 0; i < nN; i++)
		if (nodes[i].bc == 1)
			std::cout << i + 1 << "  ";

	std::cout << "\n";
}
