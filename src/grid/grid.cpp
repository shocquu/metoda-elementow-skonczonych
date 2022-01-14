#include "grid.h"

namespace plt = matplotlibcpp;
using namespace std;

struct Data {
	double alpha = 300;
	double simTime = 60;
	double simStepTime = 1;
	double conductivity = 25;
	double density = 7800;
	double specificHeat = 700;
	double initialTemp = 100;
	double ambientTemp = 1200;
} gridData;

double Grid::alpha = 300;
double Grid::simTime = 500;
double Grid::simStepTime = 1;
double Grid::conductivity = 25;
double Grid::density = 7800;
double Grid::specificHeat = 700;
double Grid::initialTemp = 100;
double Grid::ambientTemp = 1200;

Grid::Grid() {
	width = height = 0;
	nW = nH = nN = nE = p =0;
	nodes = nullptr;
	elements = nullptr;
	aggrH = aggrC = nullptr;
	aggrP = nullptr;
}
Grid::Grid(Element4_2D el4, string path) {
	width = height = 0;
	nW = nH = nN = nE = p = 0;
	nodes = nullptr;
	elements = nullptr;
	aggrH = aggrC = nullptr;
	aggrP = nullptr;
	readFromFile(path);
	initMatrices();
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
	// @TODO delete pointers

	for (int i = 0; i < nN; i++) {
		//delete aggrC[i];
		//delete aggrH[i];
		//delete[] aggrC[i];
		//delete[] aggrH[i];
	}	
	//delete[] aggrC;
	//delete[] aggrH;	

	// ERROR
	//delete aggrP;
	//delete nodes;
	//delete elements;

	//aggrC = NULL;
	//aggrH = NULL;
	//aggrP = NULL;
	//nodes = NULL;
	//elements = NULL;
}

/**
 * Dokonywanie pomiarów dla danych wczytanych z pliku.
 *
 * @param el4 - referencja do elementu typu EL4_2D
 * @param path - œcie¿ka do pliku tekstowego zawieraj¹cego dane o siatcê
 * @param showDetails - opcja umo¿liwiaj¹ca wyœwietlanie zawartoœci macierzy
 */
void Grid::launch(Element4_2D &el4, std::string path, bool showDetails) {
	readFromFile(path);
	if(nodes) initMatrices();
	launch(el4, simTime, simStepTime, conductivity, alpha, density, specificHeat, initialTemp, ambientTemp, showDetails);
}

/**
 * Dokonywanie pomiarów dla wprowadzonych danych.
 *
 * @param simTime - czas symulacji
 * @param simStepTime - krok czasowy / detla tau
 * @param conductivity - przewodnoœæ
 * @param alpha - wspó³czynnik ...
 * @param density - gêstoœæ materia³u
 * @param specificHeat - ciep³o w³aœciwe
 * @param initialTemp - temperatura pocz¹tkowa
 * @param ambientTemp - temperatura otoczenia
 * @param showDetails - opcja umo¿liwiaj¹ca wyœwietlanie zawartoœci macierzy
 */
void Grid::launch(Element4_2D &el4, double simTime, double simStepTime, double conductivity, double alpha, double density, double specificHeat, double initialTemp, double ambientTemp, bool showDetails) {
	Grid::simTime = simTime;
	Grid::simStepTime = simStepTime;
	Grid::conductivity = conductivity;
	Grid::alpha = alpha;
	Grid::density = density;
	Grid::specificHeat = specificHeat;
	Grid::initialTemp = initialTemp;
	Grid::ambientTemp = ambientTemp;

	printBanner();
	printSimulationData();		
	calcTemperature(el4, true);
}

void Grid::calcTemperature(Element4_2D el4, bool showMinMax, bool saveToFile) {
	for (int i = 0; i < nE; i++) {
		Element currEl = elements[i];
		jacobian(el4, currEl);
		calcH(el4, currEl, conductivity);
		calcC(el4, currEl, specificHeat, density);
		calcHbc(el4, currEl, alpha, ambientTemp);
		aggregate(currEl);
	}

	const int N = nN;
	double** Ccaret = divide(aggrC, simStepTime, N, N);
	double** Hcaret = add(aggrH, Ccaret, N, N);
	double* Pcaret, *vectorC, *t1;
	double* t0 = new double[N];
	ofstream file;
	Gauss g;

	if (showMinMax) cout << " Time[s]   MinTemp[C]   MaxTemp[C]\n";
	if (saveToFile) {
		file.open("tempSaveTest.txt");
		file << " Time[s]   MinTemp[C]   MaxTemp[C]\n";
	}
	
	for (int i = 0; i < N; i++) t0[i] = initialTemp;
	for (double t = simStepTime; t <= simTime; t += simStepTime) {
		vectorC = multiply(Ccaret, t0, N);
		Pcaret = add(aggrP, vectorC, N, N);
		t1 = g.elimination(Hcaret, Pcaret, N);

		double minTemp, maxTemp;
		tie(minTemp, maxTemp) = minMax(t1, N);
		t0 = t1;

		for (int i = 0; i < N; i++) {
			if (nodes[i].bc == 0)
				nodes[i].tempAtTime.push_back(t0[i]);
			else
				nodes[i].tempAtTime.push_back(t1[i]);
		}

		if (showMinMax) {
			pair<double, double> temps = minMax(t1, N);
			cout << setw(8) << t << "   " << setw(10) << temps.first << "   " << setw(10) << temps.second << "\n";
		}

		if (saveToFile) {
			pair<double, double> temps = minMax(t1, N);
			file << setw(8) << t << "   " << setw(10) << temps.first << "   " << setw(10) << temps.second << "\n";
		}
	}

	if (saveToFile) file.close();
	delete[] t0;
}

void Grid::plotHeatMap(Element4_2D el4) {
	if (!nodes) {
		cout << "File not found.\n";
		return;
	}

	if (nodes[0].tempAtTime.empty()) {
		cout << " Calculating ...";
		calcTemperature(el4);
		cout << " DONE\n";
	}

	int nW = 31, nH = 31; // @EDIT
	double xStep = width / ((double)nW - 1);
	double yStep = height / ((double)nH - 1);

	vector<float> temp(nN);	
	vector<double> xTicks, yTicks;
	const vector<string> labels = {};
	const vector<double> extent = { 0, width, 0, height};
	const map<string, string> options = { {"interpolation", "bilinear"}, {"cmap", "plasma"}, {"origin", "lower"}, {"aspect", "auto"} };
	const float* tptr = &(temp[0]);
	const int colors = 1;
	int passedTime = 0;
	
	for (size_t i = 0; i < nW; i++) xTicks.push_back(xStep * i);
	for (size_t i = 0; i < nH; i++) yTicks.push_back(yStep * i);
	for (int t = 0; t < nodes[0].tempAtTime.size(); t++) {
		int it = 0;

		for (size_t i = 0; i < nW; i++) {
			for (size_t j = 0; j < nH; j++) {
				temp.at((size_t)nW * j + i) = nodes[it].tempAtTime[t];
				it++;
			}
		}

		passedTime += simStepTime;
		PyObject* mat;
		plt::clf();
		plt::title("Rozklad temperatury po " + to_string(passedTime) + (passedTime > 1 ? " sekundach" : " sekundzie"));
		plt::imshow(tptr, nW, nH, colors, options, &mat, extent);
		//plt::xticks(xTicks, labels);
		//plt::grid(true);
		plt::colorbar(mat);
		plt::pause(0.1);
		Py_DECREF(mat);
	}
	
	plt::show();
	plt::close();
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
void Grid::readFromFile(string path) {	
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
			Grid::simTime = stod(temp);
		}
		else if (regex_match(line, std::regex("(SimulationStepTime|dt) \\d+"))) {
			temp = regex_replace(line, std::regex("SimulationStepTime|dt "), "$1");
			Grid::simStepTime = stod(temp);
		}
		else if (regex_match(line, std::regex("(Conductivity|K) \\d+"))) {
			temp = regex_replace(line, std::regex("Conductivity|K "), "$1");
			Grid::conductivity = stod(temp);
		}
		else if (regex_match(line, std::regex("Alfa \\d+"))) {
			temp = regex_replace(line, std::regex("Alfa "), "$1");
			Grid::alpha = stod(temp);
		}
		else if (regex_match(line, std::regex("Tot \\d+"))) {
			temp = regex_replace(line, std::regex("Tot "), "$1");
			Grid::ambientTemp = stod(temp);
		}
		else if (regex_match(line, std::regex("InitialTemp \\d+"))) {
			temp = regex_replace(line, std::regex("InitialTemp "), "$1");
			Grid::initialTemp = stod(temp);
		}
		else if (regex_match(line, std::regex("Density \\d+"))) {
			temp = regex_replace(line, std::regex("Density "), "$1");
			Grid::density = stod(temp);
		}
		else if (regex_match(line, std::regex("SpecificHeat \\d+"))) {
			temp = regex_replace(line, std::regex("SpecificHeat "), "$1");
			Grid::specificHeat = stod(temp);
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

	if (nodes) {
		this->width = nodes[0].x - nodes[nodesNo - 1].x;
		this->height = nodes[0].y - nodes[nodesNo - 1].y;

		double** temp = new double* [width];
		for (size_t i = 0; i < width; i++)	{
			temp[i] = new double[height];
		}


	}
}

/**
 * Ustawia w przekazanej tablicy sumê macierzy C dla ka¿dego punktu ca³kowania.
 *
 * @param globalH - tablica do zapisu zagregowanej macierzy H
 * @param globalC - tablica do zapisu zagregowanej macierzy C
 * @param globalP - tablica do zapisu zagregowanej macierzy P
 * @param currEl - element, na którym metoda ma operowaæ
 */
void Grid::aggregate(Element &currEl) {
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
			aggrH[rowIndex][colIndex] += currEl.H[i][j] + currEl.Hbc[i][j];
			aggrC[rowIndex][colIndex] += currEl.C[i][j];
		}

		int rowIndex = nodesId[i];
		aggrP[rowIndex] += currEl.P[i];
	}
}

/**
 * Ustawia w przekazanej tablicy sumê macierzy Hbc oraz wektor P dla ka¿dego punktu ca³kowania.
 *
 * @param Hbc - tablica, do której zapisaæ sumê macierzy Hbc z ka¿dego punktu ca³kowania
 * @param P - tablica, do której zapisaæ wektor P
 * @param alpha - wspó³czynnik ...
 * @param ambientTemp - temperatura otoczenia
 */
void Grid::calcHbc(Element4_2D &el4, Element &currEl, double alpha, double ambientTemp) {
	struct Side { double* pc1, * pc2; };
	double Npc1[4][4] = { 0 }, Npc2[4][4] = { 0 };

	// @TODO - zmieniæ na pobieranie z klasy Gauss w zale¿noœci od liczby pkt ca³kowania
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
 * @param el4 - Element4 2D z pochodnymi funkcji kszta³tu
 * @param H - macierz, do której zapisaæ wynik
 * @param k - wspó³czynnik przewodzenia
 */
void Grid::calcH(Element4_2D &el4, Element &currEl, double k) {
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
void Grid::calcC(Element4_2D &el4, Element &currEl, double c, double ro) {	
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
void Grid::calcJ(Element4_2D &el4, Element &currEl) {
	// i - punkt ca³kowania
	for (int i = 0; i < p; i++) {
		int nodeId = currEl.id[i] - 1;
		double x = nodes[nodeId].x;
		double y = nodes[nodeId].y;

		currEl.J[0][0] += el4.dXYdKsi(i, x);
		currEl.J[0][1] += el4.dXYdEta(i, x);
		currEl.J[1][0] += el4.dXYdKsi(i, y);
		currEl.J[1][1] += el4.dXYdEta(i, y);
	}	
}

/**
 * Uzupe³nia Jakobian dla ka¿dego punktu ca³kowania i-tego elementu.
 *
 * @param i - punkt ca³kowania
 * @param el4 - struktura z danymi o elemencie
 * @param grid - siatka, na której wykonywana jest procedura liczenia Jakobianu
 */
void Grid::jacobian(Element4_2D &el4, Element &currEl) {
	this->p = el4.p;
	this->calcJ(el4, currEl);
	currEl.detJ = detOfMatrix(currEl.J);
	inverseMatrix(currEl.invJ, currEl.J, currEl.detJ);

	// Wype³nij tablicê pochodnych funkcji N po X i Y - jakobian
	for (int j = 0; j < p; j++) {			// j - element
		for (int k = 0; k < p; k++) {		// k - punkt ca³kowania
			el4.dNdX[j][k] = currEl.invJ[0][0] * el4.dNdKsiMatrix[j][k] + currEl.invJ[0][1] * el4.dNdEtaMatrix[j][k];
			el4.dNdY[j][k] = currEl.invJ[1][0] * el4.dNdKsiMatrix[j][k] + currEl.invJ[1][1] * el4.dNdEtaMatrix[j][k];
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
 * Inicjalizacja macierzy do agregacji.
 */
void Grid::initMatrices() {
	const int N = nN;
	aggrP = new double[N];
	aggrH = new double* [N];
	aggrC = new double* [N];

	for (int i = 0; i < N; i++) {
		aggrP[i] = 0;
		aggrH[i] = new double[N];
		aggrC[i] = new double[N];

		for (int j = 0; j < N; j++) {
			aggrH[i][j] = 0;
			aggrC[i][j] = 0;
		}
	}
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

/**
 * Wyœwietla w oknie konsoli wartoœci, dla których przeprowadzana jest symulacja.
 */
void Grid::printSimulationData() {
	cout << " SimulationTime " << simTime << "\n";
	cout << " SimulationStepTime " << simStepTime << "\n";
	cout << " Conductivity " << conductivity << "\n";
	cout << " Alpha " << alpha << "\n";
	cout << " InitialTemp " << initialTemp << "\n";
	cout << " AmbientTemp " << ambientTemp << "\n";
	cout << " Density " << density << "\n";
	cout << " SpecificHeat " << specificHeat << "\n";
	cout << " Nodes number " << nN << "\n";
	cout << " Elements number " << nE << "\n\n";
	//printNodesCoords();
	//printElementsNodes();
	//printBCNodes();
}