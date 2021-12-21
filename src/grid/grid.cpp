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

	setNodesCoords();
	setElementsNodes();
	setNodesBC();
}

Grid::~Grid() {
	// ERROR
	//delete[] this->nodes;
	//delete[] this->elements;
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
