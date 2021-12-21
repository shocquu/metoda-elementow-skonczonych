#include <iostream>

#include "utils/functions.h"
#include "utils/matrix.h"
#include "utils/gauss.h"
#include "grid/grid.h"
#include "grid/element4_2d.h"

using namespace std;

#define LAB_NO 9
 
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
	delete[] globalC;
	//delete globalP;
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
		grid.readFromFile("data/MES_31_31_v2.txt", simTime, simStepTime, conductivity, alpha, density, specificHeat, initialTemp, ambientTemp);
		const int N = grid.nN;
		initGlobalMatrices(N);

		cout << "SimulationTime " << simTime << "\n";
		cout << "SimulationStepTime  " << simStepTime << "\n";
		cout << "Conductivity " << conductivity << "\n";
		cout << "Alfa " << alpha << "\n";
		cout << "Tot " << ambientTemp << "\n";
		cout << "InitialTemp " << initialTemp << "\n";
		cout << "Density " << density << "\n";
		cout << "SpecificHeat " << specificHeat << "\n";
		cout << "Nodes number " << grid.nN << "\n";
		cout << "Elements number " << grid.nE << "\n\n";
		//grid.printNodesCoords();
		//grid.printElementsNodes();
		//grid.printBCNodes();

		for (int i = 0; i < grid.nE; i++) {
			Element currEl = grid.elements[i];

			el4.jacobian(grid, currEl);
			el4.calcH(currEl.H, conductivity);
			el4.calcC(currEl.C, specificHeat, density);
			el4.calcHbc(currEl.Hbc, currEl.P, grid, currEl, i, alpha, ambientTemp);
			el4.aggregate(globalH, globalC, globalP, currEl);
		}

		double** Ccaret, ** Hcaret;
		double* Pcaret, *vectorC, * t1;
		double* t0 = new double[N];
		for (int i = 0; i < N; i++) t0[i] = initialTemp;

		cout << " Time[s]   MinTemp[s]   MaxTemp[s]\n";

		for (double dT = 0, i = 0; dT <= simTime; dT += simStepTime, i++) {
			Ccaret = dT > 0 ? divide(globalC, dT, N, N) : divide(globalC, 1, N, N);
			Hcaret = add(globalH, Ccaret, N, N);
			Hcaret = add(globalH, Ccaret, N, N);
			vectorC = multiply(Ccaret, t0, N);
			Pcaret = add(globalP, vectorC, N, N);
			t1 = gauss.elimination(Hcaret, Pcaret, N);
			t0 = t1;

			double minTemp = 9999, maxTemp = 0;

			for (int k = 0; k < N; k++) {
				if (t1[k] < minTemp) minTemp = t1[k];
				if (t1[k] > maxTemp) maxTemp = t1[k];
			}

			//minTemp = t1[N - 2];
			//maxTemp = t1[N - 1];
			cout << setw(8) << dT << "   " << setw(10) << minTemp << "   " << setw(10) << maxTemp << "\n";
		}

		delete[] t0;
		destroyGlobalMatrices(N);

	#elif LAB_NO >= 6 && LAB_NO < 10
		Grid grid(0.1, 0.1, 4, 4);
		const int N = grid.nN;
		initGlobalMatrices(N);

		for (int i = 0; i < grid.nE; i++) {
			Element currEl = grid.elements[i];

			el4.jacobian(grid, currEl);
			el4.calcH(currEl.H, 25);
			el4.calcC(currEl.C, 700, 7800, currEl, i, grid);
			el4.calcHbc(currEl.Hbc, currEl.P, grid, currEl, i, 300, 1200);
			el4.aggregate(globalH, globalC, globalP, currEl);
		}
				
		double** Ccaret;
		double** Hcaret;
		double* Pcaret, *vectorC, *t1;
		double* t0 = new double[N];
		for (int i = 0; i < N; i++) t0[i] = 100;

		cout << " Time[s]   MinTemp[s]   MaxTemp[s]\n";

		for (size_t t = 50, i = 0; t <= 500; t += 50, i++) {
			Ccaret = t > 0 ? divide(globalC, 50, N, N) : divide(globalC, 1, N, N);
			Hcaret = add(globalH, Ccaret, N, N);
			//Hcaret = add(globalH, Ccaret, N, N);
			vectorC = multiply(Ccaret, t0, N);
			Pcaret = add(globalP, vectorC, N, N);
			t1 = gauss.elimination(Hcaret, Pcaret, N);
			t0 = t1;

			/*cout << string(89, '_') << " Iteration " << i << " " << string(89, '_') << "\n";
			cout << string(85, '_') << " Matrix ([H]+[C]/dT) " << string(85, '_') << "\n";
			printMatrix(Hcaret, N, N);
			cout << string(80, '_') << " Vector ({P}+{[C]/dT}*{T0}) " << string(81, '_') << "\n";
			printMatrix(Pcaret, N);
			cout << "\n";*/

			double minTemp = 9999, maxTemp = 0;
			for (size_t k = 0; k < N; k++)	{
				if (t1[k] < minTemp) minTemp = t1[k];
				if (t1[k] > maxTemp) maxTemp = t1[k];
			}

			//cout.precision(2); // fixed << minTemp
			cout << setw(8) << t << "   " << setw(10) << minTemp << "   " << setw(10) << maxTemp << "\n";
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
		cout << setw(16) << " " << ">>> Funkcje ksztaltu dN/dKsi <<<\n";
		printMatrix4x4(el4.ksiMatrix);
		cout << setw(16) << " " << ">>> Funkcje ksztaltu dN/dEta <<<\n";
		printMatrix4x4(el4.etaMatrix);
	#elif LAB_NO == 2
		Gauss gauss;
		int n = 3;
		double interval1d = gauss.quadrature1d(-1, 1, (vFunctionCall)f, n);		// 15.3333 OK | https://www.wolframalpha.com/input/?i2d=true&i=Integrate%5B5Power%5Bx%2C2%5D%2B3x%2B6%2C%7Bx%2C-1%2C1%7D%5D
		double interval2d = gauss.quadrature2d((vFunctionCall2)f, n);		    // 26.2222 OK | https://www.wolframalpha.com/input/?i2d=true&i=Integrate%5BIntegrate%5B5Power%5Bx%2C2%5DPower%5By%2C2%5D+%2B+3xy+%2B+6%2C%7Bx%2C-1%2C1%7D%5D%2C%7By%2C-1%2C1%7D%5D
		cout << "Kwadratura 1D dla " << n << "p: " << interval1d << endl;
		cout << "Kwadratura 2D dla " << n << "p: " << interval2d << endl;
	#elif LAB_NO == 1
		Grid grid(0.2f, 0.1f, 5, 4);
		grid.printNodesCoords();
		grid.printElementsNodes();
	#endif
	system("pause");
	return 0;
}