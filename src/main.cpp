#include <iostream>

#include "utils/functions.h"
#include "utils/matrix.h"
#include "utils/gauss.h"
#include "grid/grid.h"
#include "grid/element4_2d.h"

using namespace std;

#define LAB_NO 9

pair<double, double> minMax(double* array, const int N) {
	double min = 9999, max = 0;
	for (size_t i = 0; i < N; i++) {
		if (array[i] < min) min = array[i];
		if (array[i] > max) max = array[i];
	}

	return make_pair(min, max);
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

		double** Ccaret = divide(globalC, 50, N, N);
		double** Hcaret = add(globalH, Ccaret, N, N);
		double* Pcaret, * vectorC, * t1;
		double* t0 = new double[N];
		for (int i = 0; i < N; i++) t0[i] = initialTemp;

		cout << " Time[s]   MinTemp[s]   MaxTemp[s]\n";

		for (size_t t = 50, i = 0; t <= 500; t += 50, i++) {
			vectorC = multiply(Ccaret, t0, N);
			Pcaret = add(globalP, vectorC, N, N);
			t1 = gauss.elimination(Hcaret, Pcaret, N);
			t0 = t1;
			pair<double, double> temps = minMax(t1, N);
			cout << setw(8) << t << "   " << setw(10) << minTemp << "   " << setw(10) << maxTemp << "\n";
		}

		delete[] t0;
		destroyGlobalMatrices(N);

	#elif LAB_NO >= 6 && LAB_NO < 10
		Grid grid(0.1, 0.1, 4, 4);
		const int N = grid.nN;

		for (int i = 0; i < grid.nE; i++) {
			Element currEl = grid.elements[i];

			grid.jacobian(el4, currEl);
			grid.calcH(el4, currEl, 25);
			grid.calcC(el4, currEl, 700, 7800);
			grid.calcHbc(el4, currEl, i, 300, 1200);
			grid.aggregate(currEl);

			//printMatrix(currEl.H);
			//printMatrix(currEl.C);
		}
				
		double** Ccaret = divide(grid.globalC, 50, N, N);
		double** Hcaret = add(grid.globalH, Ccaret, N, N);
		double* Pcaret, *vectorC, *t1;
		double* t0 = new double[N];
		for (int i = 0; i < N; i++) t0[i] = 100;
		
		//printMatrix(grid.globalH, N, N);

		cout << " Time[s]   MinTemp[s]   MaxTemp[s]\n";

		for (size_t t = 50, i = 0; t <= 500; t += 50, i++) {
			vectorC = multiply(Ccaret, t0, N);
			Pcaret = add(grid.globalP, vectorC, N, N);
			t1 = gauss.elimination(Hcaret, Pcaret, N);
			t0 = t1;

			/*cout << string(89, '_') << " Iteration " << i << " " << string(89, '_') << "\n";
			cout << string(85, '_') << " Matrix ([H]+[C]/dT) " << string(85, '_') << "\n";
			printMatrix(Hcaret, N, N);
			cout << string(80, '_') << " Vector ({P}+{[C]/dT}*{T0}) " << string(81, '_') << "\n";
			printMatrix(Pcaret, N);
			cout << "\n"; */

			pair<double, double> temps = minMax(t1, N);
			//cout.precision(2); // fixed << minTemp
			cout << setw(8) << t << "   " << setw(10) << temps.first << "   " << setw(10) << temps.second << "\n";
		}

		delete[] t0;
	#endif
	system("pause");
	return 0;
}