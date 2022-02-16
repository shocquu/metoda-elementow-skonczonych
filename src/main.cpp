#include "utils/functions.h"
#include "grid/grid.h"

#define MODE 1

int main() {
	SimulationData dataset;

	#if MODE == 3
		Grid("data/MES_31_31_v2.txt").heatMap();

	#elif MODE == 2
		//dataset.simStepTime = 1;
		//dataset.simTime = 20;
		//Grid().launch("data/Test4_31_31_trapez.txt");
		Grid("data/Test4_31_31_trapez.txt").start(); // .heatMap();

	#elif MODE == 1
		Grid(0.1, 0.1, 4, 4).start(dataset);
		//Grid(0.1, 0.1, 4, 4).start(dataset, true);

	#endif

	system("pause");
	return 0;
}