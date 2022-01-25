#include "utils/functions.h"
#include "grid/grid.h"

#define MODE 2

int main() {
	SimulationData dataset;
	dataset.specificHeat = 700;
	dataset.ambientTemp = 1200;
	dataset.initialTemp = 100;
	dataset.conductivity = 25;
	dataset.density = 7800;
	dataset.alpha = 300;
	dataset.simStepTime = 50;
	dataset.simTime = 500;

	#if MODE == 3
		Grid("data/MES_31_31_v2.txt").heatMap();

	#elif MODE == 2
		Grid("data/Test4_31_31_trapez.txt").start(true);

	#elif MODE == 1
		Grid grid(0.1, 0.1, 4, 4);
		grid.launch(dataset);

	#endif

	system("pause");
	return 0;
}