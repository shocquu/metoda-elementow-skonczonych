#include <iostream>
#include "utils/functions.h"
#include "grid/grid.h"

#define LAB_NO 11

int main() {	
	Element4_2D el4(2);

	#if LAB_NO == 11		
		//Grid grid(el4, "data/Test4_31_31_trapez.txt");
		Grid grid(el4, "data/MES_31_31_v2.txt");
		grid.plotHeatMap(el4);

	#elif LAB_NO == 10
		Grid grid;
		grid.launch(el4, "data/MES_31_31_v2.txt");

	#elif LAB_NO >= 6 && LAB_NO < 10
		Grid grid(0.1, 0.1, 4, 4);
		grid.launch(el4, 500, 50, 25, 300, 7800, 700, 100, 1200);

	#endif
	system("pause");
	return 0;
}