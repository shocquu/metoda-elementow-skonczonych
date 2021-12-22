#include <iostream>

#include "utils/functions.h"
#include "grid/element4_2d.h"
#include "grid/grid.h"

#include "../include/matplotlib-cpp/matplotlibcpp.h"
#include "../include/meshgrid.hpp"
//#include <vector>

using namespace std;

namespace plt = matplotlibcpp;

#define LAB_NO 11

int main() {	
	Element4_2D el4(2);

	#if LAB_NO == 11
        Grid grid(0.1, 0.1, 31, 31);
        grid.launch(el4, 500, 50, 25, 300, 7800, 700, 100, 1200);

        int it = 0;
        vector<vector<double>> x, y, z;
        for (int i = 0; i < grid.nW; i += 1) {
            vector<double> x_row, y_row, z_row;
            for (int j = 0; j < grid.nH; j += 1) {
                double t = grid.elements[it].temperature;
                x_row.push_back(i);
                y_row.push_back(j);
                z_row.push_back(t);
                it++;
            }
            x.push_back(x_row);
            y.push_back(y_row);
            z.push_back(z_row);
        }

        //plt::colorbar();
        plt::plot_surface(x, y, z);
        plt::show();

	#elif LAB_NO == 10
		Grid().launch(el4, "data/MES_31_31_v2.txt");

	#elif LAB_NO >= 6 && LAB_NO < 10
		Grid grid(0.1, 0.1, 4, 4);
		grid.launch(el4, 500, 50, 25, 300, 7800, 700, 100, 1200);
		
	#endif
	system("pause");
	return 0;
}