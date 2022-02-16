#include "helpers.h"

using namespace std;

pair<double, double> minMax(double* array, const int N) {
	double min = 9999, max = 0;
	for (int i = 0; i < N; i++) {
		if (array[i] < min) min = array[i];
		if (array[i] > max) max = array[i];
	}

	return { min, max };
}

void printBanner() {
	cout << string(79, '.') << endl;
	cout << string(22, '.') << " ##\\      ##\\ ########\\  ######\\ " << string(24, '.') << endl;
	cout << string(22, '.') << " ###\\    ### |##  _____|##  __##\\ " << string(23, '.') << endl;
	cout << string(22, '.') << " ####\\  #### |## |      ## /  \\__| " << string(22, '.') << endl;
	cout << string(22, '.') << " ##\\##\\## ## |#####\\    \\#####\\  " << string(24, '.') << endl;
	cout << string(22, '.') << " ## \\###  ## |##  __|    \\____##\\  " << string(22, '.') << endl;
	cout << string(22, '.') << " ## |\\#  /## |## |      ##\\   ## | " << string(22, '.') << endl;
	cout << string(22, '.') << " ## | \\_/ ## |########\\ \\######  | " << string(22, '.') << endl;
	cout << string(22, '.') << " \\__|     \\__|\\________| \\______/  " << string(22, '.') << endl;
	cout << string(79, '.') << "\n\n";
}