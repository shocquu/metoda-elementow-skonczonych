#pragma once
/**
 * W�ze� b�d�cy cz�ci� elementu.
 *
 * @param x, y - wsp�rz�dne w�z��w
 * @param bc - informacja o obecno�ci warunku brzegowego
 */
struct Node {
	int id = -1;
	short int bc = 0;
	double x = 0, y = 0;	
};