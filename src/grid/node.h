#pragma once
/**
 * Wêze³ bêd¹cy czêœci¹ elementu.
 *
 * @param x, y - wspó³rzêdne wêz³ów
 * @param bc - informacja o obecnoœci warunku brzegowego
 */
struct Node {
	int id = -1;
	short int bc = 0;
	double x = 0, y = 0;	
};