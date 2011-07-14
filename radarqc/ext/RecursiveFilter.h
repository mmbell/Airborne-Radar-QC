/*
 *  RecursiveFilter.h
 *
 * Copyright 2011 Michael Bell.
 *  All rights reserved.
 *
 */

#ifndef RECURSIVEF_H
#define RECURSIVEF_H
#include "precision.h"


class RecursiveFilter
{
	
public:
	RecursiveFilter(const int& fOrder, const int& fLengthScale);
	~RecursiveFilter();
	bool filterArray(real* array, const int& arrLength);
	
private:
	int order;
	int lengthScale;
	real beta;
	real alpha[5];
	real Sn[5][5];
	void getFilterCoefficients();
	real factorial(const real& max); 		
	void solveBC(real* A, real* B);	
};

#endif
