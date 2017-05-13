// InhomogenField.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <string>
#include <array>
#include <iostream>
using namespace std;


int main()
{
	double pi = 3.14159265359;
	double l = 0.5;
	int option = 2;
	double lambda = 2;
	string system = "SFS";
	double alphaL[1] = { 0 };
	double alphaR[1] = { 0 };
	double phi = 0;
	int dimy = 2;

	int nx = 50;
	int ny = 400;
	int ntheta = 100;

	double dx = 1 / nx;
	double dy = 2 * dimy / ny;
	double dv = pi*l*l / 4;
	if (dv > 4) {dv = 0.1;}


	fun = @(theta)localCurrent(x(i, j), y(i, j), theta, l, phi, option, lambda(k), system, alphaL(k_alpha), alphaR(k_alpha));


    return 0;
}

