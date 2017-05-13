#include <iostream>
#include <gsl/gsl_sf_hyperg.h>
// Dette er et helt enkelt program som du kan kopiere og bruke i
// denne oevingen.
using namespace std;
int main() {
	cout << "Hello World!";
	cin.get();
	double a = 1;
	double b = 2;
	double c = 3;
	double res = gsl_sf_hyperg_1F1(a, b, c);
	return 0;
}