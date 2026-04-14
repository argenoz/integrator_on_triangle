#pragma once
#include <fraction/Fraction.h>
#include <vector>
#include <utility>

class Integrator
{
private:
	Fraction ax;
	Fraction ay;
	Fraction bx;
	Fraction by;
	Fraction cx;
	Fraction cy;

	Fraction Kx;
	Fraction Ky;
	Fraction acx;
	Fraction acy;
	Fraction J;

public:
	std::vector<std::vector<Fraction>> *Binome;
	//Integrator();
	Integrator(std::vector<std::pair<Fraction, Fraction>> triangle);
	Fraction Integrate(int n, int m);
	Fraction Cnk(int n, int k);

};