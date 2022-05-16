#include"Blackscholesformulas.h"
#include "Gaussienne.h"
#include <cmath>
double BlackScholesCall(double Spot,double Strike,double r,double d,double Vol,double Expiry)
{
	double standardDeviation = Vol * sqrt(Expiry);
	double moneyness = log(Spot / Strike);
	double d1 = (moneyness + (r - d) * Expiry +
		0.5 * standardDeviation * standardDeviation) / standardDeviation;
	double d2 = d1 - standardDeviation;
	return Spot * exp(-d * Expiry) * NormalCDF(d1) -Strike * exp(-r * Expiry) * NormalCDF(d2);
}

double BlackScholesCallVega(double Spot,double Strike,double r,double d,double Vol,double Expiry)
{
	double standardDeviation = Vol * sqrt(Expiry);
	double moneyness = log(Spot / Strike);
	double d1 = (moneyness + (r - d) * Expiry +
		0.5 * standardDeviation * standardDeviation) / standardDeviation;
	return Spot * exp(-d * Expiry) * sqrt(Expiry) * NormalCDF(d1);
}