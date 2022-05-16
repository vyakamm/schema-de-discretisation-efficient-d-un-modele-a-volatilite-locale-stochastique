#include"Payoff.h"
#include <minmax.h>
PayOff::PayOff(const double& strike) : _strike(strike)
{
}
double PayOff::strike()
{
	return _strike;
}
PayOffCall::PayOffCall(const double& strike) : PayOff(strike)
{
}
double PayOffCall::operator()(const double Spot) const
{
	return max(Spot - _strike, 0.0);
}

PayOffCall* PayOffCall::clone() const
{
	return new PayOffCall(*this);
}

PayOffPut::PayOffPut(const double& strike): PayOff(strike)
{
}

double PayOffPut::operator()(const double Spot) const
{
	return max(_strike - Spot, 0.0);
}

PayOffPut* PayOffPut::clone() const
{
	return new PayOffPut(*this);
}
