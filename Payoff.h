#pragma once

// Designed class for options's payoffs construction
class PayOff
{
public:
	PayOff(const double & strike);
	virtual double operator()(const double Spot) const = 0;// Gives the payoff
	virtual PayOff* clone() const = 0;
	virtual ~PayOff() {}
	//Getter method
	double strike();
protected:
	double _strike;// strike of the option
};

// Vanilla call's Payoff
class PayOffCall : public PayOff
{
public:
	PayOffCall(const double& Strike_);
	double operator()(const double Spot) const override;
	virtual ~PayOffCall() {}
	virtual PayOffCall* clone() const override;
};

//Vanilla put's payoff
class PayOffPut : public PayOff
{
public:
	PayOffPut(const double& Strike_);
	double operator()(const double Spot) const override;
	virtual ~PayOffPut() {}
	PayOffPut* clone() const override;
};
