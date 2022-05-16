#pragma once
#include"DupireLocalVolatilitySurface.h"
#include "Bins.h"
#include <utility>
using Pair = std::pair<double, double>;

// Designed class for models with local volatility like Heston and SABR Local stochastic volatility Models 
class Model_with_vol 
{
protected:
	double _kappa;
	double _risk_free_rate;
	double _correlation;
	double _vol_of_vol;
	double _theta;
	Pair _init_spot_variance;
	DupireLocalVolatilitySurface _dupire_local_volatility;
	bins* bin;
public:
	Model_with_vol(const double& risk_free_rate, const double& correlation, const Pair& init_spot_variance, const double& kappa, const double& vol_of_vol, const double& theta, const DupireLocalVolatilitySurface& dupirelocalvolatility, const bins& b);
	virtual~Model_with_vol();
	Model_with_vol(const Model_with_vol& model);
	Model_with_vol& operator=(const Model_with_vol& model);
	Pair drift_pair(const double& time, const Pair& spot_variance);
	double local_volatility(const double& time, const Pair& spot_variance, vector<Pair>& spot_function_variance, const double& spot) const;
	virtual Pair diffusion_pair(const double& time, const Pair& spot_variance,vector<Pair>& spot_function_variance, const int& r) const = 0;
	// virtuals functions
	virtual double psi_function(const double& variance) const = 0;
	virtual double variance_drift(const double& time, const double& variance) const = 0;
	virtual double variance_diffusion(const double& time, const double& variance) const = 0;
	virtual Model_with_vol* clone() const = 0;
	// getter methods
	double kappa() const;
	double risk_free_rate() const;
	double correlation() const;
	double vol_of_vol() const;
	double theta() const;
	Pair init_spot_variance() const;
	bins* get_bin() const;
	DupireLocalVolatilitySurface get_dupire_vol() const;
};

// Heston local stochastic volatility Model
class Heston_local_sto_vol_Model :public Model_with_vol
{
public:
	Heston_local_sto_vol_Model(const double& risk_free_rate, const double& correlation, const Pair& init_spot_variance, const double& kappa, const double& vol_of_vol, const double& theta, const DupireLocalVolatilitySurface& dupirelocalvolatility, const bins& b);
	Pair diffusion_pair(const double& time, const Pair& spot_variance, vector<Pair>& spot_function_variance, const int& r) const;
	double psi_function(const double& variance) const override;
	double variance_drift(const double& time, const double& variance) const override;
	double variance_diffusion(const double& time, const double& variance) const override;
	Heston_local_sto_vol_Model* clone() const override;
};

class BSCall
{
public:
	BSCall(double r_, double d_,double T, double Spot_,double Strike_);
	double Price(double Vol) const;
	double Vega(double Vol) const;
private:
	double r;
	double d;
	double T;
	double Spot;
	double Strike;
};
