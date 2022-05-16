#include "Modeles.h"
#include "Blackscholesformulas.h"

Model_with_vol::Model_with_vol(const double& risk_free_rate, const double& correlation, const Pair& init_spot_variance, const double& kappa, const double& vol_of_vol, const double& theta, const DupireLocalVolatilitySurface& dupirelocalvolatility,const bins&b) :_risk_free_rate(risk_free_rate), _correlation(correlation),_init_spot_variance(init_spot_variance),_kappa(kappa),_vol_of_vol(vol_of_vol),_theta(theta),_dupire_local_volatility(dupirelocalvolatility),bin(b.clone())
{
}

Model_with_vol::~Model_with_vol()
{
	delete bin;
}

Model_with_vol::Model_with_vol(const Model_with_vol& model):_risk_free_rate(model._risk_free_rate),_correlation(model._correlation),_init_spot_variance(model._init_spot_variance),_kappa(model._kappa),_vol_of_vol(model._vol_of_vol),_theta(model._theta),_dupire_local_volatility(model._dupire_local_volatility),bin(model.bin->clone())
{
}

Model_with_vol& Model_with_vol::operator=(const Model_with_vol& model)
{
	if (this != &model)
	{
		this->_risk_free_rate = model._risk_free_rate;
		this->_correlation = model._correlation;
		this->_init_spot_variance = model._init_spot_variance;
		this->_kappa = model._kappa;
		this->_vol_of_vol = model._vol_of_vol;
		this->_theta = model._theta;
		this->_dupire_local_volatility = model._dupire_local_volatility;
		delete bin;
		this->bin = model.bin->clone();
	}
	return*this;
}

Pair Model_with_vol::drift_pair(const double& time, const Pair& spot_variance)
{
	double drift1 = _risk_free_rate * spot_variance.first;
	double drift2 = variance_drift(time, spot_variance.second);
	return Pair(drift1, drift2);
}


double Model_with_vol::local_volatility(const double& time, const Pair& spot_variance,vector<Pair>& spot_function_variance,const double &spot) const
{
	return pow(_dupire_local_volatility.local_volatility(time, spot_variance.first),2) / bin->expectation(spot_function_variance,spot);
}

bins* Model_with_vol::get_bin() const
{
	return bin;
}

DupireLocalVolatilitySurface Model_with_vol::get_dupire_vol() const
{
	return _dupire_local_volatility;
}

Heston_local_sto_vol_Model::Heston_local_sto_vol_Model(const double& risk_free_rate, const double& correlation, const Pair& init_spot_variance, const double& kappa, const double& vol_of_vol, const double& theta, const DupireLocalVolatilitySurface& dupirelocalvolatility, const bins& b):Model_with_vol(risk_free_rate,correlation,init_spot_variance,kappa,vol_of_vol,theta, dupirelocalvolatility,b)
{
}

Pair Heston_local_sto_vol_Model::diffusion_pair(const double& time, const Pair& spot_variance, vector<Pair>& spot_function_variance, const int& r) const
{
	double diffusion1 = psi_function(spot_variance.second) * spot_variance.first * local_volatility(time, spot_variance, spot_function_variance, r);
	double diffusion2 = variance_diffusion(time, spot_variance.second);
	return Pair(diffusion1, diffusion2);
}

double Heston_local_sto_vol_Model::psi_function(const double& variance) const
{
	return sqrt(variance);
}

double Heston_local_sto_vol_Model::variance_drift(const double& time, const double& variance) const
{
	return _kappa*(_theta-variance);
}

double Heston_local_sto_vol_Model::variance_diffusion(const double& time, const double& variance) const
{
	return _vol_of_vol*psi_function(variance);
}

Heston_local_sto_vol_Model* Heston_local_sto_vol_Model::clone() const
{
	return new Heston_local_sto_vol_Model(*this);
}

BSCall::BSCall(double r_, double d_,double T_, double Spot_,double Strike_):r(r_), d(d_),T(T_), Spot(Spot_),Strike(Strike_)
{}
double BSCall::Price(double Vol) const
{
	return BlackScholesCall(Spot, Strike, r, d, Vol, T);
}
double BSCall::Vega(double Vol) const
{
	return BlackScholesCallVega(Spot, Strike, r, d, Vol, T);
}

double Model_with_vol::kappa() const
{
	return _kappa;
}

double Model_with_vol::risk_free_rate() const
{
	return _risk_free_rate;
}

double Model_with_vol::correlation() const
{
	return _correlation;
}

double Model_with_vol::vol_of_vol() const
{
	return _vol_of_vol;
}

double Model_with_vol::theta() const
{
	return _theta;
}

Pair Model_with_vol::init_spot_variance() const
{
	return _init_spot_variance;
}
