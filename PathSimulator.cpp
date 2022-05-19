#include "PathSimulator.h"

PathSimulator::PathSimulator(const Model_with_vol& model, const vector<double> time_points) :model(model.clone()), _time_points(time_points)
{
}

PathSimulator::PathSimulator(const PathSimulator & path_simulator)
	: model((path_simulator.model)->clone()), _time_points(path_simulator._time_points)
{
}

PathSimulator& PathSimulator::operator=(const PathSimulator& path_simulator)
{
	if (this != &path_simulator)
	{
		_time_points = path_simulator._time_points;
		delete model;
		model = path_simulator.model->clone();
	}
	return *this;
}

PathSimulator::~PathSimulator()
{
	delete model;
}



Pair PathSimulator::path_maturity(size_t time_idx)
{
	Pair p=model->init_spot_variance();
	size_t number_of_simulations = model->get_bin()->number_of_simulations();
	vector<Pair> samples(number_of_simulations, p);
	vector<Pair> v(number_of_simulations);
	
	// for each time's step
	for (size_t i = 0; i < time_idx - 1; i++)
	{
		v = samples;
		// samples for the computation of the conditional expectation
		for (size_t j = 0; j < number_of_simulations; j++)
		{
			samples[j] = next_step(i+1, p, v);
		}
		// simulating the next step
		p = next_step(i + 1, p, v);
	}
	return p;
}

Model_with_vol* PathSimulator::get_model()
{
	return model;
}


PathSimulator* PathSimulatorEuler::clone() const
{
	return new PathSimulatorEuler(*this);
}

PathSimulatorEuler::PathSimulatorEuler(const Heston_local_sto_vol_Model& model, const vector<double> time_points):PathSimulator(model, time_points)
{
}

Pair PathSimulatorEuler::next_step(const size_t& time_idx, const Pair& spot_variance, const vector<Pair>& spot_variance_sample) const
{
	if (time_idx == 0)
		return model->init_spot_variance();
	else
	{
		double next1, next2;
		//Gives the discretization step
		double delta_t = _time_points[time_idx] - _time_points[time_idx - 1];
		// Asset's diffusion brownian
		double z1 = Gausienne();
		// Gives the asset discretization
		if (time_idx == 1)
			next1 = spot_variance.first + model->risk_free_rate() * delta_t * spot_variance.first + pow((model->get_dupire_vol()).local_volatility(_time_points[time_idx - 1], spot_variance.first), 2) / model->psi_function(spot_variance.second) * spot_variance.first * model->psi_function(spot_variance.second) * sqrt(delta_t) * z1;
		else
		{
			vector<Pair>v(spot_variance_sample.size());
			for (size_t i = 0; i < v.size(); i++)
			{
				v[i].second = pow(model->psi_function(spot_variance_sample[i].second), 2);
				v[i].first = spot_variance_sample[i].first;
			}
			next1 = spot_variance.first + model->risk_free_rate() * delta_t * spot_variance.first + model->local_volatility(_time_points[time_idx - 1], spot_variance, v, 2);
		}
		// Volatility's diffusion Brownian
		double z2 = Gausienne();
		double z3 = model->correlation() * z1 + sqrt(1 - pow(model->correlation(), 2)) * z2;
		//Gives the volatility discretization 
		next2 = spot_variance.second + model->variance_drift(_time_points[time_idx-1], spot_variance.second) * delta_t + model->variance_diffusion(_time_points[time_idx-1], spot_variance.second) * delta_t*z3;
		return Pair(next1, next2);
	}
}

PathSimulatorQE* PathSimulatorQE::clone() const
{
	return new PathSimulatorQE(*this);
}

PathSimulatorQE::PathSimulatorQE(const Model_with_vol& model, const vector<double> time_points):PathSimulator(model,time_points)
{
}

Pair PathSimulatorQE::next_step(const size_t& time_idx, const Pair& spot_variance, const vector<Pair>& spot_variance_sample) const
{
	if (time_idx == 0)
		return model->init_spot_variance();
	else
	{
		//Gives the discretization step
		double delta_t = _time_points[time_idx] - _time_points[time_idx - 1];
		// Brownian of the volatility
		double z1 = Gausienne();
		// QE algorithm to simulate de volatility of the Model
		double m = model->theta() + (spot_variance.second - model->theta()) * exp(-model->kappa() * delta_t);
		double s1= spot_variance.second*pow(model->vol_of_vol(),2)* exp(-model->kappa() * delta_t)/model->kappa();
		double s2 = model->theta() * pow(model->vol_of_vol(), 2) / (2 * model->kappa());
		double s3 = 1 - exp(-model->kappa() * delta_t);
		double s_quarred = s1 * s3 + s2 * pow(s3, 2);
		double psi = s_quarred / pow(m, 2);
		double u_v = rand() / (RAND_MAX + 1.);
		double next2{};
		if (psi <= 1.5)
		{
			double b = 2 * 1 / psi - 1 + sqrt(2 / psi) * sqrt(2 / psi - 1);
			double a = m / (1 + pow(b, 2));
			double z_v = NormalCDFInverse(u_v);
			next2 = a * pow(b + z_v, 2);
		}
		else
		{
			double p = (psi - 1) / (psi + 1);
			double beta = (1 - p) / m;
			if (u_v >= 0 && u_v <= p)
				next2 = 0;
			else if (u_v > p && u_v <= 1)
				next2=(1 / beta) * log((1 - p) / (1 - u_v));
		}
		// Bownian of the asset's diffusion
		double z2 = model->correlation() * z1 + sqrt(1 - pow(model->correlation(), 2)) * Gausienne();
		// Transforming the stock in log in other to simulate
		double current_log_spot = log(spot_variance.first);
		// Gives de quantity to condition in other to compute the conditionnal expectation
		
		vector<Pair>v(model->get_bin()->number_of_simulations());
		for (int i = 0; i < v.size(); i++)
		{
			v[i].second = pow(model->psi_function(spot_variance_sample[i].second), 2);
			v[i].first = spot_variance_sample[i].first;
		}
		double c1 = model->kappa() * delta_t - 1;
		double rho1 = sqrt(1 - pow(model->correlation(), 2));
		// QE algorithm to simulate the stock
		if (time_idx == 1)
		{
			double r1 = sqrt(pow((model->get_dupire_vol()).local_volatility(_time_points[time_idx - 1], spot_variance.first), 2) / model->psi_function(spot_variance.second)) * (next2 - model->kappa() * model->theta() * delta_t + spot_variance.second * c1);
			double r2 = sqrt(pow((model->get_dupire_vol()).local_volatility(_time_points[time_idx - 1], spot_variance.first), 2) / model->psi_function(spot_variance.second) * spot_variance.second * delta_t);
			double current_log_spot_next = current_log_spot + model->risk_free_rate()*delta_t - 0.5 * pow((model->get_dupire_vol()).local_volatility(_time_points[time_idx - 1], spot_variance.first), 2) / model->psi_function(spot_variance.second) * spot_variance.second * delta_t + (model->correlation() / model->vol_of_vol()) * r1 + rho1 * r2 * z2;
			double next1 = exp(current_log_spot_next);
			return Pair(next1, next2);
		}
		else
		{
			double r1 = sqrt(model->local_volatility(_time_points[time_idx - 1], spot_variance, v, spot_variance.first)) * (next2 - model->kappa() * model->theta() * delta_t + spot_variance.second * c1);
			double r2 = sqrt(model->local_volatility(_time_points[time_idx - 1], spot_variance, v, spot_variance.first) * spot_variance.second * delta_t);
			double current_log_spot_next = current_log_spot + model->risk_free_rate() - 0.5 * model->local_volatility(_time_points[time_idx - 1], spot_variance, v, spot_variance.first) * spot_variance.second * delta_t + (model->correlation() / model->vol_of_vol()) * r1 + rho1 * r2 * z2;
			double next1 = exp(current_log_spot_next);
			return Pair(next1, next2);
		}
	}
}
