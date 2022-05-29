#include"bins.h"
#include<vector>
#include <iostream>
#include <algorithm>

// abstract class bins's constructor with arguments
bins::bins(const size_t &number_of_simulations, const size_t &number_of_bins) :_number_of_simulations(number_of_simulations), _number_of_bins(number_of_bins)
{
}

double bins::expectation(vector< pair<double,double> >& _spot_function_variance, const double &spot) const
{
	size_t i;
	vector<double> sorted_spot;
	// Sorting the spot samples
	std::sort(_spot_function_variance.begin(), _spot_function_variance.end());
	for (i = 0; i < _spot_function_variance.size(); i++)
		sorted_spot.push_back(_spot_function_variance[i].first); 
	sorted_spot = bins_creation_manner(sorted_spot);
	double sum = 0;
	double number_of_samples = 0;
	// find the bounds
	for (i = 0; i < sorted_spot.size() && spot > sorted_spot[i]; i++);
	double bound_2 = sorted_spot[i];
	double bound_1 = sorted_spot[i-1];
	for (size_t i = 0; i < _spot_function_variance.size(); i++)
	{
		if (_spot_function_variance[i].first <= bound_2 && _spot_function_variance[i].first > bound_1)
		{
			sum += _spot_function_variance[i].second;
			number_of_samples++;
		}
	}
	return ((1 / (_number_of_simulations * (number_of_samples / _number_of_simulations))) * sum);
	
}

// Getter method
size_t bins::number_of_simulations()
{
	return _number_of_simulations;
}

equidistant::equidistant(const int& number_of_simulations, const int& number_of_bins):bins(number_of_simulations,number_of_bins)
{
}

// creation of bins with respect to an equidistant grid
vector<double> equidistant::bins_creation_manner( const vector<double>& spot) const
{
	vector<double> _bins;
	//cout << spot_variance[_number_of_simulations - 1].first<<"\n";
	_bins.push_back(0);
	for (int i =1 ; i <= _number_of_bins; i++)
	{
		_bins.push_back(spot[0] + ((static_cast<double>(i) / _number_of_bins) * (spot[_number_of_simulations-1] - spot[0])));
	}
	return _bins;
}

equidistant* equidistant::clone() const
{
	return new equidistant(*this);
}


equal_number::equal_number(const int& number_of_simulations, const int& number_of_bins):bins(number_of_simulations,number_of_bins)
{
}

//creation of bins with the same number of samples
vector<double> equal_number::bins_creation_manner(const vector<double>& spot) const
{
	vector<double> bin;
	bin.push_back(0);
	for (size_t i = 1; i <= _number_of_bins; i++)
		bin.push_back(spot[i *(_number_of_simulations / _number_of_bins)  - 1]);
	return bin;
}

equal_number* equal_number::clone() const
{
	return new equal_number(*this);
}


