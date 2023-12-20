#include"bins.h"
#include<vector>
#include <iostream>
#include <algorithm>
#include <numeric>


// abstract class bins's constructor with arguments
bins::bins(const size_t& number_of_simulations, const size_t& number_of_bins) :_number_of_simulations(number_of_simulations), _number_of_bins(number_of_bins)
{
}

double bins::expectation(const vector< pair<double, double> >& _spot_function_variance, const double& spot) const
{
    // Vector of Indexes
    std::vector<size_t> indices(_spot_function_variance.size());
    std::iota(indices.begin(), indices.end(), 0);

    // Sorting using Indexes
    std::sort(indices.begin(), indices.end(),
        [&_spot_function_variance](size_t a, size_t b) {
            return _spot_function_variance[a].first <= _spot_function_variance[b].first;
        });

    // Creating a vector containing the sorted values of the stock in order to compute bins bounds
    std::vector<double> sorted_spot;
    for (size_t i : indices) {
        sorted_spot.push_back(_spot_function_variance[i].first);
    }

    // Computing bins bounds
    sorted_spot = bins_creation_manner(sorted_spot);

    double sum = 0;
    double number_of_samples = 0;

    // Find Bounds
    auto l = std::lower_bound(sorted_spot.begin(), sorted_spot.end(), spot);
    size_t index = std::distance(sorted_spot.begin(), l);
    double bound_2 = sorted_spot[index];
    double bound_1 = sorted_spot[(index > 0) ? (index - 1) : 0];

    // Compute the expectation
    for (size_t i = 0; i < _spot_function_variance.size(); ++i)
    {
        if (_spot_function_variance[i].first > bound_1 && _spot_function_variance[i].first <= bound_2)
        {
            sum += _spot_function_variance[i].second;
            ++number_of_samples;
        }
    }

    return ((1.0 / (_number_of_simulations * (number_of_samples / _number_of_simulations))) * sum);
}

// Getter method
size_t bins::number_of_simulations()
{
	return _number_of_simulations;
}

equidistant::equidistant(const int& number_of_simulations, const int& number_of_bins) :bins(number_of_simulations, number_of_bins)
{
}

// creation of bins with respect to an equidistant grid
vector<double> equidistant::bins_creation_manner(const vector<double>& spot) const
{
	vector<double> _bins;
	//cout << spot_variance[_number_of_simulations - 1].first<<"\n";
	_bins.push_back(0);
	for (int i = 1; i <= _number_of_bins; i++)
	{
		_bins.push_back(spot[0] + ((static_cast<double>(i) / _number_of_bins) * (spot[_number_of_simulations - 1] - spot[0])));
	}
	return _bins;
}

equidistant* equidistant::clone() const
{
	return new equidistant(*this);
}


equal_number::equal_number(const int& number_of_simulations, const int& number_of_bins) :bins(number_of_simulations, number_of_bins)
{
}

//creation of bins with the same number of samples
vector<double> equal_number::bins_creation_manner(const vector<double>& spot) const
{
	vector<double> bin;
	bin.push_back(0);
	for (size_t i = 1; i <= _number_of_bins; i++)
		bin.push_back(spot[i * (_number_of_simulations / _number_of_bins) - 1]);
	return bin;
}

equal_number* equal_number::clone() const
{
	return new equal_number(*this);
}
