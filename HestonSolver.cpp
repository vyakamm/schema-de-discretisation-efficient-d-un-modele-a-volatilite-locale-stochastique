#pragma once
#include <iostream>
#include<thread>
#include <ppl.h>
#include<cmath>
#include "HestonSolver.h"
#include"Payoff.h"
#include"Gaussienne.h"
#include "PathSimulator.h"
#include "NewtonRaphson.h"
#include "Modeles.h"
using namespace std;

double HestonSolver(const double& tau, const double& kappa, const double& theta, const double& vol_of_vol, const double& rho, const double& strike, const double & init_spot ,const double& init_variance, const double& risk_free_rate)
{
    int N = 10000, N_max = 1000, i;
    double aa = theta * kappa * tau / pow(vol_of_vol, 2);
    double bb = -2 * theta * kappa / pow(vol_of_vol, 2);
    double du = static_cast<double>(N_max) / N;
    complex<double> I(0, 1), P(0, 0), u1(0, 0), kap(kappa, 0), lamb2(vol_of_vol * vol_of_vol, 0), un(1, 0), t(-tau, 0), con(log(init_spot / strike) + risk_free_rate * tau, 0),
        v0(init_variance, 0), a(theta * kappa * tau / pow(vol_of_vol, 2), 0), d(du, 0), deux(2, 0), con1(init_spot / strike - exp(-risk_free_rate * tau)), pi(3.14, 0), rho_vol(rho * vol_of_vol, 0);
    
    for (i = 1; i < N; ++i)
    {
        double u2 = i * du;
        complex<double> u1 = complex<float>(u2, -1);
        complex<double> a1 = rho_vol * u1 * I;
        complex<double> a2 = rho_vol * u2 * I;
        complex<double> d1 = sqrt((a1 - kap) * (a1 - kap) + lamb2 * (u1 * I + u1 * u1));
        complex<double> d2 = sqrt((a2 - kap) * (a2 - kap) + lamb2 * (u2 * I + u2 * u2));
        complex<double> g1 = (kap - a1 - d1) / (kap - a1 + d1);
        complex<double> g2 = (kap - a2 - d2) / (kap - a2 + d2);
        complex<double> b1 = exp(u1 * I * (con)) * pow(((un - g1 * exp(d1 * t)) / (un - g1)), bb);
        complex<double> b2 = exp(u2 * I * (con)) * pow(((un - g2 * exp(d2 * t)) / (un - g2)), bb);
        complex<double> phi1 = b1 * exp(a * (kap - a1 - d1) + v0 * (kap - a1 - d1) * (un - exp(d1 * t)) / (un - g1 * exp(d1 * t)) / lamb2);
        complex<double> phi2 = b2 * exp(a * (kap - a2 - d2) + v0 * (kap - a2 - d2) * (un - exp(d2 * t)) / (un - g2 * exp(d2 * t)) / lamb2);
        P += ((phi1 - phi2) / (u2 * I)) * d;
    }
    return strike * real((con1) / deux + P / pi);
}

double MonteCarlo(const PayOff& thePayOff,double Expiry,PathSimulator& path,double r,unsigned long NumberOfPaths)
{
    // on tire des seeds indépendantes pour chaque coeur lors de la parallélisation
    srand(time(NULL) * std::hash<std::thread::id>{}(std::this_thread::get_id()));
    double thisSpot;
    double runningSum = 0;
    for (unsigned long i = 0; i < NumberOfPaths; i++)
    {
        thisSpot = path.path_maturity().first;
        double thisPayOff = thePayOff(thisSpot);
        runningSum += thisPayOff;
    }
    double mean = runningSum / NumberOfPaths;
   
    return  mean;
}
