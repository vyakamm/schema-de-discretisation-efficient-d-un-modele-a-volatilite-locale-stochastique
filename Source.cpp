
#include "Bins.h"
#include"Pricer.h"
#include"NewtonRaphson.h"
#include"Blackscholesformulas.h"
#include "HestonSolver.h"
#include<iostream>
#include <functional>
#include "ImpliedVolatilitySurface.h"
#include<vector>
using namespace std;

double objective_function(double tau, double kappa,  double theta, double vol_of_vol, double rho,  double strike,double init_spot,double init_variance,double risk_free_rate)
{

    vector<double> s(9);
    vector<double> t(8);
    vector< vector<double> >sigma(8, vector<double>(9));
    s[0] = 20; s[1] = 40; s[2] = 60; s[3] = 80; s[4] = 100; s[5] = 120; s[6] = 140; s[7] = 160; s[8] = 180;
    t[0] = 0.25; t[1] = 0.5; t[2] = 0.75; t[3] = 1; t[4] = 2; t[5] = 3; t[6] = 4; t[7] = 5;
    sigma[0][0] = 0.39; sigma[0][1] = 0.31; sigma[0][2] = 0.24; sigma[0][3] = 0.22; sigma[0][4] = 0.16; sigma[0][5] = 0.19; sigma[0][6] = 0.23; sigma[0][7] = 0.29; sigma[0][8] = 0.38;
    sigma[1][0] = 0.44; sigma[1][1] = 0.36; sigma[1][2] = 0.27; sigma[1][3] = 0.21; sigma[1][4] = 0.17; sigma[1][5] = 0.21; sigma[1][6] = 0.27; sigma[1][7] = 0.35; sigma[1][8] = 0.4;
    sigma[2][0] = 0.45; sigma[2][1] = 0.3; sigma[2][2] = 0.25; sigma[2][3] = 0.21; sigma[2][4] = 0.18; sigma[2][5] = 0.22; sigma[2][6] = 0.29; sigma[2][7] = 0.37; sigma[2][8] = 0.45;
    sigma[3][0] = 0.48; sigma[3][1] = 0.42; sigma[3][2] = 0.34; sigma[3][3] = 0.28; sigma[3][4] = 0.2; sigma[3][5] = 0.26; sigma[3][6] = 0.31; sigma[3][7] = 0.42; sigma[3][8] = 0.5;
    sigma[4][0] = 0.52; sigma[4][1] = 0.43; sigma[4][2] = 0.34; sigma[4][3] = 0.26; sigma[4][4] = 0.21; sigma[4][5] = 0.27; sigma[4][6] = 0.38; sigma[4][7] = 0.45; sigma[4][8] = 0.55;
    sigma[5][0] = 0.54; sigma[5][1] = 0.46; sigma[5][2] = 0.34; sigma[5][3] = 0.27; sigma[5][4] = 0.23; sigma[5][5] = 0.28; sigma[5][6] = 0.36; sigma[5][7] = 0.49; sigma[5][8] = 0.58;
    sigma[6][0] = 0.57; sigma[6][1] = 0.50; sigma[6][2] = 0.46; sigma[6][3] = 0.35; sigma[6][4] = 0.25; sigma[6][5] = 0.32; sigma[6][6] = 0.45; sigma[6][7] = 0.54; sigma[6][8] = 0.6;
    sigma[7][0] = 0.60; sigma[7][1] = 0.52; sigma[7][2] = 0.41; sigma[7][3] = 0.31; sigma[7][4] = 0.26; sigma[7][5] = 0.34; sigma[7][6] = 0.4; sigma[7][7] = 0.55; sigma[7][8] = 0.62;
    ImpliedVolatilitySurface surface(t,s, sigma,risk_free_rate);
    int n = 40;
    double sum = 0;
    for(int i=0;i<n;i++)
    {
        BSCall theCall(risk_free_rate, 0, tau, init_spot, strike);
        sum += pow((surface.implied_volatility(tau, strike) - NewtonRaphson<BSCall, &BSCall::Price, &BSCall::Vega>(HestonSolver(tau, kappa, theta, vol_of_vol, rho, strike, init_spot, init_variance, risk_free_rate), 0.5, 0.1, theCall))/ surface.implied_volatility(tau, strike), 2);
            strike += 5;
    }
    return sum;
}

int main()
{
  
       /*BSCall theCall(0, 0, 5, 100, 100);
        double vol = NewtonRaphson<BSCall, &BSCall::Price, &BSCall::Vega>(21.71805807260404, 0.5, 0.1, theCall);
        double PriceTwo = BlackScholesCall(100, 100, 0, 0, vol,5);
        cout << vol;*/
    //graine aléatoire
    // srand(clock());*
     srand(1);
     Pair init_spot_variance(100, 0.0945);
     equal_number b(1000000,20);
     vector<double> s(9);
     vector<double> t(8);
     vector< vector<double> >sigma(8, vector<double>(9));
     s[0] = 20; s[1] = 40; s[2] = 60; s[3] = 80; s[4] = 100; s[5] = 120; s[6] = 140; s[7] = 160; s[8] = 180;
     t[0] = 0.25; t[1] = 0.5; t[2] = 0.75; t[3] = 1; t[4] = 2; t[5] = 3; t[6] = 4; t[7] = 5;
     sigma[0][0] = 0.39; sigma[0][1] = 0.31; sigma[0][2] = 0.24; sigma[0][3] = 0.22; sigma[0][4] = 0.16; sigma[0][5] = 0.19; sigma[0][6] = 0.23; sigma[0][7] = 0.29; sigma[0][8] = 0.38;
     sigma[1][0] = 0.44; sigma[1][1] = 0.36; sigma[1][2] = 0.27; sigma[1][3] = 0.21; sigma[1][4] = 0.17; sigma[1][5] = 0.21; sigma[1][6] = 0.27; sigma[1][7] = 0.35; sigma[1][8] = 0.4;
     sigma[2][0] = 0.45; sigma[2][1] = 0.3; sigma[2][2] = 0.25; sigma[2][3] = 0.21; sigma[2][4] = 0.18; sigma[2][5] = 0.22; sigma[2][6] = 0.29; sigma[2][7] = 0.37; sigma[2][8] = 0.45;
     sigma[3][0] = 0.48; sigma[3][1] = 0.42; sigma[3][2] = 0.34; sigma[3][3] = 0.28; sigma[3][4] = 0.2; sigma[3][5] = 0.26; sigma[3][6] = 0.31; sigma[3][7] = 0.42; sigma[3][8] = 0.5;
     sigma[4][0] = 0.52; sigma[4][1] = 0.43; sigma[4][2] = 0.34; sigma[4][3] = 0.26; sigma[4][4] = 0.21; sigma[4][5] = 0.27; sigma[4][6] = 0.38; sigma[4][7] = 0.45; sigma[4][8] = 0.55;
     sigma[5][0] = 0.54; sigma[5][1] = 0.46; sigma[5][2] = 0.34; sigma[5][3] = 0.27; sigma[5][4] = 0.23; sigma[5][5] = 0.28; sigma[5][6] = 0.36; sigma[5][7] = 0.49; sigma[5][8] = 0.58;
     sigma[6][0] = 0.57; sigma[6][1] = 0.50; sigma[6][2] = 0.46; sigma[6][3] = 0.35; sigma[6][4] = 0.25; sigma[6][5] = 0.32; sigma[6][6] = 0.45; sigma[6][7] = 0.54; sigma[6][8] = 0.6;
     sigma[7][0] = 0.60; sigma[7][1] = 0.52; sigma[7][2] = 0.41; sigma[7][3] = 0.31; sigma[7][4] = 0.26; sigma[7][5] = 0.34; sigma[7][6] = 0.4; sigma[7][7] = 0.55; sigma[7][8] = 0.62;
     ImpliedVolatilitySurface surface(t, s, sigma,0);
     DupireLocalVolatilitySurface dlvs(surface, 0.0001, 0.0001, 100);
     Heston_local_sto_vol_Model model(0,-0.315,init_spot_variance,1.05,0.95, 0.0855,dlvs,b);
     vector<double> v{0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
     PathSimulatorQE p(model, v);
     vector<Pair>ssv(1000000, init_spot_variance);
     cout << p.next_step(1, init_spot_variance, ssv).second;
     
}