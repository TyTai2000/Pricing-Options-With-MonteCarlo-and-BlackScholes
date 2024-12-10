#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#define fast_cin() ios_base::sync_with_stdio(false); cin.tie(NULL); cout.tie(NULL)


#define _USE_MATH_DEFINES
using namespace std;

// Std. Probability density function
double norm_pdf(const double& x) {
    return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);
}

double norm_cdf(const double& x) {
    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));

    if (x >= 0.0) {
        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
    } else {
        return 1.0 - norm_cdf(-x);
    }
}

double d_j(const int& j, const double& S, const double& K, const double& r, const double& sigma, const double& T) {
    return (log(S/K) + (r + (pow(-1,j-1))*0.5*sigma*sigma)*T)/(sigma*(pow(T,0.5)));
}

double call_price(const double& S, const double& K, const double& r, const double& sigma, const double& T) {
    return S * norm_cdf(d_j(1, S, K, r, sigma, T))-K*exp(-r*T) * norm_cdf(d_j(2, S, K, r, sigma, T));
}

double put_price(const double& S, const double& K, const double& r, const double& sigma, const double& T) {
    return -S*norm_cdf(-d_j(1, S, K, r, sigma, T))+K*exp(-r*T) * norm_cdf(-d_j(2, S, K, r, sigma, T));
}

int main() {
    fast_cin();

    // Option parameters
    double S = 100.0;  // Underlying price
    double K = 105.0;  // Strike price
    double r = 0.05;   // Risk-free rate
    double sigma = 0.2; // Implied Volatility
    double T = 1.0;    // Time till Expiry (1 year)

    // Calculate call and put prices
    double call = call_price(S, K, r, sigma, T);
    double put = put_price(S, K, r, sigma, T);

    // Output results
    cout << "Option Parameters:\n";
    cout << "Strike Price:     " << K << endl;
    cout << "Risk-Free Rate:   " << r << endl;
    cout << "Volatility:       " << sigma << endl;
    cout << "Maturity (Years): " << T << endl;

    cout << "\nCalculated Prices:\n";
    cout << "Call Price:       " << call << endl;
    cout << "Put Price:        " << put << endl;

    return 0;
}
