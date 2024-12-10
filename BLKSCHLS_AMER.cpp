#include <iostream>
#include <cmath>

#define _USE_MATH_DEFINES

using namespace std;

// Standard Cumulative Distribution Function (CDF)
double norm_cdf(const double& x) {
    double k = 1.0 / (1.0 + 0.2316419 * x);
    double k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));

    if (x >= 0.0) {
        return (1.0 - (1.0 / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x) * k_sum);
    } else {
        return 1.0 - norm_cdf(-x);
    }
}

// Calculating d1 and d2 terms
double d_j(const int& j, const double& S, const double& K, const double& r, const double& sigma, const double& T) {
    return (log(S / K) + (r + (pow(-1, j - 1)) * 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
}

// Black-Scholes Put Price (treated as if it were for an American option)
double put_price(const double& S, const double& K, const double& r, const double& sigma, const double& T) {
    return -S * norm_cdf(-d_j(1, S, K, r, sigma, T)) + K * exp(-r * T) * norm_cdf(-d_j(2, S, K, r, sigma, T));
}

int main() {
    // Parameters for the option
    double S = 100.0;     // Underlying price
    double K = 105.0;     // Strike price
    double r = 0.05;      // Risk-free rate
    double sigma = 0.2;   // Volatility
    double T = 1.0;       // Time to expiration in years

    // Calculate Black-Scholes price for the American put
    double american_put_price = put_price(S, K, r, sigma, T);

    // Output the result
    cout << "Attempting to price an American Put Option with Black-Scholes:\n";
    cout << "Calculated Price (using Black-Scholes): " << american_put_price << endl;
    cout << "(Note: This price may underestimate the true value due to lack of early exercise consideration.)" << endl;

    return 0;
}
