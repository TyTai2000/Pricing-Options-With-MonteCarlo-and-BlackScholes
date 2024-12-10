#include <iostream>
#include <cmath>
#include <random>
#include <vector>

#define _USE_MATH_DEFINES

using namespace std;

// Function to generate random numbers from a standard normal distribution
double generate_random_normal() {
    static random_device rd;
    static mt19937 generator(rd());
    static normal_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}

// Monte Carlo simulation for European option pricing
double monte_carlo_option_price(const double& S, const double& K, const double& r, const double& sigma, const double& T, const int& num_simulations, bool is_call) {
    double payoff_sum = 0.0;

    for (int i = 0; i < num_simulations; ++i) {
        // Simulate the end price using GBM
        double Z = generate_random_normal();  // Standard normal random variable
        double S_T = S * exp((r - 0.5 * sigma * sigma) * T + sigma * sqrt(T) * Z);

        // Calculate payoff
        double payoff = is_call ? max(S_T - K, 0.0) : max(K - S_T, 0.0);
        payoff_sum += payoff;
    }

    // Discount the average payoff to present value
    return exp(-r * T) * (payoff_sum / num_simulations);
}

int main() {
    // Option parameters
    double S = 100.0;       // Initial stock price
    double K = 105.0;       // Strike price
    double r = 0.05;        // Risk-free rate
    double sigma = 0.2;     // Volatility
    double T = 1.0;         // Time to maturity (in years)
    int num_simulations = 100000;  // Number of simulations

    // Calculate option prices
    double call_price = monte_carlo_option_price(S, K, r, sigma, T, num_simulations, true);
    double put_price = monte_carlo_option_price(S, K, r, sigma, T, num_simulations, false);

    // Output results
    cout << "Monte Carlo Simulation Results (European Options):\n";
    cout << "Call Option Price: " << call_price << endl;
    cout << "Put Option Price: " << put_price << endl;

    return 0;
}
