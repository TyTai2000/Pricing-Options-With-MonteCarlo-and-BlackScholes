#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <numeric>

#define _USE_MATH_DEFINES

using namespace std;

// Generate standard normal random numbers
double generate_random_normal() {
    static random_device rd;
    static mt19937 generator(rd());
    static normal_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}

// Simulate asset price paths
vector<vector<double>> simulate_paths(double S, double r, double sigma, double T, int num_paths, int num_steps) {
    vector<vector<double>> paths(num_paths, vector<double>(num_steps + 1, S));
    double dt = T / num_steps;

    for (int i = 0; i < num_paths; ++i) {
        for (int j = 1; j <= num_steps; ++j) {
            double Z = generate_random_normal();
            paths[i][j] = paths[i][j - 1] * exp((r - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * Z);
        }
    }
    return paths;
}

// Polynomial regression for quadratic basis
vector<double> polynomial_regression(const vector<double>& X, const vector<double>& Y) {
    int n = X.size();
    vector<vector<double>> A(3, vector<double>(3, 0.0)); // 3x3 matrix for quadratic regression
    vector<double> b(3, 0.0); // Right-hand side

    // Fill matrix A and vector b
    for (int i = 0; i < n; ++i) {
        A[0][0] += 1;
        A[0][1] += X[i];
        A[0][2] += X[i] * X[i];
        A[1][1] += X[i] * X[i];
        A[1][2] += X[i] * X[i] * X[i];
        A[2][2] += X[i] * X[i] * X[i] * X[i];
        b[0] += Y[i];
        b[1] += X[i] * Y[i];
        b[2] += X[i] * X[i] * Y[i];
    }

    A[1][0] = A[0][1];
    A[2][0] = A[0][2];
    A[2][1] = A[1][2];

    // Solve Ax = b using Gaussian elimination
    vector<double> beta(3, 0.0); // Coefficients

    // Forward elimination
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < 3; ++k) {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    // Back substitution
    for (int i = 2; i >= 0; --i) {
        beta[i] = b[i];
        for (int j = i + 1; j < 3; ++j) {
            beta[i] -= A[i][j] * beta[j];
        }
        beta[i] /= A[i][i];
    }

    return beta;
}

double longstaff_schwartz(const vector<vector<double>>& paths, double K, double r, double T, int num_steps, bool is_call) {
    int num_paths = paths.size();
    double dt = T / num_steps;
    vector<double> cashflows(num_paths, 0.0);

    // Initialize cashflows with option payoffs at maturity
    for (int i = 0; i < num_paths; ++i) {
        cashflows[i] = is_call ? max(paths[i][num_steps] - K, 0.0) : max(K - paths[i][num_steps], 0.0);
    }

    // Backward induction
    for (int step = num_steps - 1; step >= 0; --step) {
        vector<double> X, Y; // Stock prices and discounted cashflows for regression

        // Collect in-the-money paths
        for (int i = 0; i < num_paths; ++i) {
            double S = paths[i][step];
            double payoff = is_call ? max(S - K, 0.0) : max(K - S, 0.0);

            if (payoff > 0.0) { // Only consider in-the-money paths
                X.push_back(S);
                Y.push_back(cashflows[i] * exp(-r * dt)); // Discounted future cashflow
            }
        }

        // Perform regression if there are in-the-money paths
        if (!X.empty()) {
            vector<double> beta = polynomial_regression(X, Y);

            // Update cashflows based on exercise decision
            for (int i = 0; i < num_paths; ++i) {
                double S = paths[i][step];
                double payoff = is_call ? max(S - K, 0.0) : max(K - S, 0.0);
                double continuation_value = beta[0] + beta[1] * S + beta[2] * S * S;

                // Exercise if immediate payoff is better
                if (payoff > continuation_value) {
                    cashflows[i] = payoff;
                }
            }
        }
    }

    // Calculate the present value of cashflows
    double option_price = 0.0;
    for (double cashflow : cashflows) {
        option_price += cashflow * exp(-r * T);
    }
    return option_price / num_paths;
}


int main() {
    // Parameters
    double S = 100.0;  // Initial stock price
    double K = 105.0;  // Strike price
    double r = 0.05;   // Risk-free rate
    double sigma = 0.2; // Volatility
    double T = 1.0;    // Time to maturity (in years)
    int num_paths = 8000; // Number of simulation paths
    int num_steps = 5;     // Number of time steps

    // Simulate paths
    vector<vector<double>> paths = simulate_paths(S, r, sigma, T, num_paths, num_steps);

    // Calculate American option prices
    double call_price = longstaff_schwartz(paths, K, r, T, num_steps, true);
    double put_price = longstaff_schwartz(paths, K, r, T, num_steps, false);

    // Output results
    cout << "American Option Prices (Monte Carlo - Longstaff-Schwartz):\n";
    cout << "Call Option Price: " << call_price << endl;
    cout << "Put Option Price: " << put_price << endl;

    return 0;
}
