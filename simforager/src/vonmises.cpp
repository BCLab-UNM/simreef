
// vonmises.cpp — Implementation of Von Mises Sampling
// ---------------------------------------------------
// This file implements the sampling of angular values from the Von Mises
// distribution, a circular analogue of the normal distribution, used for modeling
// directional persistence in fish movement. Also includes polar-to-Cartesian
// coordinate conversion for grid navigation.

#include <cmath>
#include <vector>
#include <boost/random.hpp>
#include <memory>
#include "reef.hpp"

// Draws a sample from a Von Mises distribution with given mean angle `mu`
// and concentration parameter `kappa` using rejection sampling.
double sample_vonmises(double mu, double kappa, boost::random::mt19937 &gen) {
    if (kappa == 0) {
        // Uniform distribution on [0, 2π) if no directional preference
        boost::random::uniform_real_distribution<> uniform_dist(0.0, 2 * M_PI);
        return uniform_dist(gen);
    }

    // Constants used in the sampling method
    const double tau = 1.0 + std::sqrt(1.0 + 4.0 * kappa * kappa);
    const double rho = (tau - std::sqrt(2.0 * tau)) / (2.0 * kappa);
    const double r = (1.0 + rho * rho) / (2.0 * rho);

    boost::random::uniform_real_distribution<> uniform_dist(0.0, 1.0);

    // Rejection sampling loop
    while (true) {
        double u1 = uniform_dist(gen); 
        double z = std::cos(M_PI * u1);
        double f = (1.0 + r * z) / (r + z);
        double c = kappa * (r - f);

        double u2 = uniform_dist(gen);

        // Accept or reject the sample
        if (u2 < c * (2.0 - c) || u2 <= c * std::exp(1.0 - c)) {
            double theta = mu + (uniform_dist(gen) > 0.5 ? 1 : -1) * std::acos(f);
            return theta;
        }
    }
}

// Converts a polar coordinate (radius, angle) into Cartesian (x, y),
// useful for determining movement direction in the grid
std::pair<double, double> polar_to_cartesian(double radius, double angle, std::shared_ptr<GridCoords> grid_size) {
    double x = radius * std::cos(angle);
    double y = radius * std::sin(angle);
    return {x, y};
}
