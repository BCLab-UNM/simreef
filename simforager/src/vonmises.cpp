#include <cmath>
#include <vector>
#include <boost/random.hpp>
#include <memory>
#include "reef.hpp"


// Function to sample from a Von Mises distribution
double sample_vonmises(double mu, double kappa, boost::random::mt19937 &gen) {
    if (kappa == 0) {
        // If kappa is 0, return a random angle uniformly distributed between 0 and 2 * PI
        boost::random::uniform_real_distribution<> uniform_dist(0.0, 2 * M_PI);
        return uniform_dist(gen);
    }

    // Otherwise, use the Von Mises sampling algorithm
    const double tau = 1.0 + std::sqrt(1.0 + 4.0 * kappa * kappa);
    const double rho = (tau - std::sqrt(2.0 * tau)) / (2.0 * kappa);
    const double r = (1.0 + rho * rho) / (2.0 * rho);

    boost::random::uniform_real_distribution<> uniform_dist(0.0, 1.0);

    while (true) {
        double u1 = uniform_dist(gen); 
        double z = std::cos(M_PI * u1);
        double f = (1.0 + r * z) / (r + z);
        double c = kappa * (r - f);

        double u2 = uniform_dist(gen);

        if (u2 < c * (2.0 - c) || u2 <= c * std::exp(1.0 - c)) {
            double theta = mu + (uniform_dist(gen) > 0.5 ? 1 : -1) * std::acos(f);
            return theta;
        }
    }
}


// Function to convert polar coordinates to Cartesian
std::pair<double, double> polar_to_cartesian(double radius, double angle, std::shared_ptr<GridCoords> grid_size) {
    double x = radius * std::cos(angle);
    double y = radius * std::sin(angle);
    x = std::max(0.0, std::min(x, static_cast<double>(grid_size->x)));
    y = std::max(0.0, std::min(y, static_cast<double>(grid_size->y)));
    return {x, y};
}