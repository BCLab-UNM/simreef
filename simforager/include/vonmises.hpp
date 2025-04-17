
// vonmises.hpp — Von Mises Sampling Interface
// -------------------------------------------
// Declares functions to support fish directional movement in SimForager,
// using the Von Mises distribution for turning angle sampling.

#ifndef VONMISES_HPP
#define VONMISES_HPP

#include <utility>
#include <memory>
#include "reef.hpp"

// Returns a new turning angle sampled from the Von Mises distribution.
// mu: mean angle, kappa: concentration (higher = less deviation)
double sample_vonmises(double mu, double kappa, boost::random::mt19937 &gen);

// Converts polar (radius, angle) to Cartesian (dx, dy)
std::pair<double, double> polar_to_cartesian(double radius, double angle, std::shared_ptr<GridCoords> grid_size);

#endif
