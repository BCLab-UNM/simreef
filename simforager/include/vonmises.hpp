#ifndef VONMISES_H
#define VONMISES_H

#include <cmath>
#include <utility>
#include <boost/random.hpp>
#include <memory>

struct GridCoords;

// Function to sample from a Von Mises distribution
double sample_vonmises(double mu, double kappa, boost::random::mt19937 &gen);

// Function to convert polar coordinates to Cartesian
std::pair<double, double> polar_to_cartesian(double radius, double angle, std::shared_ptr<GridCoords> grid_size);

#endif // VONMISES_H