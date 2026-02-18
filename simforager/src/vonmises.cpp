#include <cmath>
#include <utility>
#include <boost/random.hpp>
#include <memory>

#include "reef.hpp"

// Mathematical constant pi
static constexpr double kPi = 3.141592653589793238462643383279502884;

// Wrap angle to [-pi, pi)
static inline double wrap_angle_pi(double a) {
  a = std::fmod(a + kPi, 2.0 * kPi);
  if (a < 0.0) a += 2.0 * kPi;
  return a - kPi;
}

// Sample from Von Mises(mu, kappa) with rejection sampling.
// Returns an angle wrapped to [-pi, pi).
double sample_vonmises(double mu, double kappa, boost::random::mt19937 &gen) {
  // Treat non-positive kappa as uniform noise.
  if (!(kappa > 0.0)) {
    boost::random::uniform_real_distribution<> unif(-kPi, kPi);
    return unif(gen);
  }

  const double tau = 1.0 + std::sqrt(1.0 + 4.0 * kappa * kappa);
  const double rho = (tau - std::sqrt(2.0 * tau)) / (2.0 * kappa);
  const double r   = (1.0 + rho * rho) / (2.0 * rho);

  boost::random::uniform_real_distribution<> unif01(0.0, 1.0);

  for (;;) {
    const double u1 = unif01(gen);
    const double z  = std::cos(kPi * u1);
    const double f  = (1.0 + r * z) / (r + z);
    const double c  = kappa * (r - f);

    const double u2 = unif01(gen);

    if (u2 < c * (2.0 - c) || u2 <= c * std::exp(1.0 - c)) {
      const double sign = (unif01(gen) > 0.5) ? 1.0 : -1.0;

      // Guard against tiny numerical drift outside [-1, 1]
      double f_clamped = f;
      if (f_clamped > 1.0) f_clamped = 1.0;
      if (f_clamped < -1.0) f_clamped = -1.0;

      const double theta = mu + sign * std::acos(f_clamped);
      return wrap_angle_pi(theta);
    }
  }
}

// Convert polar displacement (radius, angle) to Cartesian dx,dy.
// grid_size is unused here but retained for signature compatibility.
std::pair<double, double> polar_to_cartesian(double radius, double angle,
                                             std::shared_ptr<GridCoords> /*grid_size*/) {
  const double x = radius * std::cos(angle);
  const double y = radius * std::sin(angle);
  return {x, y};
}
