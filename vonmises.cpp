#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <boost/random.hpp>
#include <cstdlib> // For std::atoi, std::atof

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
std::pair<double, double> polar_to_cartesian(double radius, double angle) {
    double x = radius * std::cos(angle);
    double y = radius * std::sin(angle);
    return {x, y};
}

// Function to simulate a random walk with Von Mises-distributed turning angles
std::vector<std::pair<double, double>> random_walk_vonmises(int steps, double step_length, double mu, double kappa, int& total_steps_completed, int total_steps) {
    // Setup the random number generator
    boost::random::mt19937 gen;
    gen.seed(static_cast<unsigned int>(total_steps_completed));  // Seed

    // Initialize the walk path and starting position
    std::vector<std::pair<double, double>> path;
    path.emplace_back(0.0, 0.0); // Start at the origin

    // Starting direction
    double current_angle = 0.0;

    // Simulate the random walk
    for (int i = 0; i < steps; ++i) {
        // Draw turning angle from Von Mises distribution
        double turning_angle = sample_vonmises(mu, kappa, gen);
        current_angle += turning_angle;

        // Move forward in the direction given by current_angle
        auto [dx, dy] = polar_to_cartesian(step_length, current_angle);

        // Update the new position
        double new_x = path.back().first + dx;
        double new_y = path.back().second + dy;

        // Add the new position to the path
        path.emplace_back(new_x, new_y);

        // Update total steps completed and display the global progress
        total_steps_completed++;
       
    }

    double progress = (static_cast<double>(total_steps_completed) / total_steps) * 100;
    std::cout << "\rProgress: " << progress << "% complete." << std::flush;
    return path;
}

// Function to write multiple paths to a CSV file
void write_paths_to_csv(const std::vector<std::vector<std::pair<double, double>>> &paths, const std::string &filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return;
    }

    // Write header
    for (size_t i = 0; i < paths.size(); ++i) {
        file << "x" << i + 1 << ",y" << i + 1;
        if (i < paths.size() - 1) {
            file << ",";
        }
    }
    file << "\n";

    // Determine the maximum number of steps (columns)
    size_t max_steps = 0;
    for (const auto &path : paths) {
        if (path.size() > max_steps) {
            max_steps = path.size();
        }
    }

    // Write the data row by row
    for (size_t step = 0; step < max_steps; ++step) {
        for (size_t path_idx = 0; path_idx < paths.size(); ++path_idx) {
            if (step < paths[path_idx].size()) {
                file << paths[path_idx][step].first << "," << paths[path_idx][step].second;
            } else {
                file << ",";  // Leave blank if this path has fewer steps
            }

            if (path_idx < paths.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
    std::cout << "\nPaths written to " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    // Check if enough arguments are provided
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <path_length> <step_length> <number_of_walkers> <kappa_value> <output_file_path>" << std::endl;
        return 1;
    }

    // Parse the command-line arguments
    int steps = std::atoi(argv[1]);         // Path length (number of steps)
    double step_length = std::atof(argv[2]); // Step length
    int num_paths = std::atoi(argv[3]);      // Number of walkers
    double kappa = std::atof(argv[4]);       // Concentration parameter (kappa)
    std::string output_file = argv[5];       // Output file path

    // Mean direction (forward)
    double mu = 0.0;              

    std::cout << "Running with kappa " << kappa << std::endl;
    if (kappa == 0)
      {
	std::cout << "kappa=0 is a special case and reduces to a uniform random walk." << std::endl;
      }
    
    // Track total steps completed for the progress bar
    int total_steps_completed = 0;
    int total_steps = steps * num_paths;

    // Generate multiple paths
    std::vector<std::vector<std::pair<double, double>>> paths;
    for (int i = 0; i < num_paths; ++i) {
        paths.push_back(random_walk_vonmises(steps, step_length, mu, kappa, total_steps_completed, total_steps));
    }

    std::cout << "Generating paths complete. Writing to file..." << std::flush;
    
    // Write paths to CSV
    write_paths_to_csv(paths, output_file);

    std::cout << "done" << std::endl;

    return 0;
}
