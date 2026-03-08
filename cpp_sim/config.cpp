#include "config.hpp"
#include <stdexcept>

namespace galaxy {

SimulationMode parse_mode(const std::string& s) {
  if (s == "galaxy") return SimulationMode::galaxy;
  if (s == "two_body_orbit") return SimulationMode::two_body_orbit;
  if (s == "symmetric_pair") return SimulationMode::symmetric_pair;
  if (s == "small_n_conservation") return SimulationMode::small_n_conservation;
  if (s == "timestep_convergence") return SimulationMode::timestep_convergence;
  throw std::runtime_error("Unknown simulation_mode: " + s);
}

}  // namespace galaxy
