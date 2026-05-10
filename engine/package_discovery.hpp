#ifndef GALAXY_PACKAGE_DISCOVERY_HPP
#define GALAXY_PACKAGE_DISCOVERY_HPP

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace galaxy {

struct PackageMetadata {
  std::string folder_name;
  std::string name;
  std::string display_name;
  std::string version;
  std::unordered_map<std::string, std::string> defaults;
  std::unordered_set<std::string> schema_keys;
  std::unordered_set<std::string> modes;
};

std::vector<PackageMetadata> discover_packages();
const PackageMetadata* find_package_metadata(const std::string& package_name);

} // namespace galaxy

#endif
