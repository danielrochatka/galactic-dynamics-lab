#include "package_discovery.hpp"

#include <cctype>
#include <dirent.h>
#include <fstream>
#include <sys/stat.h>

namespace galaxy {
namespace {

std::string trim(const std::string& s) {
  size_t a = 0;
  while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  size_t b = s.size();
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return s.substr(a, b - a);
}

bool is_dir(const std::string& p) {
  struct stat st {};
  return stat(p.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

bool exists(const std::string& p) {
  struct stat st {};
  return stat(p.c_str(), &st) == 0;
}

std::unordered_map<std::string, std::string> parse_kv_file(const std::string& p) {
  std::unordered_map<std::string, std::string> out;
  std::ifstream f(p.c_str());
  std::string line;
  while (std::getline(f, line)) {
    auto hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    line = trim(line);
    if (line.empty()) continue;
    auto eq = line.find('=');
    if (eq == std::string::npos) continue;
    out[trim(line.substr(0, eq))] = trim(line.substr(eq + 1));
  }
  return out;
}

std::unordered_set<std::string> parse_token_file(const std::string& p) {
  std::unordered_set<std::string> out;
  std::ifstream f(p.c_str());
  std::string line;
  while (std::getline(f, line)) {
    auto hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    line = trim(line);
    if (line.empty()) continue;
    auto eq = line.find('=');
    if (eq != std::string::npos) line = trim(line.substr(0, eq));
    out.insert(line);
  }
  return out;
}

}

std::vector<PackageMetadata> discover_packages() {
  std::vector<PackageMetadata> out;
  const std::string physics_dir = "physics";
  DIR* dir = opendir(physics_dir.c_str());
  if (!dir) return out;
  for (dirent* e = readdir(dir); e != nullptr; e = readdir(dir)) {
    std::string name = e->d_name;
    if (name == "." || name == "..") continue;
    std::string pkg_dir = physics_dir + "/" + name;
    if (!is_dir(pkg_dir)) continue;
    PackageMetadata pm;
    pm.folder_name = name;
    pm.name = name;
    std::string pkg_cfg = pkg_dir + "/package.cfg";
    if (exists(pkg_cfg)) {
      auto kv = parse_kv_file(pkg_cfg);
      if (kv.count("name")) pm.name = kv["name"];
      if (kv.count("display_name")) pm.display_name = kv["display_name"];
      if (kv.count("version")) pm.version = kv["version"];
    }
    std::string defs = pkg_dir + "/defaults.cfg";
    if (exists(defs)) pm.defaults = parse_kv_file(defs);
    std::string schema = pkg_dir + "/config_schema.cfg";
    if (exists(schema)) pm.schema_keys = parse_token_file(schema);
    std::string modes = pkg_dir + "/modes.cfg";
    if (exists(modes)) pm.modes = parse_token_file(modes);
    out.push_back(pm);
  }
  closedir(dir);
  return out;
}

const PackageMetadata* find_package_metadata(const std::string& package_name) {
  static std::vector<PackageMetadata> cached = discover_packages();
  for (const auto& p : cached) {
    if (p.name == package_name || p.folder_name == package_name) return &p;
  }
  return nullptr;
}

} // namespace galaxy
