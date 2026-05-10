#ifndef GALAXY_TPF_CORE_CONFIG_HPP
#define GALAXY_TPF_CORE_CONFIG_HPP

#include <string>

namespace galaxy {
struct Config;

void sync_tpfcore_legacy_fields_from_package_options(Config& config);

bool tpfcore_owns_mode_token(const std::string& mode_token);

} // namespace galaxy

#endif
