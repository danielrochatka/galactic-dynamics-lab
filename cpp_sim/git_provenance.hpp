#ifndef GALAXY_GIT_PROVENANCE_HPP
#define GALAXY_GIT_PROVENANCE_HPP

#include <string>

namespace galaxy {

/** Git metadata captured at run time (best-effort; safe fallbacks if git is unavailable). */
struct GitProvenance {
  std::string git_commit_full;
  std::string git_commit_short;
  std::string git_branch;
  std::string git_tag;
  bool git_dirty = false;
  std::string code_version_label;
};

/** Resolve once per process; never throws. Uses git when possible; otherwise "unknown" fields. */
GitProvenance resolve_git_provenance();

}  // namespace galaxy

#endif
