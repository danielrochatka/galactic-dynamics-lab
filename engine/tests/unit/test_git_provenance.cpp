#include "doctest.h"
#include "git_provenance.hpp"

TEST_CASE("resolve_git_provenance: never returns empty required strings") {
  galaxy::GitProvenance p = galaxy::resolve_git_provenance();
  CHECK(!p.git_commit_full.empty());
  CHECK(!p.git_commit_short.empty());
  CHECK(!p.code_version_label.empty());
}
