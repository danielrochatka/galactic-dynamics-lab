#include "git_provenance.hpp"

#include <array>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

namespace galaxy {

namespace {

std::string trim_ws(const std::string& s) {
  size_t a = 0;
  while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  size_t b = s.size();
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return s.substr(a, b - a);
}

std::string sh_single_quote(const std::string& path) {
  std::string r = "'";
  for (char c : path) {
    if (c == '\'')
      r += "'\\''";
    else
      r += c;
  }
  r += '\'';
  return r;
}

std::string read_popen_line(const std::string& cmd) {
  std::array<char, 4096> buf{};
  FILE* p = ::popen(cmd.c_str(), "r");
  if (!p) return "";
  std::string out;
  while (std::fgets(buf.data(), static_cast<int>(buf.size()), p)) out += buf.data();
  int st = ::pclose(p);
  (void)st;
  return trim_ws(out);
}

std::string git_in_repo(const std::string& repo_root, const char* subcmd) {
  std::string cmd = "git -C " + sh_single_quote(repo_root) + " " + subcmd + " 2>/dev/null";
  return read_popen_line(cmd);
}

bool git_repo_has_dirty_tree(const std::string& repo_root) {
  std::string por = git_in_repo(repo_root, "status --porcelain");
  return !por.empty();
}

std::string parent_dir(const std::string& p) {
  if (p.empty()) return "";
  size_t i = p.size();
  while (i > 0 && p[i - 1] == '/') --i;
  if (i == 0) return "";
  size_t j = p.rfind('/', i - 1);
  if (j == std::string::npos) {
    if (p[0] == '/' && i == 1) return "/";
    return "";
  }
  if (j == 0) return "/";
  return p.substr(0, j);
}

bool path_has_dot_git(const std::string& dir) {
  if (dir.empty()) return false;
  std::string d = dir;
  while (!d.empty() && d.back() == '/') d.pop_back();
  std::string probe = d + "/.git";
  return ::access(probe.c_str(), F_OK) == 0;
}

std::string toplevel_from_dir(const std::string& start_dir) {
  if (start_dir.empty()) return "";
  std::string out = git_in_repo(start_dir, "rev-parse --show-toplevel");
  return out;
}

std::string find_git_repo_root() {
  const char* env = std::getenv("GALAXY_REPO_ROOT");
  if (env && env[0]) {
    std::string e(env);
    std::string t = toplevel_from_dir(e);
    if (!t.empty()) return t;
    if (path_has_dot_git(e)) return e;
  }

  char cwd_buf[4096];
  if (::getcwd(cwd_buf, sizeof(cwd_buf))) {
    std::string t = toplevel_from_dir(cwd_buf);
    if (!t.empty()) return t;
  }

  char exe_buf[4096];
  ssize_t n = ::readlink("/proc/self/exe", exe_buf, sizeof(exe_buf) - 1);
  if (n > 0) {
    exe_buf[n] = '\0';
    std::string d(exe_buf);
    for (int depth = 0; depth < 24 && !d.empty(); ++depth) {
      std::string t = toplevel_from_dir(d);
      if (!t.empty()) return t;
      if (path_has_dot_git(d)) return d;
      d = parent_dir(d);
    }
  }

  return "";
}

std::string build_code_version_label(const std::string& short_sha,
                                     const std::string& branch,
                                     const std::string& tag,
                                     bool dirty) {
  if (short_sha.empty() || short_sha == "unknown") return "unknown";
  std::string base;
  if (!tag.empty())
    base = tag;
  else if (!branch.empty())
    base = branch;
  else
    base = "HEAD";
  std::string label = base + "@" + short_sha;
  if (dirty) label += "-dirty";
  return label;
}

}  // namespace

GitProvenance resolve_git_provenance() {
  static GitProvenance cached;
  static bool done = false;
  if (done) return cached;

  GitProvenance p;
  p.git_commit_full = "unknown";
  p.git_commit_short = "unknown";
  p.git_branch = "";
  p.git_tag = "";
  p.git_dirty = false;
  p.code_version_label = "unknown";

  const std::string root = find_git_repo_root();
  if (root.empty()) {
    cached = p;
    done = true;
    return cached;
  }

  const std::string full = git_in_repo(root, "rev-parse HEAD");
  if (full.empty()) {
    cached = p;
    done = true;
    return cached;
  }

  p.git_commit_full = full;
  p.git_commit_short = git_in_repo(root, "rev-parse --short HEAD");
  if (p.git_commit_short.empty() && full.size() >= 7)
    p.git_commit_short = full.substr(0, 7);
  else if (p.git_commit_short.empty())
    p.git_commit_short = full;

  p.git_branch = git_in_repo(root, "branch --show-current");
  p.git_tag = git_in_repo(root, "describe --tags --exact-match");
  p.git_dirty = git_repo_has_dirty_tree(root);
  p.code_version_label = build_code_version_label(p.git_commit_short, p.git_branch, p.git_tag, p.git_dirty);
  cached = p;
  done = true;
  return cached;
}

}  // namespace galaxy
