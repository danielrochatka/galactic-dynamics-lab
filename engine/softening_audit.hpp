#ifndef GALAXY_SOFTENING_AUDIT_HPP
#define GALAXY_SOFTENING_AUDIT_HPP
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
namespace galaxy {
struct SofteningAuditStats {
  std::uint64_t pair_count = 0, bh_pair_count = 0, star_pair_count = 0, violation_count = 0;
  double max_ratio = 0.0, max_ratio_r = 0.0, min_r = 0.0, median_r = 0.0, eps_used = 0.0, xi_metric_max_ratio = 0.0;
  std::vector<double> distances;
};
inline void softening_audit_pair(SofteningAuditStats& s, double r, double eps, bool bh, double tol = 1e-12) {
  // Scalar-factor invariant audit only: checks Plummer softened-vs-unsoftened 1/r^3 magnitude factor.
  // This does not by itself prove full acceleration vector sign correctness.
  if (r <= 0.0) return;
  ++s.pair_count; if (bh) ++s.bh_pair_count; else ++s.star_pair_count;
  s.distances.push_back(r);
  const double sf = 1.0 / std::pow(r * r + eps * eps, 1.5);
  const double uf = 1.0 / (r * r * r);
  const double ratio = (uf > 0.0) ? (sf / uf) : 1.0;
  if (ratio > s.max_ratio) { s.max_ratio = ratio; s.max_ratio_r = r; }
  if (eps > 0.0 && sf > uf * (1.0 + tol)) ++s.violation_count;
}
inline void softening_audit_finalize(SofteningAuditStats& s) {
  if (s.distances.empty()) return;
  std::sort(s.distances.begin(), s.distances.end());
  s.min_r = s.distances.front();
  s.median_r = s.distances[s.distances.size() / 2];
}
}
#endif
