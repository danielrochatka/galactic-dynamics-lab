#include "newtonian.hpp"
#include <cmath>

namespace galaxy {
namespace {
constexpr double G_SI = 6.6743e-11;
struct NewtonianRegistrar { NewtonianRegistrar(){ register_physics_package_factory("Newtonian", [](){ return std::unique_ptr<PhysicsPackage>(new NewtonianPackage());});}};
NewtonianRegistrar s_newtonian_registrar;
}

void NewtonianPackage::compute_accelerations(const State& state,double bh_mass,double softening,bool star_star,std::vector<double>& ax,std::vector<double>& ay) const {
  const int n = state.n(); ax.assign(n,0.0); ay.assign(n,0.0); const double eps2=softening*softening;
  for(int i=0;i<n;++i){
    const double rx=state.x[i], ry=state.y[i]; const double r2=rx*rx+ry*ry+eps2; const double invr3=1.0/(r2*std::sqrt(r2));
    ax[i] -= G_SI*bh_mass*rx*invr3; ay[i] -= G_SI*bh_mass*ry*invr3;
  }
  if(!star_star) return;
  for(int i=0;i<n;++i) for(int j=0;j<n;++j){ if(i==j) continue; const double dx=state.x[j]-state.x[i], dy=state.y[j]-state.y[i]; const double r2=dx*dx+dy*dy+eps2; const double invr3=1.0/(r2*std::sqrt(r2)); ax[i]+=G_SI*state.mass[j]*dx*invr3; ay[i]+=G_SI*state.mass[j]*dy*invr3; }
}

double NewtonianPackage::compute_potential_energy(const State& state,double bh_mass,double softening,bool star_star) const {
  const int n=state.n(); double pe=0.0; const double eps2=softening*softening;
  for(int i=0;i<n;++i){ double r=std::sqrt(state.x[i]*state.x[i]+state.y[i]*state.y[i]+eps2); pe -= G_SI*bh_mass*state.mass[i]/r; }
  if(!star_star) return pe;
  for(int i=0;i<n;++i) for(int j=i+1;j<n;++j){ double dx=state.x[j]-state.x[i], dy=state.y[j]-state.y[i]; double r=std::sqrt(dx*dx+dy*dy+eps2); pe -= G_SI*state.mass[i]*state.mass[j]/r; }
  return pe;
}
}
