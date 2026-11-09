// Numerical differentiation of an immutable constitutive replay, shared by
// the 2D and 3D kernels. No OpenSees types, global scratch or material law.
#ifndef OPS_ADINA_BRANCH_TANGENT_H
#define OPS_ADINA_BRANCH_TANGENT_H
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace adina { namespace numerics {
/** Differentiate the actual active trial branch at a fixed integration grid.
 * Brittle crack/crushing transitions release finite stress: a difference that
 * crosses that jump divided by a microscopic strain is not a stiffness.
 * Smooth columns use central differences, a switching column uses the matching
 * one-sided branch. At corners, shrink the perturbation before declaring that
 * no admissible branch derivative can be resolved. The caller retains full
 * ownership of history and decides which physical states belong to its branch.
 */
template<std::size_t N, class State, class Replay, class SameBranch>
std::array<std::array<double,N>,N> replayJacobian(
    const std::array<double,N>& target, const State& trial,
    double referenceStrain, Replay replay, SameBranch sameBranch) {
    std::array<std::array<double,N>,N> result {};
    for(std::size_t j=0;j<N;++j) {
        const double scale=std::max({std::abs(referenceStrain),std::abs(target[j]),1.e-4});
        double h=1.e-8*scale;
        const double floor=64.*std::numeric_limits<double>::epsilon()*scale;
        bool resolved=false;
        for(int attempt=0;attempt<10 && h>=floor;++attempt,h*=.1) {
            auto plus=target,minus=target;plus[j]+=h;minus[j]-=h;
            const State a=replay(plus),b=replay(minus);
            const bool forward=sameBranch(a),backward=sameBranch(b);
            if(!forward && !backward) continue;
            for(std::size_t i=0;i<N;++i) {
                if(forward && backward) result[i][j]=(a.stress[i]-b.stress[i])/(2.*h);
                else if(forward) result[i][j]=(a.stress[i]-trial.stress[i])/h;
                else result[i][j]=(trial.stress[i]-b.stress[i])/h;
                if(!std::isfinite(result[i][j]))
                    throw std::runtime_error("non-finite algorithmic tangent");
            }
            resolved=true;break;
        }
        if(!resolved) throw std::runtime_error("no admissible branch for algorithmic tangent column "+std::to_string(j));
    }
    return result;
}
} }
#endif
