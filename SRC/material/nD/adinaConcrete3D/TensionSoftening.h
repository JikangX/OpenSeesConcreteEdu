// Optional Mode-I crack-band law. This extension is NOT an ADINA model-5
// equation. Keeping it in a separate pure module makes that distinction
// explicit, and leaves the legacy stress map untouched when Gf == 0.
#ifndef OPS_ADINA_TENSION_SOFTENING_H
#define OPS_ADINA_TENSION_SOFTENING_H
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace adina { namespace fracture {
/** History belonging to ONE physical crack normal, not a principal-value rank.
 * kappa is the maximum equivalent normal strain attained after formation.
 * offset accounts for the elastic transverse/Poisson contribution at onset:
 * eta = epsilon_nn - offset = elastic trial sigma_nn / E at that instant.
 * strength is frozen at the current multi-axial tensile cutoff on formation.
 * Unloading/reloading is secant to eta=0; compression is handled by ADINA.
 */
struct History {
    double kappa=0., strength=0., offset=0.;
};
inline bool valid(const History& h) {
    return std::isfinite(h.kappa) && std::isfinite(h.strength) && std::isfinite(h.offset)
        && h.kappa>=0. && h.strength>=0. && ((h.kappa==0.)==(h.strength==0.))
        && (h.strength>0. || h.offset==0.);
}
inline double onset(double E,const History& h) {return h.strength/E;}
inline double ultimate(double Gf,double length,const History& h) {return 2.*Gf/(length*h.strength);}

/** Linear cohesive traction T(w)=ft*(1-w/wc), wc=2*Gf/ft.
 * w = length*(eta-sigma/E), so eta_t=ft/E and eta_u=wc/length.
 * Eliminating w gives sigma=ft*(eta_u-eta)/(eta_u-eta_t) on softening.
 * Gf > length*ft^2/(2*E) is required for a single-valued strain law. Reject
 * an oversized crack band rather than silently changing Gf or tensile slope.
 * Units: Gf is energy/crack area (stress*length); length is crack bandwidth.
 */
inline void checkBand(double E,double Gf,double length,double strength) {
    if(!(E>0. && Gf>0. && length>0. && strength>0.) ||
       !std::isfinite(E) || !std::isfinite(Gf) || !std::isfinite(length) || !std::isfinite(strength) ||
       !(2.*Gf/length/strength > strength/E))
        throw std::runtime_error("fractureEnergy requires Gf > Lch*ft^2/(2*E); reduce crack bandwidth or supply valid Gf");
}
// Validate persisted physical history together with its material parameters.
// A checkpoint must not enable damage in the disabled/default mode, restore
// an impossible pre-onset maximum, or attach damage to a missing crack frame.
inline bool compatible(const History& h,double E,double Gf,double length,bool hasFrame) {
    if(!valid(h))return false;
    if(h.strength==0.)return true;
    if(!hasFrame || Gf<=0. || h.kappa<h.strength/E)return false;
    try {checkBand(E,Gf,length,h.strength);}catch(const std::exception&) {return false;}
    return true;
}
inline void form(History& h,double E,double Gf,double length,double strength,
                 double normalStrain,double elasticPrediction) {
    if(h.strength>0.) return; // Closure/reopening must NEVER erase damage.
    checkBand(E,Gf,length,strength);
    h.strength=strength;
    h.offset=normalStrain-elasticPrediction/E;
    h.kappa=std::max(strength/E,elasticPrediction/E);
}
inline double envelope(double E,double Gf,double length,const History& h,double kappa) {
    const double et=onset(E,h),eu=ultimate(Gf,length,h);
    if(kappa<=et) return E*kappa;
    if(kappa>=eu) return 0.;
    return h.strength*(eu-kappa)/(eu-et);
}
inline double response(double E,double Gf,double length,History& h,double normalStrain) {
    if(h.strength==0.) throw std::runtime_error("unformed fracture plane");
    const double eta=normalStrain-h.offset;
    if(eta<=0.) return 0.;
    h.kappa=std::max(h.kappa,eta);
    return envelope(E,Gf,length,h,h.kappa)*eta/h.kappa;
}
inline bool isLoading(const History& h,double normalStrain) {
    // Equality uses the loading-side tangent at a committed maximum. The
    // replay Jacobian must not average loading and secant-unloading slopes.
    const double tol=64.*std::numeric_limits<double>::epsilon()
        *std::max({1.,std::abs(normalStrain),std::abs(h.offset),std::abs(h.kappa)});
    return h.strength>0. && normalStrain-h.offset>=h.kappa-tol;
}
inline double dissipation(double E,double Gf,double length,const History& h) {
    if(h.strength==0.) return 0.;
    const double et=onset(E,h),eu=ultimate(Gf,length,h),k=std::min(h.kappa,eu);
    if(k<=et) return 0.;
    const double stress=envelope(E,Gf,length,h,k);
    // Total monotonic work minus recoverable secant elastic energy. Damage
    // is held during unloading, so this quantity depends only on kappa.
    const double work=.5*h.strength*et+.5*(h.strength+stress)*(k-et);
    const double recoverable=.5*stress*k;
    return length*(work-recoverable);
}
} } // adina::fracture
#endif
