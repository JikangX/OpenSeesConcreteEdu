#include "Concrete3D.h"
#include "BranchTangent.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

namespace adina {
namespace {
constexpr int pairIndex[3][2] = {{0,1},{1,2},{0,2}};
using Tensor = Frame;

// Compare tensile eigenvalues in stress units. Rotating an identical tensor
// perturbs a root by a few ulps; a strict comparison at exactly ft can then
// release a complete stress component one loading step earlier. Use the same
// roundoff tie rule for first cracking and later planes. Preserve CMOD3D's
// strict first-onset comparison and CRAK3D's inclusive comparisons when all
// planes are closed or the third plane forms. This is not a material-strength
// tolerance: only a 64-epsilon neighborhood is treated as equality. The factor
// is dimensionless, so Pa and MPa descriptions make the same decision.
bool exceedsTensileCutoff(double stress, double cutoff, bool inclusive=false) {
    const double roundoff=64.*std::numeric_limits<double>::epsilon()
        *std::max(std::abs(stress),std::abs(cutoff));
    return inclusive ? stress>=cutoff-roundoff : stress>cutoff+roundoff;
}

Tensor tensor(const Vector6& x, bool engineering) {
    Tensor a {};
    for (int i=0; i<3; ++i) a[i][i]=x[i];
    for (int i=0; i<3; ++i) {
        const int j=pairIndex[i][0], k=pairIndex[i][1];
        a[j][k]=a[k][j]=x[i+3]*(engineering ? .5 : 1.);
    }
    return a;
}
Vector6 voigt(const Tensor& a, bool engineering) {
    Vector6 x {};
    for (int i=0; i<3; ++i) x[i]=a[i][i];
    for (int i=0; i<3; ++i)
        x[i+3]=a[pairIndex[i][0]][pairIndex[i][1]]*(engineering ? 2. : 1.);
    return x;
}
Tensor rotate(const Tensor& a, const Frame& q, bool toLocal) {
    Tensor b {};
    for (int i=0; i<3; ++i) for (int j=0; j<3; ++j)
        for (int k=0; k<3; ++k) for (int l=0; l<3; ++l)
            b[i][j]+=(toLocal ? q[k][i]*a[k][l]*q[l][j]
                              : q[i][k]*a[k][l]*q[j][l]);
    return b;
}
Vector6 multiply(const Matrix6& a, const Vector6& x) {
    Vector6 y {};
    for (int i=0; i<6; ++i) for (int j=0; j<6; ++j) y[i]+=a[i][j]*x[j];
    return y;
}
struct Curve {
    double a,b,c;
};
// Manual p.105: rational compressive law sigma=E*eps/(1+A*x+B*x²+C*x³).
// Its derivative is the printed tangent formula. Peak and ultimate data fix
// the coefficients; confinement scales stress AND strain (manual Fig.3-9-9).
Curve coefficients(const Parameters& p, const Strength& s) {
    const double r=p.epsu/p.epsc;
    const double es=s.compression/s.peakStrain;
    const double eu=(p.fu*s.compression/p.fc)/(p.epsu*s.peakStrain/p.epsc);
    const double a=(p.E/eu+(r-2.)*r*r*p.E/es-(2.*r+1.)*(r-1.)*(r-1.))
                   /(r*(r-1.)*(r-1.));
    return {a,2.*p.E/es-3.-2.*a,2.-p.E/es+a};
}

// Smallest principal value of an engineering-shear strain vector.
double minPrincipalStrain(const Vector6& e) {
    Vector6 tensorial=e;
    for (int i=3; i<6; ++i) tensorial[i]*=.5;
    return principalStress(tensorial).values[2];
}

// Three Gauss points reproduce CMOD3D's increment-averaged tangent rather
// than substituting a different damage/plasticity model for the ADINA law.
double averageTangent(const Parameters& p, const Strength& s, double a, double b) {
    constexpr double g=.77459666924148337704;
    const double m=.5*(a+b), h=.5*(b-a);
    return (5.*compressionTangent(p,s,m-g*h)+8.*compressionTangent(p,s,m)
            +5.*compressionTangent(p,s,m+g*h))/18.;
}

/** DCRAK3's local stiffness, in OpenSees shear order.
 * At low compression ADINA uses a stress-weighted isotropic modulus. Beyond
 * kappa*fc it uses three directional tangent moduli and pairwise weighted
 * coupling moduli. A crack releases its normal stress; STIFAC is the residual
 * incremental stiffness, not a tensile softening curve or fracture energy.
 */
Matrix6 localStiffness(const Parameters& p, const Vector6& sig,
                      const std::array<double,3>& yp, bool orthotropic, int cracks,
                      std::array<double,3> retention={{-1.,-1.,-1.}}) {
    // retention[i] scales the normal stiffness of open crack plane i. The
    // stress update always uses STIFAC (DCRAK3); only the secant iteration
    // matrix supplies a larger per-plane cohesive value.
    for(double& r:retention) if(r<0.) r=p.normalRetention;
    double sum=0., E=0.;
    for (int i=0; i<3; ++i) { sum+=std::abs(sig[i]); E+=std::abs(sig[i])*yp[i]; }
    // CMOD3D uses 1e-5*SIGMAT for the TOTAL stress-weighted modulus.
    // This is a constitutive branch threshold; replacing it with a guard
    // proportional to Young's modulus changes the original stiffness near
    // zero stress. Pairwise coupling below deliberately has a different rule.
    E=sum>1.e-5*p.ft ? E/sum : p.E;
    Matrix6 c=elasticTangent(E,p.nu);
    double pair[3] {};
    for (int k=0; k<3; ++k) {
        int i=pairIndex[k][0],j=pairIndex[k][1];
        const double den=std::abs(sig[i])+std::abs(sig[j]);
        // Original CMOD3D tests DENM != 0 for each pair, without the total
        // modulus threshold. Since the weights are absolute stresses the
        // denominator is nonnegative and a positive value is well defined.
        pair[k]=den>0. ? (std::abs(sig[i])*yp[i]+std::abs(sig[j])*yp[j])/den : E;
    }
    if (orthotropic) {
        const double d=(1.+p.nu)*(1.-2.*p.nu);
        for (int i=0; i<3; ++i) c[i][i]=(1.-p.nu)*yp[i]/d;
        for (int k=0; k<3; ++k) {
            int i=pairIndex[k][0],j=pairIndex[k][1];
            c[i][j]=c[j][i]=p.nu*pair[k]/d;
            c[k+3][k+3]=pair[k]/(2.*(1.+p.nu));
        }
    }
    if (cracks==0 || cracks==4) return c;
    for (int k=3; k<6; ++k) c[k][k]*=p.shearRetention;
    if (cracks==1) {
        for (int i=0; i<3; ++i) c[0][i]*=retention[0];
        // Condense the released first normal direction (plane stress), but do
        // not condense its tiny retained diagonal again.
        c[1][1]=(orthotropic ? yp[1] : E)/(1.-p.nu*p.nu);
        c[2][2]=(orthotropic ? yp[2] : E)/(1.-p.nu*p.nu);
        c[1][2]=p.nu*(orthotropic ? pair[1] : E)/(1.-p.nu*p.nu);
        c[4][4]=(orthotropic ? pair[1] : E)/(2.*(1.+p.nu));
    } else if (cracks==2) {
        // Equal factors reproduce DCRAK3's single STIFAC scaling exactly; a
        // coupling between two open planes takes the smaller factor.
        for (int i=0; i<2; ++i) for (int j=i; j<3; ++j)
            c[i][j]*=(j<2 && i!=j ? std::min(retention[i],retention[j]) : retention[i]);
        c[2][2]=orthotropic ? yp[2] : E;
    } else {
        for (int i=0; i<3; ++i) for (int j=i; j<3; ++j)
            c[i][j]*=(i!=j ? std::min(retention[i],retention[j]) : retention[i]);
    }
    for (int i=0; i<3; ++i) for (int j=i+1; j<3; ++j) c[j][i]=c[i][j];
    return c;
}

// Preserve the first crack plane while diagonalizing the transverse 2x2 block
// when a second plane forms (CRAK3D). Storing the full frame avoids recovering
// a direction's z component by an always-positive square root.
void rotateTransverse(State& state, Vector6& sig, Vector6& eps) {
    const double angle=.5*std::atan2(2.*sig[4],sig[1]-sig[2]);
    const double c=std::cos(angle),s=std::sin(angle);
    Frame q {{{{1.,0.,0.}},{{0.,c,-s}},{{0.,s,c}}}};
    sig=stressToLocal(sig,q); eps=strainToLocal(eps,q);
    const Frame old=state.frame;
    for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) {
        state.frame[i][j]=0.;
        for (int k=0; k<3; ++k) state.frame[i][j]+=old[i][k]*q[k][j];
    }
}

// Canonicalize the local frame: currently open planes precede closed planes.
// A count alone is sufficient ONLY with this invariant. Legacy CRAK3D reduced
// its count without permuting axes when (for example) the first of three cracks
// closed, and consequently released the wrong normal stress. Transform the
// complete tensor and formation history whenever the ordering changes.
void orderOpenPlanes(State& state,Vector6& sig,Vector6& eps,
                     const std::array<bool,3>& active) {
    std::array<int,3> order {};int count=0;
    for(int i=0;i<3;++i) if(active[i]) order[count++]=i;
    const int open=count;
    for(int i=0;i<3;++i) if(!active[i]) order[count++]=i;
    Frame permutation {};
    for(int j=0;j<3;++j) permutation[order[j]][j]=1.;
    sig=stressToLocal(sig,permutation);
    eps=strainToLocal(eps,permutation);
    const Frame oldFrame=state.frame;
    const auto oldHistory=state.crackStrain;
    const auto oldTension=state.tension;
    for(int j=0;j<3;++j) {
        state.crackStrain[j]=oldHistory[order[j]];
        state.tension[j]=oldTension[order[j]];
        for(int i=0;i<3;++i) state.frame[i][j]=oldFrame[i][order[j]];
    }
    state.openCracks=open;
}

// Explicit CRAK3D state machine, with axis ordering corrected on closure.
// Only temporary flags are needed: State's frame/count/history together encode
// all active planes. No global integer/sentinel doubles encode physical state.
void updateCracks(const Parameters& p,State& state,Vector6& sig,Vector6& eps) {
    const int old=state.openCracks;
    std::array<bool,3> active {{old>=1,old>=2,old>=3}};
    if(p.regularizedTension()) {
        int formed=0;for(const auto& h:state.tension) formed+=h.strength>0.;
        if(formed) {
            // Fracture energy belongs to a physical plane. After two planes
            // have damage histories, rotating them into new principal axes
            // would rotate/erase dissipated energy and spuriously heal damage.
            const double tol=64.*std::numeric_limits<double>::epsilon();
            for(int i=0;i<3;++i) if(state.tension[i].strength>0.) {
                const double eta=eps[i]-state.tension[i].offset;
                if(eta>tol)active[i]=true;
                else if(eta<-tol)active[i]=false;
            }
            if(formed==1 && state.tension[0].strength>0.) {
                const double m=.5*(sig[1]+sig[2]),r=std::hypot(.5*(sig[1]-sig[2]),sig[4]);
                const double cutoff=m-r<0. ? p.ft*(1.-(m-r)/p.fc) : p.ft;
                // Beyond the compressive failure surface this interaction
                // cutoff can become nonpositive. That is a crushing event,
                // not a new tensile cohesive plane with negative strength.
                // Leave it to step()'s compression check; never call form()
                // with a nonphysical strength in the tensile extension.
                if(cutoff>0. && m+r>cutoff) {
                    rotateTransverse(state,sig,eps);
                    for(int i=1;i<3;++i) if(sig[i]>cutoff) {
                        active[i]=true;state.crackStrain[i]=eps[i];
                        fracture::form(state.tension[i],p.E,p.fractureEnergy,p.characteristicLength,
                                       cutoff,eps[i],sig[i]);
                    }
                }
            } else for(int i=0;i<3;++i) if(state.tension[i].strength==0. && sig[i]>p.ft) {
                active[i]=true;state.crackStrain[i]=eps[i];
                fracture::form(state.tension[i],p.E,p.fractureEnergy,p.characteristicLength,p.ft,eps[i],sig[i]);
            }
            orderOpenPlanes(state,sig,eps,active);return;
        }
    }
    auto closed=[&](int i) {
        // Roundoff in a 3D rotation must not turn exact zero normal strain
        // into a compressive closure one substep earlier than an axial path.
        const double tol=64.*std::numeric_limits<double>::epsilon()
            *std::max({1.,std::abs(eps[i]),std::abs(state.crackStrain[i])});
        // Dimensionless strain needs an absolute roundoff floor at zero.
        // Accumulated nodal arithmetic must not close equivalent crack planes
        // one full loading increment apart in differently oriented elements.
        return eps[i]<-tol && eps[i]<state.crackStrain[i]-tol;
    };
    if(old==0) {
        for(int i=0;i<3;++i) if(exceedsTensileCutoff(sig[i],p.ft,true)) {
            active[i]=true;state.crackStrain[i]=eps[i];
        }
    } else if(old==1) {
        if(closed(0)) active[0]=false;
        const double m=.5*(sig[1]+sig[2]);
        const double r=std::hypot(.5*(sig[1]-sig[2]),sig[4]);
        double cutoff=p.ft;
        if(m-r<0.) cutoff=p.ft*(1.-(m-r)/p.fc);
        if(exceedsTensileCutoff(m+r,cutoff)) {
            // Diagonalize even when shear==0: the larger transverse principal
            // stress may be in the third, rather than second, stored direction.
            rotateTransverse(state,sig,eps);
            active[1]=true;state.crackStrain[1]=eps[1];
            if(exceedsTensileCutoff(m-r,cutoff)) {active[2]=true;state.crackStrain[2]=eps[2];}
        }
    } else if(old==2) {
        for(int i=0;i<2;++i) if(closed(i)) active[i]=false;
        if(exceedsTensileCutoff(sig[2],p.ft,true)) {active[2]=true;state.crackStrain[2]=eps[2];}
    } else {
        for(int i=0;i<3;++i) if(closed(i)) active[i]=false;
    }
    orderOpenPlanes(state,sig,eps,active);
}
void releaseStress(const Parameters& p, int cracks, Vector6& sig) {
    for (int i=0; i<cracks && i<3; ++i) sig[i]=0.;
    if (p.shearRetention<=.001 && cracks>0) {
        sig[3]=sig[5]=0.;
        if (cracks>1) sig[4]=0.;
    }
}

State step(const Parameters& p, const State& old, const Vector6& target) {
    State out=old;
    out.strain=target;
    const double fOld=loadingMeasure(p,old.stress);
    // Loading: the loading function has reached its recorded maximum. In
    // the original update f never exceeds the maximum (it is raised to f at
    // every update), so this equals CMOD3D's equality test; Hybrid origins can
    // carry the smaller converged maximum, see Parameters::Update::Hybrid.
    bool loading=fOld>=old.maximumLoading-std::abs(old.maximumLoading)*1.e-8;
    // PRNCP3 is also used after crushing, even when an earlier fixed crack
    // frame exists. Its smallest principal direction drives the crushed law.
    Frame frame=(old.hasCrackFrame && !old.crushed) ? old.frame : principalStress(old.stress).frame;
    Vector6 sig=stressToLocal(old.stress,frame);
    Vector6 eps0=strainToLocal(old.strain,frame),eps1=strainToLocal(target,frame);
    if(p.regularizedTension() && old.hasCrackFrame && !old.crushed) {
        // Classify the target gap before constructing the incremental matrix,
        // then transform all old quantities into that SAME physical frame.
        updateCracks(p,out,sig,eps1);frame=out.frame;
        eps0=strainToLocal(old.strain,frame);
        for(int i=0;i<3;++i) if(out.tension[i].strength>0.) {
            const double a=eps0[i]-out.tension[i].offset,b=eps1[i]-out.tension[i].offset;
            if((a>0. && b<=0.) || (a<0. && b>=0.)) {sig[i]=0.;eps0[i]=out.tension[i].offset;}
        }
    }
    // CMOD3D obtains SIGCP/EPSCP (curve and MOD45 threshold) from the KKK=1
    // call of CRAK3D. With exactly one open crack (label 111) it sorts the
    // crack-normal stress together with the two PRINCIPAL values of the
    // transverse 2x2 block, not its normal components; with 0 or 2 cracks it
    // uses the crack-frame normal stresses (label 107), which is what `sig`
    // holds. The distinction matters whenever transverse shear is present.
    std::array<double,3> envelopeStress {{sig[0],sig[1],sig[2]}};
    if(old.hasCrackFrame && !old.crushed && old.openCracks==1) {
        const double m=.5*(sig[1]+sig[2]),r=std::hypot(.5*(sig[1]-sig[2]),sig[4]);
        envelopeStress={{sig[0],m+r,m-r}};
    }
    Strength strength=old.crushed ? Strength{old.crushingStress,old.crushingStrain,p.ft}
                                 : strengthEnvelope(p,envelopeStress);
    std::array<double,3> yp {{p.E,p.E,p.E}};
    if (loading && std::abs(sig[0])+std::abs(sig[1])+std::abs(sig[2])>1.e-5*p.ft)
        for (int i=0; i<3; ++i)
            if (sig[i]<.001*strength.compression)
                yp[i]=averageTangent(p,strength,eps0[i],eps1[i]);
    bool orthotropic=*std::min_element(sig.begin(),sig.begin()+3)<p.kappa*strength.compression;
    if (old.crushed) {
        const double ultimate=p.epsu*strength.peakStrain/p.epsc;
        Vector6 delta {};
        for (int i=0; i<6; ++i) delta[i]=target[i]-old.strain[i];
        if (p.crushedUnloading==Parameters::CrushedUnloading::Elastic && old.crushedStrain>ultimate) {
            // Optional extension (Parameters::CrushedUnloading): unload and
            // reload with E above the crushed strain history, no tension.
            const double emin=minPrincipalStrain(target);
            if (emin>old.crushedStrain) {
                Vector6 trialStress=multiply(elasticTangent(p.E,p.nu),delta);
                for (int i=0; i<6; ++i) trialStress[i]+=old.stress[i];
                const Principal pr=principalStress(trialStress);
                Vector6 local {};
                // No tension; compression bounded by the crushed curve at the
                // most compressive strain reached (the unloading point).
                const double cap=compressionStress(p,strength,old.crushedStrain);
                for (int i=0; i<3; ++i) local[i]=std::clamp(pr.values[i],cap,0.);
                out.stress=stressToGlobal(local,pr.frame);
                return out;
            }
            out.crushedStrain=emin;
        }
        // Original crushed branch uses a reduced isotropic increment. It keeps
        // PGRAV=100 irreversible; beyond scaled epsu it releases ALL stress.
        double E=p.E*p.normalRetention;
        // CMOD3D (PGRAV=100): YP(3)=E unless SIGP(3) < 0.001*SIGCP, and then
        // E=min(YP(3),0). Once the smallest principal stress is no longer
        // compressive the modulus is 0: softening stops, the stress is frozen.
        if (fOld!=0.) E=sig[2]<.001*strength.compression
                        ? std::min(0.,averageTangent(p,strength,eps0[2],eps1[2])) : 0.;
        out.stress=multiply(elasticTangent(E,p.nu),delta);
        for (int i=0; i<6; ++i) out.stress[i]+=old.stress[i];
        // CMOD3D label 60 calls PRNCP3 AFTER updating STRESS, then DCRAK3
        // IK=4 checks the target normal strains in that updated stress frame.
        // Under a nonproportional increment the axes can rotate appreciably;
        // the old frame can release all stress one increment too early/late.
        const auto releaseStrain=strainToLocal(target,principalStress(out.stress).frame);
        if (*std::min_element(releaseStrain.begin(),releaseStrain.begin()+3)<=ultimate) {
            out.stress.fill(0.);
            out.crushedStrain=std::min(out.crushedStrain,ultimate); // released for good
        }
        return out;
    }
    if (old.hasCrackFrame && !p.regularizedTension()) {
        // First classify the previous local state, as CMOD3D's KKK=1 call does.
        updateCracks(p,out,sig,eps0);
        frame=out.frame;
        eps1=strainToLocal(target,frame);
    }
    const Matrix6 c=localStiffness(p,sig,yp,orthotropic,out.openCracks);
    Vector6 de {};
    for (int i=0; i<6; ++i) de[i]=eps1[i]-eps0[i];
    const Vector6 ds=multiply(c,de);
    for (int i=0; i<6; ++i) sig[i]+=ds[i];
    bool transition=false;
    if (!old.hasCrackFrame) {
        const Vector6 prediction=stressToGlobal(sig,frame);
        const Principal pr=principalStress(prediction);
        strength=strengthEnvelope(p,pr.values);
        frame=pr.frame;
        sig=stressToLocal(prediction,frame);
        eps1=strainToLocal(target,frame);
        int cracks=0;
        // In an all-compressive state FALSTR stores the compressive envelope,
        // not a tensile cutoff. CMOD3D label 150 bypasses tensile checks then.
        if (pr.values[0]>=0.)
            for (int i=0; i<3; ++i) {
                // Use the same objective tie rule in the original and the
                // optional continuous law. Material behavior away from this
                // roundoff neighborhood is unchanged.
                if (exceedsTensileCutoff(pr.values[i],strength.tension)) ++cracks;
            }
        if (cracks>0) {
            out.hasCrackFrame=true; out.frame=frame; out.openCracks=cracks;
            std::copy_n(eps1.begin(),3,out.crackStrain.begin()); transition=true;
        }
        if (pr.values[2]<=strength.compression) {
            out.crushed=true; transition=true;
        }
    } else {
        const int n=out.openCracks;
        updateCracks(p,out,sig,eps1);
        frame=out.frame;
        transition=(n!=out.openCracks || n!=old.openCracks);
        // CRAK3D only evaluates the compression envelope with 0 or 1 open
        // cracks; two open cracks use the remaining direction's fc cutoff.
        if (out.openCracks<=1) {
            auto stresses=std::array<double,3>{{sig[0],sig[1],sig[2]}};
            if (out.openCracks==1) {
                const double m=.5*(sig[1]+sig[2]);
                const double r=std::hypot(.5*(sig[1]-sig[2]),sig[4]);
                stresses[1]=m+r; stresses[2]=m-r;
            }
            strength=strengthEnvelope(p,stresses);
            if (*std::min_element(stresses.begin(),stresses.end())<=strength.compression)
                out.crushed=true;
        } else if (out.openCracks==2 && sig[2]<=p.fc) out.crushed=true;
    }
    if (out.crushed) {
        out.crushingStress=strength.compression; out.crushingStrain=strength.peakStrain;
        // DCRAK3, IK=4 retains the predicted stress at peak failure. Stress is
        // released only at the scaled ultimate strain, not at first crushing.
        out.stress=stressToGlobal(sig,frame);
        out.crushedStrain=minPrincipalStrain(target);
        const double ultimate=p.epsu*strength.peakStrain/p.epsc;
        if (*std::min_element(eps1.begin(),eps1.begin()+3)<=ultimate) {
            out.stress.fill(0.);out.crushedStrain=std::min(out.crushedStrain,ultimate);
        }
    } else {
        if(p.regularizedTension()) {
            for(int i=0;i<out.openCracks;++i) {
                fracture::form(out.tension[i],p.E,p.fractureEnergy,p.characteristicLength,
                               old.hasCrackFrame ? p.ft : strength.tension,eps1[i],sig[i]);
                sig[i]=fracture::response(p.E,p.fractureEnergy,p.characteristicLength,out.tension[i],eps1[i]);
            }
            if(p.shearRetention<=.001 && out.openCracks>0) {
                sig[3]=sig[5]=0.;if(out.openCracks>1)sig[4]=0.;
            }
        } else releaseStress(p,out.openCracks,sig);
        out.stress=stressToGlobal(sig,frame);
        const double f=loadingMeasure(p,out.stress);
        if (loading && transition) out.maximumLoading=f;
        out.maximumLoading=std::max(out.maximumLoading,f);
    }
    return out;
}

/** Material (secant) iteration matrix of a trial state: the matrix CMOD3D
 * rebuilds after its stress update when IEQREF=0 (labels 74 -> 5 -> 100).
 * - crushed: isotropic E*STIFAC;
 * - unloading: initial modulus E;
 * - loading: rational-curve tangents at the current local strains (ITEM=2,
 *   one Gauss point), stress-weighted; orthotropic beyond kappa*sigma_c;
 * - open cracks: DCRAK3 retention; with fracture energy an open plane uses its
 *   cohesive secant s(kappa)/kappa.
 * Non-positive curve tangents are replaced by the curve secant, bounded below
 * by STIFAC*E, so the matrix is symmetric positive semidefinite. It is an
 * iteration matrix only: it never enters the stress update.
 */
Matrix6 secantStiffness(const Parameters& p,const State& s) {
    if(s.crushed) {
        // Elastic crushed unloading (optional): E in the compressive
        // unloading/reloading range, otherwise ADINA's E*STIFAC.
        const double ultimate=p.epsu*s.crushingStrain/p.epsc;
        if(p.crushedUnloading==Parameters::CrushedUnloading::Elastic && s.crushedStrain>ultimate &&
           minPrincipalStrain(s.strain)>s.crushedStrain && principalStress(s.stress).values[2]<0.)
            return elasticTangent(p.E,p.nu);
        return elasticTangent(p.E*std::max(p.normalRetention,1.e-12),p.nu);
    }
    const Frame frame=s.hasCrackFrame ? s.frame : principalStress(s.stress).frame;
    const Vector6 sig=stressToLocal(s.stress,frame),eps=strainToLocal(s.strain,frame);
    const Strength strength=strengthEnvelope(p,{{sig[0],sig[1],sig[2]}});
    const double f=loadingMeasure(p,s.stress);
    const bool loading=f>=s.maximumLoading-std::abs(s.maximumLoading)*1.e-8;
    std::array<double,3> yp {{p.E,p.E,p.E}};
    if(loading && std::abs(sig[0])+std::abs(sig[1])+std::abs(sig[2])>1.e-5*p.ft)
        for(int i=0;i<3;++i) if(sig[i]<.001*strength.compression) {
            double t=compressionTangent(p,strength,eps[i]);
            if(!(t>0.)) t=eps[i]<0. ? compressionStress(p,strength,eps[i])/eps[i] : p.E;
            yp[i]=std::max(t,p.normalRetention*p.E);
        }
    const bool orthotropic=*std::min_element(sig.begin(),sig.begin()+3)<p.kappa*strength.compression;
    std::array<double,3> retention {{p.normalRetention,p.normalRetention,p.normalRetention}};
    if(p.regularizedTension()) for(int i=0;i<s.openCracks;++i) if(s.tension[i].strength>0.) {
        const auto& h=s.tension[i];
        const double secant=fracture::envelope(p.E,p.fractureEnergy,p.characteristicLength,h,h.kappa)/h.kappa;
        retention[i]=std::max(p.normalRetention,std::min(1.,secant/p.E));
    }
    const Matrix6 local=localStiffness(p,sig,yp,orthotropic,s.openCracks,retention);
    Matrix6 out {};
    for(int j=0;j<6;++j) {
        Vector6 unit {};unit[j]=1.;
        const Vector6 column=stressToGlobal(multiply(local,strainToLocal(unit,frame)),frame);
        for(int i=0;i<6;++i) out[i][j]=column[i];
    }
    for(int i=0;i<6;++i) for(int j=i+1;j<6;++j) out[i][j]=out[j][i]=.5*(out[i][j]+out[j][i]);
    return out;
}

State replay(const Parameters& p,const State& committed,const Vector6& target,int steps) {
    State state=committed;
    for (int k=1; k<=steps; ++k) {
        Vector6 eps {};
        for (int i=0; i<6; ++i)
            eps[i]=committed.strain[i]+(target[i]-committed.strain[i])*(double(k)/steps);
        state=step(p,state,eps);
        for (double x:state.stress)
            if (!std::isfinite(x)) throw std::runtime_error("non-finite stress in ADINA concrete integration");
        if(!std::isfinite(state.maximumLoading))
            throw std::runtime_error("non-finite loading history in ADINA concrete integration");
    }
    return state;
}
} // anonymous namespace

Principal principalStress(const Vector6& stress) {
    Tensor a=tensor(stress,false);
    Frame q {{{{1.,0.,0.}},{{0.,1.,0.}},{{0.,0.,1.}}}};
    // Symmetric Jacobi eigensolver. Unlike the legacy eigenvector elimination,
    // this remains well-defined for hydrostatic and repeated eigenvalues.
    double scale=0.; for (const auto& row:a) for(double x:row) scale=std::max(scale,std::abs(x));
    for (int it=0; it<40; ++it) {
        int p=0,r=1;
        for (int i=0; i<3; ++i) for (int j=i+1; j<3; ++j)
            if (std::abs(a[i][j])>std::abs(a[p][r])) {p=i;r=j;}
        if (std::abs(a[p][r])<=1.e-14*std::max(scale,1.e-100)) break;
        const double theta=.5*std::atan2(2.*a[p][r],a[r][r]-a[p][p]);
        const double c=std::cos(theta),s=std::sin(theta);
        const double app=a[p][p],arr=a[r][r],apr=a[p][r];
        a[p][p]=c*c*app-2.*s*c*apr+s*s*arr;
        a[r][r]=s*s*app+2.*s*c*apr+c*c*arr;
        a[p][r]=a[r][p]=0.;
        for(int k=0;k<3;++k) {
            if(k!=p && k!=r) {
                double x=a[k][p],y=a[k][r];
                a[k][p]=a[p][k]=c*x-s*y; a[k][r]=a[r][k]=s*x+c*y;
            }
            double x=q[k][p],y=q[k][r]; q[k][p]=c*x-s*y; q[k][r]=s*x+c*y;
        }
    }
    std::array<int,3> order {{0,1,2}};
    std::stable_sort(order.begin(),order.end(),[&](int i,int j){return a[i][i]>a[j][j];});
    Principal out;
    for (int j=0;j<3;++j) {
        out.values[j]=a[order[j]][order[j]];
        for(int i=0;i<3;++i) out.frame[i][j]=q[i][order[j]];
    }
    return out;
}
Vector6 strainToLocal(const Vector6& v,const Frame& q) {return voigt(rotate(tensor(v,true),q,true),true);}
Vector6 strainToGlobal(const Vector6& v,const Frame& q) {return voigt(rotate(tensor(v,true),q,false),true);}
Vector6 stressToLocal(const Vector6& v,const Frame& q) {return voigt(rotate(tensor(v,false),q,true),false);}
Vector6 stressToGlobal(const Vector6& v,const Frame& q) {return voigt(rotate(tensor(v,false),q,false),false);}
Matrix6 elasticTangent(double E,double nu) {
    Matrix6 c {};
    const double G=E/(2.*(1.+nu)),lambda=E*nu/((1.+nu)*(1.-2.*nu));
    for (int i=0;i<3;++i) for(int j=0;j<3;++j) c[i][j]=lambda+(i==j ? 2.*G : 0.);
    for(int i=3;i<6;++i) c[i][i]=G;
    return c;
}
Strength strengthEnvelope(const Parameters& p,std::array<double,3> sig) {
    std::sort(sig.begin(),sig.end(),std::greater<double>());
    Strength out {p.fc,p.epsc,p.ft};
    if(sig[1]>=0.) {
        if(sig[2]<0.) out.tension=p.ft*(1.-sig[2]/p.fc);
        return out;
    }
    double t=0.,a=p.sp31[0],b=p.sp32[0],c=p.sp33[0];
    if(sig[0]<0.) {
        t=sig[0]/p.fc;
        if(t>=p.sp1.back()) throw std::runtime_error("confinement exceeds SP1 envelope range");
        const auto hi=std::upper_bound(p.sp1.begin(),p.sp1.end(),t);
        const int j=int(hi-p.sp1.begin()),i=j-1;
        const double f=(t-p.sp1[i])/(p.sp1[j]-p.sp1[i]);
        a=p.sp31[i]+f*(p.sp31[j]-p.sp31[i]);
        b=p.sp32[i]+f*(p.sp32[j]-p.sp32[i]);
        c=p.sp33[i]+f*(p.sp33[j]-p.sp33[i]);
    }
    const double r=sig[1]/p.fc;
    out.compression=(r>p.beta*b)
        ? b*p.fc+(c-b)/(c-p.beta*b)*(sig[1]-p.beta*b*p.fc)
        : a*p.fc+(b-a)/(p.beta*b-t)*(sig[1]-(sig[0]>0. ? 0. : sig[0]));
    if(out.compression<p.fc) out.peakStrain=p.gamma*p.epsc*out.compression/p.fc;
    if(sig[0]>=0.) out.tension=p.ft*(1.-sig[1]/out.compression)*(1.-sig[2]/out.compression);
    else out.tension=out.compression;
    if(out.tension<.001*p.ft && sig[0]>=0.) out.tension=p.ft;
    return out;
}
double compressionStress(const Parameters& p,const Strength& s,double eps) {
    // The scalar reference curve is piecewise: CDMOD/CMOD3D use the initial
    // modulus on the tensile side (manual pp.105,130). Keep this public helper
    // consistent with compressionTangent for eps >= 0 as well; extending the
    // rational compression polynomial into tension gives a different slope.
    // Crack stress release is applied by the integration state machine, not
    // by this pre-failure scalar curve.
    if(eps>=0.) return p.E*eps;
    const Curve c=coefficients(p,s); const double x=eps/s.peakStrain;
    return p.E*eps/(1.+x*(c.a+x*(c.b+x*c.c)));
}
double compressionTangent(const Parameters& p,const Strength& s,double eps) {
    if(eps>=0.) return p.E;
    const Curve c=coefficients(p,s); const double x=eps/s.peakStrain;
    const double d=1.+x*(c.a+x*(c.b+x*c.c));
    return p.E*(1.-c.b*x*x-2.*c.c*x*x*x)/(d*d);
}
double loadingMeasure(const Parameters& p,const Vector6& s) {
    const double j2=((s[0]-s[1])*(s[0]-s[1])+(s[1]-s[2])*(s[1]-s[2])
                     +(s[2]-s[0])*(s[2]-s[0]))/6.+s[3]*s[3]+s[4]*s[4]+s[5]*s[5];
    const double f=std::sqrt(std::max(0.,j2))+p.alpha*(s[0]+s[1]+s[2]);
    if(!std::isfinite(f)) throw std::runtime_error("stress magnitude exceeds loading-function numeric range");
    return f;
}
bool Parameters::validate(std::string& error) const {
    const double scalars[]={E,nu,ft,fc,epsc,fu,epsu,beta,gamma,kappa,alpha,
                           normalRetention,shearRetention,rho,maxStrainIncrement,fractureEnergy,characteristicLength};
    for(double x:scalars) if(!std::isfinite(x)) {error="non-finite material parameter";return false;}
    if(fractureEnergy!=0. || characteristicLength!=0.) {
        try {fracture::checkBand(E,fractureEnergy,characteristicLength,ft);}
        catch(const std::exception& ex) {error=ex.what();return false;}
    }
    if(!(E>0. && nu>-1. && nu<.5 && ft>0. && fc<0. && epsc<0. &&
         fu<0. && fu>=fc && epsu<epsc && fc/epsc<=E && beta>0. && beta<1. &&
         gamma>0. && kappa>=0. && kappa<=1. && normalRetention>=0. &&
         normalRetention<=1. && shearRetention>=0. && shearRetention<=1. &&
         rho>=0. && maxStrainIncrement>=0.)) {
        error="invalid signs, strength/strain ordering, elastic constants or retention factors";return false;
    }
    for(int i=0;i<6;++i) {
        if(!std::isfinite(sp1[i]) || !std::isfinite(sp31[i]) || !std::isfinite(sp32[i]) ||
           !std::isfinite(sp33[i]) || (i==0 && sp1[i]!=0.) ||
           (i>0 && sp1[i]<=sp1[i-1]) || sp31[i]<=sp1[i] ||
           beta*sp32[i]<=sp1[i] || sp33[i]<=beta*sp32[i]) {
            error="invalid six-point triaxial envelope or singular interpolation denominator";return false;
        }
        if(i>0 && (sp31[i]<sp31[i-1] || sp32[i]<sp32[i-1] || sp33[i]<sp33[i-1])) {
            error="triaxial failure curves must be nondecreasing";return false;
        }
    }
    const Curve c=coefficients(*this,{fc,epsc,ft});
    // A sampled denominator check can miss a narrow pole. A cubic attains
    // its minimum at an endpoint or a root of its quadratic derivative.
    const double end=epsu/epsc;
    auto positive=[&](double x) {
        const double d=1.+x*(c.a+x*(c.b+x*c.c));
        return std::isfinite(d) && d>0.;
    };
    bool valid=positive(0.) && positive(end);
    if(c.c==0.) {
        if(c.b!=0.) {const double x=-c.a/(2.*c.b);if(x>0. && x<end) valid=valid && positive(x);}
    } else {
        const double discriminant=c.b*c.b-3.*c.c*c.a;
        if(discriminant>=0.) for(double sign:{-1.,1.}) {
            const double x=(-c.b+sign*std::sqrt(discriminant))/(3.*c.c);
            if(x>0. && x<end) valid=valid && positive(x);
        }
    }
    if(!valid) {error="singular rational compression curve";return false;}
    error.clear();return true;
}
State hybridOrigin(const State& committed,const State& iterate) {
    State origin=iterate;
    if(iterate.hasCrackFrame==committed.hasCrackFrame && iterate.crushed==committed.crushed &&
       iterate.openCracks==committed.openCracks && iterate.frame==committed.frame &&
       iterate.crackStrain==committed.crackStrain) {
        origin.maximumLoading=committed.maximumLoading;
        for(int i=0;i<3;++i)
            if(iterate.tension[i].strength==committed.tension[i].strength &&
               iterate.tension[i].offset==committed.tension[i].offset)
                origin.tension[i].kappa=committed.tension[i].kappa;
    }
    return origin;
}
bool integrate(const Parameters& p,const State& committed,const Vector6& strain,
               State& trial,Matrix6* tangent,std::string& error) {
    if(!p.validate(error)) return false;
    for(double x:strain) if(!std::isfinite(x)) {error="non-finite trial strain";return false;}
    for(double x:committed.strain) if(!std::isfinite(x)) {error="non-finite committed strain";return false;}
    for(double x:committed.stress) if(!std::isfinite(x)) {error="non-finite committed stress";return false;}
    for(const auto& h:committed.tension) if(!fracture::valid(h)) {error="invalid tensile fracture history";return false;}
    if(!std::isfinite(committed.maximumLoading) || committed.openCracks<0 || committed.openCracks>3) {
        error="invalid committed loading/crack history";return false;
    }
    // Frobenius norm of the STRAIN TENSOR is invariant under rotations.
    // A largest-Voigt-component criterion would pick different partitions
    // for the same physical path in differently oriented meshes.
    double delta=0.;
    for(int i=0;i<6;++i)
        delta=std::hypot(delta,(strain[i]-committed.strain[i])*(i<3 ? 1. : std::sqrt(.5)));
    // Guard conversion and pathological work before converting ceil to int.
    // maxStrainIncrement 0: no subdivision (one update, as CMOD3D).
    const double ratio=p.maxStrainIncrement>0. ? delta/p.maxStrainIncrement : 0.;
    if(!std::isfinite(ratio) || ratio>10000.) {error="trial increment requires more than 10000 substeps";return false;}
    const double roundoff=32.*std::numeric_limits<double>::epsilon()*std::max(1.,ratio);
    const int n=std::max(1,int(std::ceil(ratio-roundoff)));
    try {
        State result=replay(p,committed,strain,n);
        Matrix6 jac {};
        // Keep the integration grid fixed while differentiating the replay.
        if(tangent && p.tangent==Parameters::Tangent::Secant) jac=secantStiffness(p,result);
        else if(tangent) jac=numerics::replayJacobian<6>(strain,result,p.epsc,
            [&](const Vector6& e) {return replay(p,committed,e,n);},
            [&](const State& state) {
                // At FIRST formation the regularized traction is continuous.
                // The multiaxial cutoff contains |transverse shear|: exactly
                // at the peak BOTH signed perturbations can form a crack.
                // There is no derivative confined to the uncracked branch in
                // that direction. Use the centered generalized derivative of
                // this continuous corner, but never cross a crushing jump or
                // apply this exception to the brittle/default stress map.
                if(p.regularizedTension() && !committed.hasCrackFrame &&
                   state.crushed==result.crushed && !result.crushed &&
                   state.hasCrackFrame!=result.hasCrackFrame) return true;
                if(state.hasCrackFrame!=result.hasCrackFrame ||
                   state.crushed!=result.crushed || state.openCracks!=result.openCracks) return false;
                if(p.regularizedTension()) {
                    const auto a=strainToLocal(state.strain,state.frame),b=strainToLocal(result.strain,result.frame);
                    for(int i=0;i<state.openCracks;++i) if(
                        fracture::isLoading(state.tension[i],a[i])!=fracture::isLoading(result.tension[i],b[i])) return false;
                }
                return true;
            });
        trial=result;if(tangent) *tangent=jac;error.clear();return true;
    } catch(const std::exception& e) {error=e.what();return false;}
}
} // namespace adina
