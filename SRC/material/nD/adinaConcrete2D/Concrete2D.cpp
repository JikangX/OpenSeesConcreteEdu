#include "Concrete2D.h"
#include "../adinaConcrete3D/BranchTangent.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace adina { namespace plane {
namespace {
bool exceeds(double stress,double cutoff) {
    // Equality at a discontinuous failure surface must classify identically
    // after rotation. This is a roundoff tolerance, not strength smoothing.
    const double tol=64.*std::numeric_limits<double>::epsilon()
                     *std::max(std::abs(stress),std::abs(cutoff));
    return stress>cutoff+tol;
}
double measure(const Parameters& p, const Vector3& s) {
    // sqrt(J2)+alpha*trace(sigma), with sigma_zz=0 (manual p.107).
    return std::sqrt((s[0]*s[0]+s[1]*s[1]-s[0]*s[1])/3.+s[2]*s[2])
           +p.alpha*(s[0]+s[1]);
}
// Smallest in-plane principal strain; e[2] is the engineering shear strain.
double minPrincipalStrain(const Vector3& e) {
    return .5*(e[0]+e[1])-std::hypot(.5*(e[0]-e[1]),.5*e[2]);
}
double average(const Parameters& p,const Strength& s,double a,double b) {
    // CDMOD's three-point Gauss integration of the rational tangent.
    constexpr double g=.77459666924148337704;
    const double m=.5*(a+b),h=.5*(b-a);
    return (5.*compressionTangent(p,s,m-g*h)+8.*compressionTangent(p,s,m)
            +5.*compressionTangent(p,s,m+g*h))/18.;
}
Vector3 product(const Matrix3& c,const Vector3& e) {
    Vector3 s{};
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) s[i]+=c[i][j]*e[j];
    return s;
}

/** DCRACK local 4-component stiffness, followed by sigma_zz condensation.
 * Local order here is [11,22,gamma12,33], exactly A06's order. Keeping the
 * fourth row until after crack degradation reproduces its plane-stress law;
 * deleting that row before degradation gives incorrect Poisson coupling.
 */
Matrix3 stiffness(const Parameters& p,const Vector3& s,
                  const std::array<double,3>& yp,bool ortho,
                  const std::array<bool,2>& open,
                  std::array<double,2> retention={{-1.,-1.}}) {
    // retention[i] scales the normal stiffness of open crack plane i. The
    // stress update always uses STIFAC for both planes (A06 DCRACK); only the
    // secant iteration matrix supplies a larger per-plane cohesive value.
    for(double& r:retention) if(r<0.) r=p.normalRetention;
    const double den=std::abs(s[0])+std::abs(s[1]);
    const double E=den>1.e-5*p.ft ? (std::abs(s[0])*yp[0]+std::abs(s[1])*yp[1])/den : p.E;
    // CDMOD distinguishes the total-modulus threshold from its pair modulus:
    // E12 is weighted whenever either in-plane normal stress is nonzero.
    // The distinction matters in independently supplied orthotropic states
    // near zero stress, even though normal integration usually stays isotropic
    // there. Reproduce the original pure subroutine for either input branch.
    const double e12=den>0. ? (std::abs(s[0])*yp[0]+std::abs(s[1])*yp[1])/den : E;
    // E14=YP1 and E24=YP2 because sigma_zz=0. At a zero pair denominator
    // CDMOD uses the weighted E as a finite fallback.
    const double e14=std::abs(s[0])>0. ? yp[0] : E;
    const double e24=std::abs(s[1])>0. ? yp[1] : E;
    const double d=(1.+p.nu)*(1.-2.*p.nu);
    std::array<std::array<double,4>,4> c{};
    c[0][0]=(1.-p.nu)*(ortho ? yp[0] : E)/d;
    c[1][1]=(1.-p.nu)*(ortho ? yp[1] : E)/d;
    c[3][3]=(1.-p.nu)*(ortho ? yp[2] : E)/d;
    c[0][1]=p.nu*(ortho ? e12 : E)/d;
    c[0][3]=p.nu*(ortho ? e14 : E)/d;
    c[1][3]=p.nu*(ortho ? e24 : E)/d;
    c[2][2]=(ortho ? e12 : E)/(2.*(1.+p.nu));
    const int n=int(open[0])+int(open[1]);
    if(n) {
        c[2][2]*=p.shearRetention;
        if(n==1) {
            const int a=open[0] ? 0 : 1,b=1-a;
            for(int j=0;j<4;++j) {
                const int i=std::min(a,j),k=std::max(a,j);
                c[i][k]*=retention[a];
            }
            c[b][b]=(ortho ? yp[b] : E)/(1.-p.nu*p.nu);
            c[b][3]=p.nu*(ortho ? (b==0 ? e14 : e24) : E)/(1.-p.nu*p.nu);
            c[3][3]=(ortho ? yp[2] : E)/(1.-p.nu*p.nu);
        } else {
            // With equal factors this is DCRACK's single STIFAC scaling.
            for(int i=0;i<2;++i) for(int j=i;j<4;++j)
                c[i][j]*=(i==0 && j==1 ? std::min(retention[0],retention[1]) : retention[i]);
            c[3][3]=ortho ? yp[2] : E;
        }
    }
    for(int i=0;i<4;++i) for(int j=i+1;j<4;++j) c[j][i]=c[i][j];
    if(!std::isfinite(c[3][3]) || c[3][3]==0.)
        throw std::runtime_error("singular out-of-plane condensation in ADINA 2D");
    Matrix3 out{};
    for(int i=0;i<3;++i) for(int j=0;j<3;++j)
        out[i][j]=c[i][j]-c[i][3]*c[3][j]/c[3][3];
    return out;
}

// CRAKID's closure/reopening, expressed with independent physical flags.
// The original two-crack branch can silently swap which plane remains active.
bool classify(const Parameters& p,State& state,const Vector3& s,const Vector3& e) {
    const auto before=state.open;
    for(int i=0;i<2;++i) {
        const double tol=64.*std::numeric_limits<double>::epsilon()
            *std::max({1.,std::abs(e[i]),std::abs(state.crackStrain[i])});
        // Strain is dimensionless. Near zero, a purely relative tolerance
        // collapses to ~1e-17 and mistakes accumulated nodal-arithmetic noise
        // for closure. Different points then close a whole increment apart.
        // The absolute roundoff floor is ~1e-14 (not a softening strain or a
        // strength change); resolved compressive strains still close the crack.
        if(p.regularizedTension() && state.tension[i].strength>0.) {
            // A damaged plane reopens as its gap becomes positive; imposing
            // the virgin ft threshold again would spuriously heal the crack.
            const double eta=e[i]-state.tension[i].offset;
            if(eta<-tol) state.open[i]=false;
            else if(eta>tol) state.open[i]=true;
        } else if(state.open[i]) {
            if(e[i]<-tol && e[i]<state.crackStrain[i]-tol) state.open[i]=false;
        } else if(exceeds(s[i],p.ft)) {
            state.open[i]=true;state.crackStrain[i]=e[i];
        }
    }
    return state.open!=before;
}
State step(const Parameters& p,const State& old,const Vector3& target) {
    State out=old;out.strain=target;
    const double f=measure(p,old.stress);
    // Loading: the loading function has reached its recorded maximum. In
    // the original update f never exceeds the maximum (it is raised to f at
    // every update), so this equals CDMOD's equality test; Hybrid origins can
    // carry the smaller converged maximum, see Parameters::Update::Hybrid.
    const bool loading=f>=old.maximumLoading-std::abs(old.maximumLoading)*1.e-8;
    double angle=old.hasCrackFrame && !old.crushed ? old.angle : principalAngle(old.stress);
    Vector3 sig=stressToLocal(old.stress,angle);
    // PRNCPL's legacy EPSL(3)=0 assumes strain and stress are coaxial. The
    // NEW strain is transformed in the OLD principal stress frame during an
    // update, so its shear generally is nonzero under nonproportional loading.
    // Keep all three components; otherwise MOD45=2 loses added shear stress.
    Vector3 e0=strainToLocal(old.strain,angle),e1=strainToLocal(target,angle);
    if(p.regularizedTension() && old.hasCrackFrame && !old.crushed) {
        // Crossing eta=0 recovers the undamaged compressive branch through
        // zero traction, not by adding an elastic increment to a damaged
        // tensile stress. Preserve kappa: closing never heals tensile damage.
        for(int i=0;i<2;++i) if(old.tension[i].strength>0.) {
            const double a=e0[i]-old.tension[i].offset,b=e1[i]-old.tension[i].offset;
            if((a>0. && b<=0.) || (a<0. && b>=0.)) {sig[i]=0.;e0[i]=old.tension[i].offset;}
        }
    }
    Strength strength=old.crushed ? Strength{old.crushingStress,old.crushingStrain,p.ft}
                                  : strengthEnvelope(p,{{sig[0],sig[1],0.}});
    if(old.crushed) {
        const double ultimate=p.epsu*strength.peakStrain/p.epsc;
        Vector3 de{};for(int i=0;i<3;++i) de[i]=target[i]-old.strain[i];
        if(p.crushedUnloading==Parameters::CrushedUnloading::Elastic && !old.released) {
            // Optional extension (see Parameters::CrushedUnloading). Strains
            // less compressive than the crushed history unload/reload with E;
            // principal stresses are capped at zero. Points already released
            // at the ultimate strain stay released (original law below).
            const double emin=minPrincipalStrain(target);
            if(emin>old.crushedStrain) {
                Vector3 s=product(elasticTangent(p.E,p.nu),de);
                for(int i=0;i<3;++i) s[i]+=old.stress[i];
                const double a=principalAngle(s);
                Vector3 local=stressToLocal(s,a);
                // No tension; compression bounded by the crushed curve at the
                // most compressive strain reached (the unloading point).
                const double cap=compressionStress(p,strength,old.crushedStrain);
                for(int i=0;i<2;++i) local[i]=std::clamp(local[i],cap,0.);
                local[2]=0.;
                out.stress=stressToGlobal(local,a);
                return out;
            }
            out.crushedStrain=emin;
        }
        // CDMOD retains the peak stress, then descends with the tangent along
        // the most compressed principal direction; beyond scaled epsu all
        // stress is released. sigma_zz stays zero throughout.
        double E=p.E*p.normalRetention;
        // CDMOD label 9, crushed: N1 selects the more compressed of P2 and
        // sigma_zz(=0), and YP keeps E (so the modulus becomes 0) unless that
        // stress is below 0.001*SIGCP. Softening therefore stops once the
        // smallest principal stress has left compression: the stress is frozen,
        // it is not driven into tension by a negative isotropic modulus.
        if(f!=0.) E=sig[1]<.001*strength.compression
                    ? std::min(0.,average(p,strength,e0[1],e1[1])) : 0.;
        const Vector3 ds=product(elasticTangent(E,p.nu),de);
        for(int i=0;i<3;++i) out.stress[i]=old.stress[i]+ds[i];
        // CDMOD label 60 recomputes EP in the principal frame of the UPDATED
        // stress and DCRACK (NUMCRK=3) keeps it only while EP(1), EP(2) and the
        // plane-stress thickness strain EP(4) all exceed the scaled epsu.
        const Vector3 e=strainToLocal(target,principalAngle(out.stress));
        const double ezz=p.nu/(p.nu-1.)*(target[0]+target[1]);
        if(std::min({e[0],e[1],ezz})<=ultimate) {
            out.stress.fill(0.);out.released=true;
            out.crushedStrain=std::min(out.crushedStrain,ultimate);
        }
        return out;
    }
    bool transition=false;
    if(old.hasCrackFrame) transition=classify(p,out,sig,e0);
    if(p.regularizedTension() && old.hasCrackFrame) {
        // Choose closure/reopening stiffness for the actual target branch.
        // Virgin plane formation still uses the elastic stress prediction.
        for(int i=0;i<2;++i) if(out.tension[i].strength>0.)
            out.open[i]=e1[i]>out.tension[i].offset;
    }
    std::array<double,3> yp{{p.E,p.E,p.E}};
    if(loading && std::abs(sig[0])+std::abs(sig[1])>1.e-5*p.ft)
        for(int i=0;i<2;++i) if(sig[i]<.001*strength.compression)
            yp[i]=average(p,strength,e0[i],e1[i]);
    const bool ortho=std::min(sig[0],sig[1])<p.kappa*strength.compression;
    const Matrix3 c=stiffness(p,sig,yp,ortho,out.open);
    Vector3 de{};for(int i=0;i<3;++i) de[i]=e1[i]-e0[i];
    const Vector3 ds=product(c,de);
    for(int i=0;i<3;++i) sig[i]+=ds[i];
    if(!old.hasCrackFrame) {
        // Determine the FIRST fixed crack normal from predicted stress. atan2
        // works at zero shear and at repeated eigenvalues; no approximate
        // degree/radian conversion or 0.0001-degree sentinel is used.
        const Vector3 prediction=stressToGlobal(sig,angle);
        angle=principalAngle(prediction);
        sig=stressToLocal(prediction,angle);e1=strainToLocal(target,angle);
        strength=strengthEnvelope(p,{{sig[0],sig[1],0.}});
        for(int i=0;i<2;++i) if(exceeds(sig[i],strength.tension)) {
            out.open[i]=true;out.crackStrain[i]=e1[i];transition=true;
        }
        if(transition) {out.hasCrackFrame=true;out.angle=angle;}
    } else {
        transition=classify(p,out,sig,e1) || transition;
        // CRAKID uses the fixed-frame normal stresses for the strength check.
        strength=strengthEnvelope(p,{{sig[0],sig[1],0.}});
    }
    if(std::min(sig[0],sig[1])<=strength.compression) {
        out.crushed=true;out.crushingStress=strength.compression;
        out.crushingStrain=strength.peakStrain;
        out.crushedStrain=minPrincipalStrain(target);
        // DCRACK also tests the plane-stress thickness strain TEP(4).
        const double ezz=p.nu/(p.nu-1.)*(target[0]+target[1]);
        const double ultimate=p.epsu*strength.peakStrain/p.epsc;
        if(std::min({e1[0],e1[1],ezz})<=ultimate) {
            sig.fill(0.);out.released=true;out.crushedStrain=std::min(out.crushedStrain,ultimate);
        }
    } else {
        for(int i=0;i<2;++i) if(out.open[i]) {
            if(p.regularizedTension()) {
                const double cutoff=old.hasCrackFrame ? p.ft : strength.tension;
                fracture::form(out.tension[i],p.E,p.fractureEnergy,p.characteristicLength,
                               cutoff,e1[i],sig[i]);
                sig[i]=fracture::response(p.E,p.fractureEnergy,p.characteristicLength,out.tension[i],e1[i]);
            } else sig[i]=0.;
        }
        if((out.open[0] || out.open[1]) && p.shearRetention<=.001) sig[2]=0.;
    }
    out.stress=stressToGlobal(sig,angle);
    if(!out.crushed) {
        const double next=measure(p,out.stress);
        if(loading && transition) out.maximumLoading=next;
        out.maximumLoading=std::max(out.maximumLoading,next);
    }
    return out;
}
State replay(const Parameters& p,const State& committed,const Vector3& target,int steps) {
    State out=committed;
    for(int k=1;k<=steps;++k) {
        Vector3 e{};
        for(int i=0;i<3;++i) e[i]=committed.strain[i]+(target[i]-committed.strain[i])*double(k)/steps;
        out=step(p,out,e);
        for(double x:out.stress) if(!std::isfinite(x)) throw std::runtime_error("non-finite 2D stress");
        if(!std::isfinite(out.maximumLoading)) throw std::runtime_error("non-finite 2D loading history");
    }
    return out;
}
/** Material (secant) iteration matrix of a converged/trial state.
 * This is the matrix CDMOD rebuilds after a stress update (labels 69-110):
 * - crushed: isotropic E*STIFAC (CDMOD: E<=0 -> EKK(1)*STIFAC);
 * - unloading (loading function below its maximum): initial modulus E;
 * - loading: rational-curve tangents YP of the local strains (ITEM=2, one
 *   point), stress-weighted as in the stress update; MOD45 orthotropy;
 * - open cracks: DCRACK retention. With fracture energy, an open plane uses
 *   its cohesive secant s(kappa)/kappa, so the matrix bounds the unloading/
 *   reloading stiffness instead of under-estimating it by ~1/STIFAC.
 * Non-positive curve tangents are replaced by the curve secant (bounded below
 * by STIFAC*E): the matrix stays symmetric positive semidefinite. It is not
 * the derivative of the stress map and does not change any converged state.
 */
Matrix3 secantStiffness(const Parameters& p,const State& s) {
    if(s.crushed) {
        // Elastic crushed unloading (optional): E while the point is in its
        // compressive unloading/reloading range, otherwise ADINA's E*STIFAC.
        if(p.crushedUnloading==Parameters::CrushedUnloading::Elastic && !s.released &&
           minPrincipalStrain(s.strain)>s.crushedStrain && std::min(s.stress[0],s.stress[1])<0.)
            return elasticTangent(p.E,p.nu);
        return elasticTangent(p.E*std::max(p.normalRetention,1.e-12),p.nu);
    }
    const double angle=s.hasCrackFrame ? s.angle : principalAngle(s.stress);
    const Vector3 sig=stressToLocal(s.stress,angle),eps=strainToLocal(s.strain,angle);
    const Strength strength=strengthEnvelope(p,{{sig[0],sig[1],0.}});
    const double f=measure(p,s.stress);
    // Same test as CDMOD's IKAS; the loading function can be negative (alpha).
    const bool loading=f>=s.maximumLoading-std::abs(s.maximumLoading)*1.e-8;
    std::array<double,3> yp{{p.E,p.E,p.E}};
    if(loading && std::abs(sig[0])+std::abs(sig[1])>1.e-5*p.ft)
        for(int i=0;i<2;++i) if(sig[i]<.001*strength.compression) {
            double t=compressionTangent(p,strength,eps[i]);
            if(!(t>0.)) t=eps[i]<0. ? compressionStress(p,strength,eps[i])/eps[i] : p.E;
            yp[i]=std::max(t,p.normalRetention*p.E);
        }
    const bool ortho=std::min(sig[0],sig[1])<p.kappa*strength.compression;
    std::array<double,2> retention{{p.normalRetention,p.normalRetention}};
    if(p.regularizedTension()) for(int i=0;i<2;++i) if(s.open[i] && s.tension[i].strength>0.) {
        const auto& h=s.tension[i];
        const double secant=fracture::envelope(p.E,p.fractureEnergy,p.characteristicLength,h,h.kappa)/h.kappa;
        retention[i]=std::max(p.normalRetention,std::min(1.,secant/p.E));
    }
    const Matrix3 local=stiffness(p,sig,yp,ortho,s.open,retention);
    // Global matrix R^T C R with R the strain rotation; stressToGlobal = R^T.
    Matrix3 out{};
    for(int j=0;j<3;++j) {
        Vector3 unit{};unit[j]=1.;
        const Vector3 column=stressToGlobal(product(local,strainToLocal(unit,angle)),angle);
        for(int i=0;i<3;++i) out[i][j]=column[i];
    }
    for(int i=0;i<3;++i) for(int j=i+1;j<3;++j) out[i][j]=out[j][i]=.5*(out[i][j]+out[j][i]);
    return out;
}
bool sameActiveBranch(const State& a,const State& b) {
    return a.hasCrackFrame==b.hasCrackFrame && a.crushed==b.crushed && a.open==b.open;
}
} // anonymous namespace

Vector3 strainToLocal(const Vector3& e,double a) {
    const double c=std::cos(a),s=std::sin(a);
    return {{c*c*e[0]+s*s*e[1]+c*s*e[2],s*s*e[0]+c*c*e[1]-c*s*e[2],
             -2.*c*s*(e[0]-e[1])+(c*c-s*s)*e[2]}};
}
Vector3 strainToGlobal(const Vector3& e,double a) {return strainToLocal(e,-a);}
Vector3 stressToLocal(const Vector3& e,double a) {
    const double c=std::cos(a),s=std::sin(a);
    return {{c*c*e[0]+s*s*e[1]+2.*c*s*e[2],s*s*e[0]+c*c*e[1]-2.*c*s*e[2],
             -c*s*(e[0]-e[1])+(c*c-s*s)*e[2]}};
}
Vector3 stressToGlobal(const Vector3& e,double a) {return stressToLocal(e,-a);}
double principalAngle(const Vector3& s) {return .5*std::atan2(2.*s[2],s[0]-s[1]);}
Matrix3 elasticTangent(double E,double nu) {
    Matrix3 c{};c[0][0]=c[1][1]=E/(1.-nu*nu);
    c[0][1]=c[1][0]=nu*c[0][0];c[2][2]=E/(2.*(1.+nu));return c;
}
State hybridOrigin(const State& committed,const State& iterate) {
    State origin=iterate;
    if(iterate.hasCrackFrame==committed.hasCrackFrame && iterate.crushed==committed.crushed &&
       iterate.released==committed.released && iterate.open==committed.open &&
       iterate.angle==committed.angle && iterate.crackStrain==committed.crackStrain) {
        origin.maximumLoading=committed.maximumLoading;
        for(int i=0;i<2;++i)
            if(iterate.tension[i].strength==committed.tension[i].strength &&
               iterate.tension[i].offset==committed.tension[i].offset)
                origin.tension[i].kappa=committed.tension[i].kappa;
    }
    return origin;
}
bool validateState(const State& s,std::string& error) {
    for(double x:s.strain) if(!std::isfinite(x)) {error="non-finite committed strain";return false;}
    for(double x:s.stress) if(!std::isfinite(x)) {error="non-finite committed stress";return false;}
    for(double x:{s.maximumLoading,s.angle,s.crackStrain[0],s.crackStrain[1],s.crushingStress,s.crushingStrain})
        if(!std::isfinite(x)) {error="non-finite committed history";return false;}
    for(const auto& h:s.tension) if(!fracture::valid(h)) {error="invalid tensile fracture history";return false;}
    if((!s.hasCrackFrame && (s.open[0] || s.open[1])) ||
       (s.crushed && (s.crushingStress>=0. || s.crushingStrain>=0.)) || (s.released && !s.crushed)) {
        error="inconsistent crack/crushing history";return false;
    }
    return true;
}
bool integrate(const Parameters& p,const State& committed,const Vector3& target,
               State& trial,Matrix3* tangent,std::string& error) {
    error.clear();
    if(!p.validate(error) || !validateState(committed,error)) return false;
    for(double x:target) if(!std::isfinite(x)) {error="non-finite trial strain";return false;}
    const double dx=target[0]-committed.strain[0],dy=target[1]-committed.strain[1],
                 dg=target[2]-committed.strain[2];
    // maxStrainIncrement 0: no subdivision (one update, as CDMOD).
    const double ratio=p.maxStrainIncrement>0. ? std::hypot(std::hypot(dx,dy),dg/std::sqrt(2.))/p.maxStrainIncrement : 0.;
    if(!std::isfinite(ratio) || ratio>10000.) {error="2D substep limit (10000) exceeded";return false;}
    const int n=std::max(1,int(std::ceil(ratio-32.*std::numeric_limits<double>::epsilon()*std::max(1.,ratio))));
    try {
        const State next=replay(p,committed,target,n);
        Matrix3 c{};
        if(tangent && p.tangent==Parameters::Tangent::Secant) c=secantStiffness(p,next);
        else if(tangent) c=numerics::replayJacobian<3>(target,next,p.epsc,
            [&](const Vector3& e) {return replay(p,committed,e,n);},
            [&](const State& state) {
                if(!sameActiveBranch(state,next)) return false;
                if(p.regularizedTension()) {
                    const auto a=strainToLocal(state.strain,state.angle),b=strainToLocal(next.strain,next.angle);
                    for(int i=0;i<2;++i) if(state.open[i] &&
                       fracture::isLoading(state.tension[i],a[i])!=fracture::isLoading(next.tension[i],b[i])) return false;
                }
                return true;
            });
        trial=next;if(tangent) *tangent=c;return true;
    } catch(const std::exception& ex) {error=ex.what();return false;}
}
} } // namespace adina::plane
