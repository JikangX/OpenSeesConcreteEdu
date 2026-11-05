#include "JiangJianjing.h"
#include <algorithm>
#include <cmath>

namespace jiangjj {

bool Parameters::validate(std::string& error) const {
    for(double x:{fc,E0,nu,epsCrush,a,b,c,d,rho})
        if(!std::isfinite(x)) {error="non-finite parameter";return false;}
    if(!(fc>0.)) {error="fc must be positive (compressive strength, as in the original program)";return false;}
    if(!(E0>0.) || !(nu>-1. && nu<.5)) {error="invalid elastic constants";return false;}
    if(!(epsCrush<0.)) {error="epsCrush must be negative";return false;}
    if(!(a>0.)) {error="criterion constant a must be positive";return false;}
    if(rho<0.) {error="negative density";return false;}
    return true;
}

Matrix6 elasticStiffness(double E,double nu) {
    // Get_D: D = E/((1+nu)(1-2nu)) [1-nu, nu, nu; ...; (1-2nu)/2 on the shear diagonal]
    // Same operation order as Get_D (D = D*E/(1+nu)/(1-2nu), element-wise),
    // so that results agree with the Fortran to the last bit where possible.
    Matrix6 D {};
    for(int i=0;i<3;++i) {
        for(int j=0;j<3;++j) D[i][j]=nu;
        D[i][i]=1.-nu;
    }
    for(int i=3;i<6;++i) D[i][i]=(1.-2.*nu)*.5;
    for(auto& row:D) for(double& x:row) x=x*E/(1.+nu)/(1.-2.*nu);
    return D;
}

Index nonlinearIndex(const Parameters& p,const Vector6& sig) {
    Index out;
    const double sm=(sig[0]+sig[1]+sig[2])/3.;
    const double s[6]={sig[0]-sm,sig[1]-sm,sig[2]-sm,sig[3],sig[4],sig[5]};
    out.I1=3.*sm;
    out.J2=-s[0]*s[1]-s[1]*s[2]-s[2]*s[0]+s[3]*s[3]+s[4]*s[4]+s[5]*s[5];
    out.J3=s[0]*s[1]*s[2]+2.*s[3]*s[4]*s[5]-s[0]*s[4]*s[4]-s[1]*s[5]*s[5]-s[2]*s[3]*s[3];
    // J2 >= 0 analytically; a roundoff-negative value would give NaN in the
    // Fortran (sqrt). Clamp at zero (only differs from it in that NaN case).
    out.J2=std::max(0.,out.J2);
    const double r=std::sqrt(4.*out.J2/3.);
    if(r!=0.) {
        // cos(3 theta) = 4 J3 / r^3, clipped to [-1, 1] as in Get_DS.
        out.theta=std::acos(std::clamp(4.*out.J3/(r*r*r),-1.,1.))/3.;
    } else out.theta=0.;
    const double A=p.a/(p.fc*p.fc);
    const double B=(p.b+p.c*std::cos(out.theta))/p.fc;
    const double C=p.d*out.I1/p.fc-1.;
    const double disc=B*B-4.*A*C;
    // The failure surface closes on the hydrostatic axis at I1 = fc/d (mean
    // stress fc/(3d) = 0.104 fc for the default constants). For I1 >= fc/d
    // (C >= 0) the quadratic has no positive root, i.e. the stress lies
    // beyond the tensile apex. Get_DS only tests disc >= 0 and then gets
    // sqrt(J2f) <= 0, beta <= 0 and Es >= E0: the point never fails and
    // stiffens (original defect D2). Treat it like the no-root case.
    // For C < 0, disc > B*B >= 0 always; the root keeps the Get_DS formula
    // (bitwise comparison with the Fortran), which rounds to 0 within ~1 ulp
    // of the apex (beta = 0/0 on the hydrostatic axis): counted as the apex.
    const double root=(C<0. && disc>=0.) ? (-B+std::sqrt(disc))/(2.*A) : 0.;
    if(root>0.) {
        out.sqrtJ2f=root;
        out.beta=std::sqrt(out.J2)/out.sqrtJ2f;
    } else {
        out.solvable=false;
        out.beta=1.;
    }
    return out;
}

bool crushes(const Parameters& p,const Vector6& e) {
    double em=(e[0]+e[1]+e[2])/3.;
    if(em>0.) em=0.;
    return std::min({e[0]-em,e[1]-em,e[2]-em})<p.epsCrush;
}

namespace {
Vector6 scaled(const Vector6& v,double s) {Vector6 r;for(int i=0;i<6;++i) r[i]=s*v[i];return r;}
double indexOf(const Parameters& p,const Vector6& sig) {
    const Index x=nonlinearIndex(p,sig);
    return x.solvable ? x.beta : 1.;
}
// gradient of beta with respect to the stress (central differences)
Vector6 indexGradient(const Parameters& p,const Vector6& sig) {
    double m=p.fc;
    for(double x:sig) m=std::max(m,std::abs(x));
    const double h=1.e-6*m;
    Vector6 g;
    for(int i=0;i<6;++i) {
        Vector6 a=sig,b=sig;a[i]+=h;b[i]-=h;
        g[i]=(indexOf(p,a)-indexOf(p,b))/(2.*h);
    }
    return g;
}
}

// BetaFrom::Current: sigma = s D(E0) eps with s = Es/E0 = (1 + sqrt(1 - beta(sigma)))/2.
// h(s) = s - (1 + sqrt(1 - min(1, beta(s sigma_e))))/2 is increasing, h(0.5) <= 0 <= h(1).
State updateCurrent(const Parameters& p,const State& committed,const Vector6& strain,Matrix6* tangent) {
    State out=committed;
    out.strain=strain;
    const Matrix6 De=elasticStiffness(p.E0,p.nu);
    Vector6 se {};
    for(int i=0;i<6;++i) for(int j=0;j<6;++j) se[i]+=De[i][j]*strain[j];
    auto bOf=[&](double s) {return indexOf(p,scaled(se,s));};
    auto hOf=[&](double s) {return s-.5*(1.+std::sqrt(1.-std::min(1.,bOf(s))));};
    double s=1.,beta=0.;
    bool surface=false;
    if(bOf(.5)>=1.) {
        // the failure surface is reached on this ray below the secant E0/2:
        // stress on the surface, beta(s sigma_e) = 1. Take the bracket end on
        // the failed side (beta >= 1), so that the hold rule of the next
        // increment applies even where beta jumps across the surface (on the
        // hydrostatic axis beta is 0 up to the tensile apex and 1 beyond it).
        double lo=0.,hi=.5;
        for(int k=0;k<200 && hi-lo>1.e-15;++k) {const double m=.5*(lo+hi);(bOf(m)<1. ? lo : hi)=m;}
        s=hi;beta=1.;surface=true;
    } else {
        double lo=.5,hi=1.,flo=hOf(lo),fhi=hOf(hi);
        if(fhi<=0.) s=1.;
        else if(flo>=0.) s=.5;
        else {
            // Illinois regula falsi
            int side=0;
            for(int k=0;k<200 && hi-lo>1.e-15;++k) {
                double m=(lo*fhi-hi*flo)/(fhi-flo);
                if(!(m>lo && m<hi)) m=.5*(lo+hi);
                const double fm=hOf(m);
                if(fm==0.) {lo=hi=m;break;}
                if(fm<0.) {lo=m;flo=fm;if(side==-1) fhi*=.5;side=-1;}
                else {hi=m;fhi=fm;if(side==1) flo*=.5;side=1;}
            }
            s=.5*(lo+hi);
        }
        beta=std::min(1.,bOf(s));
    }
    out.stress=scaled(se,s);
    out.beta=beta;
    out.betaMax=std::max(committed.betaMax,beta);
    out.indexSolvable=nonlinearIndex(p,out.stress).solvable;
    out.secantModulus=s*p.E0;
    out.tangentModulus=std::sqrt(std::max(0.,1.-beta))*p.E0;
    if(tangent) {
        if(p.tangent==Parameters::Tangent::Secant) {
            // d sigma = s De d eps + sigma_e ds ;  ds = -k (De g) . d eps
            const Vector6 g=indexGradient(p,out.stress);
            Vector6 Dg {};
            for(int i=0;i<6;++i) for(int j=0;j<6;++j) Dg[i]+=De[i][j]*g[j];
            double gs=0.;
            for(int i=0;i<6;++i) gs+=g[i]*se[i];
            double k=0.;
            if(surface) k=gs!=0. ? s/gs : 0.;
            else if(beta<1.) {const double w=.25/std::sqrt(1.-beta);k=w*s/(1.+w*gs);}
            for(int i=0;i<6;++i) for(int j=0;j<6;++j) (*tangent)[i][j]=s*De[i][j]-k*se[i]*Dg[j];
        } else {
            const Matrix6 Dt=elasticStiffness(out.tangentModulus,p.nu);
            for(int i=0;i<6;++i) for(int j=0;j<6;++j) (*tangent)[i][j]=Dt[i][j]+1.e-6*De[i][j];
        }
    }
    return out;
}

State update(const Parameters& p,const State& committed,const Vector6& strain,Matrix6* tangent) {
    State out=committed;
    out.strain=strain;
    const Matrix6 De=elasticStiffness(p.E0,p.nu);
    // Crushing is tested on the strain at the START of the increment (E in
    // Get_DS), so it is detected one increment after the strain crosses it.
    if(committed.crushed || crushes(p,committed.strain)) {
        out.crushed=true;
        out.stress.fill(0.);
        out.secantModulus=0.;out.tangentModulus=0.;
        if(tangent) for(int i=0;i<6;++i) for(int j=0;j<6;++j) (*tangent)[i][j]=1.e-6*De[i][j];
        return out;
    }
    if(p.betaFrom==Parameters::BetaFrom::Current && committed.beta<0.99999) return updateCurrent(p,committed,strain,tangent);
    const Index idx=nonlinearIndex(p,committed.stress);
    double beta=idx.beta;
    out.indexSolvable=idx.solvable;
    out.betaMax=std::max(committed.betaMax,beta);
    if(beta>1.) beta=1.;
    out.beta=beta;
    out.secantModulus=p.E0/2.*(1.+std::sqrt(1.-beta));
    out.tangentModulus=std::sqrt(1.-beta)*p.E0;
    if(beta<0.99999) {
        const Matrix6 Ds=elasticStiffness(out.secantModulus,p.nu);
        for(int i=0;i<6;++i) {
            double s=0.;
            for(int j=0;j<6;++j) s+=Ds[i][j]*strain[j];
            out.stress[i]=s;
        }
    } else out.stress=committed.stress;
    if(tangent) {
        if(p.tangent==Parameters::Tangent::Secant && beta<0.99999) {
            *tangent=elasticStiffness(out.secantModulus,p.nu);
        } else {
            const Matrix6 Dt=elasticStiffness(out.tangentModulus,p.nu);
            for(int i=0;i<6;++i) for(int j=0;j<6;++j) (*tangent)[i][j]=Dt[i][j]+1.e-6*De[i][j];
        }
    }
    return out;
}
} // namespace jiangjj
