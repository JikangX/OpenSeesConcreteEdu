// Independent material-point checks. Compile without OpenSees or Fortran.
#include "Concrete3D.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <random>
#include <stdexcept>
#include <string>
using namespace adina;
namespace {
int checks=0;
void check(bool ok,const char* message) {
    ++checks; if(!ok) throw std::runtime_error(message);
}
void near(double a,double b,double tolerance,const char* message) {
    if(std::abs(a-b)>tolerance*std::max({1.,std::abs(a),std::abs(b)}))
        std::cerr<<message<<": "<<std::setprecision(17)<<a<<" vs "<<b<<'\n';
    check(std::abs(a-b)<=tolerance*std::max({1.,std::abs(a),std::abs(b)}),message);
}
State advance(const Parameters& p,const State& c,const Vector6& e,Matrix6* d=nullptr) {
    State t;std::string error;
    if(!integrate(p,c,e,t,d,error)) throw std::runtime_error(error);
    return t;
}
void tests() {
    Parameters p;std::string error;check(p.validate(error),"default parameters");
    // Derivative checks below target the optional algorithmic matrix, not
    // ADINA's positive secant iteration matrix used by the new default.
    p.tangent=Parameters::Tangent::Algorithmic;
    // Default: committed origin and no subdivision (one update, as CMOD3D).
    check(p.maxStrainIncrement==0. && p.update==Parameters::Update::Committed,"ADINA defaults");
    // The checks below were written for subdivided increments (2e-5).
    p.maxStrainIncrement=2.e-5;
    for(int i=0;i<12;++i) {
        auto bad=p;
        switch(i) {case 0:bad.E=0;break;case 1:bad.nu=.5;break;case 2:bad.fc=30;break;
        case 3:bad.ft=0;break;case 4:bad.epsu=bad.epsc;break;case 5:bad.fu=0;break;
        case 6:bad.sp1[1]=0;break;case 7:bad.sp32[5]=1;break;
        case 8:bad.normalRetention=-1;break;case 9:bad.shearRetention=2;break;
        case 10:bad.maxStrainIncrement=-1;break;case 11:bad.E=INFINITY;break;}
        check(!bad.validate(error),"reject invalid parameters");
    }
    Strength s{p.fc,p.epsc,p.ft};
    near(compressionStress(p,s,0),0,1e-13,"curve origin");
    near(compressionStress(p,s,1.e-5),p.E*1.e-5,1e-13,"linear tensile scalar curve");
    near(compressionTangent(p,s,1.e-5),p.E,1e-13,"linear tensile scalar tangent");
    near(compressionStress(p,s,p.epsc),p.fc,1e-13,"curve peak");
    near(compressionStress(p,s,p.epsu),p.fu,1e-13,"curve ultimate");
    near(compressionTangent(p,s,p.epsc),0,1e-12,"zero peak slope");
    // Numerical derivative of the STRESS formula is independent of tangent.
    for(int i=1;i<240;++i) {
        const double e=p.epsu*i/240.,h=1e-8;
        near(compressionTangent(p,s,e),
             (compressionStress(p,s,e+h)-compressionStress(p,s,e-h))/(2*h),
             2e-6,"rational derivative");
    }
    auto b=strengthEnvelope(p,{{0,-19.5,-39}});
    near(b.compression,-39,1e-13,"biaxial envelope knot SP32");
    auto c=strengthEnvelope(p,{{-7.5,-22.5,-45}});
    near(c.compression,-45,1e-13,"triaxial envelope knot SP32");
    near(c.peakStrain,-.003,1e-13,"confinement peak strain");
    near(strengthEnvelope(p,{{0,0,-15}}).tension,1.5,1e-13,"compression tensile interaction");
    near(loadingMeasure(p,{{-30,0,0,0,0,0}}),30/std::sqrt(3.),1e-13,"J2 invariant");
    std::mt19937 generator(1701);std::uniform_real_distribution<double> distribution(-20,20);
    for(int n=0;n<500;++n) {
        Vector6 stress;for(double& x:stress)x=distribution(generator);
        if(n==0) stress={{-3,-3,-3,0,0,0}};
        if(n==1) stress={{4,4,-2,0,0,0}};
        const auto pr=principalStress(stress);
        auto local=stressToLocal(stress,pr.frame);
        auto global=stressToGlobal(local,pr.frame);
        for(int i=0;i<6;++i) near(global[i],stress[i],2e-12,"spectral reconstruction");
        for(int i=0;i<3;++i) {
            near(local[i],pr.values[i],2e-12,"sorted eigenvalues");
            near(local[i+3],0,2e-12,"principal shear zero");
            for(int j=0;j<3;++j) {
                double dot=0;for(int k=0;k<3;++k)dot+=pr.frame[k][i]*pr.frame[k][j];
                near(dot,i==j?1:0,2e-12,"orthogonal eigenvectors");
            }
        }
        // Engineering shear must preserve stress : strain under rotation.
        Vector6 e;for(double& x:e)x=distribution(generator)*1e-7;
        auto el=strainToLocal(e,pr.frame);double work=0,workLocal=0;
        for(int i=0;i<6;++i){work+=e[i]*stress[i];workLocal+=el[i]*local[i];}
        near(work,workLocal,1e-13,"stress strain work conjugacy");
    }
    const auto elastic=elasticTangent(p.E,p.nu);
    for(int j=0;j<6;++j) {
        Vector6 e{};e[j]=1e-7;Matrix6 tangent;
        auto t=advance(p,{},e,&tangent);
        for(int i=0;i<6;++i) {
            near(t.stress[i],elastic[i][j]*e[j],1e-12,"six-component elasticity");
            for(int k=0;k<6;++k)near(tangent[i][k],elastic[i][k],2e-9,"elastic algorithmic tangent");
        }
    }
    auto uncoupled=p;uncoupled.nu=0;
    for(double offset: {-1.e-13,1.e-13}) {
        Matrix6 jac;
        auto onset=advance(uncoupled,{},{{uncoupled.ft/uncoupled.E+offset,0,0,0,0,0}},&jac);
        near(jac[0][0],offset<0 ? uncoupled.E : 0.,2e-6,"one-sided brittle crack tangent");
        check(onset.hasCrackFrame==(offset>0),"onset branch preserved");
    }
    State t=advance(p,{},{{.0002,0,0,0,0,0}});
    check(t.hasCrackFrame && t.openCracks==1,"first crack formation");
    near(t.stress[0],0,1e-12,"normal tensile stress release");
    auto crack=t;
    t=advance(p,t,{{.0002,.0003,0,0,0,0}});check(t.openCracks==2,"second crack formation");
    t=advance(p,t,{{.0002,.0003,.0004,0,0,0}});check(t.openCracks==3,"third crack formation");
    const State threeCracks=t;
    for(double residual: {-3.e-17,3.e-17,-1.e-15,1.e-15}) {
        auto boundary=advance(p,threeCracks,{{residual,.0003,.0004,0,0,0}});
        check(boundary.openCracks==3,"closure insensitive to accumulated zero-strain roundoff");
    }
    for(int direction=0;direction<3;++direction) {
        auto e=threeCracks.strain;e[direction]=-.0001;
        auto closed=advance(p,threeCracks,e);
        check(closed.openCracks==2,"one of three independent planes closes");
        check(closed.stress[direction]<0.,"closed plane can carry compression");
        for(int i=0;i<3;++i) if(i!=direction)near(closed.stress[i],0,1e-12,"remaining open planes release stress");
    }
    t=advance(p,t,{{-.0001,-.0001,-.0001,0,0,0}});
    check(t.openCracks<3,"crack closure");
    // A pure shear increment in an already cracked, positive-strain frame.
    auto sh=advance(p,crack,{{.0002,0,0,1e-6,0,0}});
    near(sh.stress[3],p.E/(2*(1+p.nu))*p.shearRetention*1e-6,2e-10,"retained crack shear");
    auto history=advance(p,{},{{0,0,-.0005,0,0,0}});
    auto unload=advance(p,history,{{0,0,-.0004,0,0,0}});
    check(unload.maximumLoading>=history.maximumLoading,"monotone loading history");
    Matrix6 tangent;Vector6 target{{0,0,-.000603,0,0,0}};
    auto trial=advance(p,history,target,&tangent);
    for(int j=0;j<6;++j) {
        auto plus=target,minus=target;const double h=1e-9;plus[j]+=h;minus[j]-=h;
        auto a=advance(p,history,plus),b=advance(p,history,minus);
        for(int i=0;i<6;++i) {
            const double fd=(a.stress[i]-b.stress[i])/(2*h);
            if(std::abs(tangent[i][j]-fd)>2e-5*std::max({1.,std::abs(tangent[i][j]),std::abs(fd)}))
                std::cerr<<"component "<<i<<","<<j<<" plus "<<a.stress[i]<<" minus "<<b.stress[i]<<'\n';
            near(tangent[i][j],fd,2e-5,"nonlinear algorithmic tangent");
        }
    }
    auto repeated=advance(p,history,target);
    for(int i=0;i<6;++i)near(repeated.stress[i],trial.stress[i],1e-13,"deterministic trial replay");
    Vector6 invalid{};invalid[0]=NAN;State untouched=trial;
    check(!integrate(p,history,invalid,untouched,nullptr,error),"reject nonfinite strain");
    check(untouched.stress==trial.stress,"failed update is atomic");
    invalid[0]=1;check(!integrate(p,history,invalid,untouched,nullptr,error),"bounded substep work");
    // nu=0 removes lateral strain/stress feedback, giving an independent
    // uniaxial integration oracle from the printed rational stress formula.
    auto uniaxial=p;uniaxial.nu=0;uniaxial.maxStrainIncrement=1e-6;
    State uc;
    for(int i=1;i<=80;++i) {
        double e=-.00002*i;uc=advance(uniaxial,uc,{{0,0,e,0,0,0}});
        near(uc.stress[2],compressionStress(uniaxial,s,e),7e-4,"uniaxial stress integral");
    }
    for(int i=81;i<=310;++i)uc=advance(uniaxial,uc,{{0,0,-.00002*i,0,0,0}});
    check(uc.crushed,"irreversible compression crushing");
    near(uc.stress[2],0,1e-12,"ultimate stress release");
    // Constitutive objectivity: rotate the complete loading path, including
    // fixed crack histories, into arbitrary 3D orientations with negative
    // direction cosines. Compare stresses, never individual eigenvector signs.
    for(int rotation=0;rotation<20;++rotation) {
        Vector6 random;for(double& x:random)x=distribution(generator);
        const Frame frame=principalStress(random).frame;
        State local,global;
        for(int stage=0;stage<3;++stage) for(int k=1;k<=25;++k) {
            Vector6 e{{.0002,0.,0.,0.,0.,0.}};
            if(stage==0)e[0]=.0002*k/25.;
            if(stage>=1)e[1]=stage==1 ? .0003*k/25. : .0003;
            if(stage==2)e[2]=.0004*k/25.;
            local=advance(p,local,e);global=advance(p,global,strainToGlobal(e,frame));
            auto back=stressToLocal(global.stress,frame);
            for(int i=0;i<6;++i)near(local.stress[i],back[i],2e-7,"rotated multiaxial crack path");
            check(local.openCracks==global.openCracks,"rotated crack count");
        }
        for(int direction=0;direction<3;++direction) {
            auto e=local.strain;e[direction]=-.0001;
            auto a=advance(p,local,e),b=advance(p,global,strainToGlobal(e,frame));
            auto back=stressToLocal(b.stress,frame);
            for(int i=0;i<6;++i)near(a.stress[i],back[i],2e-7,"rotated independent crack closure");
            check(a.openCracks==b.openCracks,"rotated closure count");
        }
    }
    // Include exact ft crossings, which the preceding 25-step path misses.
    // A rotated eigenvalue differs by ulps at equality; releasing stress one
    // step early produced a 9%-of-fc discrepancy in the full-path review.
    // Test the same physical loading and closing history in random frames.
    auto threshold=p;threshold.maxStrainIncrement=1.;
    for(int rotation=0;rotation<20;++rotation) {
        Vector6 random;for(double& x:random)x=distribution(generator);
        const Frame frame=principalStress(random).frame;
        State local,global;
        for(int stage=0;stage<3;++stage) for(int k=1;k<=30;++k) {
            Vector6 e{{0.,0.,0.,0.,0.,0.}};
            e[0]=stage==0 ? k*1.e-5 : stage==1 ? (30-k)*1.e-5 : -k*1.e-5;
            local=advance(threshold,local,e);
            global=advance(threshold,global,strainToGlobal(e,frame));
            const auto back=stressToLocal(global.stress,frame);
            for(int i=0;i<6;++i)near(local.stress[i],back[i],2.e-7,"exact-threshold rotated stress");
            check(local.openCracks==global.openCracks,"exact-threshold rotated history");
        }
    }
    // Integration refinement must converge toward the independent rational
    // curve, rather than simply reproducing the old increment size.
    double previousError=1.;
    for(double h:{2e-5,1e-5,5e-6,2.5e-6}) {
        auto refined=uniaxial;refined.maxStrainIncrement=h;
        auto t=advance(refined,{},{{0,0,-.0016,0,0,0}});
        const double difference=std::abs(t.stress[2]-compressionStress(refined,s,-.0016));
        check(difference<previousError,"uniaxial integration refinement");previousError=difference;
    }
    // Update::Hybrid origin (see the 2D test): overshoot/return iterates must
    // not ratchet the compressive stress; monotonic iterates equal in-place.
    for(int kind=0;kind<2;++kind) {
        State mono,inPlace,hybrid,hybridMono;double drift=0.,hybridError=0.;
        Vector6 previous{};
        for(int k=1;k<=180;++k) {
            Vector6 e{};e[2]=-1.e-5*k;if(kind==1) {e[0]=-.8e-5*k;e[3]=2.e-6*k;}
            mono=advance(p,mono,e);
            const State committedHybrid=hybrid,committedMono=hybridMono;
            for(double f:{1.3,.9,1.05,1.}) {
                Vector6 iterate;for(int i=0;i<6;++i) iterate[i]=previous[i]+f*(e[i]-previous[i]);
                inPlace=advance(p,inPlace,iterate);
                hybrid=advance(p,hybridOrigin(committedHybrid,hybrid),iterate);
            }
            hybridMono=advance(p,hybridOrigin(committedMono,hybridMono),e);
            for(int i=0;i<6;++i) {
                drift=std::max(drift,std::abs(inPlace.stress[i]-mono.stress[i]));
                if(!mono.crushed) hybridError=std::max(hybridError,std::abs(hybrid.stress[i]-mono.stress[i]));
                check(hybridMono.stress[i]==mono.stress[i],"hybrid equals in-place for monotonic iterates");
            }
            previous=e;
        }
        check(drift>.05*std::abs(p.fc),"in-place origin shows iterate-path drift");
        near(hybridError/std::abs(p.fc),0.,.01,"hybrid origin follows the monotonic path");
    }
    std::cout<<"PASS "<<checks<<" kernel checks\n";
}
}
int main(int argc,char** argv) {
    try {
        if(argc>1 && std::string(argv[1])=="path") {
            Parameters p;p.maxStrainIncrement=1.;State c;Vector6 e;
            std::cout<<std::setprecision(17);
            while(std::cin>>e[0]>>e[1]>>e[2]>>e[3]>>e[4]>>e[5]) {
                c=advance(p,c,e);for(double x:c.stress)std::cout<<x<<' ';
                std::cout<<c.openCracks<<' '<<c.crushed<<'\n';
            }
        } else tests();
    } catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
