#include "Concrete2D.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <string>
#include <vector>
using namespace adina;
namespace ps=adina::plane;
namespace {
int checks=0;
void require(bool b,const char* message) {++checks;if(!b) {std::cerr<<message<<'\n';std::exit(1);}}
void near(double a,double b,double tol,const char* message) {
    if(std::abs(a-b)>tol*std::max({1.,std::abs(a),std::abs(b)}))
        std::cerr<<std::setprecision(17)<<"actual="<<a<<" expected="<<b<<" ";
    require(std::abs(a-b)<=tol*std::max({1.,std::abs(a),std::abs(b)}),message);
}
ps::State update(const Parameters& p,const ps::State& old,const ps::Vector3& e,ps::Matrix3* c=nullptr) {
    ps::State next;std::string error;
    require(ps::integrate(p,old,e,next,c,error),error.c_str());return next;
}
std::vector<ps::Vector3> path(int kind) {
    std::vector<ps::Vector3> v;
    for(int i=1;i<=240;++i) {
        if(kind==0) v.push_back({{0.,-i*1.e-5,0.}});
        if(kind==1) v.push_back({{-i*8.e-6,-i*1.e-5,0.}});
        if(kind==2) v.push_back({{i*2.e-6,0.,0.}});
        if(kind==3) v.push_back({{0.,0.,i*2.e-6}});
        if(kind==4 && i<=120) v.push_back({{0.,-i*1.e-5,0.}});
    }
    if(kind==4) for(int i=1;i<=120;++i) v.push_back({{0.,-.0012,i*2.e-6}});
    if(kind==2) {
        for(int i=1;i<=120;++i) v.push_back({{.00048,i*2.e-6,0.}});
        for(int i=1;i<=240;++i) v.push_back({{.00048-i*3.e-6,.00024,0.}});
        for(int i=1;i<=200;++i) v.push_back({{-.00024,.00024-i*3.e-6,0.}});
        for(int i=1;i<=180;++i) v.push_back({{-.00024+i*4.e-6,-.00036,0.}});
    }
    return v;
}
}
int main(int argc,char** argv) {
    Parameters p;
    p.tangent=Parameters::Tangent::Algorithmic;
    if(argc>2) {
        std::ifstream stream(argv[2]);
        stream>>p.E>>p.nu>>p.ft>>p.fc>>p.epsc>>p.fu>>p.epsu
              >>p.normalRetention>>p.shearRetention>>p.beta>>p.gamma>>p.kappa>>p.alpha;
        for(auto* array:{&p.sp1,&p.sp31,&p.sp32,&p.sp33}) for(double& value:*array) stream>>value;
        require(bool(stream),"read specimen constitutive parameters");
        std::string error;require(p.validate(error),error.c_str());
    }
    if(argc>1) {
        ps::State state;ps::Vector3 e;
        std::cout<<std::setprecision(17);
        while(std::cin>>e[0]>>e[1]>>e[2]) {
            state=update(p,state,e);
            for(double x:state.stress) std::cout<<x<<' ';
            std::cout<<state.open[0]<<' '<<state.open[1]<<' '<<state.crushed<<'\n';
        }
        return 0;
    }
    // Default: no subdivision, one update per increment as one CDMOD call.
    {
        Parameters d;std::string error;require(d.maxStrainIncrement==0. && d.validate(error),"default single update");
        require(d.update==Parameters::Update::Committed,"default update origin is ADINA's (committed)");
        ps::State a=update(d,{},{{0.,-.0012,0.}}),b;
        Parameters one=d;one.maxStrainIncrement=1.;b=update(one,{},{{0.,-.0012,0.}});
        for(int i=0;i<3;++i) require(a.stress[i]==b.stress[i],"maxStrainIncrement 0 equals one substep");
    }
    // The checks below were written for subdivided increments (2e-5), e.g.
    // closure part-way through a large increment; keep that control explicit.
    p.maxStrainIncrement=2.e-5;
    // Independent analytical small-strain plane stress and work conjugacy.
    for(int j=0;j<80;++j) {
        const double a=j*.037;
        ps::Vector3 e{{-1.e-7,-2.e-7,3.e-8}},s{{-7.,-13.,2.}};
        const auto el=ps::strainToLocal(e,a),sl=ps::stressToLocal(s,a);
        const auto eg=ps::strainToGlobal(el,a),sg=ps::stressToGlobal(sl,a);
        double w=0.,wl=0.;
        for(int i=0;i<3;++i) {near(e[i],eg[i],1.e-14,"strain round trip");
            near(s[i],sg[i],1.e-14,"stress round trip");w+=e[i]*s[i];wl+=el[i]*sl[i];}
        near(w,wl,1.e-14,"work conjugacy");
        const auto ep=ps::stressToLocal(s,ps::principalAngle(s));
        near(ep[2],0.,1.e-13,"principal shear zero");require(ep[0]>=ep[1],"principal ordering");
    }
    ps::Matrix3 c;auto elastic=update(p,{},{{-1.e-8,-2.e-8,3.e-9}},&c);
    const auto ce=ps::elasticTangent(p.E,p.nu);
    for(int i=0;i<3;++i) {
        double expected=0.;for(int j=0;j<3;++j) expected+=ce[i][j]*elastic.strain[j];
        near(elastic.stress[i],expected,1.e-13,"analytical plane stress");
        for(int j=0;j<3;++j) near(c[i][j],ce[i][j],1.e-7,"elastic tangent");
    }
    // nu=0 uniaxial compression isolates the manual's rational curve from
    // Poisson constraint effects. Check integral accuracy and refinement.
    Parameters zero=p;zero.nu=0.;
    // A central difference crossing the brittle stress jump used to produce
    // |Cxx| of order ft/h. At initiation require the elastic branch from the
    // uncracked side and the released branch from the cracked side. Check
    // their independent analytical derivatives without smoothing the stress.
    const double onset=zero.ft/zero.E;
    for(double offset:{-1.e-13,1.e-13}) {
        const auto at=update(zero,{},{{onset+offset,0.,0.}},&c);
        near(c[0][0],offset<0 ? zero.E : 0.,2.e-6,"one-sided crack initiation tangent");
        if(offset<0) require(!at.hasCrackFrame,"uncracked side of onset");
        else {require(at.hasCrackFrame,"cracked side of onset");near(at.stress[0],0.,1.e-14,"brittle stress release retained");}
    }
    // CDMOD intentionally uses E until |sigma| exceeds 0.001*|fc|.
    // Account analytically for that short elastic predictor, then verify
    // the integral of the manual's rational tangent on the remaining path.
    constexpr double seed=-1.001e-6;
    ps::State uniaxial=update(zero,{},{{0.,seed,0.}});
    const double offset=zero.E*seed-compressionStress(zero,{zero.fc,zero.epsc,zero.ft},seed);
    for(int i=1;i<=160;++i) {
        uniaxial=update(zero,uniaxial,{{0.,-i*1.e-5,0.}});
        near(uniaxial.stress[1],compressionStress(zero,{zero.fc,zero.epsc,zero.ft},-i*1.e-5)+offset,
             2.e-6,"rational compression curve");
        near(uniaxial.stress[0],0.,1.e-13,"zero Poisson transverse stress");
    }
    // Objectivity for monotonic compression, biaxial compression, tension,
    // shear and sequential two-plane closure/reopening. Identical tensor
    // increment norm produces identical substep count after rotation.
    for(int kind=0;kind<5;++kind) {
        const auto history=path(kind);std::vector<ps::State> baseline;
        ps::State state;
        for(const auto& e:history) {state=update(p,state,e);baseline.push_back(state);}
        if(kind==4) require(baseline[120].stress[2]>1.e-4,"orthotropic branch must respond to added shear");
        for(int rot=1;rot<=16;++rot) {
            const double a=rot*.093;state={};
            for(std::size_t k=0;k<history.size();++k) {
                state=update(p,state,ps::strainToGlobal(history[k],a));
                const auto expected=ps::stressToGlobal(baseline[k].stress,a);
                for(int i=0;i<3;++i) {
                    if(std::abs(state.stress[i]-expected[i])>2.e-8*std::max(1.,std::abs(expected[i])))
                        std::cerr<<"kind="<<kind<<" rot="<<rot<<" step="<<k+1<<" component="<<i<<" ";
                    near(state.stress[i],expected[i],2.e-8,"rotational objectivity");
                }
                require(state.crushed==baseline[k].crushed,"rotation crushing state");
                require(int(state.open[0])+int(state.open[1])==
                        int(baseline[k].open[0])+int(baseline[k].open[1]),"rotation crack count");
            }
        }
    }
    ps::State two;
    for(const auto& e:std::vector<ps::Vector3>{{{.0002,0.,0.}},{{.0002,.0002,0.}}}) two=update(p,two,e);
    require(two.open[0] && two.open[1],"two perpendicular cracks");
    for(double residual: {-3.e-17, 3.e-17, -1.e-15, 1.e-15}) {
        const auto boundary=update(p,two,{{residual,.0002,0.}});
        require(boundary.open[0] && boundary.open[1],"closure insensitive to accumulated zero-strain roundoff");
    }
    auto one=update(p,two,{{-.0001,.0002,0.}});
    require(!one.open[0] && one.open[1],"only first crack closes");
    require(one.stress[0]<-1.,"closed first plane carries compression");
    near(one.stress[1],0.,1.e-13,"second plane remains stress-free");
    auto none=update(p,one,{{-.0001,-.0001,0.}});
    require(!none.open[0] && !none.open[1],"both cracks close");
    auto reopened=update(p,none,{{.0004,-.0001,0.}});
    require(reopened.open[0] && !reopened.open[1],"first crack reopens independently");
    near(reopened.stress[0],0.,1.e-13,"reopened first plane releases tension");
    require(reopened.stress[1]<0.,"closed second plane carries compression");
    // Numerical Jacobian checked by a different perturbation size on a smooth
    // compressed state, and on an open-crack state with retained shear.
    for(auto old:{uniaxial,one}) {
        auto e=old.strain;e[0]-=1.e-6;e[1]-=2.e-6;e[2]+=3.e-7;
        update(p,old,e,&c);
        for(int j=0;j<3;++j) {
            const double h=2.e-9;auto a=e,b=e;a[j]+=h;b[j]-=h;
            const auto sa=update(p,old,a),sb=update(p,old,b);
            for(int i=0;i<3;++i) near(c[i][j],(sa.stress[i]-sb.stress[i])/(2.*h),3.e-5,"algorithmic tangent");
        }
    }
    Parameters fine=p;fine.maxStrainIncrement=1.e-5;
    Parameters finer=p;finer.maxStrainIncrement=5.e-6;
    auto a=update(p,{},{{0.,-.0012,0.}}),b=update(fine,{},{{0.,-.0012,0.}}),
         d=update(finer,{},{{0.,-.0012,0.}});
    require(std::abs(b.stress[1]-d.stress[1])<std::abs(a.stress[1]-d.stress[1]),"substep refinement");
    // Exact immutability on failed input, including tangent and committed data.
    ps::State out=one;auto saved=out.stress;ps::Matrix3 sentinel{};sentinel[0][0]=42.;std::string error;
    require(!ps::integrate(p,one,{{NAN,0.,0.}},out,&sentinel,error),"reject NaN");
    require(out.stress==saved && sentinel[0][0]==42.,"atomic failure");
    require(!ps::integrate(p,one,{{1.e10,0.,0.}},out,&sentinel,error),"reject substep overflow");
    Parameters bad=p;bad.nu=.5;require(!ps::integrate(bad,one,{},out,nullptr,error),"reject singular Poisson");
    bad=p;bad.epsu=bad.epsc;require(!ps::integrate(bad,one,{},out,nullptr,error),"reject singular curve");
    // A zero residual normal stiffness is supported, including two open cracks.
    bad=p;bad.normalRetention=0.;bad.shearRetention=0.;
    const auto released=update(bad,{},{{.0002,.0002,.0001}});
    for(double x:released.stress) near(x,0.,1.e-14,"zero-retention complete release");
    // Ultimate compression, unloading after crushing, and nonproportional
    // cyclic paths over the parameter ranges exercised in structural usage.
    auto crushed=update(p,{},{{0.,-.007,0.}});
    require(crushed.crushed,"irreversible crushing indicator");
    for(double x:crushed.stress) near(x,0.,1.e-13,"ultimate strain release");
    auto unload=update(p,crushed,{{0.,-.003,0.}});
    require(unload.crushed,"unloading cannot erase crushing");
    for(double nu:{0.,.2,.35}) for(double retention:{0.,.0001,.02}) {
        Parameters settings=p;settings.nu=nu;settings.normalRetention=retention;
        ps::State state;
        for(int i=1;i<=800;++i) {
            const double t=i*.025;
            ps::Vector3 e{{.00035*std::sin(t),-.0008*(1.-std::cos(t)),.0003*std::sin(t*.73)}};
            state=update(settings,state,e);
            require(ps::validateState(state,error),"finite cyclic history");
            for(double x:state.stress) require(std::isfinite(x),"finite cyclic stress");
        }
    }
    // Update::Hybrid origin. Each step is reached through Newton-like iterates
    // that overshoot and return. The in-place (ADINA) origin ratchets the
    // compressive stress above the curve; the hybrid origin must stay on the
    // monotonic path, and equal the in-place origin for monotonic iterates.
    for(int kind:{0,1,4}) {
        const auto strains=path(kind);
        ps::State mono,inPlace,hybrid,hybridMono;double drift=0.,hybridError=0.;
        ps::Vector3 previous{};
        for(const auto& e:strains) {
            mono=update(p,mono,e);
            const ps::State committedInPlace=inPlace,committedHybrid=hybrid,committedMono=hybridMono;
            for(double f:{1.3,.9,1.05,1.}) {
                ps::Vector3 iterate;for(int i=0;i<3;++i) iterate[i]=previous[i]+f*(e[i]-previous[i]);
                inPlace=update(p,inPlace,iterate);
                hybrid=update(p,ps::hybridOrigin(committedHybrid,hybrid),iterate);
            }
            (void)committedInPlace;
            hybridMono=update(p,ps::hybridOrigin(committedMono,hybridMono),e);
            for(int i=0;i<3;++i) {
                drift=std::max(drift,std::abs(inPlace.stress[i]-mono.stress[i]));
                if(!mono.crushed) hybridError=std::max(hybridError,std::abs(hybrid.stress[i]-mono.stress[i]));
                require(hybridMono.stress[i]==mono.stress[i],"hybrid equals in-place for monotonic iterates");
            }
            previous=e;
        }
        require(drift>.05*std::abs(p.fc),"in-place origin shows iterate-path drift");
        near(hybridError/std::abs(p.fc),0.,.01,"hybrid origin follows the monotonic path");
    }
    std::cout<<"PASS "<<checks<<" native 2D kernel checks\n";
}
