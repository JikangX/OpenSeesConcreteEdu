#include "AdinaConcrete3D.h"
#include <classTags.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <cstring>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
// A diverging global iterate is rejected at every integration point. The
// analysis receives each rejection; print the first ones and then a count,
// so that millions of identical lines do not dominate the run time.
void reportRejectedTrial(const char* name,int tag,const std::string& error) {
    static long count=0;
    ++count;
    if(count<=20 || count%100000==0)
        opserr << name << " " << tag << ": " << error.c_str() << " (rejected trials: " << count << ")" << endln;
}
}


AdinaConcrete3D::AdinaConcrete3D():NDMaterial(0,ND_TAG_AdinaConcrete3D) {revertToStart();}
AdinaConcrete3D::AdinaConcrete3D(int tag,const adina::Parameters& p)
    :NDMaterial(tag,ND_TAG_AdinaConcrete3D),parameters(p) {
    std::string error;
    if(!parameters.validate(error)) throw std::invalid_argument(error);
    revertToStart();
}
void AdinaConcrete3D::synchronize(const adina::Matrix6& c) {
    for(int i=0;i<6;++i) {
        strain(i)=trial.strain[i];stress(i)=trial.stress[i];
        for(int j=0;j<6;++j) tangent(i,j)=c[i][j];
    }
}
int AdinaConcrete3D::setTrialStrain(const Vector& v) {
    // The update starts from the previous iteration (ADINA, default) or from
    // the last commit (-update committed); see adina::Parameters::Update.
    // Do not assign output buffers until the full stress and tangent
    // calculation has succeeded (strong failure guarantee).
    if(v.Size()!=6) {opserr << "AdinaConcrete3D requires six engineering strain components\n";return -1;}
    adina::Vector6 eps;for(int i=0;i<6;++i) eps[i]=v(i);
    // A repeated request for the current trial strain must not apply a
    // second (zero) increment to the in-place iterative history.
    if(parameters.update!=adina::Parameters::Update::Committed && eps==trial.strain) return 0;
    // Hybrid: previous iteration with the converged loading reference.
    const adina::State origin=parameters.update==adina::Parameters::Update::Iterative ? trial :
        parameters.update==adina::Parameters::Update::Hybrid ? adina::hybridOrigin(committed,trial) : committed;
    adina::State next;adina::Matrix6 c;std::string error;
    if(!adina::integrate(parameters,origin,eps,next,&c,error)) {
        reportRejectedTrial("AdinaConcrete3D",getTag(),error);return -1;
    }
    trial=next;synchronize(c);return 0;
}
int AdinaConcrete3D::setTrialStrainIncr(const Vector& v) {
    // OpenSees defines this increment relative to the current trial strain.
    if(v.Size()!=6) return -1;
    Vector eps(strain);eps+=v;return setTrialStrain(eps);
}
int AdinaConcrete3D::commitState() {
    // Save the tangent as well as history: reverting must restore the exact
    // returned trial buffers, without silently doing a second material update.
    committed=trial;
    for(int i=0;i<6;++i) for(int j=0;j<6;++j) committedTangent[i][j]=tangent(i,j);
    return 0;
}
int AdinaConcrete3D::revertToLastCommit() {trial=committed;synchronize(committedTangent);return 0;}
int AdinaConcrete3D::revertToStart() {
    committed=adina::State{};trial=committed;
    committedTangent=adina::elasticTangent(parameters.E,parameters.nu);
    for(int i=0;i<6;++i) for(int j=0;j<6;++j) initialTangent(i,j)=committedTangent[i][j];
    synchronize(committedTangent);return 0;
}
NDMaterial* AdinaConcrete3D::getCopy() {
    auto* copy=new AdinaConcrete3D(*this);
    // A material point cloned from a checkpoint owns a new Channel identity.
    // Keeping the source dbTag would make sibling integration points overwrite
    // one another when an element sends their histories to a database.
    copy->setDbTag(0);
    return copy;
}
NDMaterial* AdinaConcrete3D::getCopy(const char* type) {
    return type && (!std::strcmp(type,"ThreeDimensional") || !std::strcmp(type,"3D")) ? getCopy() : nullptr;
}

namespace {
// Versioned flat serialization; no raw struct bytes, pointers, padding or local
// scratch data travel through the Channel. BOTH states and BOTH tangents are
// preserved, so restart between a trial and its commit is deterministic too.
// Version 2 appends optional fracture parameters and physical-plane history;
// a v1 checkpoint remains readable as the default original ADINA mode.
// Version 3 appends the tangent, update and crushed-unloading options and
// each state's crushed-strain history. v1/v2 messages come
// from the replay-from-commit, finite-difference-tangent formulation and are
// restored with exactly those options (crushed unloading: original).
constexpr int wireVersion=3,wireSize=198,v2WireSize=193,legacyWireSize=173;
void packState(Vector& v,int& k,const adina::State& s) {
    for(double x:s.strain) v(k++)=x;for(double x:s.stress) v(k++)=x;
    v(k++)=s.maximumLoading;v(k++)=s.hasCrackFrame;v(k++)=s.crushed;v(k++)=s.openCracks;
    for(const auto& row:s.frame) for(double x:row) v(k++)=x;
    for(double x:s.crackStrain) v(k++)=x;
    v(k++)=s.crushingStress;v(k++)=s.crushingStrain;
}
bool unpackState(const Vector& v,int& k,adina::State& s) {
    for(double& x:s.strain) x=v(k++);for(double& x:s.stress) x=v(k++);
    // Validate wire enums BEFORE converting floating-point storage to int.
    if((v(k+1)!=0. && v(k+1)!=1.) || (v(k+2)!=0. && v(k+2)!=1.) ||
       v(k+3)<0. || v(k+3)>3. || std::floor(v(k+3))!=v(k+3)) return false;
    s.maximumLoading=v(k++);s.hasCrackFrame=v(k++)!=0.;s.crushed=v(k++)!=0.;s.openCracks=int(v(k++));
    for(auto& row:s.frame) for(double& x:row) x=v(k++);
    for(double& x:s.crackStrain) x=v(k++);
    s.crushingStress=v(k++);s.crushingStrain=v(k++);
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) {
        double dot=0.;for(int a=0;a<3;++a) dot+=s.frame[a][i]*s.frame[a][j];
        if(std::abs(dot-(i==j ? 1. : 0.))>1.e-10) return false;
    }
    return !(s.crushed && (s.crushingStress>=0. || s.crushingStrain>=0.));
}
}
int AdinaConcrete3D::sendSelf(int commitTag,Channel& channel) {
    Vector v(wireSize);int k=0;v(k++)=wireVersion;v(k++)=getTag();
    const double scalars[]={parameters.E,parameters.nu,parameters.ft,parameters.fc,parameters.epsc,
        parameters.fu,parameters.epsu,parameters.beta,parameters.gamma,parameters.kappa,parameters.alpha,
        parameters.normalRetention,parameters.shearRetention,parameters.rho,parameters.maxStrainIncrement};
    for(double x:scalars) v(k++)=x;
    for(const auto* a:{&parameters.sp1,&parameters.sp31,&parameters.sp32,&parameters.sp33})
        for(double x:*a) v(k++)=x;
    packState(v,k,committed);packState(v,k,trial);
    for(const auto& row:committedTangent) for(double x:row) v(k++)=x;
    for(int i=0;i<6;++i) for(int j=0;j<6;++j) v(k++)=tangent(i,j);
    v(k++)=parameters.fractureEnergy;v(k++)=parameters.characteristicLength;
    for(const auto* state:{&committed,&trial})for(const auto& h:state->tension) {
        v(k++)=h.kappa;v(k++)=h.strength;v(k++)=h.offset;
    }
    v(k++)=parameters.tangent==adina::Parameters::Tangent::Secant ? 0. : 1.;
    v(k++)=double(int(parameters.update)); // 0 iterative, 1 committed, 2 hybrid
    v(k++)=parameters.crushedUnloading==adina::Parameters::CrushedUnloading::Original ? 0. : 1.;
    v(k++)=committed.crushedStrain;v(k++)=trial.crushedStrain;
    return k==wireSize ? channel.sendVector(getDbTag(),commitTag,v) : -1;
}
int AdinaConcrete3D::recvSelf(int commitTag,Channel& channel,FEM_ObjectBroker&) {
    Vector v(wireSize);
    // FileDatastore keys vectors by size: older versions need their length.
    if(channel.recvVector(getDbTag(),commitTag,v)<0) {
        v.resize(v2WireSize);
        if(channel.recvVector(getDbTag(),commitTag,v)<0) {
            v.resize(legacyWireSize);if(channel.recvVector(getDbTag(),commitTag,v)<0)return -1;
        }
    }
    const bool legacy=v.Size()==legacyWireSize && v(0)==1.;
    const bool v2=v.Size()==v2WireSize && v(0)==2.;
    if(!legacy && !v2 && (v.Size()!=wireSize || v(0)!=wireVersion))return -1;
    for(int i=0;i<v.Size();++i) if(!std::isfinite(v(i))) return -1;
    if(v(1)<std::numeric_limits<int>::min() || v(1)>std::numeric_limits<int>::max() ||
       std::floor(v(1))!=v(1)) return -1;
    int k=2;adina::Parameters p;
    double* scalars[]={&p.E,&p.nu,&p.ft,&p.fc,&p.epsc,&p.fu,&p.epsu,&p.beta,&p.gamma,&p.kappa,&p.alpha,
        &p.normalRetention,&p.shearRetention,&p.rho,&p.maxStrainIncrement};
    for(double* x:scalars) *x=v(k++);
    for(auto* a:{&p.sp1,&p.sp31,&p.sp32,&p.sp33}) for(double& x:*a) x=v(k++);
    std::string error;if(!p.validate(error)) return -1;
    adina::State cs,ts;adina::Matrix6 cc,tc;
    if(!unpackState(v,k,cs) || !unpackState(v,k,ts)) return -1;
    for(auto& row:cc) for(double& x:row) x=v(k++);
    for(auto& row:tc) for(double& x:row) x=v(k++);
    if(!legacy) {
        p.fractureEnergy=v(k++);p.characteristicLength=v(k++);
        for(auto* state:{&cs,&ts})for(auto& h:state->tension) {
            h.kappa=v(k++);h.strength=v(k++);h.offset=v(k++);
            if(!adina::fracture::valid(h))return -1;
        }
        if(!p.validate(error))return -1;
    }
    p.tangent=adina::Parameters::Tangent::Algorithmic;p.update=adina::Parameters::Update::Committed;
    if(!legacy && !v2) {
        if((v(k)!=0. && v(k)!=1.) || (v(k+1)!=0. && v(k+1)!=1. && v(k+1)!=2.) || (v(k+2)!=0. && v(k+2)!=1.))
            return -1;
        p.tangent=v(k++)==0. ? adina::Parameters::Tangent::Secant : adina::Parameters::Tangent::Algorithmic;
        p.update=static_cast<adina::Parameters::Update>(int(v(k++)));
        p.crushedUnloading=v(k++)==0. ? adina::Parameters::CrushedUnloading::Original
                                      : adina::Parameters::CrushedUnloading::Elastic;
        cs.crushedStrain=v(k++);ts.crushedStrain=v(k++);
    }
    if(k!=v.Size()) return -1;
    for(const auto* state:{&cs,&ts})for(const auto& h:state->tension)
        if(!adina::fracture::compatible(h,p.E,p.fractureEnergy,p.characteristicLength,state->hasCrackFrame))return -1;
    for(const auto* state:{&cs,&ts}) {
        // Crushing supersedes the cohesive law. The saved tensile maxima
        // remain historical values, not bounds on later crushed strains.
        if(state->crushed) continue;
        const auto e=adina::strainToLocal(state->strain,state->frame);
        for(int i=0;i<3;++i)if(state->tension[i].strength>0.) {
            const double eta=e[i]-state->tension[i].offset;
            const double tol=64.*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(e[i]));
            if((i<state->openCracks && eta<-tol) || eta>state->tension[i].kappa+tol)return -1;
        }
    }
    parameters=p;setTag(int(v(1)));revertToStart();
    committed=cs;trial=ts;committedTangent=cc;synchronize(tc);return 0;
}
void AdinaConcrete3D::Print(OPS_Stream& s,int) {
    s << "AdinaConcrete3D tag=" << getTag() << " E=" << parameters.E << " nu=" << parameters.nu
      << " ft=" << parameters.ft << " fc=" << parameters.fc << " cracks=" << trial.openCracks
      << " fractureEnergy=" << parameters.fractureEnergy << " Lch=" << parameters.characteristicLength
      << " crushed=" << int(trial.crushed) << endln;
}
Response* AdinaConcrete3D::setResponse(const char** argv,int argc,OPS_Stream& s) {
    // Normal stress/strain/tangent responses use NDMaterial's implementation.
    // Extra responses expose physical history rather than legacy DCA sentinels.
    if(argc<=0) return nullptr;
    if(!std::strcmp(argv[0],"fractureState")) {
        s.tag("NdMaterialOutput");s.attr("matType",getClassType());s.attr("matTag",getTag());
        const char* names[]={"Gf","Lch","kappa1","strength1","offset1","dissipation1",
                            "kappa2","strength2","offset2","dissipation2",
                            "kappa3","strength3","offset3","dissipation3"};
        for(const char* name:names)s.tag("ResponseType",name);
        s.endTag();return new MaterialResponse(this,102,Vector(14));
    }
    int id=0;
    if(!std::strcmp(argv[0],"crackState") || !std::strcmp(argv[0],"state")) id=100;
    if(!std::strcmp(argv[0],"crackFrame")) id=101;
    if(id) {
        s.tag("NdMaterialOutput");s.attr("matType",getClassType());s.attr("matTag",getTag());
        const char* stateNames[]={"openCracks","crushed","maximumLoading","hasCrackFrame",
            "crackStrain1","crackStrain2","crackStrain3","crushingStress","crushingStrain"};
        const char* frameNames[]={"q11","q12","q13","q21","q22","q23","q31","q32","q33"};
        const char** names=id==100 ? stateNames : frameNames;
        for(int i=0;i<9;++i) s.tag("ResponseType",names[i]);
        s.endTag();return new MaterialResponse(this,id,Vector(9));
    }
    return NDMaterial::setResponse(argv,argc,s);
}
int AdinaConcrete3D::getResponse(int id,Information& info) {
    if(id==102) {
        Vector v(14);int k=0;v(k++)=parameters.fractureEnergy;v(k++)=parameters.characteristicLength;
        for(const auto& h:trial.tension) {
            v(k++)=h.kappa;v(k++)=h.strength;v(k++)=h.offset;
            v(k++)=adina::fracture::dissipation(parameters.E,parameters.fractureEnergy,parameters.characteristicLength,h);
        }
        return info.setVector(v);
    }
    if(id==100) {
        Vector v(9);v(0)=trial.openCracks;v(1)=trial.crushed;v(2)=trial.maximumLoading;
        v(3)=trial.hasCrackFrame;for(int i=0;i<3;++i) v(4+i)=trial.crackStrain[i];
        v(7)=trial.crushingStress;v(8)=trial.crushingStrain;return info.setVector(v);
    }
    if(id==101) {
        Vector v(9);int k=0;for(const auto& row:trial.frame) for(double x:row) v(k++)=x;
        return info.setVector(v);
    }
    return NDMaterial::getResponse(id,info);
}
