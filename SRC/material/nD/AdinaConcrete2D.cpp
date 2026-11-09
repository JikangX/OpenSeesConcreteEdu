#include "AdinaConcrete2D.h"
#include <classTags.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <cmath>
#include <cstring>
#ifdef ADINA_TRACE
#include <cstdio>
#include <cstdlib>
#endif
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

AdinaConcrete2D::AdinaConcrete2D():NDMaterial(0,ND_TAG_AdinaConcrete2D) {revertToStart();}
AdinaConcrete2D::AdinaConcrete2D(int tag,const adina::Parameters& p,bool shell,double G)
    :NDMaterial(tag,ND_TAG_AdinaConcrete2D),parameters(p),plate(shell),
     transverseG(G==-1. ? p.E/(2.*(1.+p.nu)) : G) {
    std::string error;
    if(!p.validate(error)) throw std::invalid_argument(error);
    if(!std::isfinite(transverseG) || transverseG<=0.)
        throw std::invalid_argument("AdinaConcrete2D requires a positive finite transverse G");
    resizeBuffers();revertToStart();
}
void AdinaConcrete2D::resizeBuffers() {
    const int n=getOrder();strain.resize(n);stress.resize(n);
    tangent.resize(n,n);initialTangent.resize(n,n);
}
void AdinaConcrete2D::initializeTangent() {
    const auto c=adina::plane::elasticTangent(parameters.E,parameters.nu);
    initialTangent.Zero();
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) initialTangent(i,j)=c[i][j];
    if(plate) initialTangent(3,3)=initialTangent(4,4)=transverseG;
}
void AdinaConcrete2D::synchronize(const adina::plane::Matrix3& c) {
    tangent.Zero();
    for(int i=0;i<3;++i) {
        strain(i)=trial.strain[i];stress(i)=trial.stress[i];
        for(int j=0;j<3;++j) tangent(i,j)=c[i][j];
    }
    if(plate) for(int i=0;i<2;++i) {
        strain(i+3)=trialShear[i];stress(i+3)=transverseG*trialShear[i];
        tangent(i+3,i+3)=transverseG;
    }
}
int AdinaConcrete2D::setTrialStrain(const Vector& v) {
    if(v.Size()!=getOrder()) return -1;
    for(int i=0;i<v.Size();++i) if(!std::isfinite(v(i))) return -1;
    // In-place iteration is incremental from the current trial: requesting
    // that exact same strain is a read, even immediately after commit/restore
    // or rollback. A zero-increment replay can otherwise reclassify a brittle
    // boundary and change the history merely because Domain::update() ran.
    // Committed replay retains its former Jacobian-recomputation convention.
    // Exact comparison avoids suppressing a real, small strain increment.
    bool identical=true,differentOrigin=false;
    for(int i=0;i<v.Size();++i) {
        identical=identical && v(i)==strain(i);
        const double origin=i<3 ? committed.strain[i] : committedShear[i-3];
        differentOrigin=differentOrigin || v(i)!=origin;
    }
    if(identical && (differentOrigin || parameters.update!=adina::Parameters::Update::Committed)) return 0;
    // Compute ALL prospective outputs first. Failed substeps or tangents must
    // not replace trial history, transverse strain, or any returned buffer.
    adina::plane::Vector3 e{{v(0),v(1),v(2)}};
    // Origin: previous iteration (ADINA's in-place WA update, default) or
    // last commit (-update committed); see adina::Parameters::Update.
    // Hybrid: previous iteration with the converged loading reference.
    const adina::plane::State origin=
        parameters.update==adina::Parameters::Update::Iterative ? trial :
        parameters.update==adina::Parameters::Update::Hybrid ? adina::plane::hybridOrigin(committed,trial) : committed;
    adina::plane::State next;adina::plane::Matrix3 c;std::string error;
    if(!adina::plane::integrate(parameters,origin,e,next,&c,error)) {
        reportRejectedTrial("AdinaConcrete2D",getTag(),error);return -1;
    }
    if(plate && (!std::isfinite(transverseG*v(3)) || !std::isfinite(transverseG*v(4)))) return -1;
#ifdef ADINA_TRACE
    // Diagnostic builds only: report state changes between successive trials
    // while the environment variable ADINA_TRACE is set.
    if(std::getenv("ADINA_TRACE") && (trial.open!=next.open || trial.crushed!=next.crushed ||
       trial.released!=next.released || trial.hasCrackFrame!=next.hasCrackFrame))
        std::fprintf(stderr,"TRACE %p open %d%d->%d%d crushed %d->%d released %d->%d e=(%.3e %.3e %.3e) s=(%.3g %.3g %.3g)\n",
            (void*)this,int(trial.open[0]),int(trial.open[1]),int(next.open[0]),int(next.open[1]),
            int(trial.crushed),int(next.crushed),int(trial.released),int(next.released),e[0],e[1],e[2],
            next.stress[0],next.stress[1],next.stress[2]);
#endif
    trial=next;if(plate) trialShear={{v(3),v(4)}};
    synchronize(c);return 0;
}
int AdinaConcrete2D::setTrialStrainIncr(const Vector& v) {
    if(v.Size()!=getOrder()) return -1;
    Vector total(strain);total+=v;return setTrialStrain(total);
}
int AdinaConcrete2D::commitState() {
    committed=trial;committedShear=trialShear;
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) committedTangent[i][j]=tangent(i,j);
    return 0;
}
int AdinaConcrete2D::revertToLastCommit() {
    trial=committed;trialShear=committedShear;synchronize(committedTangent);return 0;
}
int AdinaConcrete2D::revertToStart() {
    trial=committed=adina::plane::State{};trialShear=committedShear={{0.,0.}};
    committedTangent=adina::plane::elasticTangent(parameters.E,parameters.nu);
    initializeTangent();synchronize(committedTangent);return 0;
}
NDMaterial* AdinaConcrete2D::getCopy() {
    auto* copy=new AdinaConcrete2D(*this);
    // A clone has a new Channel identity even if its source was restored from
    // a database. Sharing a nonzero dbTag would overwrite a sibling layer.
    copy->setDbTag(0);return copy;
}
NDMaterial* AdinaConcrete2D::getCopy(const char* type) {
    if(!type) return nullptr;
    const bool shell=!std::strcmp(type,"PlateFiber");
    if(!shell && std::strcmp(type,"PlaneStress") && std::strcmp(type,"PlaneStress2D")) return nullptr;
    auto* copy=new AdinaConcrete2D(*this);
    copy->setDbTag(0);
    // Changing interface size preserves the membrane trial AND commit. In
    // particular, section construction cannot alias the prototype's buffers.
    adina::plane::Matrix3 c{};
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) c[i][j]=tangent(i,j);
    copy->plate=shell;copy->resizeBuffers();copy->initializeTangent();copy->synchronize(c);
    return copy;
}

namespace {
// Explicit versioned wire format: header(4), parameters(39), histories(32),
// transverse strains(4), two 3x3 tangents(18). Never serialize struct padding.
// Version 2 appends Gf/Lch and two planes' three-value histories for BOTH
// committed and trial points. The first 97 entries retain the v1 layout.
// Version 3 appends the tangent, update and crushed-unloading options and
// each state's crushed-strain history and release flag. v1/v2 messages come
// from the replay-from-commit, finite-difference-tangent formulation and are
// restored with exactly those options (crushed unloading: original).
constexpr int wireVersion=3,wireSize=118,v2WireSize=111,legacyWireSize=97;
void pack(Vector& v,int& k,const adina::plane::State& s) {
    for(double x:s.strain) v(k++)=x;
    for(double x:s.stress) v(k++)=x;
    v(k++)=s.maximumLoading;v(k++)=s.angle;v(k++)=s.hasCrackFrame;v(k++)=s.crushed;
    for(bool x:s.open) v(k++)=x;
    for(double x:s.crackStrain) v(k++)=x;
    v(k++)=s.crushingStress;v(k++)=s.crushingStrain;
}
bool unpack(const Vector& v,int& k,adina::plane::State& s) {
    for(double& x:s.strain) x=v(k++);
    for(double& x:s.stress) x=v(k++);
    s.maximumLoading=v(k++);s.angle=v(k++);
    for(int i=0;i<4;++i) if(v(k+i)!=0. && v(k+i)!=1.) return false;
    s.hasCrackFrame=v(k++)!=0.;s.crushed=v(k++)!=0.;
    for(auto& x:s.open) x=v(k++)!=0.;
    for(double& x:s.crackStrain) x=v(k++);
    s.crushingStress=v(k++);s.crushingStrain=v(k++);
    std::string error;return adina::plane::validateState(s,error);
}
}
int AdinaConcrete2D::sendSelf(int commitTag,Channel& channel) {
    Vector v(wireSize);int k=0;v(k++)=wireVersion;v(k++)=getTag();v(k++)=plate;v(k++)=transverseG;
    const double scalars[]={parameters.E,parameters.nu,parameters.ft,parameters.fc,parameters.epsc,
        parameters.fu,parameters.epsu,parameters.beta,parameters.gamma,parameters.kappa,parameters.alpha,
        parameters.normalRetention,parameters.shearRetention,parameters.rho,parameters.maxStrainIncrement};
    for(double x:scalars) v(k++)=x;
    for(const auto* a:{&parameters.sp1,&parameters.sp31,&parameters.sp32,&parameters.sp33})
        for(double x:*a) v(k++)=x;
    pack(v,k,committed);pack(v,k,trial);
    for(double x:committedShear) v(k++)=x;
    for(double x:trialShear) v(k++)=x;
    for(const auto& row:committedTangent) for(double x:row) v(k++)=x;
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) v(k++)=tangent(i,j);
    v(k++)=parameters.fractureEnergy;v(k++)=parameters.characteristicLength;
    for(const auto* state:{&committed,&trial}) for(const auto& h:state->tension) {
        v(k++)=h.kappa;v(k++)=h.strength;v(k++)=h.offset;
    }
    v(k++)=parameters.tangent==adina::Parameters::Tangent::Secant ? 0. : 1.;
    v(k++)=double(int(parameters.update)); // 0 iterative, 1 committed, 2 hybrid
    v(k++)=parameters.crushedUnloading==adina::Parameters::CrushedUnloading::Original ? 0. : 1.;
    v(k++)=committed.crushedStrain;v(k++)=trial.crushedStrain;
    v(k++)=committed.released;v(k++)=trial.released;
    return k==wireSize ? channel.sendVector(getDbTag(),commitTag,v) : -1;
}
int AdinaConcrete2D::recvSelf(int commitTag,Channel& channel,FEM_ObjectBroker&) {
    Vector v(wireSize);
    if(channel.recvVector(getDbTag(),commitTag,v)<0) {
        // FileDatastore keys vectors by size; older versions need their length.
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
       std::floor(v(1))!=v(1) || (v(2)!=0. && v(2)!=1.) || v(3)<=0.) return -1;
    int k=4;adina::Parameters p;
    double* scalars[]={&p.E,&p.nu,&p.ft,&p.fc,&p.epsc,&p.fu,&p.epsu,&p.beta,&p.gamma,&p.kappa,&p.alpha,
        &p.normalRetention,&p.shearRetention,&p.rho,&p.maxStrainIncrement};
    for(double* x:scalars) *x=v(k++);
    for(auto* a:{&p.sp1,&p.sp31,&p.sp32,&p.sp33}) for(double& x:*a) x=v(k++);
    std::string error;if(!p.validate(error)) return -1;
    adina::plane::State cs,ts;adina::plane::Matrix3 cc{},tc{};
    std::array<double,2> sc{},st{};
    if(!unpack(v,k,cs) || !unpack(v,k,ts)) return -1;
    for(double& x:sc) x=v(k++);
    for(double& x:st) x=v(k++);
    for(double x:sc) if(!std::isfinite(v(3)*x)) return -1;
    for(double x:st) if(!std::isfinite(v(3)*x)) return -1;
    for(auto& row:cc) for(double& x:row) x=v(k++);
    for(auto& row:tc) for(double& x:row) x=v(k++);
    if(!legacy) {
        p.fractureEnergy=v(k++);p.characteristicLength=v(k++);
        for(auto* state:{&cs,&ts}) for(auto& h:state->tension) {
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
        if((v(k)!=0. && v(k)!=1.) || (v(k+1)!=0. && v(k+1)!=1.)) return -1;
        cs.released=v(k++)!=0.;ts.released=v(k++)!=0.;
        if((cs.released && !cs.crushed) || (ts.released && !ts.crushed)) return -1;
    }
    if(k!=v.Size()) return -1;
    for(const auto* state:{&cs,&ts})for(const auto& h:state->tension)
        if(!adina::fracture::compatible(h,p.E,p.fractureEnergy,p.characteristicLength,state->hasCrackFrame))return -1;
    for(const auto* state:{&cs,&ts}) {
        // After crushing the tensile law is inactive and its per-plane
        // history is frozen. Subsequent strains need not lie below that old
        // tensile maximum; reject impossible active tensile histories only.
        if(state->crushed) continue;
        const auto e=adina::plane::strainToLocal(state->strain,state->angle);
        for(int i=0;i<2;++i)if(state->tension[i].strength>0.) {
            const double eta=e[i]-state->tension[i].offset;
            const double tol=64.*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(e[i]));
            if((state->open[i] && eta<-tol) || eta>state->tension[i].kappa+tol)return -1;
        }
    }
    // Publish only after validating the entire message: corrupted checkpoints
    // cannot partially overwrite a working material or change its formulation.
    parameters=p;plate=v(2)!=0.;transverseG=v(3);setTag(int(v(1)));resizeBuffers();initializeTangent();
    committed=cs;trial=ts;committedShear=sc;trialShear=st;committedTangent=cc;synchronize(tc);return 0;
}
void AdinaConcrete2D::Print(OPS_Stream& s,int) {
    s << "AdinaConcrete2D tag=" << getTag() << " type=" << getType() << " E=" << parameters.E
      << " nu=" << parameters.nu << " transverseG=" << transverseG
      << " fractureEnergy=" << parameters.fractureEnergy << " Lch=" << parameters.characteristicLength
      << " openCracks=" << int(trial.open[0])+int(trial.open[1])
      << " crushed=" << int(trial.crushed) << endln;
}
Response* AdinaConcrete2D::setResponse(const char** argv,int argc,OPS_Stream& output) {
    if(argc>0 && !std::strcmp(argv[0],"fractureState")) {
        output.tag("NdMaterialOutput");output.attr("matType",getClassType());output.attr("matTag",getTag());
        const char* names[]={"Gf","Lch","kappa1","strength1","offset1","dissipation1",
                            "kappa2","strength2","offset2","dissipation2"};
        for(const char* name:names)output.tag("ResponseType",name);
        output.endTag();return new MaterialResponse(this,101,Vector(10));
    }
    if(argc>0 && (!std::strcmp(argv[0],"crackState") || !std::strcmp(argv[0],"state"))) {
        output.tag("NdMaterialOutput");output.attr("matType",getClassType());output.attr("matTag",getTag());
        const char* names[]={"open1","open2","crushed","maximumLoading","hasCrackFrame",
                             "angleRadians","crackStrain1","crackStrain2","crushingStress","crushingStrain"};
        for(const char* name:names) output.tag("ResponseType",name);
        output.endTag();return new MaterialResponse(this,100,Vector(10));
    }
    return NDMaterial::setResponse(argv,argc,output);
}
int AdinaConcrete2D::getResponse(int id,Information& info) {
    if(id==101) {
        Vector v(10);int k=0;v(k++)=parameters.fractureEnergy;v(k++)=parameters.characteristicLength;
        for(const auto& h:trial.tension) {
            v(k++)=h.kappa;v(k++)=h.strength;v(k++)=h.offset;
            v(k++)=adina::fracture::dissipation(parameters.E,parameters.fractureEnergy,parameters.characteristicLength,h);
        }
        return info.setVector(v);
    }
    if(id==100) {
        Vector v(10);v(0)=trial.open[0];v(1)=trial.open[1];v(2)=trial.crushed;
        v(3)=trial.maximumLoading;v(4)=trial.hasCrackFrame;v(5)=trial.angle;
        v(6)=trial.crackStrain[0];v(7)=trial.crackStrain[1];
        v(8)=trial.crushingStress;v(9)=trial.crushingStrain;return info.setVector(v);
    }
    return NDMaterial::getResponse(id,info);
}
