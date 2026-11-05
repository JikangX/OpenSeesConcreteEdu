#include "JiangJianjingConcrete.h"
#include <classTags.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <elementAPI.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>

void* OPS_JiangJianjingConcrete() {
    // nDMaterial JiangJianjingConcrete tag fc E0 nu <epsCrush>
    //     <-criterion a b c d> <-tangent original|secant> <-beta committed|current> <-rho rho>
    // fc > 0 (compressive strength), epsCrush < 0; defaults are the values of
    // the original MSC.Marc program (epsCrush -0.0033; a b c d for ft/fc 0.1).
    if(OPS_GetNumRemainingInputArgs()<4) {
        opserr << "Want: nDMaterial JiangJianjingConcrete tag fc E0 nu <epsCrush> "
                  "<-criterion a b c d> <-tangent original|secant> <-beta committed|current> <-rho rho>\n";
        return nullptr;
    }
    int tag,n=1;
    if(OPS_GetIntInput(&n,&tag)<0) return nullptr;
    double data[3];n=3;
    if(OPS_GetDoubleInput(&n,data)<0) return nullptr;
    jiangjj::Parameters p;
    p.fc=data[0];p.E0=data[1];p.nu=data[2];
    while(OPS_GetNumRemainingInputArgs()>0) {
        // Numbers and keywords are both accepted (Tcl and Python): the
        // optional epsCrush is a plain number after nu.
        char buffer[64];
        const char* option=OPS_GetStringFromAll(buffer,int(sizeof(buffer)));
        if(!option) {opserr << "JiangJianjingConcrete: invalid argument" << endln;return nullptr;}
        double values[4];
        if(!std::strcmp(option,"-tangent")) {
            char vb[32];
            const char* value=OPS_GetNumRemainingInputArgs()>0 ? OPS_GetStringFromAll(vb,int(sizeof(vb))) : nullptr;
            if(value && !std::strcmp(value,"original")) p.tangent=jiangjj::Parameters::Tangent::Original;
            else if(value && !std::strcmp(value,"secant")) p.tangent=jiangjj::Parameters::Tangent::Secant;
            else {opserr << "JiangJianjingConcrete: -tangent original|secant" << endln;return nullptr;}
        } else if(!std::strcmp(option,"-beta")) {
            char vb[32];
            const char* value=OPS_GetNumRemainingInputArgs()>0 ? OPS_GetStringFromAll(vb,int(sizeof(vb))) : nullptr;
            if(value && !std::strcmp(value,"committed")) p.betaFrom=jiangjj::Parameters::BetaFrom::Committed;
            else if(value && !std::strcmp(value,"current")) p.betaFrom=jiangjj::Parameters::BetaFrom::Current;
            else {opserr << "JiangJianjingConcrete: -beta committed|current" << endln;return nullptr;}
        } else if(!std::strcmp(option,"-criterion")) {
            n=4;
            if(OPS_GetNumRemainingInputArgs()<4 || OPS_GetDoubleInput(&n,values)<0) {
                opserr << "JiangJianjingConcrete: -criterion needs a b c d\n";return nullptr;
            }
            p.a=values[0];p.b=values[1];p.c=values[2];p.d=values[3];
        } else if(!std::strcmp(option,"-rho")) {
            n=1;
            if(OPS_GetNumRemainingInputArgs()<1 || OPS_GetDoubleInput(&n,values)<0) {
                opserr << "JiangJianjingConcrete: -rho needs a value\n";return nullptr;
            }
            p.rho=values[0];
        } else {
            // A plain number directly after nu is epsCrush.
            char* end=nullptr;
            const double x=std::strtod(option,&end);
            if(end==option || *end!='\0') {
                opserr << "JiangJianjingConcrete: unknown option " << option << endln;return nullptr;
            }
            p.epsCrush=x;
        }
    }
    std::string error;
    if(!p.validate(error)) {opserr << "JiangJianjingConcrete: " << error.c_str() << endln;return nullptr;}
    return new JiangJianjingConcrete(tag,p);
}

JiangJianjingConcrete::JiangJianjingConcrete()
    :NDMaterial(0,ND_TAG_JiangJianjingConcrete) {revertToStart();}
JiangJianjingConcrete::JiangJianjingConcrete(int tag,const jiangjj::Parameters& p)
    :NDMaterial(tag,ND_TAG_JiangJianjingConcrete),parameters(p) {revertToStart();}

void JiangJianjingConcrete::synchronize(const jiangjj::Matrix6& c) {
    for(int i=0;i<6;++i) {
        strain(i)=trial.strain[i];stress(i)=trial.stress[i];
        for(int j=0;j<6;++j) tangent(i,j)=c[i][j];
    }
}
int JiangJianjingConcrete::setTrialStrain(const Vector& v) {
    if(v.Size()!=6) {opserr << "JiangJianjingConcrete requires six engineering strain components\n";return -1;}
    jiangjj::Vector6 e;
    for(int i=0;i<6;++i) {
        e[i]=v(i);
        if(!std::isfinite(e[i])) return -1;
    }
    // A request for exactly the committed strain is not a new increment.
    // Get_DS takes beta from the stress at the START of the increment, so a
    // zero-increment re-evaluation would NOT reproduce the committed stress
    // (it applies the newer secant modulus to the same strain). OpenSees
    // re-evaluates elements at unchanged displacements (e.g. after commit,
    // database restore, element queries); keep the converged state then.
    if(e==committed.strain) {
        trial=committed;
        synchronize(committedTangent);
        return 0;
    }
    jiangjj::Matrix6 c;
    trial=jiangjj::update(parameters,committed,e,&c);
    synchronize(c);
    return 0;
}
int JiangJianjingConcrete::setTrialStrainIncr(const Vector& v) {
    if(v.Size()!=6) return -1;
    Vector e(strain);e+=v;return setTrialStrain(e);
}
int JiangJianjingConcrete::commitState() {
    committed=trial;
    for(int i=0;i<6;++i) for(int j=0;j<6;++j) committedTangent[i][j]=tangent(i,j);
    return 0;
}
int JiangJianjingConcrete::revertToLastCommit() {trial=committed;synchronize(committedTangent);return 0;}
int JiangJianjingConcrete::revertToStart() {
    committed=jiangjj::State{};
    committed.secantModulus=committed.tangentModulus=parameters.E0;
    trial=committed;
    committedTangent=jiangjj::elasticStiffness(parameters.E0,parameters.nu);
    for(int i=0;i<6;++i) for(int j=0;j<6;++j) initialTangent(i,j)=committedTangent[i][j];
    synchronize(committedTangent);
    return 0;
}
NDMaterial* JiangJianjingConcrete::getCopy() {
    auto* copy=new JiangJianjingConcrete(*this);
    copy->setDbTag(0);
    return copy;
}
NDMaterial* JiangJianjingConcrete::getCopy(const char* type) {
    return type && (!std::strcmp(type,"ThreeDimensional") || !std::strcmp(type,"3D")) ? getCopy() : nullptr;
}

namespace {
// Flat, versioned serialization of parameters, both states and both tangents.
// Version 2 adds the BetaFrom flag after the tangent flag; version-1
// records (without it) are read as BetaFrom::Committed.
constexpr int wireVersion=2;
constexpr int stateSize=6+6+6;
constexpr int wireSize=2+11+2*stateSize+2*36;
constexpr int wireSizeV1=wireSize-1;
void packState(Vector& v,int& k,const jiangjj::State& s) {
    for(double x:s.strain) v(k++)=x;
    for(double x:s.stress) v(k++)=x;
    v(k++)=s.crushed;v(k++)=s.beta;v(k++)=s.betaMax;
    v(k++)=s.secantModulus;v(k++)=s.tangentModulus;v(k++)=s.indexSolvable;
}
bool unpackState(const Vector& v,int& k,jiangjj::State& s) {
    for(double& x:s.strain) x=v(k++);
    for(double& x:s.stress) x=v(k++);
    if((v(k)!=0. && v(k)!=1.) || (v(k+5)!=0. && v(k+5)!=1.)) return false;
    s.crushed=v(k++)!=0.;s.beta=v(k++);s.betaMax=v(k++);
    s.secantModulus=v(k++);s.tangentModulus=v(k++);s.indexSolvable=v(k++)!=0.;
    return true;
}
}
int JiangJianjingConcrete::sendSelf(int commitTag,Channel& channel) {
    Vector v(wireSize);int k=0;
    v(k++)=wireVersion;v(k++)=getTag();
    for(double x:{parameters.fc,parameters.E0,parameters.nu,parameters.epsCrush,
                  parameters.a,parameters.b,parameters.c,parameters.d,parameters.rho,
                  parameters.tangent==jiangjj::Parameters::Tangent::Secant ? 1. : 0.,
                  parameters.betaFrom==jiangjj::Parameters::BetaFrom::Current ? 1. : 0.}) v(k++)=x;
    packState(v,k,committed);packState(v,k,trial);
    for(const auto& row:committedTangent) for(double x:row) v(k++)=x;
    for(int i=0;i<6;++i) for(int j=0;j<6;++j) v(k++)=tangent(i,j);
    return k==wireSize ? channel.sendVector(getDbTag(),commitTag,v) : -1;
}
int JiangJianjingConcrete::recvSelf(int commitTag,Channel& channel,FEM_ObjectBroker&) {
    Vector v(wireSize);
    if(channel.recvVector(getDbTag(),commitTag,v)<0) {
        // a version-1 record is one entry shorter
        Vector w(wireSizeV1);
        if(channel.recvVector(getDbTag(),commitTag,w)<0 || w(0)!=1.) return -1;
        for(int i=0;i<12;++i) v(i)=w(i);
        v(12)=0.;
        for(int i=12;i<wireSizeV1;++i) v(i+1)=w(i);
        v(0)=wireVersion;
    }
    if(v(0)!=wireVersion) return -1;
    for(int i=0;i<v.Size();++i) if(!std::isfinite(v(i))) return -1;
    if(v(1)<std::numeric_limits<int>::min() || v(1)>std::numeric_limits<int>::max() ||
       std::floor(v(1))!=v(1)) return -1;
    int k=2;jiangjj::Parameters p;
    for(double* x:{&p.fc,&p.E0,&p.nu,&p.epsCrush,&p.a,&p.b,&p.c,&p.d,&p.rho}) *x=v(k++);
    if(v(k)!=0. && v(k)!=1.) return -1;
    p.tangent=v(k++)==1. ? jiangjj::Parameters::Tangent::Secant : jiangjj::Parameters::Tangent::Original;
    if(v(k)!=0. && v(k)!=1.) return -1;
    p.betaFrom=v(k++)==1. ? jiangjj::Parameters::BetaFrom::Current : jiangjj::Parameters::BetaFrom::Committed;
    std::string error;if(!p.validate(error)) return -1;
    jiangjj::State cs,ts;jiangjj::Matrix6 cc,tc;
    if(!unpackState(v,k,cs) || !unpackState(v,k,ts)) return -1;
    for(auto& row:cc) for(double& x:row) x=v(k++);
    for(auto& row:tc) for(double& x:row) x=v(k++);
    if(k!=wireSize) return -1;
    parameters=p;setTag(int(v(1)));revertToStart();
    committed=cs;trial=ts;committedTangent=cc;synchronize(tc);
    return 0;
}
void JiangJianjingConcrete::Print(OPS_Stream& s,int) {
    s << "JiangJianjingConcrete tag: " << getTag() << " fc " << parameters.fc << " E0 " << parameters.E0
      << " nu " << parameters.nu << " epsCrush " << parameters.epsCrush << " criterion "
      << parameters.a << " " << parameters.b << " " << parameters.c << " " << parameters.d << endln;
    s << "  committed: crushed " << committed.crushed << " beta " << committed.beta
      << " betaMax " << committed.betaMax << endln;
}
Response* JiangJianjingConcrete::setResponse(const char** argv,int argc,OPS_Stream& output) {
    if(argc>0 && (!std::strcmp(argv[0],"state") || !std::strcmp(argv[0],"nonlinearIndex"))) {
        output.tag("NdMaterialOutput");output.attr("matType",getClassType());output.attr("matTag",getTag());
        const char* names[]={"crushed","beta","betaMax","secantModulus","tangentModulus","indexSolvable"};
        for(const char* name:names) output.tag("ResponseType",name);
        output.endTag();
        return new MaterialResponse(this,100,Vector(6));
    }
    return NDMaterial::setResponse(argv,argc,output);
}
int JiangJianjingConcrete::getResponse(int id,Information& info) {
    if(id==100) {
        Vector v(6);
        v(0)=trial.crushed;v(1)=trial.beta;v(2)=trial.betaMax;
        v(3)=trial.secantModulus;v(4)=trial.tangentModulus;v(5)=trial.indexSolvable;
        return info.setVector(v);
    }
    return NDMaterial::getResponse(id,info);
}
