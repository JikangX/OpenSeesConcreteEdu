// Tcl/Python parsing is separate from the NDMaterial lifecycle and numerics.
#include "AdinaConcrete2D.h"
#include "AdinaConcreteOptions.h"
#include <elementAPI.h>
#include <cstring>
#include <cmath>

namespace {
void* create(bool plate) {
    // This shared factory serves Tcl, Python and the newer runtime. Compression
    // signs follow OpenSees (negative); envelope arrays are dimensionless.
    if(OPS_GetNumRemainingInputArgs()<8) {
        opserr << "Want: nDMaterial AdinaConcrete2D tag E nu ft fc epsc fu epsu "
                "<-G transverseG> <-rho rho> <-retention normal shear> <-controls beta gamma kappa alpha> "
                "<-envelope 24_values> <-maxStrainIncrement value> <-fractureEnergy Gf Lch> "
                "<-tangent secant|algorithmic> <-update committed|iterative|hybrid> "
                "<-crushedUnloading original|elastic>\n";
        return nullptr;
    }
    int tag,n=1;
    if(OPS_GetIntInput(&n,&tag)<0) return nullptr;
    double data[7];n=7;
    if(OPS_GetDoubleInput(&n,data)<0) return nullptr;
    adina::Parameters p;
    p.E=data[0];p.nu=data[1];p.ft=data[2];p.fc=data[3];
    p.epsc=data[4];p.fu=data[5];p.epsu=data[6];
    double G=p.E/(2.*(1.+p.nu));
    adina::SolverOptionSelection solverOptions;
    while(OPS_GetNumRemainingInputArgs()>0) {
        const char* option=OPS_GetString();
        solverOptions.observe(option);
        int keyword=adina::parseKeywordOption(option,p);
        if(keyword<0) {opserr << "AdinaConcrete2D: invalid value for " << option << endln;return nullptr;}
        if(keyword>0) continue;
        double values[24];
        if(!std::strcmp(option,"-G") || !std::strcmp(option,"-rho") || !std::strcmp(option,"-maxStrainIncrement")) n=1;
        else if(!std::strcmp(option,"-retention") || !std::strcmp(option,"-fractureEnergy")) n=2;
        else if(!std::strcmp(option,"-controls")) n=4;
        else if(!std::strcmp(option,"-envelope")) n=24;
        else {opserr << "AdinaConcrete2D: unknown option " << option << endln;return nullptr;}
        if(OPS_GetNumRemainingInputArgs()<n || OPS_GetDoubleInput(&n,values)<0) {
            opserr << "AdinaConcrete2D: incomplete option " << option << endln;return nullptr;
        }
        if(!std::strcmp(option,"-G")) G=values[0];
        else if(!std::strcmp(option,"-rho")) p.rho=values[0];
        else if(!std::strcmp(option,"-maxStrainIncrement")) p.maxStrainIncrement=values[0];
        else if(!std::strcmp(option,"-fractureEnergy")) {p.fractureEnergy=values[0];p.characteristicLength=values[1];}
        else if(!std::strcmp(option,"-retention")) {p.normalRetention=values[0];p.shearRetention=values[1];}
        else if(!std::strcmp(option,"-controls")) {
            p.beta=values[0];p.gamma=values[1];p.kappa=values[2];p.alpha=values[3];
        } else for(int i=0;i<6;++i) {
            p.sp1[i]=values[i];p.sp31[i]=values[i+6];p.sp32[i]=values[i+12];p.sp33[i]=values[i+18];
        }
    }
    solverOptions.applyDefaults(p);
    std::string error;
    if(!p.validate(error)) {opserr << "AdinaConcrete2D: " << error.c_str() << endln;return nullptr;}
    if(!std::isfinite(G) || G<=0.) {opserr << "AdinaConcrete2D: -G must be positive and finite\n";return nullptr;}
    return new AdinaConcrete2D(tag,p,plate,G);
}

} // parsing is independent of constitutive integration
void* OPS_AdinaConcrete2D() {return create(false);}
void* OPS_AdinaConcrete2DPlateFiber() {return create(true);}
