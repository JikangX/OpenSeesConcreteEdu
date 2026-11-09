// Keyword options shared by the AdinaConcrete2D and AdinaConcrete3D factories.
#ifndef OPS_ADINA_CONCRETE_OPTIONS_H
#define OPS_ADINA_CONCRETE_OPTIONS_H
#include <elementAPI.h>
#include <cstring>
#include "adinaConcrete3D/Concrete3D.h"

namespace adina {
/** Record explicit solver options so order-independent defaults can be used
 * for the optional continuous cohesive law. Brittle original mode uses
 * secant/committed (as ADINA); regularized mode defaults to committed replay
 * and its algorithmic Jacobian. Newton trial overshoots must not become irreversible
 * physical cohesive damage before an increment has actually converged.
 * Explicit user options always take precedence, including research runs of
 * the path-dependent iterative variant; those require separate verification.
 */
struct SolverOptionSelection {
    bool tangent=false, update=false;
    void observe(const char* option) {
        tangent=tangent || !std::strcmp(option,"-tangent");
        update=update || !std::strcmp(option,"-update");
    }
    void applyDefaults(Parameters& p) const {
        if(p.regularizedTension()) {
            if(!tangent) p.tangent=Parameters::Tangent::Algorithmic;
            if(!update) p.update=Parameters::Update::Committed;
        }
    }
};
/** Parse "-tangent", "-update" or "-crushedUnloading" and its keyword value.
 * Returns 1 when `option` was one of them and has been applied, 0 when it is
 * some other option, -1 when its value is missing or not recognised.
 */
inline int parseKeywordOption(const char* option,Parameters& p) {
    const bool tangent=!std::strcmp(option,"-tangent"),update=!std::strcmp(option,"-update"),
               crushed=!std::strcmp(option,"-crushedUnloading");
    if(!tangent && !update && !crushed) return 0;
    const char* value=OPS_GetNumRemainingInputArgs()>0 ? OPS_GetString() : "";
    auto is=[&](const char* x) {return !std::strcmp(value,x);};
    if(tangent && is("secant")) p.tangent=Parameters::Tangent::Secant;
    else if(tangent && is("algorithmic")) p.tangent=Parameters::Tangent::Algorithmic;
    else if(update && is("iterative")) p.update=Parameters::Update::Iterative;
    else if(update && is("committed")) p.update=Parameters::Update::Committed;
    else if(update && is("hybrid")) p.update=Parameters::Update::Hybrid;
    else if(crushed && is("original")) p.crushedUnloading=Parameters::CrushedUnloading::Original;
    else if(crushed && is("elastic")) p.crushedUnloading=Parameters::CrushedUnloading::Elastic;
    else return -1;
    return 1;
}
}
#endif
