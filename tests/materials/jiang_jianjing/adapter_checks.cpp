// Adapter checks: JiangJianjingConcrete linked against the real OpenSees NDMaterial,
// PlaneStress and PlateFiber classes (parser, uniaxial strengths, commit/revert,
// sendSelf/recvSelf restart, state response, crushing, -beta current).
// Build (from the repository root; INC = -ISRC plus every directory under SRC and OTHER):
//   c++ -std=c++17 $INC tests/materials/jiang_jianjing/adapter_checks.cpp \
//     SRC/material/nD/JiangJianjingConcrete/*.cpp SRC/material/nD/{NDMaterial,PlaneStressMaterial,PlateFiberMaterial, \
//     BeamFiberMaterial,BeamFiberMaterial2d,BeamFiberMaterial2dPS}.cpp SRC/material/Material.cpp \
//     SRC/matrix/{Vector,Matrix,ID}.cpp SRC/handler/{OPS_Stream,StandardStream,DummyStream}.cpp \
//     SRC/actor/actor/MovableObject.cpp SRC/actor/channel/Channel.cpp SRC/tagged/TaggedObject.cpp \
//     SRC/tagged/storage/{MapOfTaggedObjects,MapOfTaggedObjectsIter}.cpp \
//     SRC/recorder/response/{MaterialResponse,Response}.cpp SRC/element/Information.cpp -llapack -lblas
// Driver: JiangJianjingConcrete adapter linked against real OpenSees classes.
#include <elementAPI.h>
#include <JiangJianjingConcrete.h>
#include <PlaneStressMaterial.h>
#include <PlateFiberMaterial.h>
#include <StandardStream.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Response.h>
#include <DummyStream.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

StandardStream sserr;
OPS_Stream* opserrPtr = &sserr;

// --- minimal OPS_ argument API over a string vector -------------------------
static std::vector<std::string> args; static size_t cur = 0;
extern "C" int OPS_GetNumRemainingInputArgs() { return int(args.size() - cur); }
extern "C" int OPS_GetIntInput(int* n, int* d) {
    for (int i = 0; i < *n; ++i) { if (cur >= args.size()) return -1; char* e; long v = strtol(args[cur].c_str(), &e, 10); if (*e) return -1; d[i] = int(v); ++cur; }
    return 0;
}
extern "C" int OPS_GetDoubleInput(int* n, double* d) {
    for (int i = 0; i < *n; ++i) { if (cur >= args.size()) return -1; char* e; double v = strtod(args[cur].c_str(), &e); if (*e) return -1; d[i] = v; ++cur; }
    return 0;
}
extern "C" const char* OPS_GetStringFromAll(char* buf, int len) {
    if (cur >= args.size()) return nullptr; strncpy(buf, args[cur++].c_str(), len - 1); buf[len - 1] = 0; return buf;
}
void* OPS_JiangJianjingConcrete();
static NDMaterial* make(const std::string& cmd) {
    args.clear(); cur = 0; char* s = strdup(cmd.c_str());
    for (char* t = strtok(s, " "); t; t = strtok(nullptr, " ")) args.push_back(t);
    free(s); return static_cast<NDMaterial*>(OPS_JiangJianjingConcrete());
}

// --- in-memory channel -------------------------------------------------------
struct MemChannel : Channel {
    std::map<std::pair<int,int>, Vector> store;
    char* addToProgram() override { return nullptr; }
    int setUpConnection() override { return 0; }
    int setNextAddress(const ChannelAddress&) override { return 0; }
    ChannelAddress* getLastSendersAddress() override { return nullptr; }
    int sendObj(int, MovableObject&, ChannelAddress*) override { return -1; }
    int recvObj(int, MovableObject&, FEM_ObjectBroker&, ChannelAddress*) override { return -1; }
    int sendMsg(int, int, const Message&, ChannelAddress*) override { return -1; }
    int recvMsg(int, int, Message&, ChannelAddress*) override { return -1; }
    int sendMatrix(int, int, const Matrix&, ChannelAddress*) override { return -1; }
    int recvMatrix(int, int, Matrix&, ChannelAddress*) override { return -1; }
    int sendVector(int db, int ct, const Vector& v, ChannelAddress*) override { store[{db, ct}] = v; return 0; }
    int recvVector(int db, int ct, Vector& v, ChannelAddress*) override {
        auto it = store.find({db, ct}); if (it == store.end() || it->second.Size() != v.Size()) return -1; v = it->second; return 0;
    }
    int sendID(int, int, const ID&, ChannelAddress*) override { return -1; }
    int recvID(int, int, ID&, ChannelAddress*) override { return -1; }
};

static int checks = 0, fails = 0;
static void check(bool c, const char* msg) { ++checks; if (!c) { ++fails; printf("FAIL %s\n", msg); } }

// uniaxial stress in yy via PlaneStress wrapper: solve sigma_xx = 0 for eps_xx
static double uniaxialPeak(NDMaterial& ps, double epsEnd, int steps, bool& ok) {
    Vector e(3); double peak = 0.; ok = true;
    for (int k = 1; k <= steps; ++k) {
        e(1) = epsEnd * k / steps;
        int it = 0;
        for (; it < 50; ++it) {
            ps.setTrialStrain(e);
            const Vector& s = ps.getStress(); const Matrix& C = ps.getTangent();
            if (std::fabs(s(0)) < 1e-9) break;
            e(0) -= s(0) / C(0, 0);
        }
        if (it == 50) { ok = false; break; }
        ps.commitState();
        double syy = ps.getStress()(1);
        if (std::fabs(syy) > std::fabs(peak)) peak = syy;
    }
    return peak;
}

int main() {
    // 1. parser
    check(make("1 30 30000 0.18") != nullptr, "default arguments");
    check(make("2 30 30000 0.18 -0.0035 -criterion 1.2856 1.4268 10.2551 3.2128 -rho 2.4e-9") != nullptr, "all options");
    check(make("3 -30 30000 0.18") == nullptr, "negative fc rejected");
    check(make("4 30 30000 0.18 0.0033") == nullptr, "positive epsCrush rejected");
    check(make("5 30 30000 0.18 -unknown 1") == nullptr, "unknown option rejected");
    check(make("6 30 30000 0.18 -beta current -tangent secant") != nullptr, "-beta current");
    check(make("7 30 30000 0.18 -beta later") == nullptr, "bad -beta rejected");
    check(make("8 30 30000") == nullptr, "too few args rejected");

    // 2. uniaxial compression / tension, PlaneStress wrapper, secant tangent
    NDMaterial* m = make("10 30 30000 0.18 -tangent secant");
    check(m && !strcmp(m->getType(), "ThreeDimensional") && m->getOrder() == 6, "3D type");
    {
        PlaneStressMaterial ps(11, *m); bool ok;
        double pk = uniaxialPeak(ps, -0.0030, 300, ok);
        printf("INFO PlaneStress uniaxial compression peak %.3f MPa (fc 30)\n", pk);
        check(ok && std::fabs(-pk - 30.) < 0.5, "compressive peak ~ fc");
    }
    {
        PlaneStressMaterial ps(12, *m); bool ok;
        double pk = uniaxialPeak(ps, 0.0003, 600, ok);
        printf("INFO PlaneStress uniaxial tension peak (committed beta, 600 steps) %.4f MPa (0.1 fc = 3)\n", pk);
        PlaneStressMaterial ps2(13, *m); bool ok2;
        double pk2 = uniaxialPeak(ps2, 0.0003, 6000, ok2);
        printf("INFO PlaneStress uniaxial tension peak (committed beta, 6000 steps) %.4f MPa\n", pk2);
        check(ok && ok2 && std::fabs(pk2 - 3.) < std::fabs(pk - 3.), "explicit-beta overshoot shrinks with step size");
        NDMaterial* mc = make("14 30 30000 0.18 -beta current -tangent secant");
        PlaneStressMaterial ps3(15, *mc); bool ok3;
        double pk3 = uniaxialPeak(ps3, 0.0003, 600, ok3);
        printf("INFO PlaneStress uniaxial tension peak (-beta current) %.4f MPa\n", pk3);
        check(ok3 && std::fabs(pk3 - 3.) < 0.01, "tensile peak = 0.1 fc with -beta current");
    }
    // first step is linear elastic: sigma = D(E0) eps
    {
        NDMaterial* a = m->getCopy(); Vector e(6); e(2) = -1e-5;
        a->setTrialStrain(e);
        double lam = 30000 * 0.18 / (1.18 * (1 - 0.36)), mu = 30000 / 2.36;
        check(std::fabs(a->getStress()(2) - (lam + 2 * mu) * -1e-5) < 1e-9, "first step elastic");
        check(std::fabs(a->getInitialTangent()(0, 0) - (lam + 2 * mu)) < 1e-6, "initial tangent");
        delete a;
    }
    // 3. trial / revert / commit independence and restart
    {
        NDMaterial* a = make("20 30 30000 0.18");
        Vector e(6);
        for (int k = 1; k <= 40; ++k) { e(0) = 0.0002 * k / 40; e(1) = e(0); e(2) = -0.0016 * k / 40; e(3) = 0.0002 * k / 40; a->setTrialStrain(e); a->commitState(); }
        Vector committed = a->getStress();
        Vector e2(e); e2(2) -= 1e-4;
        a->setTrialStrain(e2); Vector s1 = a->getStress();
        a->setTrialStrain(e2); a->setTrialStrain(e2); Vector s2 = a->getStress();
        check(s1 == s2, "repeated trial identical");
        a->revertToLastCommit();
        check(a->getStress() == committed, "revert restores committed stress");
        a->setTrialStrain(e);
        check(a->getStress() == committed, "re-evaluation at committed strain keeps state");
        // response 'state'
        DummyStream ds; const char* av[] = {"state"};
        Response* r = a->setResponse(av, 1, ds);
        check(r != nullptr, "state response");
        if (r) { r->getResponse(); const Vector& v = *r->getInformation().theVector;
                 printf("INFO state: crushed %g beta %.4f betaMax %.4f Es %.1f Et %.1f solvable %g\n", v(0), v(1), v(2), v(3), v(4), v(5));
                 check(v.Size() == 6 && v(1) > 0 && v(1) < 1, "beta in (0,1)"); delete r; }
        // send / receive
        MemChannel ch; a->setDbTag(5);
        check(a->sendSelf(1, ch) == 0, "sendSelf");
        JiangJianjingConcrete b; b.setDbTag(5); FEM_ObjectBroker* broker = nullptr;
        check(b.recvSelf(1, ch, *broker) == 0, "recvSelf");
        check(b.getTag() == 20, "restored tag");
        check(b.getStress() == committed, "restored stress");
        a->setTrialStrain(e2); b.setTrialStrain(e2);
        bool same = a->getStress() == b.getStress();
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) same = same && a->getTangent()(i, j) == b.getTangent()(i, j);
        check(same, "restart continues bitwise");
        delete a;
    }
    // 4. crushing (irreversible)
    {
        NDMaterial* a = make("30 30 30000 0.18 -tangent secant"); Vector e(6);
        for (int k = 1; k <= 80; ++k) { e(2) = -0.006 * k / 80; e(0) = e(1) = 0.0004 * k / 80; a->setTrialStrain(e); a->commitState(); }
        DummyStream ds; const char* av[] = {"state"}; Response* r = a->setResponse(av, 1, ds); r->getResponse();
        check((*r->getInformation().theVector)(0) == 1., "crushed beyond epsCrush");
        check(a->getStress().Norm() < 1e-12, "stress zero after crushing");
        delete r; delete a;
    }
    // 5. PlateFiber wrapper (layered shell)
    {
        PlateFiberMaterial pf(40, *m); Vector e(5); bool finite = true;
        for (int k = 1; k <= 50; ++k) { e(0) = -0.0015 * k / 50; e(2) = 0.0003 * k / 50; if (pf.setTrialStrain(e)) finite = false; pf.commitState();
            for (int i = 0; i < 5; ++i) finite = finite && std::isfinite(pf.getStress()(i)); }
        printf("INFO PlateFiber sigma11 %.3f MPa\n", pf.getStress()(0));
        check(finite && pf.getStress()(0) < -10., "PlateFiber wrapper works");
    }
    // 6. -beta current option
    {
        NDMaterial* a = make("50 30 30000 0.18 -beta current -tangent secant"); NDMaterial* b = make("51 30 30000 0.18 -beta current -tangent secant");
        Vector e(6); e(2) = -0.0012; e(0) = e(1) = 0.00015;
        for (int k = 1; k <= 10; ++k) { Vector x(e); x *= k / 10.; a->setTrialStrain(x); a->commitState(); }
        b->setTrialStrain(e); b->commitState();
        Vector d = a->getStress(); d -= b->getStress();
        printf("INFO -beta current: 10 steps vs 1 step stress difference %.2e MPa\n", d.Norm());
        check(d.Norm() < 1e-8, "-beta current is path independent");
        delete a; delete b;
    }
    printf("%s %d adapter checks (%d failed)\n", fails ? "FAIL" : "PASS", checks, fails);
    return fails != 0;
}
