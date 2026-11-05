// OpenSees adapter of the Jiang Jianjing nonlinear-elastic concrete model
// (teaching model, see JiangJianjing.h). Three-dimensional material; use the
// standard PlaneStress / PlateFiber wrappers for 2D and shell elements.
#ifndef OPS_JIANG_JIANJING_CONCRETE_H
#define OPS_JIANG_JIANJING_CONCRETE_H
#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include "JiangJianjing.h"

class JiangJianjingConcrete : public NDMaterial {
public:
    JiangJianjingConcrete(); // object broker / restart constructor
    JiangJianjingConcrete(int tag, const jiangjj::Parameters& parameters);
    const char* getClassType() const override { return "JiangJianjingConcrete"; }
    double getRho() override { return parameters.rho; }
    int setTrialStrain(const Vector&) override;
    int setTrialStrain(const Vector& v, const Vector&) override { return setTrialStrain(v); }
    int setTrialStrainIncr(const Vector&) override;
    int setTrialStrainIncr(const Vector& v, const Vector&) override { return setTrialStrainIncr(v); }
    const Vector& getStress() override { return stress; }
    const Vector& getStrain() override { return strain; }
    const Matrix& getTangent() override { return tangent; }
    const Matrix& getInitialTangent() override { return initialTangent; }
    int commitState() override;
    int revertToLastCommit() override;
    int revertToStart() override;
    NDMaterial* getCopy() override;
    NDMaterial* getCopy(const char*) override;
    const char* getType() const override { return "ThreeDimensional"; }
    int getOrder() const override { return 6; }
    int sendSelf(int, Channel&) override;
    int recvSelf(int, Channel&, FEM_ObjectBroker&) override;
    void Print(OPS_Stream&, int flag = 0) override;
    Response* setResponse(const char**, int, OPS_Stream&) override;
    int getResponse(int, Information&) override;
private:
    // Every trial is evaluated from the committed state (as Marc calls HYPELA2
    // with the stress and strain at the start of the increment), so repeated
    // trials, rejected iterations and reverts never change the history.
    jiangjj::Parameters parameters;
    jiangjj::State committed, trial;
    jiangjj::Matrix6 committedTangent {};
    Vector strain = Vector(6), stress = Vector(6);
    Matrix tangent = Matrix(6, 6), initialTangent = Matrix(6, 6);
    void synchronize(const jiangjj::Matrix6&);
};
#endif
