// OpenSees adapter for the separately testable ADINA material-point kernel.
#ifndef OPS_ADINA_CONCRETE3D_H
#define OPS_ADINA_CONCRETE3D_H
#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include "adinaConcrete3D/Concrete3D.h"

class AdinaConcrete3D : public NDMaterial {
public:
    AdinaConcrete3D(); // object broker / restart constructor
    AdinaConcrete3D(int tag, const adina::Parameters& parameters);
    const char* getClassType() const override { return "AdinaConcrete3D"; }
    double getRho() override { return parameters.rho; }
    int setTrialStrain(const Vector&) override;
    int setTrialStrain(const Vector& v,const Vector&) override {return setTrialStrain(v);}
    int setTrialStrainIncr(const Vector&) override;
    int setTrialStrainIncr(const Vector& v,const Vector&) override {return setTrialStrainIncr(v);}
    const Vector& getStress() override { return stress; }
    const Vector& getStrain() override { return strain; }
    const Matrix& getTangent() override { return tangent; }
    const Matrix& getInitialTangent() override { return initialTangent; }
    int commitState() override;
    int revertToLastCommit() override;
    int revertToStart() override;
    NDMaterial* getCopy() override;
    NDMaterial* getCopy(const char*) override;
    const char* getType() const override {return "ThreeDimensional";}
    int getOrder() const override {return 6;}
    int sendSelf(int,Channel&) override;
    int recvSelf(int,Channel&,FEM_ObjectBroker&) override;
    void Print(OPS_Stream&,int flag=0) override;
    Response* setResponse(const char**,int,OPS_Stream&) override;
    int getResponse(int,Information&) override;
private:
    // Each object owns both histories and all returned OpenSees buffers.
    // A rejected Newton trial never mutates committed state or sibling copies.
    adina::Parameters parameters;
    adina::State committed,trial;
    adina::Matrix6 committedTangent {};
    Vector strain=Vector(6),stress=Vector(6);
    Matrix tangent=Matrix(6,6),initialTangent=Matrix(6,6);
    void synchronize(const adina::Matrix6&);
};
#endif
