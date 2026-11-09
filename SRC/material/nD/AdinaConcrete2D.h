#ifndef OPS_ADINA_CONCRETE2D_H
#define OPS_ADINA_CONCRETE2D_H
#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include "adinaConcrete2D/Concrete2D.h"

/** OpenSees lifecycle adapter for the native two-dimensional model 5 kernel.
 * PlaneStress has three components [xx,yy,xy]. A PlateFiber copy has five
 * [xx,yy,xy,yz,xz], following PlateFromPlaneStressMaterial in this repository.
 * LayeredShell requests getCopy("PlateFiber") automatically. Each returned
 * copy owns its histories and buffers, including both transverse strains.
 * Transverse shear is an explicit elastic shell extension, not an ADINA 2D
 * failure law. G defaults to E/(2(1+nu)); -G overrides it in either factory.
 */
class AdinaConcrete2D : public NDMaterial {
public:
    AdinaConcrete2D(); // broker constructor; recvSelf restores the formulation
    AdinaConcrete2D(int tag,const adina::Parameters&,bool plate=false,double G=-1.);
    const char* getClassType() const override {return "AdinaConcrete2D";}
    const char* getType() const override {return plate ? "PlateFiber" : "PlaneStress";}
    int getOrder() const override {return plate ? 5 : 3;}
    double getRho() override {return parameters.rho;}
    int setTrialStrain(const Vector&) override;
    int setTrialStrain(const Vector& v,const Vector&) override {return setTrialStrain(v);}
    int setTrialStrainIncr(const Vector&) override;
    int setTrialStrainIncr(const Vector& v,const Vector&) override {return setTrialStrainIncr(v);}
    const Vector& getStrain() override {return strain;}
    const Vector& getStress() override {return stress;}
    const Matrix& getTangent() override {return tangent;}
    const Matrix& getInitialTangent() override {return initialTangent;}
    int commitState() override;
    int revertToLastCommit() override;
    int revertToStart() override;
    NDMaterial* getCopy() override;
    NDMaterial* getCopy(const char*) override;
    int sendSelf(int,Channel&) override;
    int recvSelf(int,Channel&,FEM_ObjectBroker&) override;
    void Print(OPS_Stream&,int flag=0) override;
    Response* setResponse(const char**,int,OPS_Stream&) override;
    int getResponse(int,Information&) override;
private:
    adina::Parameters parameters;
    adina::plane::State committed,trial;
    adina::plane::Matrix3 committedTangent {};
    std::array<double,2> committedShear {},trialShear {};
    bool plate=false;
    double transverseG=12500.;
    Vector strain=Vector(3),stress=Vector(3);
    Matrix tangent=Matrix(3,3),initialTangent=Matrix(3,3);
    void resizeBuffers();
    void synchronize(const adina::plane::Matrix3&);
    void initializeTangent();
};
#endif
