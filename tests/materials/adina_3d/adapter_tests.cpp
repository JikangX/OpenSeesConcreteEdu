// Exercise the real OpenSees NDMaterial adapter and Channel serialization.
// No mock constitutive model or substitute OpenSees interface is used.
#include <AdinaConcrete3D.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <StandardStream.h>
#include <OPS_Globals.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>


class MemoryChannel final : public Channel {
public:
    Vector data;
    char* addToProgram() override {return nullptr;}
    int setUpConnection() override {return 0;}
    int setNextAddress(const ChannelAddress&) override {return 0;}
    ChannelAddress* getLastSendersAddress() override {return nullptr;}
    int sendObj(int,MovableObject&,ChannelAddress*) override {return -1;}
    int recvObj(int,MovableObject&,FEM_ObjectBroker&,ChannelAddress*) override {return -1;}
    int sendMsg(int,int,const Message&,ChannelAddress*) override {return -1;}
    int recvMsg(int,int,Message&,ChannelAddress*) override {return -1;}
    int sendMatrix(int,int,const Matrix&,ChannelAddress*) override {return -1;}
    int recvMatrix(int,int,Matrix&,ChannelAddress*) override {return -1;}
    int sendVector(int,int,const Vector& v,ChannelAddress*) override {data=v;return 0;}
    int recvVector(int,int,Vector& v,ChannelAddress*) override {v=data;return 0;}
    int sendID(int,int,const ID&,ChannelAddress*) override {return -1;}
    int recvID(int,int,ID&,ChannelAddress*) override {return -1;}
};
namespace {
int checks=0;
void check(bool condition,const char* message) {
    ++checks;if(!condition)throw std::runtime_error(message);
}
void equal(const Vector& a,const Vector& b) {
    check(a.Size()==b.Size(),"vector sizes");
    for(int i=0;i<a.Size();++i)check(a(i)==b(i),"exact state/buffer restoration");
}
void equal(const Matrix& a,const Matrix& b) {
    for(int i=0;i<6;++i)for(int j=0;j<6;++j)check(a(i,j)==b(i,j),"exact tangent restoration");
}
}
int main() {
    try {
        adina::Parameters p;p.rho=2.4;p.normalRetention=.0002;p.shearRetention=.6;
        p.tangent=adina::Parameters::Tangent::Algorithmic;
        p.update=adina::Parameters::Update::Committed;
        AdinaConcrete3D material(42,p);Vector e(6);
        check(material.getRho()==2.4 && material.getOrder()==6,"density/order");
        check(material.getCopy("PlaneStress")==nullptr,"reject unsupported formulation");
        e(2)=-.0005;check(material.setTrialStrain(e)==0,"compression trial");
        material.commitState();const Vector committedStress=material.getStress();
        const Vector committedStrain=material.getStrain();const Matrix committedTangent=material.getTangent();
        e(2)=-.0007;check(material.setTrialStrain(e)==0,"second trial");
        const Vector trialStress=material.getStress();const Matrix trialTangent=material.getTangent();
        material.revertToLastCommit();equal(material.getStress(),committedStress);
        equal(material.getStrain(),committedStrain);equal(material.getTangent(),committedTangent);
        check(material.setTrialStrain(e)==0,"repeat second trial");equal(material.getStress(),trialStress);
        equal(material.getTangent(),trialTangent);
        std::unique_ptr<NDMaterial> copy(material.getCopy("ThreeDimensional"));
        check(copy!=nullptr,"three dimensional clone");const Vector& borrowed=material.getStress();
        Vector ownStress=borrowed;e(0)=.0002;copy->setTrialStrain(e);
        equal(borrowed,ownStress);equal(material.getStress(),ownStress);
        MemoryChannel channel;FEM_ObjectBroker broker;AdinaConcrete3D restored;
        check(material.sendSelf(8,channel)==0,"sendSelf");check(channel.data.Size()==198,"v3 versioned wire size");
        check(restored.recvSelf(8,channel,broker)==0,"recvSelf");
        check(restored.getTag()==42 && restored.getRho()==p.rho,"serialized inputs and tag");
        equal(restored.getStress(),material.getStress());equal(restored.getStrain(),material.getStrain());
        equal(restored.getTangent(),material.getTangent());equal(restored.getInitialTangent(),material.getInitialTangent());
        const Vector currentWire=channel.data;channel.data.resize(173);channel.data(0)=1.;
        check(restored.recvSelf(8,channel,broker)==0,"backward compatible v1 checkpoint");
        equal(restored.getStress(),material.getStress());equal(restored.getTangent(),material.getTangent());
        channel.data=currentWire;channel.data.resize(193);channel.data(0)=2.;
        check(restored.recvSelf(8,channel,broker)==0,"backward compatible v2 checkpoint");
        equal(restored.getStress(),material.getStress());equal(restored.getTangent(),material.getTangent());
        channel.data=currentWire;
        restored.revertToLastCommit();equal(restored.getStress(),committedStress);equal(restored.getTangent(),committedTangent);
        restored.setTrialStrain(material.getStrain());equal(restored.getStress(),material.getStress());
        const Vector validWire=channel.data;
        for(int field:{0,1,54,56,57,173,174,175,176,177,193,194,195}) {
            channel.data=validWire;channel.data(field)=field==0?99.:1e100;
            Vector before=restored.getStress();
            check(restored.recvSelf(8,channel,broker)!=0,"reject corrupt wire data");equal(restored.getStress(),before);
        }
        channel.data=validWire;channel.data(105)=NAN;
        check(restored.recvSelf(8,channel,broker)!=0,"reject nonfinite serialized tangent");
        Vector invalid(6);invalid(0)=NAN;Vector before=material.getStress();Matrix beforeTangent=material.getTangent();
        check(material.setTrialStrain(invalid)!=0,"failed strain trial");equal(material.getStress(),before);equal(material.getTangent(),beforeTangent);
        check(material.setTrialStrain(Vector(3))!=0,"reject incorrect strain dimension");
        material.revertToStart();Vector zero(6);equal(material.getStress(),zero);equal(material.getStrain(),zero);
        equal(material.getTangent(),material.getInitialTangent());
        e.Zero();e(0)=1e-7;material.setTrialStrain(e);Vector increment(6);increment(0)=1e-7;
        material.setTrialStrainIncr(increment);check(material.getStrain()(0)==2e-7,"increment relative to current trial");
        material.commitState();copy.reset(material.getCopy());material.revertToStart();
        check(copy->getStrain()(0)==2e-7,"committed clone independent of original reset");
        // Round-trip a genuinely rotated, multiply cracked material. Committed
        // and trial frames differ after closing one plane; all nine direction
        // cosines and both histories must survive the Channel unchanged.
        AdinaConcrete3D rotated(55,p),rotatedRestored;
        auto frame=adina::principalStress({{3.,-2.,5.,1.,-3.,2.}}).frame;
        adina::Vector6 local{{.0002,.0003,.0004,0,0,0}};
        auto global=adina::strainToGlobal(local,frame);Vector strain3D(6);
        for(int i=0;i<6;++i)strain3D(i)=global[i];
        check(rotated.setTrialStrain(strain3D)==0,"rotated cracking trial");rotated.commitState();
        local[1]=-.0001;global=adina::strainToGlobal(local,frame);
        for(int i=0;i<6;++i)strain3D(i)=global[i];
        check(rotated.setTrialStrain(strain3D)==0,"rotated closure trial");
        check(rotated.sendSelf(9,channel)==0,"send rotated history");Vector rotatedWire=channel.data;
        check(rotatedRestored.recvSelf(9,channel,broker)==0,"receive rotated history");
        equal(rotatedRestored.getStress(),rotated.getStress());equal(rotatedRestored.getTangent(),rotated.getTangent());
        check(rotatedRestored.sendSelf(9,channel)==0,"repack rotated history");equal(channel.data,rotatedWire);
        rotated.revertToLastCommit();rotatedRestored.revertToLastCommit();
        equal(rotatedRestored.getStress(),rotated.getStress());
        rotated.setTrialStrain(strain3D);rotatedRestored.setTrialStrain(strain3D);
        equal(rotatedRestored.getStress(),rotated.getStress());equal(rotatedRestored.getTangent(),rotated.getTangent());
        auto gf=p;gf.nu=0.;gf.fractureEnergy=.0002;gf.characteristicLength=.1;
        AdinaConcrete3D regularized(181,gf),recovered;
        Vector de(6);de(0)=.0005;
        check(regularized.setTrialStrain(de)==0,"regularized damage trial");regularized.commitState();
        de(0)=.0007;check(regularized.setTrialStrain(de)==0,"regularized trial damage");
        check(regularized.sendSelf(10,channel)==0,"regularized send");const Vector gfWire=channel.data;
        check(recovered.recvSelf(10,channel,broker)==0,"regularized receive");
        equal(regularized.getStress(),recovered.getStress());equal(regularized.getTangent(),recovered.getTangent());
        recovered.sendSelf(10,channel);equal(gfWire,channel.data);
        copy.reset(regularized.getCopy());
        regularized.revertToLastCommit();recovered.revertToLastCommit();copy->revertToLastCommit();
        for(double strain:{.0002,0.,-.00002,0.,.0002,.0005,.0008}) {
            de(0)=strain;
            check(regularized.setTrialStrain(de)==0 && recovered.setTrialStrain(de)==0 && copy->setTrialStrain(de)==0,"regularized continuation");
            equal(regularized.getStress(),recovered.getStress());equal(regularized.getStress(),copy->getStress());
            equal(regularized.getTangent(),recovered.getTangent());
            regularized.commitState();recovered.commitState();copy->commitState();
        }
        for(int field:{173,174,175,176,177}) {
            channel.data=gfWire;channel.data(field)=field==173?0.:(field==175?-1.:1e100);
            Vector before=recovered.getStress();
            check(recovered.recvSelf(10,channel,broker)!=0,"reject inconsistent fracture checkpoint");equal(before,recovered.getStress());
        }
        regularized.setDbTag(73);copy.reset(regularized.getCopy());
        check(copy->getDbTag()==0 && regularized.getDbTag()==73,"restored 3D clone gets independent database identity");
        for(auto unloading:{adina::Parameters::CrushedUnloading::Original,adina::Parameters::CrushedUnloading::Elastic}) {
            auto mixed=gf;mixed.crushedUnloading=unloading;
            AdinaConcrete3D source(198,mixed),restoredMixed;Vector e(6);e(0)=.0005;
            check(source.setTrialStrain(e)==0,"form tensile history before crushing");source.commitState();
            e(2)=-.003;check(source.setTrialStrain(e)==0,"crush damaged 3D point");source.commitState();
            e(0)=.0015;check(source.setTrialStrain(e)==0,"advance beyond frozen tensile maximum");
            check(source.sendSelf(11,channel)==0,"send mixed tension / crush");const Vector wire=channel.data;
            check(wire(55)==1.,"mixed path actually reached crushed state");
            check(restoredMixed.recvSelf(11,channel,broker)==0,"receive frozen tensile history after crushing");
            restoredMixed.sendSelf(11,channel);equal(channel.data,wire);
            equal(restoredMixed.getStress(),source.getStress());equal(restoredMixed.getTangent(),source.getTangent());
        }
        std::cout<<"PASS "<<checks<<" adapter checks\n";
    } catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
