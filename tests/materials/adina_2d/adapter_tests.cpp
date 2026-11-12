// Exercise the real OpenSees NDMaterial adapter and Channel serialization.
// No mock constitutive model or substitute OpenSees interface is used.
#include <AdinaConcrete2D.h>
#include <AdinaConcreteOptions.h>
#include <Channel.h>
#include <LayeredShellFiberSection.h>
#include <FEM_ObjectBrokerAllClasses.h>
#include <classTags.h>
#include <algorithm>
#include <FEM_ObjectBroker.h>
#include <StandardStream.h>
#include <OPS_Globals.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <deque>
#include <vector>
#include <PlateFromPlaneStressMaterial.h>
#include <PlateRebarMaterial.h>
#include <ElasticMaterial.h>
#include <stdexcept>


class MemoryChannel : public Channel {
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
// Keep message sequences for nested section -> plate -> material packets.
class PacketChannel final : public MemoryChannel {
public:
    std::deque<Vector> vectors;std::deque<Matrix> matrices;std::deque<ID> ids;
    int sendVector(int,int,const Vector& v,ChannelAddress*) override {vectors.push_back(v);return 0;}
    int recvVector(int,int,Vector& v,ChannelAddress*) override {if(vectors.empty())return -1;v=vectors.front();vectors.pop_front();return 0;}
    int sendMatrix(int,int,const Matrix& v,ChannelAddress*) override {matrices.push_back(v);return 0;}
    int recvMatrix(int,int,Matrix& v,ChannelAddress*) override {if(matrices.empty())return -1;v=matrices.front();matrices.pop_front();return 0;}
    int sendID(int,int,const ID& v,ChannelAddress*) override {ids.push_back(v);return 0;}
    int recvID(int,int,ID& v,ChannelAddress*) override {if(ids.empty())return -1;v=ids.front();ids.pop_front();return 0;}
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
    check(a.noRows()==b.noRows() && a.noCols()==b.noCols(),"matrix sizes");
    for(int i=0;i<a.noRows();++i)for(int j=0;j<a.noCols();++j)check(a(i,j)==b(i,j),"exact tangent restoration");
}
}

void near(double a,double b,double tol,const char* msg) {
    if(std::abs(a-b)>tol*std::max({1.,std::abs(a),std::abs(b)}))
        std::cerr<<msg<<": "<<a<<" vs "<<b<<'\n';
    check(std::abs(a-b)<=tol*std::max({1.,std::abs(a),std::abs(b)}),msg);
}
int main() {
  try {
    adina::Parameters p;p.rho=2.4;
    for(bool regularized:{false,true}) for(bool explicitTangent:{false,true}) for(bool explicitUpdate:{false,true}) {
      adina::Parameters options;
      if(regularized) {options.fractureEnergy=.0002;options.characteristicLength=.1;}
      adina::SolverOptionSelection selected;
      if(explicitTangent) selected.observe("-tangent");
      // An explicit "-update iterative" must survive the fracture defaults.
      if(explicitUpdate) {selected.observe("-update");options.update=adina::Parameters::Update::Iterative;}
      selected.applyDefaults(options);
      check(options.tangent==(regularized&&!explicitTangent ? adina::Parameters::Tangent::Algorithmic : adina::Parameters::Tangent::Secant),"fracture default / explicit tangent precedence");
      check(options.update==(explicitUpdate ? adina::Parameters::Update::Iterative : adina::Parameters::Update::Committed),"fracture default / explicit update precedence");
    }
    p.tangent=adina::Parameters::Tangent::Algorithmic;
    p.update=adina::Parameters::Update::Committed;
    for(bool plate:{false,true}) {
      AdinaConcrete2D m(42,p,plate,9000.);const int n=m.getOrder();Vector e(n);
      check(n==(plate ? 5 : 3),"formulation order");check(m.getRho()==2.4,"density");
      check(m.getCopy("ThreeDimensional")==nullptr,"reject 3D request");
      e(1)=-.0004;if(plate) {e(3)=.0003;e(4)=-.0002;}
      check(m.setTrialStrain(e)==0,"initial trial");m.commitState();
      const Vector cs=m.getStress(),ce=m.getStrain();const Matrix cc=m.getTangent();
      e(0)=.0002;e(1)=-.0006;if(plate) {e(3)=-.0008;e(4)=.0005;}
      check(m.setTrialStrain(e)==0,"cracked trial");
      const Vector ts=m.getStress(),te=m.getStrain();const Matrix tc=m.getTangent();
      if(plate) {near(ts(3),-7.2,1e-14,"yz shear");near(ts(4),4.5,1e-14,"xz shear");}
      m.revertToLastCommit();equal(m.getStress(),cs);equal(m.getStrain(),ce);equal(m.getTangent(),cc);
      m.setTrialStrain(e);equal(m.getStress(),ts);equal(m.getTangent(),tc);
      m.setDbTag(123);
      std::unique_ptr<NDMaterial> copy(m.getCopy());const Vector& borrowed=m.getStress();
      check(copy->getDbTag()==0,"clone Channel identity is independent");
      copy->revertToStart();equal(borrowed,ts);equal(m.getStress(),ts);
      MemoryChannel ch;FEM_ObjectBrokerAllClasses broker;
      check(m.sendSelf(1,ch)==0 && ch.data.Size()==118,"v3 serialized size");const Vector wire=ch.data;
      std::unique_ptr<NDMaterial> restored(broker.getNewNDMaterial(ND_TAG_AdinaConcrete2D));
      check(restored!=nullptr,"actual object broker registration");
      check(restored->recvSelf(1,ch,broker)==0,"checkpoint receive");
      equal(restored->getStress(),ts);equal(restored->getStrain(),te);equal(restored->getTangent(),tc);
      equal(restored->getInitialTangent(),m.getInitialTangent());
      check(restored->getOrder()==n && restored->getTag()==42,"serialized interface and tag");
      restored->sendSelf(1,ch);equal(wire,ch.data);
      ch.data=wire;ch.data.resize(97);ch.data(0)=1.;
      check(restored->recvSelf(1,ch,broker)==0,"backward compatible v1 checkpoint");
      equal(restored->getStress(),ts);equal(restored->getTangent(),tc);
      restored->revertToLastCommit();equal(restored->getStress(),cs);equal(restored->getStrain(),ce);
      restored->setTrialStrain(e);equal(restored->getStress(),ts);
      ch.data=wire;ch.data.resize(111);ch.data(0)=2.;
      check(restored->recvSelf(1,ch,broker)==0,"backward compatible v2 checkpoint");
      equal(restored->getStress(),ts);equal(restored->getTangent(),tc);
      for(int field:{0,1,2,3,5,51,52,53,54,79,97,98,99,100,101,111,112,113,116,117}) {
        ch.data=wire;ch.data(field)=field==79 ? NAN : (field==0 ? 9. : (field==3 ? -1. : 1.e100));
        check(restored->recvSelf(1,ch,broker)!=0,"reject corrupted wire field");equal(restored->getStress(),ts);
        equal(restored->getStrain(),te);equal(restored->getTangent(),tc);
      }
      ch.data=wire;ch.data(96)=NAN;check(restored->recvSelf(1,ch,broker)!=0,"reject wire NaN");
      for(int size:{0,96}) {
        ch.data=wire;ch.data.resize(size);
        check(restored->recvSelf(1,ch,broker)!=0,"reject incomplete wire vector");
        equal(restored->getStress(),ts);equal(restored->getStrain(),te);
      }
      Vector bad=e;bad(0)=NAN;check(m.setTrialStrain(bad)!=0,"reject strain NaN");
      equal(m.getStress(),ts);equal(m.getStrain(),te);equal(m.getTangent(),tc);
      check(m.setTrialStrain(Vector(n+1))!=0,"reject wrong order");
      if(plate) {bad=e;bad(3)=1e308;check(m.setTrialStrain(bad)!=0,"reject shear stress overflow");equal(m.getStress(),ts);}
      m.revertToStart();Vector zero(n);equal(m.getStress(),zero);equal(m.getStrain(),zero);
      Vector inc(n);inc(0)=1e-8;check(m.setTrialStrainIncr(inc)==0,"increment trial");
      check(m.setTrialStrainIncr(inc)==0,"second increment");near(m.getStrain()(0),2e-8,1e-16,"increment relative to trial");
    }
    // A native PlaneStress prototype can be assigned directly to LayeredShell.
    // getCopy("PlateFiber") must preserve its membrane histories, not initialize
    // a fresh concrete point or borrow the prototype's buffers.
    AdinaConcrete2D membrane(1,p);Vector e(3);e(0)=.0002;membrane.setTrialStrain(e);membrane.commitState();
    std::unique_ptr<NDMaterial> plate(membrane.getCopy("PlateFiber"));
    check(plate && plate->getOrder()==5,"direct native -> PlateFiber conversion");
    for(int i=0;i<3;++i) near(plate->getStress()(i),membrane.getStress()(i),1e-14,"copy membrane history");
    plate->revertToLastCommit();near(plate->getStrain()(0),.0002,1e-16,"copy committed history");

    // Real LayeredShellFiberSection, asymmetric thickness and heterogeneous
    // elastic properties. Verify every column of the eight-component section
    // map against independent integration of N, M, V and the layer centroids.
    adina::Parameters q=p;q.E=40000.;q.ft=4.;q.fc=-40.;q.fu=-8.;
    AdinaConcrete2D a(10,p),b(11,q,false,8000.);
    NDMaterial* mats[]={&a,&b,&a};double h[]={.02,.04,.03};const double total=.09;
    const double z[]={-.035,-.005,.030};
    LayeredShellFiberSection section(7,3,h,mats);near(section.getRho(),total*2.4,1e-14,"section areal density");
    const Matrix initialSection=section.getInitialTangent();
    for(int col=0;col<8;++col) {
      Vector sec(8);sec(col)=1e-8;check(section.setTrialSectionDeformation(sec)==0,"elastic section trial");
      const Vector actual=section.getStressResultant();const Matrix tangent=section.getSectionTangent();
      Vector expected(8);
      for(int l=0;l<3;++l) {
        const auto c=adina::plane::elasticTangent(l==1 ? q.E : p.E,p.nu);
        double eps[3],sig[3]={};for(int i=0;i<3;++i) eps[i]=sec(i)-z[l]*sec(i+3);
        for(int i=0;i<3;++i) for(int j=0;j<3;++j) sig[i]+=c[i][j]*eps[j];
        for(int i=0;i<3;++i) {expected(i)+=h[l]*sig[i];expected(i+3)+=h[l]*z[l]*sig[i];}
        const double G=l==1 ? 8000. : 12500.;
        expected(6)+=G*h[l]*sec(6);expected(7)+=G*h[l]*sec(7);
      }
      for(int i=0;i<8;++i) {
        near(actual(i),expected(i),1e-13,"elastic layered resultant");
        near(tangent(i,col),expected(i)/1e-8,2e-7,"section tangent mapping");
      }
    }
    section.revertToStart();std::array<adina::plane::State,3> history{};
    // Nonlinear membrane + curvature + independent transverse shear history.
    // Compare the real section against three independently advanced points;
    // commit only after all layers/resultants are checked.
    for(int step=1;step<=240;++step) {
      Vector sec(8);sec(0)=step*1e-6;sec(1)=-step*4e-6;sec(3)=step*2e-5;
      sec(6)=std::sin(step*.03)*1e-4;sec(7)=std::cos(step*.02)*2e-4;
      check(section.setTrialSectionDeformation(sec)==0,"nonlinear layered trial");Vector expected(8);
      for(int l=0;l<3;++l) {
        adina::plane::Vector3 eps{{sec(0)-z[l]*sec(3),sec(1),0.}};
        adina::plane::State next;std::string error;
        check(adina::plane::integrate(l==1?q:p,history[l],eps,next,nullptr,error),"independent layer update");history[l]=next;
        for(int i=0;i<3;++i) {expected(i)+=h[l]*next.stress[i];expected(i+3)+=h[l]*z[l]*next.stress[i];}
        const double G=l==1 ? 8000. : 12500.;expected(6)+=G*h[l]*sec(6);expected(7)+=G*h[l]*sec(7);
      }
      const Vector actual=section.getStressResultant();
      for(int i=0;i<8;++i) near(actual(i),expected(i),1e-10,"nonlinear section resultant");
      check(section.commitState()==0,"section commit");
    }
    const Vector committed=section.getStressResultant(),committedDeformation=section.getSectionDeformation();
    const Matrix kt=section.getSectionTangent();
    equal(section.getInitialTangent(),initialSection);
    equal(section.getStressResultant(),committed);
    equal(section.getSectionTangent(),kt);
    Vector large(8);large(0)=.001;large(6)=.0007;section.setTrialSectionDeformation(large);
    section.revertToLastCommit();equal(section.getStressResultant(),committed);equal(section.getSectionTangent(),kt);
    equal(section.getSectionDeformation(),committedDeformation);
    std::unique_ptr<SectionForceDeformation> clone(section.getCopy());equal(clone->getStressResultant(),committed);
    equal(clone->getSectionDeformation(),committedDeformation);
    section.revertToStart();equal(clone->getStressResultant(),committed);
    // Middle layer rejects a trial AFTER the first layer has accepted it.
    // The section must report failure and restore ALL committed quantities.
    auto guarded=p;guarded.maxStrainIncrement=1.e-8;
    AdinaConcrete2D guardedMaterial(12,guarded);
    NDMaterial* heterogeneous[]={&a,&guardedMaterial,&a};
    LayeredShellFiberSection guardedSection(8,3,h,heterogeneous);
    Vector small(8);small(1)=-5.e-8;
    check(guardedSection.setTrialSectionDeformation(small)==0,"guarded section seed");
    guardedSection.commitState();const Vector seedForce=guardedSection.getStressResultant();
    const Matrix seedTangent=guardedSection.getSectionTangent();
    Vector rejected(8);rejected(1)=-.0002;
    check(guardedSection.setTrialSectionDeformation(rejected)!=0,"middle layer rejection propagates");
    equal(guardedSection.getSectionDeformation(),small);
    equal(guardedSection.getStressResultant(),seedForce);equal(guardedSection.getSectionTangent(),seedTangent);
    // Five-component wrapper kinematics must survive commit, failed steps,
    // cloning and nested checkpoints, including both transverse shear strains.
    AdinaConcrete2D virgin(90,p);
    PlateFromPlaneStressMaterial wrapped(91,virgin,8000.);
    ElasticMaterial elasticBar(92,200000.);
    PlateRebarMaterial rebar(93,elasticBar,37.);
    for(NDMaterial* object:std::vector<NDMaterial*>{&wrapped,&rebar}) {
      Vector seed(5);seed(0)=-2.e-5;seed(1)=-3.e-5;seed(2)=1.e-5;seed(3)=.002;seed(4)=-.003;
      check(object->setTrialStrain(seed)==0,"wrapper seed");object->commitState();
      const Vector force=object->getStress();const Matrix stiffness=object->getTangent();
      Vector later(seed);later*=2;object->setTrialStrain(later);object->revertToLastCommit();
      equal(object->getStrain(),seed);equal(object->getStress(),force);equal(object->getTangent(),stiffness);
      std::unique_ptr<NDMaterial> duplicated(object->getCopy());equal(duplicated->getStrain(),seed);equal(duplicated->getStress(),force);
      PacketChannel packets;FEM_ObjectBrokerAllClasses actualBroker;
      check(object->sendSelf(10,packets)==0,"nested wrapper packet send");
      std::unique_ptr<NDMaterial> recovered(actualBroker.getNewNDMaterial(object->getClassTag()));
      check(recovered && recovered->recvSelf(10,packets,actualBroker)==0,"nested wrapper packet restore");
      equal(recovered->getStrain(),seed);equal(recovered->getStress(),force);equal(recovered->getTangent(),stiffness);
      check(recovered->getTag()==object->getTag(),"wrapper tag survives checkpoint");
      recovered->setTrialStrain(later);recovered->revertToLastCommit();equal(recovered->getStrain(),seed);
      // Previous CUDA wrappers used a positive nested class tag, four ID
      // entries and a six-value angle/G + strain payload. Preserve that wire
      // contract while checking the real nested material and broker.
      PacketChannel legacy;
      check(object->sendSelf(10,legacy)==0,"legacy wrapper fixture send");
      check(legacy.ids.front().Size()==4,"CUDA wrapper header layout");
      legacy.ids.front()(1)=-legacy.ids.front()(1);
      Vector oldPayload(6);
      for(int i=0;i<6;++i) oldPayload(i)=legacy.vectors.front()(i);
      legacy.vectors.front()=oldPayload;
      std::unique_ptr<NDMaterial> oldRecovered(actualBroker.getNewNDMaterial(object->getClassTag()));
      check(oldRecovered && oldRecovered->recvSelf(10,legacy,actualBroker)==0,"read previous CUDA wrapper packet");
      equal(oldRecovered->getStrain(),seed);equal(oldRecovered->getStress(),force);
      oldRecovered->setTrialStrain(later);oldRecovered->revertToLastCommit();equal(oldRecovered->getStrain(),seed);
      object->revertToStart();Vector zero(5);equal(object->getStrain(),zero);
    }
    // Optional histories must round-trip BOTH committed and uncommitted
    // damage, then produce the same future unloading/reopening response.
    auto gf=p;gf.nu=0.;gf.fractureEnergy=.0002;gf.characteristicLength=.1;
    for(bool plate:{false,true}) {
      AdinaConcrete2D source(181,gf,plate),recovered;
      Vector e(source.getOrder());e(0)=.0005;
      check(source.setTrialStrain(e)==0,"regularized damage trial");source.commitState();
      e(0)=.0007;check(source.setTrialStrain(e)==0,"regularized uncommitted damage");
      MemoryChannel channel;FEM_ObjectBrokerAllClasses broker;
      check(source.sendSelf(1,channel)==0,"regularized send");const Vector wire=channel.data;
      check(recovered.recvSelf(1,channel,broker)==0,"regularized receive");
      equal(source.getStress(),recovered.getStress());equal(source.getTangent(),recovered.getTangent());
      recovered.sendSelf(1,channel);equal(channel.data,wire);
      std::unique_ptr<NDMaterial> copy(source.getCopy());
      source.revertToLastCommit();recovered.revertToLastCommit();copy->revertToLastCommit();
      for(double strain:{.0002,0.,-.00002,0.,.0002,.0005,.0008}) {
        e(0)=strain;
        check(source.setTrialStrain(e)==0 && recovered.setTrialStrain(e)==0 && copy->setTrialStrain(e)==0,"regularized continuation");
        equal(source.getStress(),recovered.getStress());equal(source.getStress(),copy->getStress());
        equal(source.getTangent(),recovered.getTangent());
        source.commitState();recovered.commitState();copy->commitState();
      }
      for(int field:{97,98,99,100,101}) {
        channel.data=wire;channel.data(field)=field==97?0.:(field==99?-1.:1e100);
        Vector before=recovered.getStress();
        check(recovered.recvSelf(1,channel,broker)!=0,"reject inconsistent fracture checkpoint");equal(before,recovered.getStress());
      }
    }
    // Frozen tensile histories survive later crushing, even if a subsequent
    // crushed strain exceeds the tensile maximum. They must not make valid
    // v3 checkpoints unreadable. Exercise original and optional unloading.
    for(auto unloading:{adina::Parameters::CrushedUnloading::Original,adina::Parameters::CrushedUnloading::Elastic}) {
      auto mixed=p;mixed.nu=0.;mixed.fractureEnergy=.0002;mixed.characteristicLength=.1;mixed.crushedUnloading=unloading;
      AdinaConcrete2D source(198,mixed),recovered;Vector e(3);e(0)=.0005;
      check(source.setTrialStrain(e)==0,"form tensile history before crushing");source.commitState();
      e(1)=-.003;check(source.setTrialStrain(e)==0,"crush damaged 2D point");source.commitState();
      e(0)=.0015;check(source.setTrialStrain(e)==0,"advance beyond frozen tensile maximum");
      MemoryChannel channel;FEM_ObjectBrokerAllClasses broker;
      check(source.sendSelf(11,channel)==0,"send mixed tension / crush");const Vector wire=channel.data;
      check(wire(52)==1.,"mixed path actually reached crushed state");
      check(recovered.recvSelf(11,channel,broker)==0,"receive valid frozen tensile history after crushing");
      recovered.sendSelf(11,channel);equal(channel.data,wire);
      equal(recovered.getStress(),source.getStress());equal(recovered.getTangent(),source.getTangent());
    }
    std::cout<<"PASS "<<checks<<" adapter / layered section checks\n";
  } catch(const std::exception& ex) {std::cerr<<"FAIL: "<<ex.what()<<'\n';return 1;}
}
