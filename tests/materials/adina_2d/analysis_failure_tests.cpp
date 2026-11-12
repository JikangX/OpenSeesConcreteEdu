// Integration regression for a rejected Newton line-search trial.
// The domain, yielding truss, assembly, solver and analysis are production
// OpenSees objects. Only the search strategy injects a deterministic rejection.
// Before the return-value repair the analysis ignored it and committed the load.
#include <Domain.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <LoadPattern.h>
#include <LinearSeries.h>
#include <NodalLoad.h>
#include <Truss.h>
#include <Steel01.h>
#include <PlainHandler.h>
#include <PlainNumberer.h>
#include <AnalysisModel.h>
#include <FullGenLinSOE.h>
#include <FullGenLinLapackSolver.h>
#include <LoadControl.h>
#include <CTestNormUnbalance.h>
#include <NewtonLineSearch.h>
#include <LineSearch.h>
#include <StaticAnalysis.h>
#include <Vector.h>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
class RejectingSearch final : public LineSearch {
public:
    int calls = 0;
    RejectingSearch(): LineSearch(0) {}
    int newStep(LinearSOE&) override { return 0; }
    int search(double, double, LinearSOE&, IncrementalIntegrator&) override {
        ++calls;
        return -1;
    }
    void Print(OPS_Stream&, int = 0) override {}
    int sendSelf(int, Channel&) override { return -1; }
    int recvSelf(int, Channel&, FEM_ObjectBroker&) override { return -1; }
};
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
}

int main() {
    try {
        Domain domain;
        require(domain.addNode(new Node(1, 2, 0., 0.)), "node 1");
        require(domain.addNode(new Node(2, 2, 1., 0.)), "node 2");
        require(domain.addSP_Constraint(new SP_Constraint(1, 0, 0., true)), "support x");
        require(domain.addSP_Constraint(new SP_Constraint(1, 1, 0., true)), "support y");
        require(domain.addSP_Constraint(new SP_Constraint(2, 1, 0., true)), "guide y");
        Steel01 steel(1, 10., 200000., .01);
        require(domain.addElement(new Truss(1, 2, 1, 2, steel, 1.)), "truss");
        auto* pattern = new LoadPattern(1);
        pattern->setTimeSeries(new LinearSeries(1));
        require(domain.addLoadPattern(pattern), "load pattern");
        Vector load(2); load(0) = 50.;
        require(domain.addNodalLoad(new NodalLoad(1, 2, load), 1), "nodal load");

        auto* handler = new PlainHandler;
        auto* numberer = new PlainNumberer;
        auto* model = new AnalysisModel;
        auto* solver = new FullGenLinLapackSolver;
        auto* system = new FullGenLinSOE(*solver);
        auto* integrator = new LoadControl(1., 1, 1., 1.);
        auto* test = new CTestNormUnbalance(1.e-10, 20, 0);
        auto* search = new RejectingSearch;
        auto* algorithm = new NewtonLineSearch(*test, search);
        StaticAnalysis analysis(domain, *handler, *numberer, *model,
                                *algorithm, *system, *integrator, test);
        const int result = analysis.analyze(1);
        // StaticAnalysis must reject and revert the entire load step. Checking
        // the time alone is insufficient: the nodal displacement and material
        // history must also return to their last committed (zero) state.
        require(search->calls == 1, "line search was not exercised exactly once");
        require(result < 0, "rejected line search incorrectly accepted");
        require(domain.getCurrentTime() == 0., "failed step committed its time");
        require(domain.getNode(2)->getDisp()(0) == 0., "failed step committed displacement");
        auto* truss = static_cast<Truss*>(domain.getElement(1));
        require(truss->getMaterialStrain() == 0., "failed step committed material strain");
        require(truss->getMaterialStress() == 0., "failed step committed material stress");
        analysis.clearAll(); // Analysis explicitly owns these strategy objects.
        std::cout << "PASS production line-search rejection / full-domain rollback\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
