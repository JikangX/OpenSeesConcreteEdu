// ADINA concrete model 5: independent, reentrant material-point kernel.
// Theory: ADINA/ADINAT Chinese manual, printed pp. 105-108 and 129-132.
// Source correspondence: A16.FOR, CMOD3D / PRNCP3 / CRAK3D / DCRAK3.
// No OpenSees types, global variables, I/O, or shared scratch buffers occur here.
#ifndef OPS_ADINA_CONCRETE3D_KERNEL_H
#define OPS_ADINA_CONCRETE3D_KERNEL_H

#include <array>
#include <string>
#include "TensionSoftening.h"

namespace adina {
using Vector6 = std::array<double, 6>;
using Matrix6 = std::array<std::array<double, 6>, 6>;
using Frame = std::array<std::array<double, 3>, 3>;

/** Material input, in a consistent stress/strain/unit system.
 * Tensile stress is positive; fc, epsc, fu, epsu MUST be negative.
 * The four six-point arrays are the manual's SP1, SP31, SP32, SP33.
 * They are dimensionless, normalized by the negative uniaxial strength fc.
 * Keeping the envelope in the input makes confinement a constitutive property,
 * rather than hard-coding a particular concrete strength into the integrator.
 */
struct Parameters {
    double E = 30000., nu = 0.2, ft = 3., fc = -30.;
    double epsc = -0.002, fu = -6., epsu = -0.006;
    double beta = 0.5, gamma = 1., kappa = 0.7, alpha = 0.;
    double normalRetention = 1.e-4, shearRetention = 0.5, rho = 0.;
    // Optional subdivision of a trial increment: maximum Frobenius norm of the
    // strain tensor increment per substep (rotation invariant; a numerical
    // control, not a material parameter). 0 (default) = no subdivision: one
    // update per increment, as CDMOD/CMOD3D are called once per element
    // integration point with the total increment from the converged state.
    double maxStrainIncrement = 0.;
    // Optional extension, disabled by default. Gf has units stress*length;
    // Lch must be supplied for this material point's actual crack bandwidth.
    double fractureEnergy = 0., characteristicLength = 0.;
    bool regularizedTension() const {return fractureEnergy>0.;}
    /** Matrix returned to the global solver; it never changes the stress map.
     * Secant: CDMOD/CMOD3D-type material stiffness of the trial state (crack
     * retention, unloading modulus, STIFAC after crushing), symmetric and
     * positive semidefinite, with the cohesive secant of a softening crack.
     * Algorithmic: finite-difference Jacobian of the complete trial replay.
     * It is exact on smooth branches but costs one replay per strain column
     * and has no meaning across the brittle stress jumps of the original law.
     */
    enum class Tangent {Secant, Algorithmic};
    Tangent tangent = Tangent::Secant;
    /** Origin of each trial update; used by the OpenSees adapters only.
     * Committed (default, ADINA): every trial integrates the total increment
     * from the last converged state. This is what the original program does:
     * in UNBLD each equilibrium iteration reads the element group data
     * (stresses, strains, history WA) from unit IT2 and never writes them back;
     * the history is stored only when the stiffness of the next step is formed
     * from the converged increment. Verified on the compiled all_adina.for with
     * WRITEs in CDMOD (B07 deck: all iterations of a step start from the same
     * converged stresses and EVMAX) and by running the same finite-element
     * models in ADINA and OpenSees (tests/adinaConcreteReview/ADINA84). A brittle
     * crack can appear and disappear between iterations; when the equilibrium
     * point falls inside the stress jump at ft, Newton may cycle. ADINA does
     * the same (it stops with the iteration limit).
     * Iterative (OPTIONAL, not ADINA): each trial starts from the previous
     * iteration, so crack and crushing events met during the iterations of a
     * step persist within that step. The result depends on the iteration path.
     * Hybrid (OPTIONAL, not ADINA): as Iterative, but while a point has no new
     * event in the step its loading/unloading reference (maximum of the
     * loading function, cohesive damage maximum) is the converged value. In
     * Iterative mode overshoot/return iterates unload and reload with E and
     * ratchet the compressive stress above the curve (up to +36% at half the
     * peak strain for overshoot/return iterates of 30%); Hybrid removes this,
     * and with monotonic iterates it equals Iterative.
     * All modes restart from the committed state after revertToLastCommit.
     */
    enum class Update {Iterative, Committed, Hybrid};
    Update update = Update::Committed;
    /** OPTIONAL extension, not an ADINA equation (default Original).
     * Original (CDMOD/CMOD3D, PGRAV=100): a crushed point always follows the
     * softening-curve tangent, also when its strain decreases. Unloading thus
     * climbs back up the curve and, past the peak strain, the modulus becomes
     * 0: the point keeps the full peak compressive stress with zero stiffness,
     * even at tensile strain. Under cyclic loading the locked stresses of a
     * crushed wall toe dominate the response.
     * Elastic: below the most compressive principal strain reached since
     * crushing, the point unloads/reloads with the initial modulus E (ADINA's
     * unloading rule for uncrushed concrete) and its principal stresses are
     * capped at zero (no tension). Loading beyond that strain, and the release
     * of all stress at the scaled ultimate strain, follow the original law.
     */
    enum class CrushedUnloading {Original, Elastic};
    CrushedUnloading crushedUnloading = CrushedUnloading::Original;
    std::array<double, 6> sp1 {{0., .25, .5, .75, 1., 1.2}};
    std::array<double, 6> sp31 {{1., 1.4, 1.7, 2.2, 2.5, 2.8}};
    std::array<double, 6> sp32 {{1.3, 1.5, 2., 2.3, 2.7, 3.2}};
    std::array<double, 6> sp33 {{1.25, 1.45, 1.95, 2.25, 2.65, 3.15}};
    bool validate(std::string& error) const;
};

/** Complete committed OR trial state, with value semantics.
 * Voigt order is OpenSees [xx, yy, zz, xy, yz, xz]; shear strains are engineering
 * gamma_ij, whereas shear stresses are tensor components tau_ij.
 * Frame columns are the crack axes in global coordinates. All nine entries are
 * retained: ADINA's four direction cosines plus +/-51/100 sentinels lose signs
 * for general 3D rotations. closed cracks retain their frame/history.
 */
struct State {
    Vector6 strain {}, stress {};
    double maximumLoading = 0.;
    bool hasCrackFrame = false, crushed = false;
    int openCracks = 0;
    Frame frame {{{{1.,0.,0.}}, {{0.,1.,0.}}, {{0.,0.,1.}}}};
    std::array<double, 3> crackStrain {};
    std::array<fracture::History,3> tension {};
    double crushingStress = 0., crushingStrain = 0.;
    // Most compressive principal strain reached since crushing (0 before).
    // Used only by CrushedUnloading::Elastic.
    double crushedStrain = 0.;
};

struct Principal {
    std::array<double, 3> values {};
    Frame frame {};
};
struct Strength {
    double compression = 0., peakStrain = 0., tension = 0.;
};

// Small pure numerical kernels exposed for independent verification.
Principal principalStress(const Vector6& stress);
Vector6 strainToLocal(const Vector6& strain, const Frame& frame);
Vector6 strainToGlobal(const Vector6& strain, const Frame& frame);
Vector6 stressToLocal(const Vector6& stress, const Frame& frame);
Vector6 stressToGlobal(const Vector6& stress, const Frame& frame);
Matrix6 elasticTangent(double E, double nu);
Strength strengthEnvelope(const Parameters&, std::array<double, 3> stresses);
// Rational compression / initial linear tensile curve before crack release.
// Stress and tangent agree on both sides of zero strain.
double compressionStress(const Parameters&, const Strength&, double strain);
double compressionTangent(const Parameters&, const Strength&, double strain);
double loadingMeasure(const Parameters&, const Vector6& stress);

/** Replay one prescribed total-strain trial from an immutable committed state.
 * A failure leaves the caller's output untouched. Finite differences replay the
 * SAME committed state, including crack formation, closure and all substeps.
 * This is the algorithmic tangent of the implemented stress map, not the old
 * solver's positive replacement stiffness on a descending branch.
 */
bool integrate(const Parameters&, const State& committed, const Vector6& strain,
               State& trial, Matrix6* tangent, std::string& error);
/** Origin of the next trial update for Parameters::Update::Hybrid: the
 * previous iterate, with the converged loading-function maximum and cohesive
 * damage maxima when no crack/crushing event occurred since the commit.
 */
State hybridOrigin(const State& committed, const State& iterate);
} // namespace adina
#endif
