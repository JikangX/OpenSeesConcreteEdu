// Native ADINA model 5 plane-stress kernel, independent of OpenSees.
// Reference: Chinese ADINA/ADINAT manual, printed pp.105-108;
// A06.FOR: CDMOD, PRNCPL, CRAKID, DCRACK (ITYP2D >= 2).
#ifndef OPS_ADINA_CONCRETE2D_KERNEL_H
#define OPS_ADINA_CONCRETE2D_KERNEL_H

// Reuse the validated input and scalar compression/envelope laws of model 5.
// This is shared numerical infrastructure; the 3D integrator is never called.
#include "../adinaConcrete3D/Concrete3D.h"

namespace adina { namespace plane {
using Vector3 = std::array<double, 3>;
using Matrix3 = std::array<std::array<double, 3>, 3>;

/** Complete history of one integration point. All quantities are owned values.
 * Strain: [epsilon_xx, epsilon_yy, gamma_xy]; stress: [sigma_xx,sigma_yy,tau_xy].
 * angle (radians) fixes two perpendicular crack normals. Each normal has its
 * OWN active flag and formation strain: closing crack 1 must not close crack 2.
 * A closed crack retains its frame and can reopen. Crushing is irreversible.
 * No legacy angle offsets (+180, -361, -722) encode physical state here.
 */
struct State {
    Vector3 strain {}, stress {};
    double maximumLoading = 0., angle = 0.;
    bool hasCrackFrame = false, crushed = false;
    std::array<bool, 2> open {{false, false}};
    std::array<double, 2> crackStrain {};
    // Separate histories for the optional fracture-energy mode. The original
    // formation strains retain their ADINA meaning and are not repurposed.
    std::array<fracture::History,2> tension {};
    double crushingStress = 0., crushingStrain = 0.;
    // Most compressive in-plane principal strain reached since crushing (0
    // before). Used only by Parameters::CrushedUnloading::Elastic.
    double crushedStrain = 0.;
    // All stress released at the scaled ultimate strain (DCRACK, NUMCRK=3).
    bool released = false;
};

// Work-conjugate transformations. Engineering shear requires DIFFERENT
// stress and strain transformation matrices; using one for both loses energy.
Vector3 strainToLocal(const Vector3&, double angle);
Vector3 strainToGlobal(const Vector3&, double angle);
Vector3 stressToLocal(const Vector3&, double angle);
Vector3 stressToGlobal(const Vector3&, double angle);
double principalAngle(const Vector3&);
Matrix3 elasticTangent(double E, double nu);

/** Immutable committed -> trial update, with rotation-invariant substepping.
 * Optional tangent differentiates the COMPLETE replay at fixed substep count.
 * A failed update leaves trial/tangent unchanged. No shared scratch or I/O.
 */
bool integrate(const Parameters&, const State& committed, const Vector3& target,
               State& trial, Matrix3* tangent, std::string& error);
bool validateState(const State&, std::string& error);
// Parameters::Update::Hybrid origin; see adina::hybridOrigin.
State hybridOrigin(const State& committed, const State& iterate);
} } // namespace adina::plane
#endif
