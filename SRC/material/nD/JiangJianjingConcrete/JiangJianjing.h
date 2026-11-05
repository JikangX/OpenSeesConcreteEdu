// Jiang Jianjing concrete: nonlinear-elastic total (secant) E-nu model with
// the four-parameter Jiang failure criterion and the sqrt(J2) nonlinear index.
// Teaching model of the Tsinghua graduate course "Finite element analysis of
// reinforced concrete" (Lu Xinzheng). Source: JiangJJ2007.for (MSC.Marc
// HYPELA2), routine Get_DS; theory: lecture 7, slides 25-39 and 66-77.
//
// Independent, reentrant material-point kernel: no OpenSees types, no global
// state, no I/O. The OpenSees adapter is JiangJianjingConcrete.{h,cpp}.
#ifndef OPS_JIANG_JIANJING_KERNEL_H
#define OPS_JIANG_JIANJING_KERNEL_H

#include <array>
#include <string>

namespace jiangjj {
using Vector6 = std::array<double, 6>;
using Matrix6 = std::array<std::array<double, 6>, 6>;

/** Material input. Voigt order [xx, yy, zz, xy, yz, zx]; shear STRAINS are
 * engineering (gamma), shear STRESSES are tensor components (as in Marc).
 * Compression is negative for stresses and strains; fc is entered POSITIVE
 * (C%fc=30 in the original program), epsCrush NEGATIVE (C%EPS_Crush=-0.0033).
 */
struct Parameters {
    double fc = 30., E0 = 30000., nu = 0.18, epsCrush = -0.0033;
    // Jiang four-parameter criterion
    //   a J2/fc^2 + (b + c cos(theta)) sqrt(J2)/fc + d I1/fc - 1 = 0,
    // the constants of the original program (ft/fc = 0.1, lecture slide 35).
    double a = 1.2856, b = 1.4268, c = 10.2551, d = 3.2128;
    double rho = 0.;
    /** Matrix returned to the solver (never changes the stress update).
     * Original (default, as Get_DS): D(Et, nu) + 1e-6 D(E0, nu), with the
     * tangent modulus Et = E0 sqrt(1-beta) of the equivalent uniaxial curve.
     * Secant: D(Es, nu), the exact derivative of this update (beta is taken
     * from the converged stress, so the trial stress is linear in the trial
     * strain). With Original, full Newton diverges once Es/Et > 2
     * (beta > 8/9); use ModifiedNewton -initial or the Secant matrix.
     */
    enum class Tangent {Original, Secant};
    Tangent tangent = Tangent::Original;
    /** Stress from which the nonlinear index beta is evaluated.
     * Committed (default, as Get_DS): beta of the stress at the START of the
     * increment (explicit; the result depends on the step size and a
     * re-evaluation at an unchanged strain changes the stress).
     * Current (option, not in the original program): beta of the stress
     * being computed, i.e. the total secant law of the lecture
     * sigma = D(Es(beta(sigma)), nu) eps solved exactly at every call
     * (scalar equation for Es/E0 on the ray of D(E0) eps); path independent
     * before the ultimate stress. A point whose trial reaches the failure
     * surface gets the stress on the surface on that ray (beta = 1); from the
     * next increment on it is held exactly as in the original program, and
     * crushing is unchanged. With Tangent::Secant the matrix is the exact
     * (unsymmetric) derivative of this update.
     */
    enum class BetaFrom {Committed, Current};
    BetaFrom betaFrom = BetaFrom::Committed;
    bool validate(std::string& error) const;
};

/** Converged or trial state. beta is the nonlinear index used in the last
 * update; betaMax its largest value so far (recorded, but - exactly as in the
 * original program - not used: the secant modulus follows the current beta).
 */
struct State {
    Vector6 strain {}, stress {};
    bool crushed = false;
    double beta = 0., betaMax = 0.;
    double secantModulus = 0., tangentModulus = 0.;
    bool indexSolvable = true;
};

/** Nonlinear index of a stress state: beta = sqrt(J2)/sqrt(J2f), where J2f is
 * the J2 of the failure surface at the same I1 and Lode angle theta
 * (Jiang's sqrt(J2) method, slide 25). solvable=false when the quadratic for
 * sqrt(J2f) has no positive root, beta = 1 then: no real root (Get_DS writes
 * "Wrong" and uses beta = 1), or tension beyond the apex of the failure
 * surface, I1 >= fc/d (Get_DS returns beta <= 0 there, original defect D2).
 */
struct Index {
    double I1 = 0., J2 = 0., J3 = 0., theta = 0., sqrtJ2f = 0., beta = 0.;
    bool solvable = true;
};
Index nonlinearIndex(const Parameters&, const Vector6& stress);

/** Crushing test of Get_DS on a strain state: with em = min(mean strain, 0),
 * crushed when min_i(eps_ii - em) < epsCrush (normal components only).
 */
bool crushes(const Parameters&, const Vector6& strain);

/** Isotropic stiffness in engineering-shear Voigt form (Get_D). */
Matrix6 elasticStiffness(double E, double nu);

/** One update exactly as Get_DS, from the converged state `committed` to the
 * total strain `strain` (any number of trials may be evaluated from the same
 * committed state; this is what Marc does within an increment):
 *  1. crushed (history, or the crushing test on the COMMITTED strain):
 *     stress 0, tangent 1e-6 E0 elastic matrix;
 *  2. beta from the COMMITTED stress (beta = 1 if not solvable), cut at 1;
 *  3. Es = E0/2 (1 + sqrt(1-beta)), Et = E0 sqrt(1-beta), nu_s = nu;
 *  4. stress = D(Es, nu) * strain if beta < 0.99999, else the committed
 *     stress (ultimate stress reached: stress kept);
 *  5. tangent = D(Et, nu) + 1e-6 D(E0, nu) (Tangent::Original), or D(Es, nu)
 *     before the ultimate stress with Tangent::Secant.
 * The update never fails. Returns the tangent in `tangent` if non-null.
 */
State update(const Parameters&, const State& committed, const Vector6& strain, Matrix6* tangent);
/** The BetaFrom::Current update (called by update(); public for tests). */
State updateCurrent(const Parameters&, const State& committed, const Vector6& strain, Matrix6* tangent);
} // namespace jiangjj
#endif
