// Independent checks of the Jiang Jianjing kernel (no Fortran, no OpenSees).
#include "JiangJianjing.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <iostream>

namespace {
int checks = 0;
void check(bool ok, const char* msg) {
    ++checks;
    if(!ok) {std::cerr << "FAIL: " << msg << '\n';std::exit(1);}
}
void near(double a, double b, double tol, const char* msg) {
    if(std::abs(a - b) > tol) std::cerr << msg << ": " << a << " vs " << b << '\n';
    check(std::abs(a - b) <= tol, msg);
}
}

int main() {
    using namespace jiangjj;
    // Lecture 7, slides 33-39: fc 20, ft 2 (ft/fc 0.1), E0 20 GPa, nu 0.18,
    // stress {-6 -6 -12 2 2 1}, strain {-1.62 -1.62 -5.80 2.78 2.78 1.39}e-4,
    // strain increment {-2 -1 -3 1 2 1}e-5. Slide values: sqrt(J2) 4.583,
    // J3 -2, theta 31.03 deg, sqrt(J2f) 8.998, beta 0.5093, Es 17.00 GPa,
    // stress {-6.53 -6.38 -12.69 2.07 2.07 1.07} (rounded slide arithmetic).
    Parameters p;p.fc = 20.;p.E0 = 20000.;p.nu = 0.18;
    const Vector6 s0 {{-6, -6, -12, 2, 2, 1}};
    const Index idx = nonlinearIndex(p, s0);
    near(idx.I1, -24., 1e-12, "lecture I1");
    near(std::sqrt(idx.J2), 4.583, 5e-4, "lecture sqrt(J2)");
    near(idx.J3, -2., 1e-12, "lecture J3");
    near(idx.theta * 45. / std::atan(1.), 31.03, 5e-3, "lecture theta");
    near(idx.sqrtJ2f, 8.998, 1e-3, "lecture sqrt(J2f)");
    near(idx.beta, 0.5093, 1e-4, "lecture beta");
    State st;st.stress = s0;st.strain = {{-1.62e-4, -1.62e-4, -5.80e-4, 2.78e-4, 2.78e-4, 1.39e-4}};
    Vector6 e1;const double de[6] = {-2e-5, -1e-5, -3e-5, 1e-5, 2e-5, 1e-5};
    for(int i = 0;i < 6;++i) e1[i] = st.strain[i] + de[i];
    Matrix6 c;
    const State s1 = update(p, st, e1, &c);
    near(s1.secantModulus, 17000., 5., "lecture Es");
    // Slide 39 prints 2.07 for the yz component, but its own data give
    // G (2.78e-4 + 2e-5) = 7.20e3 x 2.98e-4 = 2.15: the slide used 1e-5 for the
    // yz strain increment (arithmetic slip on the slide). All other components
    // agree with the slide.
    const double slide[6] = {-6.53, -6.38, -12.69, 2.07, 2.147, 1.07};
    for(int i = 0;i < 6;++i) near(s1.stress[i], slide[i], 0.02, "lecture stress");
    std::printf("INFO lecture example: beta %.4f sqrtJ2f %.4f Es %.1f stress %.3f %.3f %.3f %.3f %.3f %.3f\n",
                idx.beta, idx.sqrtJ2f, s1.secantModulus, s1.stress[0], s1.stress[1], s1.stress[2],
                s1.stress[3], s1.stress[4], s1.stress[5]);

    // Failure surface: uniaxial compression and tension strengths implied by
    // the criterion constants (beta = 1 exactly at sigma = -fc, and at ft).
    Parameters q; // fc 30
    near(nonlinearIndex(q, {{0, 0, -30, 0, 0, 0}}).beta, 1.0, 2e-4, "uniaxial compressive strength");
    {
        // find the uniaxial tensile strength by bisection
        double lo = 0., hi = 10.;
        for(int k = 0;k < 200;++k) {
            const double m = .5 * (lo + hi);
            (nonlinearIndex(q, {{0, 0, m, 0, 0, 0}}).beta < 1. ? lo : hi) = m;
        }
        near(lo / q.fc, 0.1, 2e-3, "tensile strength ft/fc = 0.1");
        std::printf("INFO uniaxial tensile strength %.4f MPa (fc 30)\n", lo);
    }
    // Objectivity: beta is invariant under rotations of the stress state.
    {
        const Vector6 s {{-12, -4, 3, 2.5, -1.5, 4}};
        const double b0 = nonlinearIndex(q, s).beta;
        for(int k = 0;k < 50;++k) {
            const double a = 0.13 * k, ca = std::cos(a), sa = std::sin(a);
            // rotation about z
            double S[3][3] = {{s[0], s[3], s[5]}, {s[3], s[1], s[4]}, {s[5], s[4], s[2]}};
            double R[3][3] = {{ca, -sa, 0}, {sa, ca, 0}, {0, 0, 1}}, T[3][3] = {}, U[3][3] = {};
            for(int i = 0;i < 3;++i) for(int j = 0;j < 3;++j) for(int m = 0;m < 3;++m) T[i][j] += R[i][m] * S[m][j];
            for(int i = 0;i < 3;++i) for(int j = 0;j < 3;++j) for(int m = 0;m < 3;++m) U[i][j] += T[i][m] * R[j][m];
            const Vector6 r {{U[0][0], U[1][1], U[2][2], U[0][1], U[1][2], U[0][2]}};
            near(nonlinearIndex(q, r).beta, b0, 1e-12, "beta objectivity");
        }
    }
    // Total (secant) model: from a zero state the response is linear elastic;
    // a crushed point carries no stress and keeps a tiny stiffness.
    {
        Matrix6 t;
        const State a = update(q, State{}, {{1e-5, -2e-5, 3e-5, 1e-5, 0, 0}}, &t);
        const Matrix6 D = elasticStiffness(q.E0, q.nu);
        for(int i = 0;i < 6;++i) {
            double s = 0.;for(int j = 0;j < 6;++j) s += D[i][j] * a.strain[j];
            near(a.stress[i], s, 1e-12, "first step linear elastic");
        }
        State cr;cr.strain = {{0, 0, -0.006, 0, 0, 0}}; // eps_zz - eps_m = -0.004 < -0.0033
        const State z = update(q, cr, {{0, 0, -0.0061, 0, 0, 0}}, &t);
        check(z.crushed, "crushing test on start-of-increment strain");
        for(double x:z.stress) check(x == 0., "crushed stress zero");
        near(t[0][0], 1e-6 * D[0][0], 1e-12, "crushed stiffness 1e-6 E0");
    }
    // Repeated trials from the same committed state are identical (reentrant).
    {
        State c0;c0.stress = {{-8, -2, -1, 1, 0, 0}};c0.strain = {{-3e-4, -1e-4, 0, 1e-4, 0, 0}};
        const Vector6 e {{-3.5e-4, -1e-4, 0, 1e-4, 0, 0}};
        const State x = update(q, c0, e, nullptr), y = update(q, c0, e, nullptr);
        for(int i = 0;i < 6;++i) check(x.stress[i] == y.stress[i], "reentrant update");
    }
    // -tangent secant returns the exact derivative of the update (central
    // differences from the same committed state); the original tangent does not.
    {
        Parameters s = q;s.tangent = Parameters::Tangent::Secant;
        State c0;c0.stress = {{-14, -3, -2, 2, -1, 1.5}};c0.strain = {{-5e-4, -1e-4, -5e-5, 1.5e-4, -7e-5, 1e-4}};
        const Vector6 e {{-5.5e-4, -1.1e-4, -6e-5, 1.6e-4, -8e-5, 1.2e-4}};
        Matrix6 ts, to;
        update(s, c0, e, &ts);update(q, c0, e, &to);
        double maxSecant = 0., maxOriginal = 0.;
        for(int j = 0;j < 6;++j) {
            Vector6 ep = e, em = e;const double h = 1e-9;ep[j] += h;em[j] -= h;
            const State a = update(s, c0, ep, nullptr), b = update(s, c0, em, nullptr);
            for(int i = 0;i < 6;++i) {
                const double fd = (a.stress[i] - b.stress[i]) / (2 * h);
                maxSecant = std::max(maxSecant, std::abs(fd - ts[i][j]) / q.E0);
                maxOriginal = std::max(maxOriginal, std::abs(fd - to[i][j]) / q.E0);
            }
        }
        check(maxSecant < 1e-6, "secant tangent is the exact derivative of the update");
        check(maxOriginal > 1e-2, "original tangent differs from the derivative");
        std::printf("INFO tangent vs finite difference: secant %.1e, original %.3f (relative to E0)\n",
                    maxSecant, maxOriginal);
    }
    // -beta current (option): the secant law sigma = D(Es(beta(sigma))) eps is
    // satisfied exactly, independent of the path and of the step size before the
    // ultimate stress; the explicit original update converges to it as the steps
    // shrink; the secant tangent is the exact derivative.
    {
        Parameters ic = q;ic.betaFrom = Parameters::BetaFrom::Current;
        Parameters ics = ic;ics.tangent = Parameters::Tangent::Secant;
        const Vector6 e {{-5.0e-4, -1.2e-4, 0.6e-4, 1.4e-4, -0.8e-4, 1.1e-4}};
        const State one = update(ic, State{}, e, nullptr);
        const double beta = nonlinearIndex(ic, one.stress).beta;
        near(one.secantModulus / ic.E0, .5 * (1. + std::sqrt(1. - beta)), 1e-12, "current beta: secant law satisfied");
        near(one.beta, beta, 1e-12, "current beta: beta of the returned stress");
        check(beta > 0.5 && beta < 0.99, "current beta: test point in the nonlinear range");
        // two different paths in 37 steps
        State x, y;
        for(int k = 1;k <= 37;++k) {
            Vector6 a, b;
            for(int i = 0;i < 6;++i) {
                a[i] = e[i] * k / 37.;
                b[i] = e[i] * std::pow(k / 37., i % 2 ? 2. : .5);
            }
            x = update(ic, x, a, nullptr);y = update(ic, y, b, nullptr);
        }
        for(int i = 0;i < 6;++i) {
            near(x.stress[i], one.stress[i], 1e-9, "current beta: step-size independent");
            near(y.stress[i], one.stress[i], 1e-9, "current beta: path independent");
        }
        // zero increment reproduces the converged stress
        const State z = update(ic, one, e, nullptr);
        for(int i = 0;i < 6;++i) near(z.stress[i], one.stress[i], 1e-10, "current beta: zero increment");
        // explicit original update with n steps -> current-beta result
        double prevErr = 1e300;
        for(int n:{10, 100, 1000}) {
            State w;
            for(int k = 1;k <= n;++k) {Vector6 a;for(int i = 0;i < 6;++i) a[i] = e[i] * k / n;w = update(q, w, a, nullptr);}
            double er = 0.;
            for(int i = 0;i < 6;++i) er = std::max(er, std::abs(w.stress[i] - one.stress[i]));
            check(er < prevErr, "explicit update converges to the current-beta solution");
            std::printf("INFO explicit (committed beta) with %d steps vs current beta: max stress difference %.3e MPa\n", n, er);
            prevErr = er;
        }
        check(prevErr < 0.02, "explicit update with 1000 steps close to current beta");
        // exact tangent (finite differences), before and on the failure surface
        const Vector6 et {{2.0e-4, -0.3e-4, -0.2e-4, 0.5e-4, 0, 0}}; // beyond ft
        for(const Vector6& es:{e, et}) {
            Matrix6 t;
            const State st = update(ics, State{}, es, &t);
            double mx = 0.;
            for(int j = 0;j < 6;++j) {
                Vector6 ep = es, em = es;const double h = 1e-9;ep[j] += h;em[j] -= h;
                const State a = update(ics, State{}, ep, nullptr), b = update(ics, State{}, em, nullptr);
                for(int i = 0;i < 6;++i) mx = std::max(mx, std::abs((a.stress[i] - b.stress[i]) / (2 * h) - t[i][j]) / q.E0);
            }
            std::printf("INFO current beta: tangent vs finite difference %.1e (beta %.4f)\n", mx, st.beta);
            check(mx < 2e-5, "current beta: secant tangent is the derivative of the update");
        }
        // on the failure surface: beta = 1 and held from the next increment on
        {
            const Vector6 es = et;
            const State f = update(ic, State{}, es, nullptr);
            near(nonlinearIndex(ic, f.stress).beta, 1.0, 1e-9, "current beta: stress on the failure surface");
            Vector6 e2 = es;e2[0] *= 1.1;
            const State g = update(ic, f, e2, nullptr);
            for(int i = 0;i < 6;++i) check(g.stress[i] == f.stress[i], "current beta: ultimate stress held");
        }
        // crushing unchanged
        State cr;cr.strain = {{0, 0, -0.006, 0, 0, 0}};
        check(update(ic, cr, {{0, 0, -0.0061, 0, 0, 0}}, nullptr).crushed, "current beta: crushing unchanged");
    }
    // Tension beyond the apex of the failure surface (I1 >= fc/d, mean stress
    // >= fc/(3d) = 0.104 fc). Original defect D2: the quadratic for sqrt(J2f)
    // has only a non-positive root, Get_DS gets beta <= 0 and Es >= E0, and
    // (near-)hydrostatic tension never fails (45 MPa for fc 30 in 4 coarse
    // steps). Here such states count as failed (beta = 1).
    {
        const double apex = q.fc / (3. * q.d); // hydrostatic tensile strength (mean stress)
        const Vector6 beyond[] = {{{4, 4, 4, 0, 0, 0}}, {{4, 4, 4.5, 0, 0, 0}}, {{3, 3, 6, 0, 0, 0}}, {{6, 6, 6, 1, 0, 0}}};
        for(const Vector6& s:beyond) {
            const Index x = nonlinearIndex(q, s);
            check(!x.solvable && x.beta == 1., "beyond the tensile apex: beta = 1");
        }
        check(nonlinearIndex(q, {{2, 2, 2, 0, 0, 0}}).beta == 0., "hydrostatic tension below the apex: beta = 0");
        check(nonlinearIndex(q, {{2.9, 2.9, 0, 0, 0, 0}}).solvable, "biaxial tension below failure: solvable");
        double worst = 0.;
        for(const double dz:{1.0, 1.1}) // exactly hydrostatic, and near-hydrostatic
        for(const int n:{60, 4}) // fine and coarse steps
        for(const Parameters::BetaFrom bf:{Parameters::BetaFrom::Committed, Parameters::BetaFrom::Current}) {
            Parameters t = q;t.betaFrom = bf;
            State pt;
            double held = 0.;
            bool secantOk = true, heldOk = true;
            for(int k = 1;k <= 3 * n;++k) {
                const double x = 3e-4 * k / n;
                pt = update(t, pt, {{x, x, dz * x, 0, 0, 0}}, nullptr);
                secantOk = secantOk && pt.secantModulus <= t.E0 * (1. + 1e-12);
                const double sm = (pt.stress[0] + pt.stress[1] + pt.stress[2]) / 3.;
                if(k == 3 * n / 2) held = sm;
                if(k > 3 * n / 2) heldOk = heldOk && sm == held;
            }
            check(secantOk, "tension beyond the apex: Es <= E0 (finite)");
            check(heldOk, "tension beyond the apex: stress held");
            const double sm = (pt.stress[0] + pt.stress[1] + pt.stress[2]) / 3.;
            worst = std::max(worst, sm);
            // committed beta: at most one increment beyond the surface (mean
            // stress = 3 K x mean strain, K the bulk modulus); current beta:
            // on the surface
            const double bulk = t.E0 / (3. * (1. - 2. * t.nu));
            check(sm < apex + 3. * bulk * 3.3e-4 / n + 1e-9, "tension beyond the apex: mean stress bounded");
        }
        std::printf("INFO tension beyond the apex: largest held mean stress %.3f MPa (apex %.3f MPa, fc 30)\n", worst, apex);
    }
    std::string err;
    Parameters bad = q;bad.fc = -30.;check(!bad.validate(err), "reject negative fc");
    bad = q;bad.epsCrush = 0.0033;check(!bad.validate(err), "reject positive epsCrush");
    bad = q;bad.nu = 0.5;check(!bad.validate(err), "reject nu 0.5");
    std::printf("PASS %d Jiang Jianjing kernel checks\n", checks);
    return 0;
}
