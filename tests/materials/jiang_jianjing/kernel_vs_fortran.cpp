// Material-point driver of the C++ kernel with the same protocol as
// oracle_driver.f90 (one converged update per increment from the previous
// converged state). stdin: fc E0 nu epsCrush, then one total strain row
// (6 components, engineering shear) per increment.
// stdout per increment: stress(6) crushed betaMax D11 D12 D44 0
// argv[1] == "fortran-literals": use the criterion constants as the Fortran
// program actually evaluates them. JiangJJ2007.for writes 1.2856, 1.4268,
// 10.2551, 3.2128 as default-kind (REAL*4) literals, so they enter the
// REAL*8 arithmetic rounded to single precision (relative error ~1e-8).
#include "JiangJianjing.h"
#include <cstdio>
#include <iostream>

#include <cstring>
int main(int argc, char** argv) {
    jiangjj::Parameters p;
    if(argc > 1 && !std::strcmp(argv[1], "fortran-literals")) {
        p.a = double(1.2856f); p.b = double(1.4268f); p.c = double(10.2551f); p.d = double(3.2128f);
    }
    if(!(std::cin >> p.fc >> p.E0 >> p.nu >> p.epsCrush)) return 2;
    jiangjj::State s;
    jiangjj::Vector6 e;
    while(std::cin >> e[0] >> e[1] >> e[2] >> e[3] >> e[4] >> e[5]) {
        jiangjj::Matrix6 c;
        // Marc passes E and DE = enew - E; Get_DS forms E + DE, which is not
        // bitwise equal to enew. Feed the kernel the same total strain.
        jiangjj::Vector6 total;
        for(int i = 0; i < 6; ++i) total[i] = s.strain[i] + (e[i] - s.strain[i]);
        s = jiangjj::update(p, s, total, &c);
        for(double x : s.stress) std::printf("%25.16e", x);
        std::printf("%25.16e%25.16e%25.16e%25.16e%25.16e%25.16e\n",
                    double(s.crushed ? 100 : 0), s.betaMax, c[0][0], c[0][1], c[3][3], 0.);
    }
    return 0;
}
