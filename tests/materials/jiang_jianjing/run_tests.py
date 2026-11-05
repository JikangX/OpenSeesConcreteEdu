"""Verification of the JiangJianjingConcrete migration.

1. C++ kernel vs the original MSC.Marc Fortran (JiangJJ2007.for, compiled
   unchanged, routine Get_DS driven like HYPELA2) on many strain paths.
2. Lecture example (lecture 7, slides 33-39): beta, Es, stress.
3. OpenSees Tcl checks (single stdBrick, PlaneStress wrapper, restart,
   parser) with the new interpreter.

usage: python run_tests.py --build <new dir> [--opensees bm/OpenSees.exe]
       [--fortran SRC/material/nD/JiangJianjingConcrete/JiangJJ2007.for]
"""
import argparse
import json
import math
import os
import random
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parents[0]
ROOT = HERE.parents[2]
SRC = ROOT / 'SRC/material/nD/JiangJianjingConcrete'


def run(cmd, **kw):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)
    if r.returncode:
        raise RuntimeError(f'{cmd[0]} failed ({r.returncode})\n{r.stdout}\n{r.stderr}')
    return r.stdout + r.stderr


def paths():
    out = {}
    lin = lambda a, n: [[x * k / n for x in a] for k in range(1, n + 1)]
    out['uniaxial_compression'] = lin([0.00052, 0.00052, -0.0040, 0, 0, 0], 100)
    out['uniaxial_stress_like'] = lin([0.0006, 0.0006, -0.0030, 0, 0, 0], 60)
    out['biaxial_compression'] = lin([0.0008, -0.0025, -0.0025, 0, 0, 0], 80)
    out['triaxial_hydrostatic'] = lin([-0.002, -0.002, -0.002, 0, 0, 0], 40)
    out['confined_then_axial'] = lin([-0.0002, -0.0002, -0.0002, 0, 0, 0], 10) + \
        [[-0.0002 + 0.0001 * k / 60, -0.0002 + 0.0001 * k / 60, -0.0002 - 0.0035 * k / 60, 0, 0, 0] for k in range(1, 61)]
    out['uniaxial_tension'] = lin([-0.00004, -0.00004, 0.0003, 0, 0, 0], 60)
    out['pure_shear'] = lin([0, 0, 0, 0.0015, 0, 0], 60)
    out['general_shear'] = lin([-0.0004, 0.0001, -0.0012, 0.0006, -0.0004, 0.0003], 80)
    up = lin([0.0003, 0.0003, -0.0018, 0, 0, 0], 40)
    out['unload_reload'] = up + [[x * (40 - k) / 40 for x in up[-1]] for k in range(1, 31)] + \
        [[x * (10 + k) / 40 for x in up[-1]] for k in range(1, 41)]
    rng = random.Random(20261007)
    for r in range(6):
        e = [0.0] * 6
        rows = []
        for k in range(120):
            e = [x + rng.uniform(-1, 1) * 4e-5 - (2e-5 if i < 3 else 0) for i, x in enumerate(e)]
            rows.append(list(e))
        out[f'random{r}'] = rows
    return out


def compare(oracle, kernel, params, rows, mode='fortran-literals'):
    data = ' '.join(f'{x:.17g}' for x in params) + '\n' + '\n'.join(' '.join(f'{x:.17g}' for x in r) for r in rows) + '\n'
    a = [[float(x) for x in l.split()] for l in run([oracle], input=data, cwd=oracle.parent).splitlines() if l.strip()]
    b = [[float(x) for x in l.split()] for l in run([kernel, mode], input=data).splitlines() if l.strip()]
    assert len(a) == len(b) == len(rows), (len(a), len(b), len(rows))
    fc = params[0]
    err = [max(abs(x - y) for x, y in zip(ra[:6], rb[:6])) / fc for ra, rb in zip(a, b)]
    crush = sum(ra[6] != rb[6] for ra, rb in zip(a, b))
    beta = max(abs(ra[7] - rb[7]) for ra, rb in zip(a, b))
    tang = max(abs(x - y) for ra, rb in zip(a, b) for x, y in zip(ra[8:11], rb[8:11])) / params[1]
    first = next((k + 1 for k, e in enumerate(err) if e > 1e-10), None)
    # Same algorithm; 1-ulp differences between the C and Fortran math
    # libraries (sqrt, acos) are amplified only near failure:
    # d(theta)/d(cos 3theta) is infinite on the meridians (cos 3theta = +-1),
    # Es ~ sqrt(1 - beta) and the stress-hold switch at beta = 0.99999.
    # Pre-failure branch: steps while the Fortran betaMax < 0.98.
    pre = [k for k, ra in enumerate(a) if ra[7] < 0.98]
    return dict(steps=len(rows), max_stress_error_over_fc=max(err), crush_flag_mismatches=crush,
                pre_failure_steps=len(pre),
                pre_failure_max_stress_error_over_fc=max((err[k] for k in pre), default=0.),
                bitwise_identical_steps=sum(e == 0. for e in err),
                max_betaMax_error=beta, max_tangent_error_over_E0=tang, first_step_over_1e10fc=first,
                crushed_at_end=bool(a[-1][6]), betaMax_end=a[-1][7], oracle_betaMax=[ra[7] for ra in a])


def lecture_example(kernel_lib_exe):
    out = run([kernel_lib_exe, 'lecture'])
    return json.loads(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--build', type=Path, required=True)
    ap.add_argument('--opensees', type=Path, default=None)
    ap.add_argument('--fortran', type=Path, default=SRC / 'JiangJJ2007.for')
    ap.add_argument('--cxx', default='g++')
    ap.add_argument('--fc', default='gfortran')
    a = ap.parse_args()
    b = a.build.resolve()
    if b.exists():
        raise FileExistsError(f'choose a fresh build directory: {b}')
    b.mkdir(parents=True)
    report = {}
    # 1. Fortran oracle (the original file, unchanged) and C++ driver
    shutil.copyfile(a.fortran, b / 'JiangJJ2007.for')
    run([a.fc, '-std=legacy', '-ffixed-line-length-none', '-O0', '-c', 'JiangJJ2007.for'], cwd=b)
    run([a.fc, '-O0', HERE / 'oracle_driver.f90', 'JiangJJ2007.o', '-o', 'oracle.exe'], cwd=b)
    # Diagnostic copy with ONE change: the intrinsic matmul(C%Ds,C%Strain) is
    # replaced by the plain sequential sum. libgfortran's matmul accumulates in
    # a different order (1-ulp different stresses; another Fortran runtime,
    # e.g. Marc's, differs again). With it the C++ must agree bit for bit.
    src = (b / 'JiangJJ2007.for').read_bytes().decode('latin1')
    key = 'C%Stress=matmul(C%Ds,C%Strain);'
    assert src.count(key) == 1, 'unexpected source layout'
    src = src.replace(key, 'call SeqMatVec(C%Ds,C%Strain,C%Stress)')
    src += (
        '\n'
        '\tsubroutine SeqMatVec(D,E,S)\n'
        '\t\timplicit none\n'
        '\t\treal*8 D(6,6),E(6),S(6)\n'
        '\t\tinteger I,J\n'
        '\t\tdo I=1,6\n'
        '\t\t\tS(I)=0.d0\n'
        '\t\t\tdo J=1,6\n'
        '\t\t\t\tS(I)=S(I)+D(I,J)*E(J)\n'
        '\t\t\tend do\n'
        '\t\tend do\n'
        '\tend subroutine\n'
        '\n'
    )
    (b / 'seq').mkdir()
    (b / 'seq' / 'JiangJJ2007_seq.for').write_bytes(src.encode('latin1'))
    run([a.fc, '-std=legacy', '-ffixed-line-length-none', '-O0', '-c', 'JiangJJ2007_seq.for'], cwd=b / 'seq')
    run([a.fc, '-O0', HERE / 'oracle_driver.f90', 'JiangJJ2007_seq.o', '-o', 'oracle_seq.exe'], cwd=b / 'seq')
    # Corrected reference: additionally J2 = max(J2, 0). For (near-)hydrostatic
    # stresses the deviatoric J2 can be -1e-30 by roundoff; Get_DS then takes
    # sqrt of a negative number (NaN), the NaN test fails, the routine writes
    # "Wrong", sets beta = 1 and the stress is frozen from then on (original
    # defect D1). The C++ clamps J2 at 0 (the exact value).
    key2 = 'J2=-s(1)*s(2)-s(2)*s(3)-s(3)*s(1)+s(4)**2+s(5)**2+s(6)**2'
    assert src.count(key2) == 1, 'unexpected source layout (J2)'
    src2 = src.replace(key2, key2 + '; if(J2<0.d0) J2=0.d0')
    # And: for I1 >= fc/d (TempC >= 0, tension beyond the apex of the failure
    # surface) the quadratic has only a non-positive root; Get_DS then gets
    # beta <= 0 and Es >= E0 and the point never fails (original defect D2).
    # The C++ treats it like the no-root branch (beta = 1).
    key3 = 'if(TempB**2-4.*TempA*TempC>=0.) then'
    assert src2.count(key3) == 1, 'unexpected source layout (discriminant)'
    src2 = src2.replace(key3, 'if(TempB**2-4.*TempA*TempC>=0. .and. TempC<0.) then')
    (b / 'corr').mkdir()
    (b / 'corr' / 'JiangJJ2007_corr.for').write_bytes(src2.encode('latin1'))
    run([a.fc, '-std=legacy', '-ffixed-line-length-none', '-O0', '-c', 'JiangJJ2007_corr.for'], cwd=b / 'corr')
    run([a.fc, '-O0', HERE / 'oracle_driver.f90', 'JiangJJ2007_corr.o', '-o', 'oracle_corr.exe'], cwd=b / 'corr')
    run([a.cxx, '-std=c++17', '-O2', '-Wall', '-Wextra', '-I' + str(SRC), HERE / 'kernel_vs_fortran.cpp',
         SRC / 'JiangJianjing.cpp', '-o', b / 'kernel.exe'])
    run([a.cxx, '-std=c++17', '-O2', '-Wall', '-Wextra', '-I' + str(SRC), HERE / 'kernel_checks.cpp',
         SRC / 'JiangJianjing.cpp', '-o', b / 'kernel_checks.exe'])
    sets = {'fc30_E30000': [30., 30000., 0.18, -0.0033], 'fc20_E20000': [20., 20000., 0.18, -0.0033],
            'fc50_E34500_nu0.2': [50., 34500., 0.2, -0.0033]}
    report['fortran_comparison'] = {}
    worst = 0.
    for sname, prm in sets.items():
        for pname, rows in paths().items():
            # Primary: corrected Fortran (sequential mat-vec sum, J2 >= 0),
            # same REAL*4 literals as the Fortran evaluates.
            r = compare(b / 'corr' / 'oracle_corr.exe', b / 'kernel.exe', prm, rows)
            # Information: Fortran with only the summation order changed, the
            # unmodified Fortran, and the production default constants.
            r['sequential_sum_fortran'] = compare(b / 'seq' / 'oracle_seq.exe', b / 'kernel.exe', prm, rows)
            r['unmodified_fortran'] = compare(b / 'oracle.exe', b / 'kernel.exe', prm, rows)
            r['exact_constants'] = compare(b / 'oracle.exe', b / 'kernel.exe', prm, rows, 'exact')
            # D1: the unpatched routine jumps to beta = 1 ("Wrong") while the
            # corrected one is far from failure.
            r['defect_D1_exercised'] = any(x == 1.0 and y < 0.98 for x, y in zip(
                r['sequential_sum_fortran']['oracle_betaMax'], r['oracle_betaMax']))
            for v in (r, r['sequential_sum_fortran'], r['unmodified_fortran'], r['exact_constants']):
                v.pop('oracle_betaMax')
            report['fortran_comparison'][f'{sname}:{pname}'] = r
            worst = max(worst, r['max_stress_error_over_fc'])
            print(f'{sname}:{pname:22s} steps={r["steps"]:4d} bitwise={r["bitwise_identical_steps"]:4d} '
                  f'pre-failure({r["pre_failure_steps"]:3d})={r["pre_failure_max_stress_error_over_fc"]:.1e} '
                  f'all={r["max_stress_error_over_fc"]:.1e} '
                  f'unmodified={r["unmodified_fortran"]["max_stress_error_over_fc"]:.1e} '
                  f'crush-mismatch={r["crush_flag_mismatches"]} crushed={r["crushed_at_end"]} '
                  f'betaMax={r["betaMax_end"]:.4f}' + (' [original defect D1]' if r['defect_D1_exercised'] else ''))
            # On the meridians (axisymmetric paths) a 1-ulp change of cos(3 theta)
            # moves theta by ~sqrt(eps): the problem itself is conditioned at the
            # sqrt(eps) ~ 1.5e-8 level there. Require 1e-6 fc before failure.
            assert r['pre_failure_max_stress_error_over_fc'] < 1e-6, 'Fortran mismatch before failure'
            assert r['crush_flag_mismatches'] == 0, 'crushing flag mismatch'
    report['worst_stress_error_over_fc'] = worst
    report['worst_unmodified_fortran'] = max(v['unmodified_fortran']['max_stress_error_over_fc']
                                             for v in report['fortran_comparison'].values())
    report['worst_with_exact_constants'] = max(v['exact_constants']['max_stress_error_over_fc']
                                               for v in report['fortran_comparison'].values())
    print('worst stress error / fc (all steps, incl. post-failure): vs corrected Fortran %.2e; vs unmodified Fortran %.2e; '
          'exact decimal constants %.2e' % (worst, report['worst_unmodified_fortran'],
                                            report['worst_with_exact_constants']))
    # 2. kernel checks incl. the lecture example
    out = run([b / 'kernel_checks.exe'])
    print(out, end='')
    report['kernel_checks'] = out.strip().splitlines()[-1]
    # 3. Tcl checks
    if a.opensees:
        env = dict(os.environ)
        if os.name == 'nt':
            env['PATH'] = 'C:/msys64/mingw64/bin' + os.pathsep + env.get('PATH', '')
        out = run([a.opensees.resolve(), HERE / 'tcl_checks.tcl', b / 'tcl'], env=env)
        if 'PASS JiangJianjingConcrete Tcl checks' not in out:
            raise AssertionError('Tcl checks did not finish:\n' + out)
        print([l for l in out.splitlines() if l.startswith('PASS')][-1])
        report['tcl'] = [l for l in out.splitlines() if l.startswith(('PASS', 'INFO'))]
    (b / 'report.json').write_text(json.dumps(report, indent=1))
    print('PASS all JiangJianjingConcrete tests')


if __name__ == '__main__':
    main()
