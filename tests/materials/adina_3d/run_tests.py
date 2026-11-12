"""Reproducible kernel and original-source comparison; Python standard library.

Example: python run_tests.py --cxx g++ --fc gfortran --adina C:/Work/Adina/all_adina.for
The ADINA source is optional for kernel checks, mandatory for --adina comparisons.
Original source is read-only. All generated files go to --build (default tmp).
"""
import argparse
from pathlib import Path
import subprocess
import re
import json
import shlex

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[0]

def run(command, **kwargs):
    combine_output = kwargs.pop('combine_output', False)
    result = subprocess.run([str(x) for x in command], capture_output=True, text=True, **kwargs)
    if result.returncode:
        raise RuntimeError(f"{command[0]} failed ({result.returncode})\n{result.stdout}\n{result.stderr}")
    # The Tcl interpreter routes its console channel to OPS' error stream on
    # some builds. Preserve that channel when checking its PASS sentinel.
    return result.stdout + result.stderr if combine_output else result.stdout

def extract(source, name):
    # Fixed-format legacy source: retain exactly the routine's original text.
    match = re.search(r"^      SUBROUTINE " + name + r"\b.*?^      END\s*$",
                      source, re.M | re.S)
    if not match:
        raise RuntimeError(f"Cannot find complete ADINA routine {name}")
    return match.group(0) + "\n"

def paths():
    elastic = [[-i*1e-6, -i*.5e-6, -i*2e-6, i*1e-7, i*2e-7, -i*3e-7] for i in range(1, 21)]
    compression = [[0, 0, -i*1e-5, 0, 0, 0] for i in range(1, 551)]
    biaxial = [[0, -i*5e-6, -i*1e-5, 0, 0, 0] for i in range(1, 401)]
    triaxial = [[-i*1e-6, -i*3e-6, -i*8e-6, 0, 0, 0] for i in range(1, 301)]
    cyclic = [[0, 0, -i*1e-5, 0, 0, 0] for i in list(range(1, 81))+list(range(79, 19, -1))+list(range(21, 101))]
    tension = [[i*2e-6, 0, 0, 0, 0, 0] for i in range(1, 101)]
    tension += [[.0002, i*2e-6, 0, 0, 0, 0] for i in range(1, 101)]
    tension += [[.0002, .0002, i*2e-6, 0, 0, 0] for i in range(1, 101)]
    closure = [[i*2e-6, 0, 0, 0, 0, 0] for i in range(1, 101)]
    closure += [[i*2e-6, 0, 0, 0, 0, 0] for i in range(99, -101, -1)]
    return dict(elastic=elastic, compression=compression, biaxial=biaxial,
                triaxial=triaxial, cyclic=cyclic, tension=tension, closure=closure)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cxx', default='g++'); ap.add_argument('--fc', default='gfortran')
    ap.add_argument('--adina', type=Path)
    ap.add_argument('--opensees-build', type=Path,
                    help='CMake build tree containing libOpenSeesLIB.a (adapter tests)')
    ap.add_argument('--opensees', type=Path, help='Newly built Tcl interpreter (brick tests)')
    ap.add_argument('--build', type=Path, default=ROOT/'tmp'/'adinaConcrete3D-tests')
    args = ap.parse_args(); args.build.mkdir(parents=True, exist_ok=True)
    exe = args.build/'kernel_tests.exe'
    run([args.cxx, '-std=c++17', '-O2', '-Wall', '-Wextra', '-pedantic',
         '-I'+str(ROOT / 'SRC/material/nD/adinaConcrete3D'), HERE / 'kernel_tests.cpp',
         ROOT / 'SRC/material/nD/adinaConcrete3D/Concrete3D.cpp', '-o', exe])
    print(run([exe]), end='')
    report = {}
    if args.adina:
        original = args.adina.read_text(encoding='latin1')
        routines = args.build/'original_routines.f'
        routines.write_text(''.join(extract(original, n) for n in ['CMOD3D','PRNCP3','CRAK3D','DCRAK3']), encoding='latin1')
        oracle = args.build/'oracle.exe'
        run([args.fc, '-std=legacy', '-fno-automatic', '-ffixed-line-length-none', '-O0',
             routines, HERE / 'oracle_driver.f90', '-o', oracle])
        for name, strains in paths().items():
            data = '\n'.join(' '.join(f'{x:.17g}' for x in strain) for strain in strains)+'\n'
            expected = [[float(x) for x in line.split()] for line in run([oracle], input=data, cwd=args.build).splitlines()]
            actual = [[float(x) for x in line.split()[:6]] for line in run([exe, 'path'], input=data).splitlines()]
            if len(expected) != len(strains) or len(actual) != len(strains):
                raise AssertionError(f'{name}: incomplete original-source run')
            error = max(abs(a-b) for aa,bb in zip(actual,expected) for a,b in zip(aa,bb))
            scale = max(1., max(abs(x) for row in expected for x in row))
            report[name] = dict(steps=len(strains), max_absolute_error=error, relative_to_peak=error/scale)
            print(name, json.dumps(report[name]))
            if name == 'tension':
                # CRAK3D skips reordering transverse axes when SIGP(6)==0.
                # Under sequential axis-aligned tension it consequently zeros
                # the wrong normal component at the second crack. The port
                # diagonalizes the transverse block even for zero shear.
                # This is an intentional, independently tested correction.
                assert abs(actual[135][1]) < 1e-12 and expected[135][1] > 3.
                assert max(abs(a-b) for aa,bb in zip(actual[:135],expected[:135])
                           for a,b in zip(aa,bb)) <= 2e-6*scale
                assert all(max(row[:3]) <= 3.+1e-12 for row in actual)
                assert max(abs(x) for x in actual[-1]) < 1e-12
                report[name]['intentional_legacy_axis_correction'] = True
                report[name]['matching_prefix_steps'] = 135
                print('PASS documented zero-shear transverse-axis correction')
            elif error > 2e-6*scale:
                for i,(aa,bb) in enumerate(zip(actual,expected)):
                    if max(abs(a-b) for a,b in zip(aa,bb)) > 2e-6*scale:
                        print('first mismatch',i+1,'strain',strains[i],'C++',aa,'ADINA',bb);break
        (args.build/'oracle_results.json').write_text(json.dumps(report, indent=2))
        assert all(x.get('intentional_legacy_axis_correction') or x['relative_to_peak'] <= 2e-6
                   for x in report.values()), 'Original-source comparison failed'
    if args.opensees_build:
        flags = (args.opensees_build/'CMakeFiles/OPS_Material.dir/flags.make').read_text()
        defines = shlex.split(re.search(r'^CXX_DEFINES = (.*)$',flags,re.M)[1])
        includes = shlex.split(re.search(r'^CXX_INCLUDES = (.*)$',flags,re.M)[1])
        includes = [re.sub(r'^-I/([a-zA-Z])/',r'-I\1:/',x) for x in includes]
        adapter = args.build/'adapter_tests.exe'
        # Reuse this MinGW CMake target's actual dependency list. Linking only
        # libOpenSeesLIB omits the broker, Fortran and interpreter libraries.
        # Keep the executable objects as an archive: our test supplies main(),
        # so the interpreter's main object will not be extracted.
        build_make = (args.opensees_build/'CMakeFiles/OpenSees.dir/build.make').read_text()
        link_line = next(line for line in build_make.splitlines()
                         if ' -o OpenSees.exe ' in line and 'libOPS_InterpTcl.a' in line)
        dependency_tokens = shlex.split(link_line[link_line.index('libOPS_InterpTcl.a'):])
        dependencies = [args.opensees_build/'CMakeFiles/OpenSees.dir/objects.a']
        for token in dependency_tokens:
            token = re.sub(r'^/([a-zA-Z])/',r'\1:/',token)
            dependencies.append(args.opensees_build/token if token.endswith('.a') and ':' not in token else token)
        run([args.cxx,'-std=c++17','-O2','-ffunction-sections','-fdata-sections',*defines,*includes,
             HERE / 'adapter_tests.cpp',ROOT / 'SRC/material/nD/AdinaConcrete3D.cpp',
             ROOT / 'SRC/material/nD/AdinaConcrete3DFactory.cpp',
             ROOT / 'SRC/material/nD/adinaConcrete3D/Concrete3D.cpp',
             '-Wl,--gc-sections','-Wl,--start-group',*dependencies,'-Wl,--end-group',
             '-o',adapter])
        print(run([adapter]),end='')
    if args.opensees:
        output = run([args.opensees.resolve(),HERE / 'brick_tests.tcl',args.build.resolve()/'brick'],combine_output=True)
        if 'OpenSees brick / restart / parser checks' not in output:
            raise AssertionError('Brick test did not finish: '+output)
        print(output,end='')
    print('PASS all requested tests')

if __name__ == '__main__':
    main()
