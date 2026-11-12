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
import math
import csv

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
    elastic = [[-i*1e-7, -i*2e-7, i*1e-8] for i in range(1,21)]
    compression = [[0, -i*1e-5, 0] for i in range(1,551)]
    biaxial = [[-i*8e-6, -i*1e-5, 0] for i in range(1,401)]
    cyclic = [[0,-i*1e-5,0] for i in list(range(1,81))+list(range(79,19,-1))+list(range(21,101))]
    tension = [[i*2e-6,0,0] for i in range(1,101)]
    sequential = tension+[[.0002,i*2e-6,0] for i in range(1,101)]
    closure = tension+[[i*2e-6,0,0] for i in range(99,-101,-1)]
    shear = [[0,0,i*2e-6] for i in range(1,151)]
    nonproportional = compression[:120]+[[0,-.0012,i*2e-6] for i in range(1,121)]
    # A nonzero angle avoids PRNCPL's 180-degree sentinel collision and
    # validates tensile/closure mechanics independently of that legacy defect.
    c,s=math.cos(.3),math.sin(.3)
    oblique = [[c*c*e[0]+s*s*e[1]-c*s*e[2],
                s*s*e[0]+c*c*e[1]+c*s*e[2],
                2*c*s*(e[0]-e[1])+(c*c-s*s)*e[2]] for e in closure]
    return dict(elastic=elastic, compression=compression, biaxial=biaxial,
                cyclic=cyclic, tension=tension, sequential=sequential, closure=closure,
                oblique=oblique, shear=shear,nonproportional=nonproportional)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cxx', default='g++'); ap.add_argument('--fc', default='gfortran')
    ap.add_argument('--adina', type=Path)
    ap.add_argument('--opensees-build', type=Path,
                    help='CMake build tree containing libOpenSeesLIB.a (adapter tests)')
    ap.add_argument('--opensees', type=Path, help='Newly built Tcl interpreter (quad and layered-shell tests)')
    ap.add_argument('--build', type=Path, default=ROOT/'tmp'/'adinaConcrete2D-tests')
    args = ap.parse_args(); args.build.mkdir(parents=True, exist_ok=True)
    exe = args.build/'kernel_tests.exe'
    run([args.cxx, '-std=c++17', '-O2', '-Wall', '-Wextra', '-pedantic',
         '-I'+str(ROOT / 'SRC/material/nD/adinaConcrete2D'), HERE / 'kernel_tests.cpp',
         ROOT / 'SRC/material/nD/adinaConcrete2D/Concrete2D.cpp',
         ROOT / 'SRC/material/nD/adinaConcrete3D/Concrete3D.cpp', '-o', exe])
    print(run([exe]), end='')
    report = {}
    if args.adina:
        original = args.adina.read_text(encoding='latin1')
        routines = args.build/'original_routines.f'
        routines.write_text(''.join(extract(original, n) for n in ['CDMOD','PRNCPL','CRAKID','DCRACK']), encoding='latin1')
        oracle = args.build/'oracle.exe'
        run([args.fc, '-std=legacy', '-fno-automatic', '-ffixed-line-length-none', '-O0',
             routines, HERE / 'oracle_driver.f90', '-o', oracle])
        # Diagnostic reference only: localized mathematical corrections.
        # The original routines remain available and are always compared too.
        fixed_text=routines.read_text(encoding='latin1')
        assert fixed_text.count('      ANG=ANG/2.')==1, 'unexpected PRNCPL source layout'
        assert fixed_text.count('EPSL(3)=-D12*SG + EPS(3)*CG')==1, 'unexpected CRAKID source layout'
        fixed_text=fixed_text.replace('      ANG=ANG/2.',
            '      ANG=ANG/2.\n      IF (ANG.GE.180.D0) ANG=ANG-180.D0\n'
            '      IF (ANG.EQ.0.D0) ANG=1.D-10')
        fixed_text=fixed_text.replace('EPSL(3)=-D12*SG + EPS(3)*CG',
                                     'EPSL(3)=-2.D0*D12*SG + EPS(3)*CG')
        # PRNCPL is reused for the NEW strain in an OLD principal stress frame.
        # Stress shear vanishes in that frame, strain shear generally does not
        # during nonproportional loading. Restore the complete strain tensor.
        fixed_text=fixed_text.replace('      EPSL(3)=0.D0',
            '      EPSL(3)=-2.D0*SG*CG*(EPS(1)-EPS(2))\n'
            '     1         +(CG*CG-SG*SG)*EPS(3)',1)
        # Use exact conversion so the corrected reference is objective.
        fixed_text=fixed_text.replace('57.296*ATAN(DUM)',
                                     '(180.D0/(4.D0*ATAN(1.D0)))*ATAN(DUM)')
        corrected_routines=args.build/'corrected_diagnostic_routines.f'
        corrected_routines.write_text(fixed_text,encoding='latin1')
        corrected_oracle=args.build/'corrected_oracle.exe'
        run([args.fc,'-std=legacy','-fno-automatic','-ffixed-line-length-none','-O0',
             corrected_routines,HERE / 'oracle_driver.f90','-o',corrected_oracle])
        for name, strains in paths().items():
            data = '\n'.join(' '.join(f'{x:.17g}' for x in strain) for strain in strains)+'\n'
            expected = [[float(x) for x in line.split()] for line in run([oracle], input=data, cwd=args.build).splitlines()]
            actual = [[float(x) for x in line.split()[:3]] for line in run([exe, 'path'], input=data).splitlines()]
            if len(expected) != len(strains) or len(actual) != len(strains):
                raise AssertionError(f'{name}: incomplete original-source run')
            error = max(abs(a-b) for aa,bb in zip(actual,expected) for a,b in zip(aa,bb))
            scale = max(1., max(abs(x) for row in expected for x in row))
            report[name] = dict(steps=len(strains), max_absolute_error=error, relative_to_peak=error/scale)
            print(name, json.dumps(report[name]))
            fixed=[[float(x) for x in line.split()] for line in run([corrected_oracle],input=data,cwd=args.build).splitlines()]
            assert len(fixed)==len(strains), f'{name}: incomplete corrected reference'
            fixed_error=max(abs(a-b) for aa,bb in zip(actual,fixed) for a,b in zip(aa,bb))
            report[name]['corrected_reference_error']=fixed_error
            print('corrected reference error',fixed_error)
            # Preserve every data point for independent inspection and plotting.
            with (args.build/f'{name}.csv').open('w',newline='') as stream:
                writer=csv.writer(stream)
                writer.writerow(['step','eps_xx','eps_yy','gamma_xy',
                                 'native_xx','native_yy','native_xy',
                                 'original_xx','original_yy','original_xy',
                                 'diagnostic_xx','diagnostic_yy','diagnostic_xy'])
                for k,(strain,native,raw,diagnostic) in enumerate(zip(strains,actual,expected,fixed),1):
                    writer.writerow([k,*strain,*native,*raw,*diagnostic])
            assert fixed_error<=2e-6*scale, f'{name}: corrected-reference mismatch'
            if name in ('tension','sequential','closure','oblique','nonproportional'):
                # Discontinuities in the raw legacy routine are expected ONLY
                # for these documented angle/shear defects. The first crack
                # predictor must still match before the affected path begins.
                prefix=120 if name=='nonproportional' else 49
                prefix_error=max(abs(a-b) for aa,bb in zip(actual[:prefix],expected[:prefix])
                                 for a,b in zip(aa,bb))
                assert prefix_error<=2e-6*scale, f'{name}: initial predictor mismatch'
                report[name].update(intentional_legacy_correction=True,
                                    matching_prefix_steps=prefix,prefix_error=prefix_error)
            else:
                assert error<=2e-6*scale, f'{name}: original-source mismatch'
            if error > 2e-6*scale:
                for i,(aa,bb) in enumerate(zip(actual,expected)):
                    if max(abs(a-b) for a,b in zip(aa,bb)) > 2e-6*scale:
                        print('first mismatch',i+1,'strain',strains[i],'C++',aa,'ADINA',bb);break
        (args.build/'oracle_results.json').write_text(json.dumps(report, indent=2))

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
             HERE / 'adapter_tests.cpp',ROOT / 'SRC/material/nD/AdinaConcrete2D.cpp',
             ROOT / 'SRC/material/nD/AdinaConcrete2DFactory.cpp',
             ROOT / 'SRC/material/nD/adinaConcrete2D/Concrete2D.cpp',
             ROOT / 'SRC/material/nD/adinaConcrete3D/Concrete3D.cpp',
             '-Wl,--gc-sections','-Wl,--start-group',*dependencies,'-Wl,--end-group',
             '-o',adapter])
        print(run([adapter]),end='')
        failure = args.build/'analysis_failure_tests.exe'
        # The Tcl build exports the legacy element API in a direct executable
        # object. Force that same object ahead of the archives: otherwise the
        # unified-interpreter API is extracted first and duplicates Tcl symbols.
        # Do not suppress duplicate-symbol diagnostics or substitute API stubs.
        tcl_api=args.opensees_build/'CMakeFiles/OPS_InterpTcl.dir/SRC/api/elementAPI_TCL.cpp.obj'
        run([args.cxx,'-std=c++17','-O2','-ffunction-sections','-fdata-sections',*defines,*includes,
             HERE / 'analysis_failure_tests.cpp',tcl_api,'-Wl,--gc-sections','-Wl,--start-group',
             *dependencies,'-Wl,--end-group','-o',failure])
        print(run([failure],combine_output=True),end='')
    if args.opensees:
        output = run([args.opensees.resolve(),HERE / 'shell_tests.tcl',args.build.resolve()/'shell'],combine_output=True)
        if 'OpenSees layered shell / restart / parser checks' not in output:
            raise AssertionError('Layered-shell test did not finish: '+output)
        print(output,end='')
        for name,token,arguments in [
            ('failure_tests.tcl','PASS shell constitutive failure propagation',[]),
            # The coupon gates require an EXACTLY uniform state (free transverse
            # displacement 0 at nu=0, mesh invariance to 2e-8). Only the
            # replay-from-commit formulation guarantees that. The default
            # in-place (ADINA) update makes the state depend on the iteration
            # path; roundoff between Gauss points then selects different crack
            # events (7e-8 m transverse displacement, see REVIEW). Test the
            # committed formulation strictly; do not loosen the gates.
            ('layered_lifecycle_tests.tcl','native RC membrane free-DOF',
             ['native',args.build.resolve()/'lifecycle','-tangent algorithmic -update committed'])]:
            output=run([args.opensees.resolve(),HERE/name,*arguments],combine_output=True)
            if token not in output:
                raise AssertionError(name+' did not finish: '+output)
            print(output,end='')
    print('PASS all requested tests')

if __name__ == '__main__':
    main()
