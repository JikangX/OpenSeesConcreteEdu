# Native plane stress -> LayeredShell -> ShellMITC4 integration verification.
# Run: OpenSees.exe shell_tests.tcl <temporary-output-directory>
set out [lindex $argv 0]
if {$out eq ""} {set out "adinaConcrete2D-output"}
file mkdir $out
set checks 0
proc check {condition message} {
    global checks
    incr checks
    if {!$condition} {error $message}
}
proc near {a b tol message} {
    check [expr {abs($a-$b)<$tol*max(1.0,abs($a),abs($b))}] "$message: $a vs $b"
}
proc configure {increment} {
    constraints Transformation
    numberer RCM
    system FullGeneral
    test NormDispIncr 1e-12 30
    algorithm Newton
    integrator LoadControl $increment
    analysis Static
}
proc adinaShellPatch {deformation steps {type AdinaConcrete2D} {layers 3}} {
    wipe
    model BasicBuilder -ndm 3 -ndf 6
    nDMaterial $type 1 30000 .2 3 -30 -.002 -6 -.006 -G 9000 -rho 2.4
    # Assign the native three-component material DIRECTLY to the section.
    # Its PlateFiber copy is produced by the section's ordinary API request.
    section LayeredShell 1 $layers 1 .12
    set coords {{0 0 0} {1 0 0} {1 1 0} {0 1 0}}
    set tag 0
    foreach xyz $coords {incr tag;node $tag {*}$xyz}
    element ShellMITC4 1 1 2 3 4 1
    node 101 2 0 0;node 102 2 0 0
    fix 101 1 1 1 1 1 1;fix 102 0 1 1 1 1 1
    uniaxialMaterial Elastic 2 1000
    element zeroLength 2 101 102 -mat 2 -dir 1
    lassign $deformation ex ey gxy kx ky kxy gx gy
    timeSeries Linear 1
    pattern Plain 1 1 {
        set tag 0
        foreach xyz $coords {
            incr tag;lassign $xyz x y z
            sp $tag 1 [expr {$ex*$x+.5*$gxy*$y}]
            sp $tag 2 [expr {$ey*$y+.5*$gxy*$x}]
            sp $tag 3 [expr {.5*$kx*$x*$x+.5*$ky*$y*$y+.5*$kxy*$x*$y+$gx*$x+$gy*$y}]
            sp $tag 4 [expr {$ky*$y+.5*$kxy*$x}]
            sp $tag 5 [expr {-$kx*$x-.5*$kxy*$y}]
            sp $tag 6 0
        }
        load 102 1 0 0 0 0 0
    }
    configure [expr {1.0/$steps}]
    check [expr {[analyze $steps]==0}] "prescribed shell patch failed"
}
# A plane-stress quad exercises the original three-component formulation.
# Constrained transverse displacement leaves an exactly uniform stress field,
# while the two loaded edge nodes remain free for nonlinear Newton iterations.
wipe
model BasicBuilder -ndm 2 -ndf 2
nDMaterial AdinaConcrete2D 1 30000 .2 3 -30 -.002 -6 -.006
node 1 0 0;node 2 1 0;node 3 1 1;node 4 0 1
fix 1 1 1;fix 2 0 1;fix 3 0 1;fix 4 1 1
element quad 1 1 2 3 4 1.0 PlaneStress 1
timeSeries Linear 1
pattern Plain 1 1 {load 2 -5 0;load 3 -5 0}
configure .025
for {set step 1} {$step<=40} {incr step} {
    check [expr {[analyze 1]==0}] "plane-stress quad Newton step"
    for {set gp 1} {$gp<=4} {incr gp} {
        set stress [eleResponse 1 material $gp stress]
        check [expr {[llength $stress]==3}] "native quad material order"
        near [lindex $stress 0] [expr {-.25*$step}] 2e-8 "quad compression equilibrium"
    }
}
reactions
near [expr {[nodeReaction 1 1]+[nodeReaction 4 1]}] 10 1e-9 "plane-stress quad reaction balance"
# All eight section components: membrane, curvature, transverse shear.
foreach type {AdinaConcrete2D AdinaConcrete2DPlateFiber} {
  for {set component 0} {$component<8} {incr component} {
    set deformation {0 0 0 0 0 0 0 0}
    lset deformation $component 1e-7
    adinaShellPatch $deformation 1 $type
    for {set gp 1} {$gp<=4} {incr gp} {
      set sec [eleResponse 1 section $gp deformation]
      check [expr {[llength $sec]==8}] "section deformation order"
      for {set i 0} {$i<8} {incr i} {near [lindex $sec $i] [lindex $deformation $i] 1e-12 "shell patch deformation $i"}
      set resultant [eleResponse 1 section $gp force]
      set expected {0 0 0 0 0 0 0 0}
      for {set layer 1} {$layer<=3} {incr layer} {
        set z [expr {-.06+($layer-.5)*.04}]
        set eps {}
        for {set i 0} {$i<3} {incr i} {lappend eps [expr {[lindex $sec $i]-$z*[lindex $sec [expr {$i+3}]]}]}
        # PlateFiber order [xx,yy,xy,yz,xz], section order shear [xz,yz].
        lappend eps [lindex $sec 7] [lindex $sec 6]
        set sig [list [expr {31250*[lindex $eps 0]+6250*[lindex $eps 1]}] \
                      [expr {6250*[lindex $eps 0]+31250*[lindex $eps 1]}] \
                      [expr {12500*[lindex $eps 2]}] \
                      [expr {9000*[lindex $eps 3]}] [expr {9000*[lindex $eps 4]}]]
        set actualE [eleResponse 1 section $gp fiber $layer strain]
        set actualS [eleResponse 1 section $gp fiber $layer stress]
        check [expr {[llength $actualS]==5}] "plate-fiber stress order"
        for {set i 0} {$i<5} {incr i} {
          near [lindex $actualE $i] [lindex $eps $i] 1e-13 "layer strain $i"
          near [lindex $actualS $i] [lindex $sig $i] 1e-11 "layer stress $i"
        }
        for {set i 0} {$i<3} {incr i} {
          lset expected $i [expr {[lindex $expected $i]+.04*[lindex $sig $i]}]
          set j [expr {$i+3}]
          lset expected $j [expr {[lindex $expected $j]+.04*$z*[lindex $sig $i]}]
        }
        lset expected 6 [expr {[lindex $expected 6]+.04*[lindex $sig 4]}]
        lset expected 7 [expr {[lindex $expected 7]+.04*[lindex $sig 3]}]
      }
      for {set i 0} {$i<8} {incr i} {near [lindex $resultant $i] [lindex $expected $i] 1e-11 "integrated N/M/V $i"}
    }
  }
}
# Midpoint layer integration converges to the analytical h^3/12 bending
# stiffness. This check distinguishes section quadrature from material error.
set previousError 1e100
foreach n {3 6 12 24} {
    adinaShellPatch {0 0 0 1e-7 0 0 0 0} 1 AdinaConcrete2D $n
    set moment [lindex [eleResponse 1 section 1 force] 3]
    set exact [expr {-31250*pow(.12,3)/12*1e-7}]
    set predicted [expr {$exact*(1.0-1.0/($n*$n))}]
    near $moment $predicted 1e-12 "midpoint bending stiffness for $n layers"
    set err [expr {abs($moment-$exact)}]
    check [expr {$err<$previousError}] "thickness quadrature did not converge"
    set previousError $err
}
# Force-controlled nonlinear shell compression: real free structural DOFs,
# constitutive Jacobian, reaction balance, equal stresses at all shell points.
wipe
model BasicBuilder -ndm 3 -ndf 6
nDMaterial AdinaConcrete2D 1 30000 .2 3 -30 -.002 -6 -.006
section LayeredShell 1 6 1 .12
node 1 0 0 0;node 2 1 0 0;node 3 1 1 0;node 4 0 1 0
fix 1 1 1 1 1 1 1;fix 2 0 1 1 1 1 1;fix 3 0 1 1 1 1 1;fix 4 1 1 1 1 1 1
element ShellMITC4 1 1 2 3 4 1
timeSeries Linear 1
pattern Plain 1 1 {load 2 -.6 0 0 0 0 0;load 3 -.6 0 0 0 0 0}
configure .025
for {set step 1} {$step<=40} {incr step} {
    check [expr {[analyze 1]==0}] "force-controlled shell Newton step"
    for {set gp 1} {$gp<=4} {incr gp} {
      for {set layer 1} {$layer<=6} {incr layer} {
        near [lindex [eleResponse 1 section $gp fiber $layer stress] 0] [expr {-.25*$step}] 2e-8 "shell compression equilibrium"
      }
    }
}
reactions
near [expr {[nodeReaction 1 1]+[nodeReaction 4 1]}] 1.2 1e-9 "shell compression reaction balance"
check [expr {[nodeDisp 2 1]<0}] "shell compressed displacement sign"
# Bending puts the upper and lower layers on different constitutive paths.
# Check actual cracking and independent layer histories before restart.
adinaShellPatch {0 0 0 .006 0 0 .0001 -.0002} 180 AdinaConcrete2D 6
for {set gp 1} {$gp<=4} {incr gp} {
    set bottom [eleResponse 1 section $gp fiber 1 crackState]
    set top [eleResponse 1 section $gp fiber 6 crackState]
    check [expr {[llength $bottom]==10 && [llength $top]==10}] "crack history order"
    check [expr {[lindex $bottom 0]+[lindex $bottom 1]==1}] "bottom tensile layer cracked"
    check [expr {[lindex $top 0]+[lindex $top 1]==0}] "top compressive layer remains uncracked"
    near [lindex [eleResponse 1 section $gp fiber 1 stress] 0] 0 1e-12 "tensile layer releases normal stress"
    check [expr {[lindex [eleResponse 1 section $gp fiber 6 stress] 0]<0}] "upper layer in compression"
}
database File [file join $out restart]
save 8
integrator LoadControl .02
check [expr {[analyze 4]==0}] "post-checkpoint shell loading"
set expected {}
for {set gp 1} {$gp<=4} {incr gp} {
  for {set layer 1} {$layer<=6} {incr layer} {
    lappend expected [eleResponse 1 section $gp fiber $layer stress] [eleResponse 1 section $gp fiber $layer crackState]
  }
}
set expectedResultant [eleResponse 1 stresses]
wipe
model BasicBuilder -ndm 3 -ndf 6
database File [file join $out restart]
restore 8
wipeAnalysis
configure .02
check [expr {[analyze 4]==0}] "shell loading after empty-domain restart"
set k 0
for {set gp 1} {$gp<=4} {incr gp} {
  for {set layer 1} {$layer<=6} {incr layer} {
    foreach query {stress crackState} {
      set actual [eleResponse 1 section $gp fiber $layer $query]
      set reference [lindex $expected $k];incr k
      for {set i 0} {$i<[llength $reference]} {incr i} {near [lindex $actual $i] [lindex $reference $i] 1e-11 "restarted layer $query"}
    }
  }
}
set actualResultant [eleResponse 1 stresses]
for {set i 0} {$i<32} {incr i} {near [lindex $actualResultant $i] [lindex $expectedResultant $i] 1e-11 "restarted shell resultant"}
# Shared Tcl/Python parser: invalid and incomplete controls must fail.
# Reinforced layered shell under displacement-controlled tension -> closure
# -> reopening. PlateRebar remains a separate material layer; the concrete
# kernel receives no artificial steel stress or MCFT-specific state feedback.
wipe
model BasicBuilder -ndm 3 -ndf 6
nDMaterial AdinaConcrete2D 1 30000 .2 3 -30 -.002 -6 -.006
uniaxialMaterial Steel02 2 400 200000 .01 18 .925 .15
nDMaterial PlateRebar 3 2 0
section LayeredShell 1 5 1 .04 3 .001 1 .038 3 .001 1 .04
node 1 0 0 0;node 2 1 0 0;node 3 1 1 0;node 4 0 1 0
fix 1 1 1 1 1 1 1;fix 2 0 1 1 1 1 1;fix 3 0 1 1 1 1 1;fix 4 1 1 1 1 1 1
equalDOF 2 3 1
element ShellMITC4 1 1 2 3 4 1
timeSeries Linear 1
pattern Plain 1 1 {load 2 .5 0 0 0 0 0;load 3 .5 0 0 0 0 0}
configure .01
foreach history {{250 .000002} {400 -.000002} {400 .000002}} {
    lassign $history steps increment
    integrator DisplacementControl 2 1 $increment
    for {set step 1} {$step<=$steps} {incr step} {
        check [expr {[analyze 1]==0}] "reinforced shell cyclic displacement step"
        check [expr {abs([getLoadFactor 1])<100}] "nonfinite or unbounded reinforced shell reaction"
    }
    reactions
    near [expr {[nodeReaction 1 1]+[nodeReaction 4 1]+[getLoadFactor 1]}] 0 1e-8 "reinforced shell cyclic equilibrium"
    near [nodeDisp 2 1] [nodeDisp 3 1] 1e-12 "reinforced shell uniform edge displacement"
}
for {set gp 1} {$gp<=4} {incr gp} {
    for {set layer 1} {$layer<=5} {incr layer 2} {
        near [lindex [eleResponse 1 section $gp fiber $layer stress] 0] 0 1e-12 "reopened concrete tensile release"
    }
    set steel [lindex [eleResponse 1 section $gp fiber 2 stress] 0]
    near [lindex [eleResponse 1 section $gp force] 0] [expr {.002*$steel}] 1e-9 "steel carries reopened shell tension"
}
# Shared Tcl/Python parser: invalid and incomplete controls must fail.
foreach options {{-G 0} {-G -1} {-rho -2} {-retention 1.1 .5} {-maxStrainIncrement -1} {-unknown 1} {-envelope 0}} {
    check [catch {nDMaterial AdinaConcrete2D 99 30000 .2 3 -30 -.002 -6 -.006 {*}$options}] "invalid controls accepted: $options"
}
check [catch {nDMaterial AdinaConcrete2D 99 30000 .5 3 -30 -.002 -6 -.006}] "singular Poisson ratio accepted"
puts "PASS $checks OpenSees layered shell / restart / parser checks"
wipe
