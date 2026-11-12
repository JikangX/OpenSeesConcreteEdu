# Run with the newly built OpenSees interpreter, from any working directory:
# OpenSees.exe brick_tests.tcl <output-directory>
# A fully prescribed affine cube isolates the material's engineering shear
# convention from equilibrium errors. A separate elastic spring supplies one
# free equation to exercise the normal OpenSees analysis lifecycle.
set out [lindex $argv 0]
if {$out eq ""} {set out "adinaConcrete3D-output"}
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
proc configureAnalysis {increment} {
    constraints Transformation
    numberer RCM
    system FullGeneral
    test NormDispIncr 1e-12 10
    algorithm Newton
    integrator LoadControl $increment
    analysis Static
}
proc cube {strain steps} {
    wipe
    model BasicBuilder -ndm 3 -ndf 3
    nDMaterial AdinaConcrete3D 1 30000 0.2 3 -30 -0.002 -6 -0.006
    set coords {{0 0 0} {1 0 0} {1 1 0} {0 1 0} {0 0 1} {1 0 1} {1 1 1} {0 1 1}}
    set tag 0
    foreach xyz $coords {incr tag;node $tag {*}$xyz}
    element stdBrick 1 1 2 3 4 5 6 7 8 1
    node 101 2 0 0;node 102 2 0 0;node 103 2 0 0
    fix 101 1 1 1;fix 102 0 1 1;fix 103 0 1 1
    equalDOF 102 103 1
    uniaxialMaterial Elastic 2 1000
    element zeroLength 2 101 102 -mat 2 -dir 1
    lassign $strain ex ey ez gxy gyz gxz
    timeSeries Linear 1
    pattern Plain 1 1 {
        set tag 0
        foreach xyz $coords {
            incr tag;lassign $xyz x y z
            sp $tag 1 [expr {$ex*$x+0.5*$gxy*$y+0.5*$gxz*$z}]
            sp $tag 2 [expr {$ey*$y+0.5*$gxy*$x+0.5*$gyz*$z}]
            sp $tag 3 [expr {$ez*$z+0.5*$gxz*$x+0.5*$gyz*$y}]
        }
        load 102 1 0 0
    }
    configureAnalysis [expr {1.0/$steps}]
    check [expr {[analyze $steps]==0}] "affine brick analysis failed"
}
# All six individual components, then a mixed tensor: independent Hooke law.
foreach strain {{1e-7 0 0 0 0 0} {0 1e-7 0 0 0 0} {0 0 1e-7 0 0 0}
                {0 0 0 1e-7 0 0} {0 0 0 0 1e-7 0} {0 0 0 0 0 1e-7}
                {-2e-7 -1e-7 -3e-7 2e-7 -3e-7 4e-7}} {
    cube $strain 1
    set trace [expr {[lindex $strain 0]+[lindex $strain 1]+[lindex $strain 2]}]
    for {set gp 1} {$gp<=8} {incr gp} {
        set stress [eleResponse 1 material $gp stress]
        set actualStrain [eleResponse 1 material $gp strain]
        check [expr {[llength $stress]==6}] "material stress response size"
        for {set i 0} {$i<6} {incr i} {
            set expected [expr {$i<3 ? 25000*[lindex $strain $i]+8333.333333333333*$trace : 12500*[lindex $strain $i]}]
            near [lindex $stress $i] $expected 1e-11 "elastic brick stress $gp,$i"
            near [lindex $actualStrain $i] [lindex $strain $i] 1e-13 "engineering strain $gp,$i"
        }
    }
}
# Load-controlled compression exercises the constitutive tangent in actual
# Newton equilibrium iterations (the affine tests prescribe every cube DOF).
wipe
model BasicBuilder -ndm 3 -ndf 3
nDMaterial AdinaConcrete3D 1 30000 .2 3 -30 -.002 -6 -.006
set coords {{0 0 0} {1 0 0} {1 1 0} {0 1 0} {0 0 1} {1 0 1} {1 1 1} {0 1 1}}
set tag 0
foreach xyz $coords {
    incr tag;node $tag {*}$xyz
    fix $tag 1 1 [expr {$tag<=4 ? 1 : 0}]
}
element stdBrick 1 1 2 3 4 5 6 7 8 1
timeSeries Linear 1
pattern Plain 1 1 {foreach node {5 6 7 8} {load $node 0 0 -2.5}}
configureAnalysis .05
for {set step 1} {$step<=20} {incr step} {
    check [expr {[analyze 1]==0}] "nonlinear force-controlled brick Newton step"
    for {set gp 1} {$gp<=8} {incr gp} {
        near [lindex [eleResponse 1 material $gp stress] 2] [expr {-.5*$step}] 2e-8 "force-controlled compression stress"
    }
}
reactions
set balance 0
foreach node {1 2 3 4} {set balance [expr {$balance+[nodeReaction $node 3]}]}
near $balance 10 1e-9 "compression reaction balance"
check [expr {[nodeDisp 5 3]<0}] "negative compression displacement"
cube {0.0002 0 0 0 0 0} 100
for {set gp 1} {$gp<=8} {incr gp} {
    set state [eleResponse 1 material $gp crackState]
    check [expr {[llength $state]==9 && [lindex $state 0]==1 && [lindex $state 3]==1}] "brick fixed crack state"
    near [lindex [eleResponse 1 material $gp stress] 0] 0 1e-12 "brick tensile normal release"
}
# Database reconstruction must use the registered object broker and sendSelf /
# recvSelf, not a retained pointer. Compare subsequent nonlinear history too.
database File [file join $out restart]
save 7
integrator LoadControl 0.05
check [expr {[analyze 2]==0}] "post-save analysis"
set expected [eleResponse 1 material 1 stress]
set expectedState [eleResponse 1 material 1 crackState]
# Start from an EMPTY domain, requiring actual broker reconstruction. Reopen
# the store after wipe; the helper's equalDOF also tests distinct MP/EQ fields.
wipe
model BasicBuilder -ndm 3 -ndf 3
database File [file join $out restart]
restore 7
wipeAnalysis
configureAnalysis 0.05
check [expr {[analyze 2]==0}] "post-restore analysis"
set actual [eleResponse 1 material 1 stress]
set actualState [eleResponse 1 material 1 crackState]
for {set i 0} {$i<6} {incr i} {near [lindex $actual $i] [lindex $expected $i] 1e-12 "restart stress $i"}
for {set i 0} {$i<9} {incr i} {near [lindex $actualState $i] [lindex $expectedState $i] 1e-12 "restart history $i"}
near [nodeDisp 103 1] [nodeDisp 102 1] 1e-12 "restart MP constraint"
# Parser rejects physically invalid inputs and incomplete options.
check [catch {nDMaterial AdinaConcrete3D 99 30000 .5 3 -30 -.002 -6 -.006}] "invalid nu accepted"
check [catch {nDMaterial AdinaConcrete3D 99 30000 .2 3 -30 -.002 -6 -.006 -envelope 0}] "short envelope accepted"
puts "PASS $checks OpenSees brick / restart / parser checks"
wipe
