# JiangJianjingConcrete OpenSees checks.
# usage: OpenSees tcl_checks.tcl <output-directory>
set out [lindex $argv 0]
if {$out eq ""} {set out jiangjj-tcl}
file mkdir $out
set checks 0
proc check {cond msg} {global checks; incr checks; if {!$cond} {error "FAIL: $msg"}}
proc near {a b tol msg} {check [expr {abs($a-$b) <= $tol*max(1.0,abs($b))}] "$msg: $a vs $b"}

# --- 1. parser ---------------------------------------------------------------
wipe
model basic -ndm 3 -ndf 3
check [expr {![catch {nDMaterial JiangJianjingConcrete 1 30.0 30000.0 0.18}]}] "default arguments"
check [expr {![catch {nDMaterial JiangJianjingConcrete 2 30.0 30000.0 0.18 -0.0035 -criterion 1.2856 1.4268 10.2551 3.2128 -rho 2.4e-9}]}] "all options"
check [catch {nDMaterial JiangJianjingConcrete 3 -30.0 30000.0 0.18}] "negative fc rejected"
check [catch {nDMaterial JiangJianjingConcrete 4 30.0 30000.0 0.18 0.0033}] "positive epsCrush rejected"
check [catch {nDMaterial JiangJianjingConcrete 5 30.0 30000.0 0.18 -unknown 1}] "unknown option rejected"
check [expr {![catch {nDMaterial JiangJianjingConcrete 6 30.0 30000.0 0.18 -beta current -tangent secant}]}] "-beta current"
check [expr {![catch {nDMaterial JiangJianjingConcrete 7 30.0 30000.0 0.18 -beta committed}]}] "-beta committed"
check [catch {nDMaterial JiangJianjingConcrete 8 30.0 30000.0 0.18 -beta later}] "bad -beta rejected"

# --- 2. single brick, uniaxial compression by prescribed top displacement ----
# lateral faces free; symmetric supports. Peak stress must be fc (beta = 1)
# and the stress is held at the ultimate value afterwards (total model).
proc brick {mat confinement} {
    wipe
    model basic -ndm 3 -ndf 3
    eval nDMaterial $mat
    set k 0
    foreach z {0.0 1.0} {foreach {x y} {0 0 1 0 1 1 0 1} {node [incr k] $x $y $z}}
    fix 1 1 1 1; fix 2 0 1 1; fix 3 0 0 1; fix 4 1 0 1
    element stdBrick 1 1 2 3 4 5 6 7 8 1
    equalDOF 5 6 3; equalDOF 5 7 3; equalDOF 5 8 3
    if {$confinement != 0.0} {
        # lateral pressure (stress controlled) on the x=1 and y=1 faces
        timeSeries Linear 9
        set f [expr {-$confinement/4.0}]
        pattern Plain 9 9 {
            foreach n {2 3 6 7} {load $n $f 0 0}
            foreach n {3 4 7 8} {load $n 0 $f 0}
        }
        constraints Plain; numberer Plain; system FullGeneral
        test NormUnbalance 1e-8 50 0; algorithm Newton
        integrator LoadControl 0.1; analysis Static
        for {set i 0} {$i < 10} {incr i} {if {[analyze 1] != 0} {error "confinement stage failed"}}
        loadConst -time 0.0
        wipeAnalysis
    }
    timeSeries Linear 1
    pattern Plain 1 1 {load 5 0 0 -1.0}
    constraints Plain; numberer Plain; system FullGeneral
    test NormUnbalance 1e-6 400 0; algorithm ModifiedNewton -initial
    integrator DisplacementControl 5 3 -2.5e-5
    analysis Static
}
proc axialPeak {nsteps fname} {
    set f [open $fname w]
    set peak 0.0
    for {set i 1} {$i <= $nsteps} {incr i} {
        if {[analyze 1] != 0} {break}
        set s [expr {-[getLoadFactor 1]}]
        puts $f "[expr {-[nodeDisp 5 3]}] $s [eleResponse 1 material 1 state]"
        if {-$s > $peak} {set peak [expr {-$s}]}
    }
    close $f
    return $peak
}
brick {JiangJianjingConcrete 1 30.0 30000.0 0.18} 0.0
set p [axialPeak 200 $out/uniaxial.txt]
# The secant modulus of a step uses beta of the PREVIOUS converged stress
# (explicit in beta, as Get_DS), so the computed peak depends on the step
# size and approaches fc from below; Newton stalls once beta -> 1
# (tangent sqrt(1-beta) E0 -> 0). Because of the lag the stress can also
# exceed fc by up to one increment before it is held.
check [expr {$p > 0.95*30.0 && $p < 30.0*1.02}] "uniaxial strength close to fc: $p"
puts "INFO uniaxial peak $p MPa"

# --- 3./4. stress-controlled cube: biaxial and confined compression ---------
# Unit cube, symmetry supports, uniform face loads (each node 1/4 of the face
# force). Stage 1: lateral confinement p (if any); stage 2: proportional
# increase of the "active" stresses in steps of 0.05 MPa until Newton fails.
# The last converged stress is the computed strength (the course reports, for
# fc 30 with MSC.Marc: uniaxial 30, biaxial 36/36, confinement 2 MPa -> 44,
# 6 MPa -> 66.6).
proc cube {opts} {
    wipe
    model basic -ndm 3 -ndf 3
    eval nDMaterial JiangJianjingConcrete 1 30.0 30000.0 0.18 $opts
    set k 0
    foreach z {0.0 1.0} {foreach {x y} {0 0 1 0 1 1 0 1} {node [incr k] $x $y $z}}
    fix 1 1 1 1; fix 2 0 1 1; fix 3 0 0 1; fix 4 1 0 1; fix 5 1 1 0; fix 6 0 1 0; fix 8 1 0 0
    element stdBrick 1 1 2 3 4 5 6 7 8 1
}
proc faceLoads {sx sy sz} {
    foreach n {2 3 6 7} {load $n [expr {$sx/4.0}] 0 0}
    foreach n {3 4 7 8} {load $n 0 [expr {$sy/4.0}] 0}
    foreach n {5 6 7 8} {load $n 0 0 [expr {$sz/4.0}]}
}
proc strength {lateral active {opts ""} {alg {ModifiedNewton -initial}}} {
    cube $opts
    constraints Plain; numberer Plain; system FullGeneral
    test NormUnbalance 1e-9 400 0; algorithm {*}$alg
    if {$lateral != 0.0} {
        timeSeries Linear 1
        pattern Plain 1 1 "faceLoads [expr {-$lateral}] [expr {-$lateral}] 0.0"
        integrator LoadControl 0.05; analysis Static
        for {set i 0} {$i < 20} {incr i} {if {[analyze 1] != 0} {error "confinement stage failed"}}
        loadConst -time 0.0
        wipeAnalysis
        constraints Plain; numberer Plain; system FullGeneral
        test NormUnbalance 1e-9 400 0; algorithm {*}$alg
    }
    timeSeries Linear 2
    pattern Plain 2 2 "faceLoads {*}$active"
    integrator LoadControl 0.05; analysis Static
    set last 0.0
    for {set i 0} {$i < 3000} {incr i} {
        if {[analyze 1] != 0} break
        set last [getLoadFactor 2]
    }
    return $last
}
# Original tangent (as Get_DS) with full Newton diverges once Es/Et > 2
# (beta > 8/9): record how far it gets, then use the recommended settings.
set uN [strength 0.0 {0.0 0.0 -1.0} "" Newton]
puts "INFO uniaxial, original tangent + full Newton: last converged $uN MPa"
set uS [strength 0.0 {0.0 0.0 -1.0} "-tangent secant" Newton]
puts "INFO uniaxial, -tangent secant + full Newton: last converged $uS MPa"
near $uS 30.0 0.025 "secant tangent + Newton reaches fc"
set u [strength 0.0 {0.0 0.0 -1.0}]
set bi [strength 0.0 {0.0 -1.0 -1.0}]
set c2 [strength 2.0 {0.0 0.0 -1.0}]
set c6 [strength 6.0 {0.0 0.0 -1.0}]
puts "INFO stress-controlled cube (original tangent, ModifiedNewton -initial): uniaxial $u, biaxial $bi, confined 2 MPa: [expr {$c2}], confined 6 MPa: [expr {$c6}] (MPa, last converged)"
puts "INFO course (Marc) results: uniaxial 30, biaxial 36, confined 2 MPa 44, confined 6 MPa 66.6"
near $u 30.0 0.025 "stress-controlled uniaxial strength"
near $bi 36.0 0.025 "biaxial strength"
near $c2 44.0 0.025 "confined 2 MPa strength"
near $c6 66.6 0.025 "confined 6 MPa strength"

# --- 5. PlaneStress wrapper with a quad, and PlateFiber with a layered shell -
# The wrappers condense sigma_33 = 0 (and the transverse shears) by an inner
# Newton on the material tangent; with the original tangent Et it diverges
# like the global Newton once Es/Et > 2. Use -tangent secant with wrappers.
proc planeStressPeak {opts} {
    wipe
    model basic -ndm 2 -ndf 2
    eval nDMaterial JiangJianjingConcrete 1 30.0 30000.0 0.18 $opts
    nDMaterial PlaneStress 2 1
    node 1 0 0; node 2 1 0; node 3 1 1; node 4 0 1
    fix 1 1 1; fix 2 0 1; fix 4 1 0
    element quad 1 1 2 3 4 1.0 PlaneStress 2
    equalDOF 3 4 2
    timeSeries Linear 1
    pattern Plain 1 1 {load 3 0 -1.0}
    constraints Transformation; numberer Plain; system FullGeneral
    test NormUnbalance 1e-6 400 0; algorithm ModifiedNewton -initial
    integrator DisplacementControl 3 2 -2.5e-5; analysis Static
    set peak 0.0
    for {set i 1} {$i <= 200} {incr i} {
        if {[analyze 1] != 0} {break}
        if {[getLoadFactor 1] > $peak} {set peak [getLoadFactor 1]}
    }
    return $peak
}
set p0 [planeStressPeak ""]
set p1 [planeStressPeak "-tangent secant"]
puts "INFO plane stress quad (PlaneStress wrapper): peak $p0 MPa with the original tangent, $p1 MPa with -tangent secant"
check [expr {$p1 > 0.95*30.0 && $p1 < 30.0*1.02}] "plane-stress wrapper peak near fc: $p1"

wipe
model basic -ndm 3 -ndf 6
nDMaterial JiangJianjingConcrete 1 30.0 30000.0 0.18 -tangent secant
nDMaterial PlateFiber 2 1
section LayeredShell 1 4 2 50 2 50 2 50 2 50
node 1 0 0 0; node 2 1000 0 0; node 3 1000 1000 0; node 4 0 1000 0
fix 1 1 1 1 1 1 1; fix 2 0 1 1 1 1 1; fix 3 0 0 1 1 1 1; fix 4 1 0 1 1 1 1
element ShellNLDKGQ 1 1 2 3 4 1
equalDOF 3 4 2
timeSeries Linear 1
pattern Plain 1 1 {load 3 0 -1.0 0 0 0 0}
constraints Transformation; numberer Plain; system FullGeneral
test NormUnbalance 1e-3 400 0; algorithm ModifiedNewton -initial
integrator DisplacementControl 3 2 -0.025; analysis Static
set ok 0
for {set i 1} {$i <= 40} {incr i} {set ok [analyze 1]; if {$ok != 0} break}
puts "INFO layered shell (PlateFiber wrapper): load factor [getLoadFactor 1] N at [nodeDisp 3 2] mm, status $ok"
check [expr {$ok == 0}] "layered shell with PlateFiber wrapper runs"

# --- 6. restart (database), both beta options ----------------------------------
foreach {opts db} {{} jj_db {-beta current -tangent secant} jj_db_current} {
wipe
model basic -ndm 3 -ndf 3
eval nDMaterial JiangJianjingConcrete 1 30.0 30000.0 0.18 $opts
set k 0
foreach z {0.0 1.0} {foreach {x y} {0 0 1 0 1 1 0 1} {node [incr k] $x $y $z}}
fix 1 1 1 1; fix 2 0 1 1; fix 3 0 0 1; fix 4 1 0 1
element stdBrick 1 1 2 3 4 5 6 7 8 1
equalDOF 5 6 3; equalDOF 5 7 3; equalDOF 5 8 3
timeSeries Linear 1
pattern Plain 1 1 {load 5 0 0 -1.0}
constraints Plain; numberer Plain; system FullGeneral
test NormUnbalance 1e-6 100 0; algorithm Newton
integrator DisplacementControl 5 3 -5e-5; analysis Static
for {set i 0} {$i < 20} {incr i} {analyze 1}
database File $out/$db
save 20
set s20 [eleResponse 1 material 1 stress]
set st20 [eleResponse 1 material 1 state]
for {set i 0} {$i < 10} {incr i} {analyze 1}
set s30 [eleResponse 1 material 1 stress]
restore 20
set r20 [eleResponse 1 material 1 stress]
set rt20 [eleResponse 1 material 1 state]
foreach a $s20 b $r20 {near $a $b 1e-12 "restart stress"}
foreach a $st20 b $rt20 {near $a $b 1e-12 "restart state"}
for {set i 0} {$i < 10} {incr i} {analyze 1}
foreach a $s30 b [eleResponse 1 material 1 stress] {near $a $b 1e-12 "restart continuation"}
}

# --- 7. -beta current: step-size independent before the ultimate stress -------
proc brickStress {opts du n} {
    wipe
    model basic -ndm 3 -ndf 3
    eval nDMaterial JiangJianjingConcrete 1 30.0 30000.0 0.18 $opts
    set k 0
    foreach z {0.0 1.0} {foreach {x y} {0 0 1 0 1 1 0 1} {node [incr k] $x $y $z}}
    fix 1 1 1 1; fix 2 0 1 1; fix 3 0 0 1; fix 4 1 0 1
    element stdBrick 1 1 2 3 4 5 6 7 8 1
    equalDOF 5 6 3; equalDOF 5 7 3; equalDOF 5 8 3
    timeSeries Linear 1
    pattern Plain 1 1 {load 5 0 0 -1.0}
    constraints Plain; numberer Plain; system FullGeneral
    test NormUnbalance 1e-8 100 0; algorithm Newton
    integrator DisplacementControl 5 3 $du; analysis Static
    for {set i 0} {$i < $n} {incr i} {if {[analyze 1] != 0} {error "brick analysis failed"}}
    return [lindex [eleResponse 1 material 1 stress] 2]
}
set a1 [brickStress {-beta current -tangent secant} -1.2e-4 10]
set a2 [brickStress {-beta current -tangent secant} -1.2e-5 100]
set c1 [brickStress {} -1.2e-4 10]
set c2 [brickStress {} -1.2e-5 100]
puts "INFO uniaxial brick at strain -0.0012: -beta current $a1 / $a2 MPa (10 / 100 steps), committed (original) $c1 / $c2 MPa"
near $a1 $a2 1e-9 "-beta current independent of the step size"
check [expr {abs($c2-$a2) < abs($c1-$a2)}] "original update approaches -beta current with smaller steps"

puts "PASS JiangJianjingConcrete Tcl checks: $checks"
wipe
